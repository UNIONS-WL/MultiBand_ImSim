# UNIONS twin image simulations

`unions/` simulates a real UNIONS ShapePipe campaign with known shear, for
measuring the multiplicative shear bias m. One config file (`config.yaml`)
describes a run completely; `twin.py` builds it.

## What is simulated

For each real exposure of the campaign and each of its 40 CCDs:

- **Geometry and photometry are real.** The CCD keeps its real header: WCS
  (TAN + PV distortion), `PHOTZP`, `FSCALE`, `GAIN`, `SATURATE`. The overscan
  is trimmed (2048 × 4612 pixels, `CRPIX1` shifted by −32).
- **Weights and flags are real.** The real weight (0/1) and flag maps,
  trimmed, so the twin inherits the real bad pixels, columns and masks.
- **Noise:** Gaussian, zero background, at the rms measured on the real CCD
  (mesh-median sky removed, 3σ-clipped std over good pixels). Seeded by
  exposure and CCD only, so every shear branch has the same noise.
- **PSF:** an elliptical Moffat, constant over the CCD in its pixel frame,
  whose HSM size T and ellipticity (g1, g2) at the CCD centre, measured in the
  sky frame as ShapePipe measures them, equal the median of that CCD's stars
  in the v2.0.0 star catalogue (`HSM_*_STAR`, `HSM_FLAG_STAR == 0`). Being
  constant in pixels, the one dictionary stamp is the PSF at every position on
  the CCD; on the sky it varies with the WCS's local distortion (≤ 0.8 % in
  scale, ≤ 0.004 in shear between a CCD's centre and corners). One Moffat β for the campaign, chosen
  so the Moffat's HSM radial fourth moment ρ4 equals the stars' median
  `HSM_RHO4_STAR`. A CCD with no v2 stars (about 2 % of them) had no PSF
  model in the data: it is drawn with its exposure's median PSF and left out
  of the PSF dictionary, so ShapePipe drops its epochs as it did in the data.
- **Galaxies:** a square lattice on the sky (tangent plane at the campaign
  centre) over the tiles' footprint, each point a SKiLLS galaxy drawn at random
  from those with `mag_column` in `mag_range`: Sérsic n, half-light radius
  √q·Re, axis ratio q from the catalogue (SKiLLS conventions and limits from
  `modules/ImSimObject.py`), a random position angle. Flux on a CCD is
  10^(−0.4 (mag − PHOTZP)) ADU. No stars.
- **Tiles:** the branch's exposures are coadded with SWarp onto each real
  tile's pixel grid (its WCS and size), as the real tiles were made:
  `COMBINE_TYPE WEIGHTED`, `RESAMPLING_TYPE LANCZOS3`, fluxes scaled to zero
  point 30 by `FSCALE`, weights `MAP_WEIGHT`. No background subtraction (the
  sims have none). The tile's HISTORY cards name the same exposures as the
  real tile's, as `input image simu_image-NNNNNNN.fits ...`.

## Output

```
{output_root}/{branch}_{layout}_{sim_num}/       one per shear branch
    images/SP_tiles/CFIS_simu_image-XXX-YYY.fits  coadd, primary HDU, ZP 30
    images/SP_tiles/CFIS_simu_weight-XXX-YYY.fits SWarp weight (inverse variance)
    images/SP_exp/simu_image-NNNNNNN.fits         empty primary + 40 CCD HDUs, float32
    images/SP_exp/simu_weight-NNNNNNN.fits        real weight, int16 (hard link)
    images/SP_exp/simu_flag-NNNNNNN.fits          real flag, uint8 (hard link)
    catalogues/input/truth.fits                   input galaxies (hard link)
    basic_info.txt                                `g_cosmic = g1 g2`, branch, inputs
    tiles.txt                                     tiles built so far, `XXX.YYY`
{output_root}/{layout}_{sim_num}_shared/          branch-independent inputs
    truth.fits, SP_exp/ (weights, flags), psf/, ccds.parquet (per-CCD rms and PSF)
{psf_dict}                                        "<exp>-<ccd>" -> {VIGNET, SHAPES}
```

The PSF dictionary is what ShapePipe's `psf_model: fake` reads: for each
exposure and CCD (0-based, as split by ShapePipe), a 51 × 51 stamp of the true
PSF on that CCD's pixel grid, unit sum, centred on the central pixel, and its
`HSM_{G1,G2,T,FLAG,M4_1,M4_2,RHO4}_PSF` measured by ShapePipe's own
`psfex_interp` functions.

## Running

Inside the sp_validation container (GalSim, astropy, pyarrow, ShapePipe), from
the repository root, on a compute allocation; SWarp from conda-forge
(`micromamba create -p /scratch/cdaley/imsim/env/swarp -c conda-forge
astromatic-swarp`):

```bash
module load apptainer
SIF=/project/rrg-kyi/cdaley/unions/containers/sp_validation-develop.sif
srun --jobid=<allocation> --overlap -n 1 -c 64 \
    apptainer exec --bind /scratch,/project $SIF \
    env OMP_NUM_THREADS=1 python -m unions.twin unions/config.yaml -j 64
```

`--tiles 205.239 206.239` builds only those tiles and the exposures they
coadd. A run skips exposures and tiles that already exist, so a subset run
followed by a full run builds each product once.

Known-answer checks of a built run:

```bash
srun ... python -m unions.checks unions/config.yaml
```

## Rungs: one ingredient changed at a time

Optional config keys each change one ingredient from the twin above; left
out, the twin is as described. The October-ingredient rungs of run 2
(`config_grid2_<rung>.yaml`, sims `*_grid_2-<rung>`) use them on 25 tiles of
run 2 (`tiles:`), with run 2's galaxies and noise seeds, so objects pair
across rungs.

| key | values (default first) | effect |
|---|---|---|
| `tiles` | all the campaign's tiles, or a list | build only these tiles and the exposures they coadd; the lattice still spans the campaign, so the truth table is the full run's |
| `camera` | `real`, `tan` | `tan`: each CCD's WCS reduced to a pure TAN projection, the TPV polynomial's affine part folded into CD and CRPIX and its higher orders (the field distortion) dropped: on run 2 a CCD moves on the sky by up to 14″ at its centre and 27″ at a corner, and its pixel scale changes by up to 1.5 %. Galaxies, PSF fit and dictionary stamp, coadd and ShapePipe all read the written header, so they stay consistent; the PSF is then constant on the sky as well as in pixels |
| `noise_rms` | `ccd`, `exposure` | `exposure`: every CCD gets the median of its exposure's 40 CCD sky rms (same noise seeds, so the noise is each CCD's rescaled) |
| `masks` | `real`, `none` | `none`: unit weight and zero flag on every pixel (the noise rms is still measured on the real good pixels) |
| `coadd` | `weighted`, `median` | SWarp `COMBINE_TYPE` |
| `exposures_from` | unset, or a sim tag | hard-link that sim's exposure images and PSFs instead of drawing them (it must have been drawn with the same galaxies, noise and camera); its weight and flag too unless `masks: none` |

`ccds.parquet` gives each CCD's measured sky rms (`CCD_RMS`) and the noise
rms drawn (`RMS`).

## What an m from these sims covers

A lattice of isolated galaxies brighter than `mag_range`'s faint limit, no
stars, the true PSF handed to the shape measurement, and a PSF constant over
each CCD. The m measured on them is the bias of the shape measurement and of
the selection on isolated galaxies. It leaves out blending, galaxies below the
detection limit, and PSF-model error.
