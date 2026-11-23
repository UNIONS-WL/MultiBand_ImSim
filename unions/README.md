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
- **PSF** (`psf` in the config), drawn into the images and handed to
  ShapePipe as the truth:
  - `model: moffat`: an elliptical Moffat, constant over the CCD in its pixel
    frame, whose HSM size T and ellipticity (g1, g2) at the CCD centre,
    measured in the sky frame as ShapePipe measures them, equal the median of
    that CCD's stars in the v2.0.0 star catalogue (`HSM_*_STAR`,
    `HSM_FLAG_STAR == 0`). Being constant in pixels, the one dictionary stamp
    is the PSF at every position on the CCD; on the sky it varies with the
    WCS's local distortion (≤ 0.8 % in scale, ≤ 0.004 in shear between a CCD's
    centre and corners). One Moffat β for the campaign, chosen so the Moffat's
    HSM radial fourth moment ρ4 equals the stars' median `HSM_RHO4_STAR`.
    Drawn with the pixel response.
  - `model: psfex`: the PSFEx model ShapePipe fitted to the real exposure's
    stars on that CCD (`psfex_models/<exp>.tar`, the exposure PSF tar
    ShapePipe's `exp_persist` keeps), moved from the real CCD's pixel frame
    (overscan included) to the trimmed one and written to
    `{shared}/psfex/star_split_ratio_80-<exp>-<ccd>.psf`. It is evaluated by
    ShapePipe's own `PSFExInterpolator.interpsfex`, so at the v2.0.0 stars it
    reproduces the star catalogue's `HSM_*_PSF` exactly. Its basis holds the
    pixel response (`PSF_SAMP = 1`), so galaxies are convolved with it and
    drawn without another pixel. `vary: false` takes the model at the CCD
    centre as the PSF everywhere on the CCD, so the dictionary stamp is exact;
    `vary: true` evaluates it at each galaxy (0-based pixel, as ShapePipe's
    multi-epoch interpolation does), and ShapePipe is then given the models in
    `{shared}/psfex/` instead of the dictionary. Models failing ShapePipe's v2
    gate (fewer than 22 accepted stars or χ² > 2) are still drawn; the CCD
    table records `PSFEX_GATE`.
- **Galaxies** (`layout`), each built by MultiBand_ImSim's own
  `_build_galaxy_model` (`modules/ImSimObject.py`): Sérsic n, half-light
  radius √q·Re, axis ratio q, truncated at 5 half-light radii with the
  catalogue flux assigned to the untruncated profile. Flux on a CCD is
  10^(−0.4 (mag − PHOTZP)) ADU.
  - `layout: grid`: a square lattice on the sky (tangent plane at the campaign
    centre) over the tiles' footprint, each point a SKiLLS galaxy drawn at
    random from those with `mag_column` in `mag_range`, at a random position
    angle. Isolated galaxies.
  - `layout: skills`: every SKiLLS galaxy in `mag_range` at its own position
    and position angle, the SKiLLS patch moved onto the campaign by keeping
    each galaxy's offset in the tangent plane, so the blending is SKiLLS's;
    kept wherever the exposures reach, not only in the tiles.
- **Tiles** (`tiles:`, optional): simulate only these tiles of the campaign and
  the exposures they coadd.
- **Stars** (`star_population`, optional): a stellar population catalogue without
  positions (SKiLLS's TRILEGAL files), stars placed uniformly at random over
  the exposures at the catalogue's density, each a catalogue row drawn at random
  in `mag_range`, drawn as the PSF and identical in every shear branch. With
  stars, ShapePipe can fit PSFEx to the simulation as it does to the data.
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
    truth.fits, stars.fits (with `star_population`), SP_exp/ (weights, flags), psf/, ccds.parquet (per-CCD rms and PSF),
    psfex/ (psf.model psfex: each CCD's PSFEx model in the trimmed frame)
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

## What an m from these sims covers

The configs are rungs of one ladder on the same campaign, noise and shears,
each changing one ingredient:

| config | galaxies | PSF | what its m covers |
|---|---|---|---|
| `config.yaml` | lattice, r 18.5–26 | Moffat per CCD, true PSF given | shape measurement and selection on isolated galaxies |
| `config_psfex.yaml` | same galaxies | real PSFEx model at the CCD centre, given | + the real PSF's profile |
| `config_psfex_vary.yaml` | same galaxies | real PSFEx model at each galaxy, given | + its variation across the CCD |
| `config_blended.yaml` | SKiLLS positions, r 15–27, + TRILEGAL stars | as above, but ShapePipe fits PSFEx to the simulated stars | + blending, undetected galaxies, PSF-model error |
| `config_grid2_blended.yaml` | as `config_blended.yaml`, on run 2's 25 rung tiles | as above | as above, at run 2's depth and seeing |

The three lattice rungs share galaxies and noise object by object, so their
differences in m can be measured on matched objects with little noise.
