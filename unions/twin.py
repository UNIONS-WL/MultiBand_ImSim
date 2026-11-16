"""UNIONS twin image simulations.

Re-simulate a real ShapePipe campaign with known shear. Every CCD of every real
exposure keeps its WCS (overscan trimmed), zero point, weight and flag maps,
and gets SKiLLS galaxies on a sky lattice, built by MultiBand_ImSim's own
galaxy construction, sheared, convolved with the CCD's PSF (a Moffat matched to
the CCD's v2.0.0 stars, or the PSFEx model ShapePipe fitted to the real
exposure; see the PSF section), plus Gaussian noise at the CCD's measured sky
rms on zero background. The noise depends on exposure and CCD only, so all shear branches
share it. Exposures are coadded onto the real tiles' pixel grids with SWarp,
as the real tiles are. The output is ShapePipe's image-simulation input layout
(see README.md).

    python -m unions.twin unions/config.yaml [-j N] [--tiles 205.239 206.239]
        [--shard I N]

`--shard I N` draws only every N-th exposure (offset I) and stops; shards can
run at once on separate machines, and a final run without it coadds the tiles
and writes the PSF dictionary.
"""

import argparse
import glob
import os
import pickle
import re
import shutil
import subprocess
import tarfile
import warnings
from multiprocessing import Pool
from pathlib import Path
from types import SimpleNamespace

import galsim
import numpy as np
import pandas as pd
import pyarrow.dataset as pads
import pyarrow.feather as feather
import yaml
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning
from shapepipe.modules.psfex_interp_package.psfex_interp import (
    PSFExInterpolator,
    _hsm_columns,
    _hsm_row,
)

from modules.ImSimObject import (
    Q_MAX,
    Q_MIN,
    SERSIC_N_MAX,
    SERSIC_N_MIN,
    TRUNC_FACTOR,
    _build_galaxy_model,
)

warnings.simplefilter("ignore", FITSFixedWarning)

N_CCD = 40
CCD_NX, CCD_NY = 2048, 4612
OVERSCAN_X = 32  # DATASEC [33:2080,1:4612] of the raw 2112x4644 CCDs
STAMP_MAX = 256  # galaxy stamp size cap, pixels


# --- configuration and layout ---------------------------------------------

def branch_shear(name, shear):
    """(g1, g2) of a branch named 1{z,p,m}2{z,p,m}."""
    sign = {"z": 0.0, "p": 1.0, "m": -1.0}
    return sign[name[1]] * shear, sign[name[3]] * shear


class Run:
    """Paths and campaign inputs of one run."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.campaign = Path(cfg["campaign_dir"])
        self.root = Path(cfg["output_root"])
        tag = f"{cfg['layout']}_{cfg['sim_num']}"
        self.shared = self.root / f"{tag}_shared"
        self.branches = {b: branch_shear(b, cfg["shear"]) for b in cfg["branches"]}
        self.sim = {b: self.root / f"{b}_{tag}" for b in cfg["branches"]}
        # tile ID -> (real header, exposures it coadds)
        self.tiles = {}
        for path in sorted(glob.glob(str(self.campaign / "tiles/CFIS.*.r.fits"))):
            tile = re.search(r"CFIS\.(\d{3}\.\d{3})\.r\.fits", path).group(1)
            header = fits.getheader(path)
            history = [str(h).strip() for h in header["HISTORY"]]
            exps = [re.search(r"(\d{7})p\.fits", h).group(1) for h in history]
            self.tiles[tile] = (header, exps)

    def exp_path(self, exp, kind=""):
        return self.campaign / "exposures" / f"{exp}p{kind}.fits.fz"


def tile_name(tile):
    return tile.replace(".", "-")


# --- truth catalogue -----------------------------------------------------------

def footprint(run, scale):
    """Tangent-plane WCS at the campaign centre with pixel `scale` (deg), each
    tile's (WCS, nx, ny), and the footprint's bounding box in that plane."""
    centres = []
    walls = {}
    for tile, (header, _) in run.tiles.items():
        w = WCS(header)
        walls[tile] = (w, header["NAXIS1"], header["NAXIS2"])
        centres.append(w.all_pix2world([[header["NAXIS1"] / 2, header["NAXIS2"] / 2]], 1)[0])
    centres = np.array(centres)
    lat = WCS(naxis=2)
    lat.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    lat.wcs.crval = [np.median(centres[:, 0]), np.median(centres[:, 1])]
    lat.wcs.crpix = [0.0, 0.0]
    lat.wcs.cdelt = [-scale, scale]
    corners = []
    for w, nx, ny in walls.values():
        c = w.all_pix2world([[0.5, 0.5], [nx + 0.5, 0.5], [0.5, ny + 0.5], [nx + 0.5, ny + 0.5]], 1)
        corners.append(lat.all_world2pix(c, 1))
    corners = np.concatenate(corners)
    return lat, walls, np.floor(corners.min(0)), np.ceil(corners.max(0))


def in_footprint(walls, ra, dec):
    keep = np.zeros(np.size(ra), bool)
    for w, nx, ny in walls.values():
        x, y = w.all_world2pix(ra, dec, 1)
        keep |= (x > 0.5) & (x < nx + 0.5) & (y > 0.5) & (y < ny + 0.5)
    return keep


def lattice(run):
    """Sky lattice over the union of the tiles, as (ra, dec) in degrees."""
    lat, walls, lo, hi = footprint(run, run.cfg["grid_spacing_arcsec"] / 3600)
    i, j = np.meshgrid(np.arange(lo[0], hi[0] + 1), np.arange(lo[1], hi[1] + 1))
    ra, dec = lat.all_pix2world(i.ravel(), j.ravel(), 1)
    keep = in_footprint(walls, ra, dec)
    return ra[keep], dec[keep]


def transplant(run, ra, dec):
    """Move sky positions from their own patch onto the campaign: the offsets
    in the tangent plane at the patch centre become offsets in the tangent
    plane at the campaign centre, so separations are kept. Returns the moved
    (ra, dec) and which of them fall in the tiles."""
    lat, walls, _, _ = footprint(run, 1 / 3600)
    src = lat.deepcopy()
    src.wcs.crval = [0.5 * (ra.min() + ra.max()), 0.5 * (dec.min() + dec.max())]
    x, y = src.all_world2pix(ra, dec, 1)
    ra2, dec2 = lat.all_pix2world(x, y, 1)
    return ra2, dec2, in_footprint(walls, ra2, dec2)


def read_skills(cfg, cols):
    parts = []
    for path in sorted(glob.glob(cfg["skills"])):
        t = feather.read_table(path, columns=cols).to_pandas()
        mag = t[cfg["mag_column"]]
        parts.append(t[(mag >= cfg["mag_range"][0]) & (mag < cfg["mag_range"][1])])
    return pd.concat(parts, ignore_index=True)


def make_truth(run):
    """Input galaxies. `layout: grid`: lattice positions with SKiLLS galaxies
    drawn at random and random position angles. `layout: skills`: every SKiLLS
    galaxy in `mag_range` at its own position and position angle, its patch
    moved onto the campaign (`transplant`), so blending is SKiLLS's."""
    cfg = run.cfg
    cols = ["index", "Re_arcsec", "BA", "shape/sersic_n", cfg["mag_column"]]
    rng = np.random.default_rng([cfg["seed"], 0])
    if cfg["layout"] == "grid":
        ra, dec = lattice(run)
        skills = read_skills(cfg, cols)
        pick = skills.iloc[rng.choice(len(skills), size=ra.size, replace=False)]
        pa = rng.uniform(0, 180, ra.size)
    else:
        skills = read_skills(cfg, cols + ["ra", "dec", "PA_random"])
        ra, dec, on = transplant(run, skills["ra"].to_numpy(), skills["dec"].to_numpy())
        pick, ra, dec = skills[on], ra[on], dec[on]
        pa = pick["PA_random"].to_numpy()
    q = np.clip(pick["BA"].to_numpy(), Q_MIN, Q_MAX)
    g = (1 - q) / (1 + q)
    mag = pick[cfg["mag_column"]].to_numpy()
    return Table(
        {
            "ID": np.arange(ra.size),
            "RA": ra,
            "DEC": dec,
            "SKILLS_INDEX": pick["index"].to_numpy().astype(np.int64),
            "MAG_R": mag,
            "FLUX_ZP30": 10 ** (-0.4 * (mag - 30)),
            "SERSIC_N": np.clip(pick["shape/sersic_n"].to_numpy(), SERSIC_N_MIN, SERSIC_N_MAX),
            "RE_ARCSEC": pick["Re_arcsec"].to_numpy(),
            "HLR_ARCSEC": pick["Re_arcsec"].to_numpy() * np.sqrt(q),
            "Q": q,
            "PA_DEG": pa,
            "G1_INT": g * np.cos(np.radians(2 * pa)),
            "G2_INT": g * np.sin(np.radians(2 * pa)),
        },
        meta={
            "LAYOUT": cfg["layout"],
            "MAGCOL": cfg["mag_column"],
            "SEED": cfg["seed"],
            "COMMENT": [
                "Shapes are reduced shears g in GalSim's local sky frame:"
                " +u West, +v North; PA_DEG is measured from +u toward +v.",
                "Galaxies are built by MultiBand_ImSim's _build_galaxy_model:"
                " HLR_ARCSEC = RE_ARCSEC * sqrt(Q) is the circularised Sersic"
                f" half-light radius, truncated at {TRUNC_FACTOR} HLR; FLUX_ZP30 is"
                " the untruncated profile's flux, so the drawn flux is the part"
                " inside the truncation.",
            ],
        },
    )


def make_stars(run):
    """Input stars (`star_population:` in the config): a stellar population catalogue
    without positions (SKiLLS's TRILEGAL files, `..._<N>deg2_...`), placed
    uniformly at random over the tiles at the catalogue's density, each star
    a catalogue row drawn at random, in `mag_range` of the catalogue's r."""
    cfg = run.cfg["star_population"]
    cat = feather.read_table(cfg["catalogue"]).to_pandas()
    area = float(re.search(r"_([\d.]+)deg2_", Path(cfg["catalogue"]).name).group(1))
    cat = cat[(cat["r"] >= cfg["mag_range"][0]) & (cat["r"] < cfg["mag_range"][1])]
    lat, walls, lo, hi = footprint(run, 1 / 3600)
    box = np.prod(hi - lo) / 3600 ** 2  # deg^2, tangent plane
    rng = np.random.default_rng([run.cfg["seed"], 1])
    n = rng.poisson(len(cat) / area * box)
    x, y = rng.uniform(lo[0], hi[0], n), rng.uniform(lo[1], hi[1], n)
    ra, dec = lat.all_pix2world(x, y, 1)
    on = in_footprint(walls, ra, dec)
    pick = cat.iloc[rng.choice(len(cat), size=on.sum())]
    mag = pick["r"].to_numpy()
    return Table({"ID": np.arange(on.sum()), "RA": ra[on], "DEC": dec[on],
                  "CATALOGUE_INDEX": pick["index"].to_numpy().astype(np.int64),
                  "MAG_R": mag, "FLUX_ZP30": 10 ** (-0.4 * (mag - 30))},
                 meta={"CATALOG": Path(cfg["catalogue"]).name, "SEED": run.cfg["seed"]})


# --- PSF --------------------------------------------------------------------
#
# Each CCD's PSF is a CcdPsf: the profile to convolve a galaxy with at a CCD
# position, the stamp ShapePipe is handed as the truth (`psf_dict`), and the
# method that draws the convolution. Two models (`psf.model` in the config):
#
# - moffat: an elliptical Moffat, constant in the CCD's pixel frame, matched at
#   the CCD centre to the median HSM T, g1, g2 of the CCD's v2.0.0 stars, one
#   beta for the campaign from the stars' median rho4. Drawn with the pixel.
# - psfex: the PSFEx model ShapePipe fitted to the real exposure's stars on that
#   CCD. Its basis already holds the pixel response (PSF_SAMP = 1), so it is
#   drawn without one. `vary: false` takes the model at the CCD centre as the
#   PSF everywhere on the CCD (the dictionary stamp is then exact); `vary: true`
#   evaluates it at each galaxy, as ShapePipe's multi-epoch psfex_interp does,
#   and ShapePipe is handed the models themselves (`{shared}/psfex/`).

PSFEX_PATTERN = "star_split_ratio_80"  # ShapePipe's ME_DOT_PSF_PATTERN
PSFEX_STAR_THRESH, PSFEX_CHI2_THRESH = 22, 2  # ShapePipe's v2 psfex_interp gate


def star_targets(run, exps):
    """Per-(exposure, CCD) median star T, g1, g2, and the campaign median rho4."""
    d = pads.dataset(run.cfg["star_population"], format="parquet")
    cols = ["EXPOSURE", "CCD_NB", "HSM_T_STAR", "HSM_G1_STAR", "HSM_G2_STAR",
            "HSM_RHO4_STAR", "HSM_FLAG_STAR"]
    s = d.to_table(columns=cols, filter=pads.field("EXPOSURE").isin([int(e) for e in exps]))
    s = s.to_pandas()
    s = s[s.HSM_FLAG_STAR == 0]
    med = s.groupby(["EXPOSURE", "CCD_NB"])[["HSM_T_STAR", "HSM_G1_STAR", "HSM_G2_STAR"]].median()
    med["N_STAR"] = s.groupby(["EXPOSURE", "CCD_NB"]).size()
    return med, float(s.HSM_RHO4_STAR.median())


def moffat_beta(rho4, scale=0.187):
    """Moffat beta whose HSM rho4 equals the stars' median rho4."""
    betas = np.arange(1.8, 8.0, 0.05)
    r = [galsim.hsm.FindAdaptiveMom(
        galsim.Moffat(beta=b, fwhm=0.7).drawImage(nx=51, ny=51, scale=scale)).moments_rho4
        for b in betas]
    return float(np.interp(rho4, r[::-1], betas[::-1]))


def psf_profile(beta, fwhm, g1, g2):
    return galsim.Moffat(beta=beta, fwhm=fwhm).shear(g1=g1, g2=g2)


def draw_psf(psf, local, size):
    """PSF stamp on the CCD pixel grid, unit sum, centred on the central pixel."""
    stamp = psf.drawImage(nx=size, ny=size, wcs=local).array
    return stamp / stamp.sum()


def hsm_sky(stamp, local):
    """HSM moments in the sky frame, as ShapePipe measures PSFs and stars."""
    return galsim.hsm.FindAdaptiveMom(
        galsim.Image(stamp, wcs=local), strict=False, use_sky_coords=True)


def fit_psf(target, beta, local, size):
    """Moffat FWHM and shear whose HSM T, g1, g2 match the stars' medians."""
    T, g1t, g2t = target
    fwhm, g1, g2 = 0.7, g1t, g2t
    for _ in range(8):
        m = hsm_sky(draw_psf(psf_profile(beta, fwhm, g1, g2), local, size), local)
        fwhm *= np.sqrt(T / (2 * m.moments_sigma ** 2))
        g1 += g1t - m.observed_shape.g1
        g2 += g2t - m.observed_shape.g2
    return fwhm, g1, g2


def psfex_name(exp, ccd):
    return f"{PSFEX_PATTERN}-{int(exp)}-{ccd}.psf"


def read_psfex(source, exp, ccd):
    """A real exposure's PSFEx model for one CCD, from the exposure's persisted
    PSF tar (`{source}/{exp}.tar`, as ShapePipe's exp_persist writes it)."""
    with tarfile.open(Path(source) / f"{int(exp)}.tar") as t:
        return fits.HDUList.fromstring(t.extractfile(psfex_name(exp, ccd)).read())


def psfex_in_sim_frame(hdul):
    """The model with its polynomial origin moved from the real CCD's pixel
    frame (overscan included, where ShapePipe fitted it) to the trimmed CCD's."""
    hdul[1].header["POLZERO1"] -= OVERSCAN_X
    return hdul


def psfex_stamps(path, x, y):
    """The PSFEx model in `path` at pixels (x, y) of its frame, by ShapePipe's
    own interpolation (PSFExInterpolator.interpsfex), ungated, unit sum."""
    gate = SimpleNamespace(_star_thresh=0, _chi2_thresh=np.inf, _w_log=None)
    stamps = PSFExInterpolator.interpsfex(gate, str(path), np.column_stack([x, y]))
    return stamps / stamps.sum(axis=(1, 2), keepdims=True)


def pixel_profile(stamp):
    """A stamp as a profile in pixel units, centred on its central pixel."""
    return galsim.InterpolatedImage(galsim.Image(stamp, scale=1.0), normalization="flux")


class CcdPsf:
    """The PSF of one CCD: the dictionary stamp, the drawing method, metadata
    for the CCD table, and `profiles`, the PSF at image positions in sky units."""

    def __init__(self, stamp, method, row, pixel=None, model=None):
        self.stamp, self.method, self.row = stamp, method, row
        self._pixel, self._model = pixel, model

    def profiles(self, x, y, locals_):
        """PSF profiles at 1-based image positions (x, y), through each
        position's local WCS."""
        if self._model is None:
            return [lw.toWorld(self._pixel) for lw in locals_]
        # ShapePipe's multi-epoch interpolation evaluates at 0-based pixels
        stamps = psfex_stamps(self._model, np.asarray(x) - 1, np.asarray(y) - 1)
        return [lw.toWorld(pixel_profile(s)) for s, lw in zip(stamps, locals_)]


def ccd_psf(run, exp, ccd, inputs, gwcs, centre):
    """This CCD's PSF under the run's `psf` config."""
    cfg = run.cfg["psf"]
    local = gwcs.local(centre)
    size = run.cfg["psf_stamp"]
    if cfg["model"] == "moffat":
        target, beta = inputs
        fwhm, pg1, pg2 = fit_psf(target, beta, local, size)
        psf = psf_profile(beta, fwhm, pg1, pg2)
        row = {"PSF_BETA": beta, "PSF_FWHM": fwhm, "PSF_G1_IN": pg1, "PSF_G2_IN": pg2,
               "STAR_T": target[0], "STAR_G1": target[1], "STAR_G2": target[2]}
        return CcdPsf(draw_psf(psf, local, size), "auto", row, pixel=local.toImage(psf))
    hdul = psfex_in_sim_frame(read_psfex(cfg["psfex_models"], exp, ccd))
    path = run.shared / "psfex" / psfex_name(exp, ccd)
    path.parent.mkdir(exist_ok=True)
    hdul.writeto(path, overwrite=True)
    h = hdul[1].header
    row = {"PSFEX_ACCEPTED": h["ACCEPTED"], "PSFEX_CHI2": h["CHI2"],
           "PSFEX_GATE": h["ACCEPTED"] >= PSFEX_STAR_THRESH and h["CHI2"] <= PSFEX_CHI2_THRESH}
    stamp = psfex_stamps(path, [centre.x - 1], [centre.y - 1])[0]
    if cfg["vary"]:
        return CcdPsf(stamp, "no_pixel", row, model=path)
    return CcdPsf(stamp, "no_pixel", row, pixel=pixel_profile(stamp))


def psf_inputs(run, exps):
    """Per-(exposure, CCD) inputs that ccd_psf needs beyond the CCD itself."""
    if run.cfg["psf"]["model"] != "moffat":
        return {(e, c): None for e in exps for c in range(N_CCD)}
    targets, rho4 = star_targets(run, sorted({e for _, es in run.tiles.values() for e in es}))
    beta = moffat_beta(rho4)
    print(f"stars: median rho4 {rho4:.4f} -> Moffat beta {beta:.3f}", flush=True)
    return {(e, c): (tuple(targets.loc[(int(e), c)].iloc[:3]), beta)
            for e in exps for c in range(N_CCD)}


# --- exposures --------------------------------------------------------------

def ccd_header(real):
    """Real CCD header with the overscan trimmed off."""
    h = real.copy()
    for key in ("BZERO", "BSCALE", "BLANK"):
        h.remove(key, ignore_missing=True)
    h["CRPIX1"] = h["CRPIX1"] - OVERSCAN_X
    h["DATASEC"] = f"[1:{CCD_NX},1:{CCD_NY}]"
    h["TRIMSEC"] = f"[1:{CCD_NX},1:{CCD_NY}]"
    return h


def trim(a):
    return a[:CCD_NY, OVERSCAN_X:OVERSCAN_X + CCD_NX]


def sky_rms(data, good, mesh=128):
    """Sky noise rms: mesh-median background removed, 3-sigma clipped std."""
    ny, nx = (s // mesh * mesh for s in data.shape)
    d = np.where(good, data, np.nan)[:ny, :nx]
    bkg = np.nanmedian(d.reshape(ny // mesh, mesh, nx // mesh, mesh), axis=(1, 3))
    r = d - np.repeat(np.repeat(bkg, mesh, 0), mesh, 1)
    r = r[np.isfinite(r)]
    sigma = 1.4826 * np.median(np.abs(r - np.median(r)))
    for _ in range(5):
        c = r[np.abs(r - np.median(r)) < 3 * sigma]
        sigma = c.std() / 0.986572  # std of a Gaussian clipped at 3 sigma
    return float(sigma)


def gal_info(row, flux):
    """A truth row as MultiBand_ImSim's galaxy record (`_build_galaxy_model`)."""
    return {"r": flux, "position_angle": row["PA_DEG"], "sersic_n": row["SERSIC_N"],
            "Re": row["RE_ARCSEC"], "axis_ratio": row["Q"]}


TRUTH = STARS = None  # set before the worker pool forks


def draw_ccd(task):
    """Draw one CCD in every branch; write image, weight and flag arrays."""
    run, exp, ccd, inputs = task
    work = run.shared / "work" / exp
    with fits.open(run.exp_path(exp)) as f:
        real = f[ccd + 1]
        header = ccd_header(real.header)
        data = trim(real.data).astype(np.float32)
    weight = trim(fits.getdata(run.exp_path(exp, ".weight"), ccd + 1))
    flag = trim(fits.getdata(run.exp_path(exp, ".flag"), ccd + 1))
    rms = sky_rms(data, (weight > 0) & (flag == 0))

    wcs = WCS(header)
    gwcs = galsim.GSFitsWCS(header=header)
    zp = header["PHOTZP"]

    centre = galsim.PositionD((CCD_NX + 1) / 2, (CCD_NY + 1) / 2)
    local = gwcs.local(centre)
    psf = ccd_psf(run, exp, ccd, inputs, gwcs, centre)
    shapes = _hsm_columns(_hsm_row(psf.stamp, hsm_sky(psf.stamp, local), local), "PSF")

    def on_ccd(table):
        """Rows whose stamp may touch the CCD, their pixel positions and PSFs."""
        ra0, dec0 = wcs.all_pix2world([[centre.x, centre.y]], 1)[0]
        near = (np.abs(table["DEC"] - dec0) < 0.16) & (
            np.abs((table["RA"] - ra0) * np.cos(np.radians(dec0))) < 0.16)
        t = table[near]
        x, y = wcs.all_world2pix(t["RA"], t["DEC"], 1)
        m = STAMP_MAX / 2
        on = (x > -m) & (x < CCD_NX + m) & (y > -m) & (y < CCD_NY + m)
        t, x, y = t[on], x[on], y[on]
        locals_ = [gwcs.local(galsim.PositionD(xi, yi)) for xi, yi in zip(x, y)]
        return t, x, y, locals_, psf.profiles(x, y, locals_)

    cand, x, y, locals_, psfs = on_ccd(TRUTH)
    canvases = {b: galsim.ImageF(CCD_NX, CCD_NY) for b in run.branches}
    for row, xi, yi, lw, psf_here in zip(cand, x, y, locals_, psfs):
        pos = galsim.PositionD(xi, yi)
        info = gal_info(row, 10 ** (-0.4 * (row["MAG_R"] - zp)))
        n = None
        for b, g in run.branches.items():
            obj = _build_galaxy_model(info, "r", 0.0, g, True, psf_here)
            n = n or min(obj.getGoodImageSize(lw.minLinearScale()), STAMP_MAX)
            st = obj.drawImage(nx=n, ny=n, wcs=lw, center=pos, method=psf.method)
            overlap = st.bounds & canvases[b].bounds
            if overlap.isDefined():
                canvases[b][overlap] += st[overlap]

    stars = galsim.ImageF(CCD_NX, CCD_NY)
    n_star = 0
    if STARS is not None:
        cand_s, xs, ys, locals_s, psfs_s = on_ccd(STARS)
        n_star = len(cand_s)
        for row, xi, yi, lw, psf_here in zip(cand_s, xs, ys, locals_s, psfs_s):
            star = psf_here.withFlux(10 ** (-0.4 * (row["MAG_R"] - zp)))
            n = min(star.getGoodImageSize(lw.minLinearScale()), STAMP_MAX)
            st = star.drawImage(nx=n, ny=n, wcs=lw, center=galsim.PositionD(xi, yi),
                                method=psf.method)
            overlap = st.bounds & stars.bounds
            if overlap.isDefined():
                stars[overlap] += st[overlap]

    # the same in every branch: noise, and the stars
    common = np.random.default_rng([run.cfg["seed"], int(exp), ccd]).normal(
        0.0, rms, (CCD_NY, CCD_NX)).astype(np.float32) + stars.array
    work.mkdir(parents=True, exist_ok=True)
    for b in run.branches:
        np.save(work / f"{ccd:02d}-{b}.npy", canvases[b].array + common)
    np.save(work / f"{ccd:02d}-weight.npy", weight)
    np.save(work / f"{ccd:02d}-flag.npy", flag)
    return {"EXPOSURE": int(exp), "CCD": ccd, "RMS": rms, "N_GAL": len(cand), "N_STAR": n_star,
            "PSF_MODEL": run.cfg["psf"]["model"], **psf.row,
            "HEADER": header.tostring(), "VIGNET": psf.stamp, "SHAPES": shapes}


def exp_file(kind, exp):
    return f"simu_{kind}-{exp}.fits"


def assemble(task):
    """Write one exposure's per-branch images and shared weight and flag."""
    run, exp, rows = task
    work = run.shared / "work" / exp
    primary = fits.getheader(run.exp_path(exp), 0)
    headers = [fits.Header.fromstring(r["HEADER"]) for r in sorted(rows, key=lambda r: r["CCD"])]

    def write(path, part):
        hdus = [fits.PrimaryHDU(header=primary)]
        for ccd, h in enumerate(headers):
            hdus.append(fits.ImageHDU(np.load(work / f"{ccd:02d}-{part}.npy"), header=h))
        tmp = path.with_suffix(".tmp")
        fits.HDUList(hdus).writeto(tmp, overwrite=True)
        os.replace(tmp, path)

    shared_exp = run.shared / "SP_exp"
    shared_exp.mkdir(parents=True, exist_ok=True)
    for kind in ("weight", "flag"):
        write(shared_exp / exp_file(kind, exp), kind)
    for b, sim in run.sim.items():
        d = sim / "images" / "SP_exp"
        d.mkdir(parents=True, exist_ok=True)
        write(d / exp_file("image", exp), b)
        for kind in ("weight", "flag"):
            link(shared_exp / exp_file(kind, exp), d / exp_file(kind, exp))
    with open(run.shared / "psf" / f"{exp}.pickle", "wb") as f:
        pickle.dump([{k: v for k, v in r.items() if k != "HEADER"} for r in rows], f)
    shutil.rmtree(work)


def link(src, dst):
    if not dst.exists():
        os.link(src, dst)


# --- tiles ------------------------------------------------------------------

def coadd(task):
    """SWarp one branch's exposures onto a real tile's pixel grid."""
    run, tile, branch, nthreads = task
    header, exps = run.tiles[tile]
    sim = run.sim[branch]
    out = sim / "images" / "SP_tiles"
    out.mkdir(parents=True, exist_ok=True)
    tmp = run.shared / "work" / "swarp" / f"{tile}-{branch}"
    tmp.mkdir(parents=True, exist_ok=True)
    img = tmp / f"CFIS_simu_image-{tile_name(tile)}.fits"
    wht = tmp / f"CFIS_simu_weight-{tile_name(tile)}.fits"
    geometry = ["NAXIS", "NAXIS1", "NAXIS2", "EQUINOX", "RADESYS", "CTYPE1", "CUNIT1",
                "CRVAL1", "CRPIX1", "CD1_1", "CD1_2", "CTYPE2", "CUNIT2", "CRVAL2",
                "CRPIX2", "CD2_1", "CD2_2"]
    head = fits.Header([(k, header[k]) for k in geometry])
    img.with_suffix(".head").write_text(head.tostring(sep="\n", endcard=True) + "\n")
    sp_exp = sim / "images" / "SP_exp"
    cmd = [
        run.cfg["swarp"],
        ",".join(str(sp_exp / exp_file("image", e)) for e in exps),
        "-IMAGEOUT_NAME", str(img), "-WEIGHTOUT_NAME", str(wht),
        "-WEIGHT_TYPE", "MAP_WEIGHT",
        "-WEIGHT_IMAGE", ",".join(str(sp_exp / exp_file("weight", e)) for e in exps),
        "-COMBINE_TYPE", "WEIGHTED", "-RESAMPLING_TYPE", "LANCZOS3",
        "-SUBTRACT_BACK", "N", "-FSCALE_KEYWORD", "FSCALE", "-GAIN_KEYWORD", "GAIN",
        "-COPY_KEYWORDS", "OBJECT,FILTER",
        "-RESAMPLE_DIR", str(tmp), "-VMEM_DIR", str(tmp),
        "-MEM_MAX", "4096", "-COMBINE_BUFSIZE", "4096",
        "-NTHREADS", str(nthreads), "-WRITE_XML", "N", "-VERBOSE_TYPE", "QUIET",
    ]
    subprocess.run(cmd, check=True, cwd=tmp)
    with fits.open(img, mode="update") as f:
        h = f[0].header
        h.remove("HISTORY", ignore_missing=True, remove_all=True)
        for line in header["HISTORY"]:
            h.add_history(re.sub(r"(\d{7})p\.fits", r"simu_image-\1.fits", str(line).strip()))
    os.replace(wht, out / wht.name)
    os.replace(img, out / img.name)
    shutil.rmtree(tmp)


# --- driver -----------------------------------------------------------------

def done_exposure(run, exp):
    return all((s / "images" / "SP_exp" / exp_file("image", exp)).exists()
               for s in run.sim.values())


def done_tile(run, tile):
    return all((s / "images" / "SP_tiles" / f"CFIS_simu_weight-{tile_name(tile)}.fits").exists()
               for s in run.sim.values())


def finish(run):
    """PSF dictionary, CCD table, tile lists and per-sim metadata from what exists."""
    rows = []
    for path in sorted((run.shared / "psf").glob("*.pickle")):
        with open(path, "rb") as f:
            rows += pickle.load(f)
    psf_dict = {f"{r['EXPOSURE']}-{r['CCD']}": {"VIGNET": r["VIGNET"], "SHAPES": r["SHAPES"]}
                for r in rows}
    tmp = Path(run.cfg["psf_dict"] + ".tmp")
    with open(tmp, "wb") as f:
        pickle.dump(psf_dict, f)
    os.replace(tmp, run.cfg["psf_dict"])
    table = pd.DataFrame([{k: v for k, v in r.items() if k not in ("VIGNET", "SHAPES")}
                          | {k: r["SHAPES"][k] for k in ("HSM_T_PSF", "HSM_G1_PSF", "HSM_G2_PSF")}
                          for r in rows])
    table.to_parquet(run.shared / "ccds.parquet")

    tiles = [t for t in run.tiles if done_tile(run, t)]
    for b, sim in run.sim.items():
        sim.mkdir(parents=True, exist_ok=True)
        (sim / "tiles.txt").write_text("".join(f"{t}\n" for t in tiles))
        g1, g2 = run.branches[b]
        (sim / "basic_info.txt").write_text(
            f"g_cosmic = {g1} {g2}\n"
            f"branch = {b}\n"
            f"campaign_dir = {run.campaign}\n"
            f"psf_dict = {run.cfg['psf_dict']}\n")
        cat = sim / "catalogues" / "input"
        cat.mkdir(parents=True, exist_ok=True)
        link(run.shared / "truth.fits", cat / "truth.fits")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("config")
    p.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    p.add_argument("--tiles", nargs="+", help="build only these tiles and their exposures")
    p.add_argument("--shard", nargs=2, type=int, metavar=("I", "N"),
                   help="draw only exposures I, I+N, ... and stop before the tiles")
    args = p.parse_args()

    with open(args.config) as f:
        run = Run(yaml.safe_load(f))
    tiles = args.tiles or list(run.tiles)
    exps = sorted({e for t in tiles for e in run.tiles[t][1]})
    run.shared.mkdir(parents=True, exist_ok=True)
    (run.shared / "psf").mkdir(exist_ok=True)

    global TRUTH, STARS
    truth_path = run.shared / "truth.fits"
    if not truth_path.exists():
        make_truth(run).write(str(truth_path) + ".tmp", format="fits")
        os.replace(str(truth_path) + ".tmp", truth_path)
    TRUTH = Table.read(truth_path)
    print(f"truth: {len(TRUTH)} galaxies", flush=True)
    if "star_population" in run.cfg:
        stars_path = run.shared / "stars.fits"
        if not stars_path.exists():
            make_stars(run).write(str(stars_path) + ".tmp", format="fits")
            os.replace(str(stars_path) + ".tmp", stars_path)
        STARS = Table.read(stars_path)
        print(f"stars: {len(STARS)}", flush=True)

    if args.shard:
        exps = exps[args.shard[0]::args.shard[1]]
    todo = [e for e in exps if not done_exposure(run, e)]
    if todo:
        inputs = psf_inputs(run, todo)
        tasks = [(run, e, c, inputs[(e, c)]) for e in todo for c in range(N_CCD)]
        rows = {e: [] for e in todo}
        with Pool(args.jobs) as pool:
            for r in pool.imap_unordered(draw_ccd, tasks):
                rows[f"{r['EXPOSURE']}"].append(r)
                print(f"drew {r['EXPOSURE']}-{r['CCD']}: {r['N_GAL']} galaxies, "
                      f"rms {r['RMS']:.2f}", flush=True)
            pool.map(assemble, [(run, e, rows[e]) for e in todo])
        print(f"exposures written: {todo}", flush=True)
    if args.shard:
        return

    todo = [(t, b) for t in tiles if not done_tile(run, t) for b in run.branches]
    if todo:
        nthreads = 4
        with Pool(max(1, args.jobs // nthreads)) as pool:
            pool.map(coadd, [(run, t, b, nthreads) for t, b in todo], chunksize=1)
        print(f"tiles written: {sorted({t for t, _ in todo})}", flush=True)

    finish(run)


if __name__ == "__main__":
    main()
