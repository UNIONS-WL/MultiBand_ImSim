"""UNIONS twin image simulations.

Re-simulate a real ShapePipe campaign with known shear. Every CCD of every real
exposure keeps its WCS (overscan trimmed), zero point, weight and flag maps,
and gets SKiLLS galaxies on a sky lattice, sheared, convolved with an
elliptical Moffat PSF, constant in the CCD's pixel frame and matched at the
CCD centre to that CCD's stars in the v2.0.0 star catalogue, plus Gaussian noise at the CCD's measured sky rms on zero
background. The noise depends on exposure and CCD only, so all shear branches
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
import warnings
from multiprocessing import Pool
from pathlib import Path

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
    _hsm_columns,
    _hsm_row,
)

from modules.ImSimObject import (
    Q_MAX,
    Q_MIN,
    SERSIC_N_MAX,
    SERSIC_N_MIN,
    TRUNC_FACTOR,
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

def lattice(run):
    """Sky lattice over the union of the tiles, as (ra, dec) in degrees."""
    spacing = run.cfg["grid_spacing_arcsec"] / 3600
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
    lat.wcs.cdelt = [-spacing, spacing]
    corners = []
    for w, nx, ny in walls.values():
        c = w.all_pix2world([[0.5, 0.5], [nx + 0.5, 0.5], [0.5, ny + 0.5], [nx + 0.5, ny + 0.5]], 1)
        corners.append(lat.all_world2pix(c, 1))
    corners = np.concatenate(corners)
    lo, hi = np.floor(corners.min(0)), np.ceil(corners.max(0))
    i, j = np.meshgrid(np.arange(lo[0], hi[0] + 1), np.arange(lo[1], hi[1] + 1))
    ra, dec = lat.all_pix2world(i.ravel(), j.ravel(), 1)
    keep = np.zeros(ra.size, bool)
    for w, nx, ny in walls.values():
        x, y = w.all_world2pix(ra, dec, 1)
        keep |= (x > 0.5) & (x < nx + 0.5) & (y > 0.5) & (y < ny + 0.5)
    return ra[keep], dec[keep]


def make_truth(run):
    """Lattice positions with SKiLLS galaxies drawn at random."""
    cfg = run.cfg
    ra, dec = lattice(run)
    cols = ["index", "Re_arcsec", "BA", "shape/sersic_n", cfg["mag_column"]]
    parts = []
    for path in sorted(glob.glob(cfg["skills"])):
        t = feather.read_table(path, columns=cols).to_pandas()
        mag = t[cfg["mag_column"]]
        parts.append(t[(mag >= cfg["mag_range"][0]) & (mag < cfg["mag_range"][1])])
    skills = pd.concat(parts, ignore_index=True)
    rng = np.random.default_rng([cfg["seed"], 0])
    pick = skills.iloc[rng.choice(len(skills), size=ra.size, replace=False)]
    q = np.clip(pick["BA"].to_numpy(), Q_MIN, Q_MAX)
    pa = rng.uniform(0, 180, ra.size)
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
            "SPACING": cfg["grid_spacing_arcsec"],
            "MAGCOL": cfg["mag_column"],
            "SEED": cfg["seed"],
            "COMMENT": [
                "Shapes are reduced shears g in GalSim's local sky frame:"
                " +u West, +v North; PA_DEG is measured from +u toward +v.",
                "HLR_ARCSEC = RE_ARCSEC * sqrt(Q), the circularised Sersic"
                " half-light radius drawn (SKiLLS convention), truncated at"
                f" {TRUNC_FACTOR} HLR; FLUX_ZP30 is the flux drawn inside the truncation.",
            ],
        },
    )


# --- PSF --------------------------------------------------------------------

def star_targets(run, exps):
    """Per-(exposure, CCD) median star T, g1, g2, and the campaign median rho4."""
    d = pads.dataset(run.cfg["stars"], format="parquet")
    cols = ["EXPOSURE", "CCD_NB", "HSM_T_STAR", "HSM_G1_STAR", "HSM_G2_STAR",
            "HSM_RHO4_STAR", "HSM_FLAG_STAR"]
    s = d.to_table(columns=cols, filter=pads.field("EXPOSURE").isin([int(e) for e in exps]))
    s = s.to_pandas()
    s = s[s.HSM_FLAG_STAR == 0]
    med = s.groupby(["EXPOSURE", "CCD_NB"])[["HSM_T_STAR", "HSM_G1_STAR", "HSM_G2_STAR"]].median()
    med["N_STAR"] = s.groupby(["EXPOSURE", "CCD_NB"]).size()
    return med, float(s.HSM_RHO4_STAR.median())


def ccd_target(targets, exp, ccd):
    """(T, g1, g2, has_stars) for one CCD. A CCD with no v2 stars had no PSF
    model in the data and so contributes no epoch to the shape measurement;
    it is drawn with its exposure's median PSF and left out of the PSF
    dictionary, which drops its epochs in ShapePipe as in the data."""
    if (exp, ccd) in targets.index:
        return (*targets.loc[(exp, ccd)].iloc[:3], True)
    return (*targets.loc[exp].iloc[:, :3].median(), False)


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


def galaxy(row):
    """Unsheared SKiLLS Sersic galaxy, in the sky frame."""
    hlr = row["HLR_ARCSEC"]
    gal = galsim.Sersic(n=row["SERSIC_N"], half_light_radius=hlr,
                        trunc=TRUNC_FACTOR * hlr, flux_untruncated=True)
    return gal.shear(q=row["Q"], beta=row["PA_DEG"] * galsim.degrees)


TRUTH = None  # set before the worker pool forks


def draw_ccd(task):
    """Draw one CCD in every branch; write image, weight and flag arrays."""
    run, exp, ccd, target, beta = task
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
    size = run.cfg["psf_stamp"]

    centre = galsim.PositionD((CCD_NX + 1) / 2, (CCD_NY + 1) / 2)
    local = gwcs.local(centre)
    fwhm, pg1, pg2 = fit_psf(target[:3], beta, local, size)
    psf = psf_profile(beta, fwhm, pg1, pg2)
    stamp = draw_psf(psf, local, size)
    # The PSF is constant in the CCD's pixel frame, so the one stamp is the
    # PSF at every position; on the sky it follows the WCS's local distortion.
    psf_pix = local.toImage(psf)
    shapes = _hsm_columns(_hsm_row(stamp, hsm_sky(stamp, local), local), "PSF")

    # galaxies whose stamp may touch the CCD
    ra0, dec0 = wcs.all_pix2world([[centre.x, centre.y]], 1)[0]
    near = (np.abs(TRUTH["DEC"] - dec0) < 0.16) & (
        np.abs((TRUTH["RA"] - ra0) * np.cos(np.radians(dec0))) < 0.16)
    cand = TRUTH[near]
    x, y = wcs.all_world2pix(cand["RA"], cand["DEC"], 1)
    margin = STAMP_MAX / 2
    on = (x > -margin) & (x < CCD_NX + margin) & (y > -margin) & (y < CCD_NY + margin)
    cand, x, y = cand[on], x[on], y[on]

    canvases = {b: galsim.ImageF(CCD_NX, CCD_NY) for b in run.branches}
    for row, xi, yi in zip(cand, x, y):
        pos = galsim.PositionD(xi, yi)
        lw = gwcs.local(pos)
        flux = 10 ** (-0.4 * (row["MAG_R"] - zp))
        gal = galaxy(row).withFlux(flux)
        psf_here = lw.toWorld(psf_pix).withFlux(1.0)
        n = None
        for b, (g1, g2) in run.branches.items():
            obj = galsim.Convolve(gal.shear(g1=g1, g2=g2), psf_here)
            n = n or min(obj.getGoodImageSize(lw.minLinearScale()), STAMP_MAX)
            st = obj.drawImage(nx=n, ny=n, wcs=lw, center=pos)
            overlap = st.bounds & canvases[b].bounds
            if overlap.isDefined():
                canvases[b][overlap] += st[overlap]

    noise = np.random.default_rng([run.cfg["seed"], int(exp), ccd]).normal(
        0.0, rms, (CCD_NY, CCD_NX)).astype(np.float32)
    work.mkdir(parents=True, exist_ok=True)
    for b in run.branches:
        np.save(work / f"{ccd:02d}-{b}.npy", canvases[b].array + noise)
    np.save(work / f"{ccd:02d}-weight.npy", weight)
    np.save(work / f"{ccd:02d}-flag.npy", flag)
    return {
        "EXPOSURE": int(exp), "CCD": ccd, "RMS": rms, "N_GAL": len(cand),
        "PSF_BETA": beta, "PSF_FWHM": fwhm, "PSF_G1_IN": pg1, "PSF_G2_IN": pg2,
        "STAR_T": target[0], "STAR_G1": target[1], "STAR_G2": target[2], "HAS_STARS": target[3],
        "HEADER": header.tostring(), "VIGNET": stamp, "SHAPES": shapes,
    }


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
    """PSF dictionary (CCDs with stars), CCD table, tile lists and per-sim metadata from what exists."""
    rows = []
    for path in sorted((run.shared / "psf").glob("*.pickle")):
        with open(path, "rb") as f:
            rows += pickle.load(f)
    psf_dict = {f"{r['EXPOSURE']}-{r['CCD']}": {"VIGNET": r["VIGNET"], "SHAPES": r["SHAPES"]}
                for r in rows if r.get("HAS_STARS", True)}
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

    global TRUTH
    truth_path = run.shared / "truth.fits"
    if not truth_path.exists():
        make_truth(run).write(str(truth_path) + ".tmp", format="fits")
        os.replace(str(truth_path) + ".tmp", truth_path)
    TRUTH = Table.read(truth_path)
    print(f"truth: {len(TRUTH)} galaxies", flush=True)

    if args.shard:
        exps = exps[args.shard[0]::args.shard[1]]
    todo = [e for e in exps if not done_exposure(run, e)]
    if todo:
        targets, rho4 = star_targets(run, sorted({e for _, es in run.tiles.values() for e in es}))
        beta = moffat_beta(rho4)
        print(f"stars: median rho4 {rho4:.4f} -> Moffat beta {beta:.3f}", flush=True)
        tasks = [(run, e, c, ccd_target(targets, int(e), c), beta)
                 for e in todo for c in range(N_CCD)]
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
