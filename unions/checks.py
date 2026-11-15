"""Known-answer checks of the UNIONS twin.

    python -m unions.checks unions/config.yaml [--tile 205.239] [-j N]

1. Astrometry and photometry: the production code (twin.draw_ccd, assemble,
   coadd) is run noise-free on one tile's exposures with a sparse truth of
   bright round galaxies and point sources, into {output_root}/kat/. Centroids
   and aperture fluxes are compared with the truth positions (through the
   written headers' WCS) and the expected fluxes, in the exposures (at PHOTZP)
   and in the coadd (at zero point 30).
2. PSF dictionary: the point sources at each CCD's centre and corners are
   compared with that CCD's dictionary stamp.
3. Noise: the rms of the built sims' CCDs, by the same estimator, against the
   real CCDs' rms; the coadd's rms against the real tile's and its weight map.
4. Branches: difference images between branches of the built run, away from
   the galaxies.
"""

import argparse
import os
import pickle
from multiprocessing import Pool

import galsim
import numpy as np
import pandas as pd
import yaml
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from unions import twin


def measure(image, x, y, half=10, rap=30):
    """First-moment centroid (1-based) and aperture flux around (x, y)."""
    ix, iy = int(round(x)) - 1, int(round(y)) - 1
    ny, nx = image.shape
    if not (rap <= ix < nx - rap and rap <= iy < ny - rap):
        return np.nan, np.nan, np.nan
    box = image[iy - half:iy + half + 1, ix - half:ix + half + 1]
    yy, xx = np.mgrid[iy - half:iy + half + 1, ix - half:ix + half + 1] + 1.0
    cx, cy = (box * xx).sum() / box.sum(), (box * yy).sum() / box.sum()
    ap = image[iy - rap:iy + rap + 1, ix - rap:ix + rap + 1]
    yy, xx = np.mgrid[-rap:rap + 1, -rap:rap + 1]
    return cx, cy, ap[(xx ** 2 + yy ** 2) <= rap ** 2].sum()


STAR_PIXELS = [(1025, 2307), (101, 101), (1948, 4512)]  # centre and two corners


def kat_truth(run, tile, spacing=60.0, mag=20.0):
    """Round bright galaxies on a lattice in the tile; point sources at the
    centre and two corners of every CCD of the tile's exposures, on pixel
    centres, keeping those more than 20 arcsec from one another."""
    header, exps = run.tiles[tile]
    w = WCS(header)
    n = header["NAXIS1"]
    step = spacing / 3600 / abs(header["CD1_1"])
    g = np.arange(step / 2, n, step)
    gx, gy = np.meshgrid(g, g)
    gra, gdec = w.all_pix2world(gx.ravel(), gy.ravel(), 1)
    stars = []
    for e in exps:
        with fits.open(run.exp_path(e)) as f:
            for c in range(twin.N_CCD):
                cw = WCS(twin.ccd_header(f[c + 1].header))
                for x, y in STAR_PIXELS:
                    r, d = cw.all_pix2world([[x, y]], 1)[0]
                    stars.append((r, d, int(e), c, x, y))

    def sep(r, d, ra, dec):
        return np.hypot((np.asarray(ra) - r) * np.cos(np.radians(d)), np.asarray(dec) - d)

    kept = []
    for s in stars:
        if not kept or sep(s[0], s[1], [k[0] for k in kept], [k[1] for k in kept]).min() > 20 / 3600:
            kept.append(s)
    kept = np.array(kept)
    far = np.array([sep(r, d, kept[:, 0], kept[:, 1]).min() > 20 / 3600 for r, d in zip(gra, gdec)])
    ns, ng = len(kept), far.sum()
    star = np.r_[np.zeros(ng, bool), np.ones(ns, bool)]
    return Table({
        "ID": np.arange(ng + ns), "RA": np.r_[gra[far], kept[:, 0]],
        "DEC": np.r_[gdec[far], kept[:, 1]], "STAR": star,
        "STAR_EXP": np.r_[np.full(ng, -1), kept[:, 2]].astype(int),
        "STAR_CCD": np.r_[np.full(ng, -1), kept[:, 3]].astype(int),
        "STAR_X": np.r_[np.zeros(ng), kept[:, 4]].astype(int),
        "STAR_Y": np.r_[np.zeros(ng), kept[:, 5]].astype(int),
        "MAG_R": np.full(ng + ns, mag), "SERSIC_N": np.ones(ng + ns),
        "HLR_ARCSEC": np.where(star, 1e-3, 0.3), "Q": np.ones(ng + ns),
        "PA_DEG": np.zeros(ng + ns),
    })


def kat(cfg, tile, jobs):
    """Run the production code noise-free on one tile and measure it."""
    cfg = dict(cfg, output_root=os.path.join(cfg["output_root"], "kat"), branches=["1z2z"])
    cfg["psf_dict"] = os.path.join(cfg["output_root"], "psf_dict.pickle")
    run = twin.Run(cfg)
    exps = run.tiles[tile][1]
    (run.shared / "psf").mkdir(parents=True, exist_ok=True)
    twin.TRUTH = kat_truth(run, tile)
    twin.TRUTH.write(run.shared / "truth.fits", overwrite=True)
    twin.sky_rms = lambda data, good: 0.0
    targets, rho4 = twin.star_targets(run, sorted({e for _, es in run.tiles.values() for e in es}))
    beta = twin.moffat_beta(rho4)
    tasks = [(run, e, c, tuple(targets.loc[(int(e), c)].iloc[:3]), beta)
             for e in exps for c in range(twin.N_CCD)]
    with Pool(jobs) as pool:
        rows = pool.map(twin.draw_ccd, tasks)
        pool.map(twin.assemble, [(run, e, [r for r in rows if r["EXPOSURE"] == int(e)])
                                 for e in exps])
    twin.coadd((run, tile, "1z2z", jobs))
    twin.finish(run)
    sim = run.sim["1z2z"]
    truth = twin.TRUTH
    zero = 10 ** (-0.4 * truth["MAG_R"])
    trunc = galsim.Sersic(n=1, half_light_radius=0.3, trunc=1.5, flux_untruncated=True).flux

    out = []
    for e in exps:
        with fits.open(sim / "images" / "SP_exp" / twin.exp_file("image", e)) as f:
            for c in range(twin.N_CCD):
                h = f[c + 1].header
                x, y = WCS(h).all_world2pix(truth["RA"], truth["DEC"], 1)
                on = (x > 0) & (x < twin.CCD_NX) & (y > 0) & (y < twin.CCD_NY)
                for i in np.flatnonzero(on):
                    cx, cy, fl = measure(f[c + 1].data, x[i], y[i])
                    expect = zero[i] * 10 ** (0.4 * h["PHOTZP"]) * trunc
                    out.append(("exp", e, c, truth["STAR"][i], cx - x[i], cy - y[i], fl / expect))
    image = fits.getdata(sim / "images" / "SP_tiles" / f"CFIS_simu_image-{twin.tile_name(tile)}.fits")
    th = fits.getheader(sim / "images" / "SP_tiles" / f"CFIS_simu_image-{twin.tile_name(tile)}.fits")
    x, y = WCS(th).all_world2pix(truth["RA"], truth["DEC"], 1)
    for i in range(len(truth)):
        cx, cy, fl = measure(image, x[i], y[i])
        expect = zero[i] * 1e12 * trunc
        out.append(("tile", tile, -1, truth["STAR"][i], cx - x[i], cy - y[i], fl / expect))
    res = pd.DataFrame(out, columns=["where", "exp", "ccd", "star", "dx", "dy", "flux_ratio"]).dropna()
    res.to_parquet(run.shared / "kat_astrometry.parquet")
    print("1. astrometry and photometry (noise-free, measured - truth)")
    for where, r in res.groupby("where"):
        print(f"  {where}: {len(r)} sources; dx median {r.dx.median():+.4f} max|dx| "
              f"{r.dx.abs().max():.4f} px; dy median {r.dy.median():+.4f} max|dy| "
              f"{r.dy.abs().max():.4f} px; flux ratio median {r.flux_ratio.median():.5f} "
              f"[{r.flux_ratio.quantile(0.05):.5f}, {r.flux_ratio.quantile(0.95):.5f}] (5-95%)")

    with open(cfg["psf_dict"], "rb") as fh:
        psf = pickle.load(fh)
    dev = []
    for st in truth[truth["STAR"]]:
        e, c, x, y = st["STAR_EXP"], st["STAR_CCD"], st["STAR_X"], st["STAR_Y"]
        with fits.open(sim / "images" / "SP_exp" / twin.exp_file("image", e)) as f:
            cut = f[c + 1].data[y - 26:y + 25, x - 26:x + 25].astype(float)
            local = galsim.GSFitsWCS(header=f[c + 1].header).local(galsim.PositionD(x, y))
        cut /= cut.sum()
        ref = psf[f"{e}-{c}"]["VIGNET"]
        m = twin.hsm_sky(cut, local)
        r = twin.hsm_sky(ref, local)
        dev.append((np.abs(cut - ref).max() / ref.max(),
                    m.moments_sigma ** 2 / r.moments_sigma ** 2 - 1,
                    m.observed_shape.g1 - r.observed_shape.g1,
                    m.observed_shape.g2 - r.observed_shape.g2))
    dev = np.array(dev)
    print(f"2. PSF dictionary stamp vs point source drawn at a CCD's centre or corner "
          f"({len(dev)} sources): max |pixel diff|/peak {dev[:, 0].max():.2e}; "
          f"max |dT/T| {np.abs(dev[:, 1]).max():.2e}; max |dg1| {np.abs(dev[:, 2]).max():.2e}; "
          f"max |dg2| {np.abs(dev[:, 3]).max():.2e} (HSM, sky frame at the source)")


def noise(cfg, tile):
    """Sim CCD rms vs real; coadd rms vs real tile and weight map."""
    run = twin.Run(cfg)
    ccds = pd.read_parquet(run.shared / "ccds.parquet")
    sim = run.sim["1z2z"]
    rows = []
    for e in run.tiles[tile][1]:
        with fits.open(sim / "images" / "SP_exp" / twin.exp_file("image", e)) as f, \
                fits.open(run.shared / "SP_exp" / twin.exp_file("weight", e)) as w, \
                fits.open(run.shared / "SP_exp" / twin.exp_file("flag", e)) as fl:
            for c in range(0, twin.N_CCD, 3):
                good = (w[c + 1].data > 0) & (fl[c + 1].data == 0)
                real = ccds[(ccds.EXPOSURE == int(e)) & (ccds.CCD == c)].RMS.item()
                rows.append(twin.sky_rms(f[c + 1].data, good) / real)
    rows = np.array(rows)
    print(f"3a. sim CCD rms / real CCD rms ({len(rows)} CCDs of the tile's exposures): "
          f"median {np.median(rows):.4f}, range [{rows.min():.4f}, {rows.max():.4f}]")
    name = twin.tile_name(tile)
    simg = fits.getdata(sim / "images" / "SP_tiles" / f"CFIS_simu_image-{name}.fits")
    swt = fits.getdata(sim / "images" / "SP_tiles" / f"CFIS_simu_weight-{name}.fits")
    rimg = fits.getdata(run.campaign / "tiles" / f"CFIS.{tile}.r.fits")
    rwt = fits.getdata(run.campaign / "tiles" / f"CFIS.{tile}.r.weight.fits.fz", 1)
    for label, img, wt in (("sim", simg, swt), ("real", rimg, rwt)):
        good = wt > 0
        rms = twin.sky_rms(img, good)
        med_w = np.median(wt[good])
        print(f"3b. {label} tile {tile}: sky rms {rms:.4f} (ZP 30 units); median weight "
              f"{med_w:.4g}; rms * sqrt(median weight) {rms * np.sqrt(med_w):.3f}")


def branches(cfg, tile):
    """Branch differences are the galaxies alone."""
    run = twin.Run(cfg)
    e = run.tiles[tile][1][0]
    names = list(run.sim)
    with fits.open(run.sim[names[1]] / "images" / "SP_exp" / twin.exp_file("image", e)) as a, \
            fits.open(run.sim[names[2]] / "images" / "SP_exp" / twin.exp_file("image", e)) as b, \
            fits.open(run.sim[names[0]] / "images" / "SP_exp" / twin.exp_file("image", e)) as z:
        c = 20
        d = a[c + 1].data - b[c + 1].data
        n = a[c + 1].data - z[c + 1].data
        print(f"4a. exposure {e} CCD {c}, {names[1]} - {names[2]}: {np.mean(d == 0):.3f} of pixels "
              f"exactly 0; rms of difference {d.std():.4f} vs image rms {z[c + 1].data.std():.3f}; "
              f"max |diff| {np.abs(d).max():.2f}; {names[1]} - {names[0]}: "
              f"{np.mean(n == 0):.3f} exactly 0")
    name = twin.tile_name(tile)
    ta = fits.getdata(run.sim[names[1]] / "images" / "SP_tiles" / f"CFIS_simu_image-{name}.fits")
    tb = fits.getdata(run.sim[names[2]] / "images" / "SP_tiles" / f"CFIS_simu_image-{name}.fits")
    d = ta - tb
    q = np.abs(d)
    print(f"4b. tile {tile}, {names[1]} - {names[2]}: median |diff| {np.median(q):.2e}, "
          f"99th pct {np.percentile(q, 99):.3f}, tile rms {ta.std():.3f}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("config")
    p.add_argument("--tile", default="205.239")
    p.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    args = p.parse_args()
    with open(args.config) as f:
        cfg = yaml.safe_load(f)
    noise(cfg, args.tile)
    branches(cfg, args.tile)
    kat(cfg, args.tile, args.jobs)


if __name__ == "__main__":
    main()
