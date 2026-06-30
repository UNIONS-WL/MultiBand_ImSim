# scratch/ — orchestration scaffolding (work-in-progress)

Surrounding scaffolding for running the UNIONS image simulations and the
downstream ShapePipe-on-the-sims step. Published here for **visibility, all in
one place** — not yet cleaned up or generalised. Paths, cluster specifics, and
hard-coded grids reflect Fabian's candide working tree.

## SP_simu_fab/ — ShapePipe run on the SKiLLS sims (downstream half)

Originally `/n09data/hervas/SP_simu_fab/` on candide. Scripts only; the multi-TB
SP_* run outputs, Ms_tiles_* link farms, and per-tile output dirs are data and
are not included.

- `hack_get_images.sh`, `link_*.sh` — stage sim images/exposures into a
  ShapePipe-style input tree per shear grid.
- `job_per_tile.job`, `job_per_tile_newversion.job`, `job_multi_tiles.job`,
  `tile_launcher.job` — SLURM drivers that run ShapePipe per tile / per batch.
- `get_tile_numbers.sh`, `copy_grids.sh`, `finalize.sh`, `fix_links_and_log.sh`,
  `update*.sh`, `update_tiles*.py`, `log_merge_final_cat.py` — bookkeeping over
  the per-grid runs.
- `star_selection.setools`, `*_tiles*.txt`, `numbers_run.txt` — selection config
  and the tile lists each grid was run over.

Shear grids: `1m2z / 1p2z / 1z2m / 1z2p / 1z2z` = the g1/g2 ± shear combinations.
