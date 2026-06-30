for f in Ms_tiles_1m2z_good_offset/final_cat*.fits; do
    ln -s "$(realpath "$f")" \
          SP_1z2m_grid/output/run_sp_combined_final/make_catalog_runner/output/
done
