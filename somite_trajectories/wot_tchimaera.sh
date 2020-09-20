#!/bin/bash

cd /nfs/research1/marioni/iimaz/embryo_integ/atlas/WOTatlas

wot optimal_transport --matrix data/pca_wot.mtx --cell_days data/time_wot.txt \
  --growth_iters 3 --local_pca 0 --out tmaps_pca/embryo --verbose

wot trajectory --tmap tmaps_pca/embryo \
  --cell_set data/cell_sets_paraxialmes_subclustering_wot.gmt --day 8.5 \
  --out trajectories_pca/paraxial_mes_E85 --verbose

wot trajectory --tmap tmaps_pca/embryo \
  --cell_set data/cell_sets_somiticmes_subclustering_wot.gmt --day 8.5 \
  --out trajectories_pca/somitic_mes_E85 --verbose
