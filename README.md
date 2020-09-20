# TChimeras2020

This repo contains the scripts used for analysis in Guibentif et al. 2020.
Please note that these scripts are not well configured for directly re-running the code.
Rather, this is a record of what was done.
I do not recommend that you try to re-run this cofe directly.
If you do want to reproduce some analyses, please use the data from the R package (below), from which you should be able to reproduce fairly easily the main parts of our work.
A singularity container with the software versions used in these analyses can be accessed with the `get_singularity.bash` file in this repo.

The data is available in the [MouseGastrulationData](https://www.bioconductor.org/packages/devel/data/experiment/html/MouseGastrulationData.html) Bioconductor package, from version 3.12 onwards.
The function `TChimeraData` will get the processed single-cell RNAseq data from the chimeras, and `EmbryoAtlasData` the data from the atlas.
For accessory files related to this paper (masses, NMP orderings, trajectories etc.) please see the function `GuibentifData` (which is to-be-implemented at the time of writing).

An interactive app for the data will be available shortly.
