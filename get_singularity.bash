#!/bin/bash
download="https://content.cruk.cam.ac.uk/jmlab/chimera_t_data/accessory/chimera2018.simg"
output="chimera-singularity.simg"
curl $download > $output
 