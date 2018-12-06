#!/bin/bash
#docker run --rm -it \
#    --security-opt seccomp=unconfined \
# where are the LINKS to the meshes
#    -v /home/polona/src/nextsim/nextsim:/nextsim \
# ???
#    -v /home/polona/data/sim/data/:/Data/sim/data/ \
# where are the LINKS to the ocean forecast forcing data?
#    -v /home/polona/data/nextsimf/data/:/Data/nextsimf/data/ \
# where are the LINKS to the forcing data?
#    -v /home/polona/data/nextsimf/src/nextsim-cxx/nextsim-develop-new/data/:/simdata/data \
#where to save outputs
#    -v /home/antonk/nextsim/nextsimf_forecasts/:/simforecast/ \  
#    nextsim /nextsim/test.cfg 7
docker run --rm -it \
    --security-opt seccomp=unconfined \
    -v /home/polona/data/sim/data/:/home/polona/data/sim/data/ \
    -v /home/polona/data/nextsimf/data/:/home/polona/data/nextsimf/data/ \
    -v /mnt/10.11.12.231/sim/data/:/mnt/10.11.12.231/sim/data/ \
    -v /home/polona/src/nextsim/mesh/:/mesh/ \
    -v /home/polona/src/nextsim/data_links/:/data/ \
    -v /home/polona/data/nextsim_output/:/docker_io/ \
    nextsim /docker_io/nextsim_model_options_parallel.cfg 7
#   nextsim
