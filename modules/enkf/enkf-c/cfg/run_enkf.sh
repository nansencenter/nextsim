#!/bin/bash
# perform data assimilation using enkf-c
# this script should be called after neXtSIM-ensemble

if [[ $# == 0 ]]; then
echo Argument needed:
echo \'help\': to see help menu
echo \'run\' : to execute the script

elif [[ $@ == "execute" ]]; then
enkf_dir=/Group/da/aliayd/enkf-c/enkf
ensdir=Moorings
ens=( /Group/da/sukeng/IO_nextsim/neXtSIM_test19_03_airdrag/1/*/Moorings.nc )

for (( mem = 0; mem < ${#ens[@]}; mem++ )); do
#  cp ${ens[${mem}]} ${ensdir}/mem$(printf "%03d" $[mem+1]).nc
  echo mem$(printf "%03d" $[mem+1]).nc
done

cd conf
#ln -sf ../Moorings/mem001.nc grid.nc;
#ln -sf ../Moorings/mem001.nc tmp.nc;
cd ../
#cd ../;  ln -sf ${enkf_dir}/bin/{enkf_prep,enkf_calc,enkf_update} .
make enkf

elif [[ $@ == "help" ]]; then
  H=()
  H+=("neXtSIM will provide mem%3d.nc in which all state variables will be on a curvilinear regular grid")
  H+=("mem%3d.nc will be linked into the FILTER directory")
  H+=("all \*.prm files will be modified by a shell script and linked into the FILTER directory")
  H+=("enkf_prep enkf_calc, enkf_update will be linked into the FILTER directory")
  H+=("observations in the assimilation cycle will be linked into the FILTER/obs directory")
  H+=("mem%3d.nc.analysis will be written by enkf-c")
  for (( com = 0; com <= ${#H[@]}; com++ )); do
      echo $[com+1]- ${H[${com}]}
  done
exit 0
fi
