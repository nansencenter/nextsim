#!/bin/bash

inml=../nml/pseudo2D.nml
onml=../bin/pseudo2D.nml
srun=run_pseudo2D.sh
fout=../POST/mean_wind_speed.txt

ifrt=../POST/frt.slp_vs_wnd_tmp.jnl
ofrt=../POST/frt.slp_vs_wnd.jnl

iplt=../POST/frt.plt.slp_vs_wnd_tmp.jnl
oplt=../POST/frt.plt.slp_vs_wnd.jnl

rm ${fout}; touch ${fout}

var_slp=( 10 20 30 40 50 60 100)
var_wnd=( 0 0.01 0.04 0.09 0.16 0.25 0.36 0.49 0.64 0.81 1 )
echo 'slp variance: '${var_slp[@]}' '${#var_slp[@]}
echo 'wind speed variance: '${var_wnd[@]}' '${#var_wnd[@]}

wnd_index=0
while [ ${wnd_index} -lt ${#var_wnd[@]} ]; do
slp_index=0
  while [ ${slp_index} -lt ${#var_slp[@]} ]; do
      rm ../IO/{randfld,synforc}.{00,01}
      sed -e "s;^vslp     =.*$;vslp     =  "${var_slp[${slp_index}]}";g"\
          -e "s;^vwndspd  =.*$;vwndspd  =  "${var_wnd[${wnd_index}]}";g"\
      ${inml} > ${onml}
      sed -e "s;^let vslp=.*$;let vslp="${var_slp[${slp_index}]}";g"\
          -e "s;^let vwndspd=.*$;let vwndspd="${var_wnd[${wnd_index}]}";g"\
      ${ifrt} > ${ofrt}
      ./${srun} run slp
      (( slp_index++ ))
  done
  (( wnd_index++ ))
done
sed -e "s;^let slp_size=.*$;let slp_size="${#var_slp[@]}";g"\
    -e "s;^let wnd_size=.*$;let wnd_size="${#var_wnd[@]}";g"\
    ${iplt} > ${oplt}
cd ..
source activate FERRET
pyferret -png -script POST/frt.plt.slp_vs_wnd.jnl
source deactivate FERRET
