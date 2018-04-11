#!/bin/bash
# example script to rename nextsim variables
# use with care
# - eg only run on a clean model so changes can be checked easily and also be reversed with 'git checkout'

oldnames=()

# DYNAMICS
#oldnames+=("simul.alea_factor")
#oldnames+=("simul.young")
#oldnames+=("simul.cfix")
#oldnames+=("simul.nu0")
#oldnames+=("simul.tan_phi")
#oldnames+=("simul.tract_coef")
#oldnames+=("simul.compr_strength")
#oldnames+=("simul.ridging_exponent")
#oldnames+=("simul.min_h")
#oldnames+=("simul.min_c")
#oldnames+=("simul.ridge_to_normal_cohesion_ratio")
#oldnames+=("simul.cohesion_thickness_normalisation")
#oldnames+=("simul.cohesion_thickness_exponent")
#oldnames+=("simul.scale_coef")
#oldnames+=("simul.use_temperature_dependent_healing")
#oldnames+=("simul.time_relaxation_damage")
#oldnames+=("simul.deltaT_relaxation_damage")
#oldnames+=("simul.undamaged_time_relaxation_sigma")
#oldnames+=("simul.exponent_relaxation_sigma")
#oldnames+=("simul.ERAi_quad_drag_coef_air")
#oldnames+=("simul.ECMWF_quad_drag_coef_air")
#oldnames+=("simul.ASR_quad_drag_coef_air")
#oldnames+=("simul.CFSR_quad_drag_coef_air")
#oldnames+=("simul.lin_drag_coef_air")
#oldnames+=("simul.quad_drag_coef_water")
#oldnames+=("simul.lin_drag_coef_water")
#oldnames+=("simul.use_coriolis")
#oldnames+=("simul.oceanic_turning_angle")
#oldnames+=("simul.Lemieux_basal_k1")
#oldnames+=("simul.Lemieux_basal_k2")
#oldnames+=("simul.Lemieux_basal_Cb")
#oldnames+=("simul.Lemieux_basal_u_0")
#oldnames+=("simul.Lemieux_basal_u_crit")

# THERMODYNAMICS
#oldnames+=("simul.Qio-type")
#oldnames+=("simul.use_thermo_forcing")
#oldnames+=("simul.albedoW")
#oldnames+=("simul.alb_scheme")
#oldnames+=("simul.flooding")
#oldnames+=("simul.alb_ice")
#oldnames+=("simul.alb_sn")
#oldnames+=("simul.I_0")
#oldnames+=("simul.Qdw")
#oldnames+=("simul.Fdw")
#oldnames+=("simul.newice_type")
#oldnames+=("simul.melt_type")
#oldnames+=("simul.hnull")
#oldnames+=("simul.PhiF")
#oldnames+=("simul.PhiM")
#oldnames+=("simul.h_thin_max")
#oldnames+=("simul.h_thin_min")
#oldnames+=("simul.drag_ice_t")
#oldnames+=("simul.drag_ocean_u")
#oldnames+=("simul.drag_ocean_t")
#oldnames+=("simul.drag_ocean_q")
#oldnames+=("simul.diffusivity_sss")
#oldnames+=("simul.diffusivity_sst")
#oldnames+=("simul.ocean_nudge_timeT")
#oldnames+=("simul.ocean_nudge_timeS")




newnames=()
#newnames+=("dynamics.alea_factor")
#newnames+=("dynamics.young")
#newnames+=("dynamics.cfix")
#newnames+=("dynamics.nu0")
#newnames+=("dynamics.tan_phi")
#newnames+=("dynamics.tract_coef")
#newnames+=("dynamics.compr_strength")
#newnames+=("dynamics.ridging_exponent")
#newnames+=("dynamics.min_h")
#newnames+=("dynamics.min_c")
#newnames+=("dynamics.ridge_to_normal_cohesion_ratio")
#newnames+=("dynamics.cohesion_thickness_normalisation")
#newnames+=("dynamics.cohesion_thickness_exponent")
#newnames+=("dynamics.scale_coef")
#newnames+=("dynamics.use_temperature_dependent_healing")
#newnames+=("dynamics.time_relaxation_damage")
#newnames+=("dynamics.deltaT_relaxation_damage")
#newnames+=("dynamics.undamaged_time_relaxation_sigma")
#newnames+=("dynamics.exponent_relaxation_sigma")
#newnames+=("dynamics.ERAi_quad_drag_coef_air")
#newnames+=("dynamics.ECMWF_quad_drag_coef_air")
#newnames+=("dynamics.ASR_quad_drag_coef_air")
#newnames+=("dynamics.CFSR_quad_drag_coef_air")
#newnames+=("dynamics.lin_drag_coef_air")
#newnames+=("dynamics.quad_drag_coef_water")
#newnames+=("dynamics.lin_drag_coef_water")
#newnames+=("dynamics.use_coriolis")
#newnames+=("dynamics.oceanic_turning_angle")
#newnames+=("dynamics.Lemieux_basal_k1")
#newnames+=("dynamics.Lemieux_basal_k2")
#newnames+=("dynamics.Lemieux_basal_Cb")
#newnames+=("dynamics.Lemieux_basal_u_0")
#newnames+=("dynamics.Lemieux_basal_u_crit")

#newnames+=("thermo.Qio-type")
#newnames+=("thermo.use_thermo_forcing")
#newnames+=("thermo.albedoW")
#newnames+=("thermo.alb_scheme")
#newnames+=("thermo.flooding")
#newnames+=("thermo.alb_ice")
#newnames+=("thermo.alb_sn")
#newnames+=("thermo.I_0")
#newnames+=("thermo.Qdw")
#newnames+=("thermo.Fdw")
#newnames+=("thermo.newice_type")
#newnames+=("thermo.melt_type")
#newnames+=("thermo.hnull")
#newnames+=("thermo.PhiF")
#newnames+=("thermo.PhiM")
#newnames+=("thermo.h_thin_max")
#newnames+=("thermo.h_thin_min")
#newnames+=("thermo.drag_ice_t")
#newnames+=("thermo.drag_ocean_u")
#newnames+=("thermo.drag_ocean_t")
#newnames+=("thermo.drag_ocean_q")
#newnames+=("thermo.diffusivity_sss")
#newnames+=("thermo.diffusivity_sst")
#newnames+=("thermo.ocean_nudge_timeT")
#newnames+=("thermo.ocean_nudge_timeS")

N=${#oldnames[@]}
for n in `seq 0 $((N-1))`
do
   v1=${oldnames[$n]}
   v2=${newnames[$n]}
   echo "$v1 -> $v2"
   lst1=`grep -l $v1 model/*.*`
   lst2=`grep -l -R $v1 core contrib research scripts`
   for f in ${lst1[@]} ${lst2[@]}
   do
      echo sed -i -e 's/'"$v1"'/'"$v2"'/g' $f
      sed -i -e 's/'"$v1"'/'"$v2"'/g' $f
   done
   echo
done
