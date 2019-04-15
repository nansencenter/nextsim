#!/bin/bash

LBREAK="/------------------------------------------------------------/"
KERNEL=$(uname -s)
MACHINE=$(uname -n)

[[ ${KERNEL} == 'Linux' ]] || DOCKER=1
[ !$# ] || [[ $@ == "EnKF" ]] && UPDATE=1

if [ ${DOCKER} ]; then
  echo ${LBREAK}
  echo " Running on a non-Linux machine  "
  echo ${LBREAK}
else
  echo ${LBREAK}
  echo " Running on a Linux machine  "
  echo ${LBREAK}
fi

[ ${UPDATE} ] && echo '/--State WILL be updated--/' || echo '/--NO state update--/'

. neXtSIM_io_paths.env      # paths to relevant directories
. ensemble_run_function.sh  # functions used in this script
. ${DYNFILE}

checkpath ${ENSPATH}        # check if ensemble root directory/create

# initialize and run ensemble
ENSEMBLE+=()                # copy files to ensemble members' paths
for (( mem=1; mem<=${ESIZE}; mem++ )); do
    MEMBER=$(leadingzero 3 ${mem}); MEMNAME=mem${MEMBER}
    ENSEMBLE+=([${mem}]=${MEMNAME})
    MEMPATH=${ENSPATH}/${MEMNAME}; checkpath ${MEMPATH}
    SRUN=run_${EXPNAME}_${MEMNAME}.sh

    ${COPY} ${RUNFILE} ${MEMPATH}/${SRUN}
    [ ${DOCKER} ] && nml_path=${RUNPATH} ||  nml_path=${MEMPATH}
    [ ${DOCKER} ] && exporter_path="/docker_io" || exporter_path=${MEMPATH}

    sed "s;^iopath.*$;iopath = '${exporter_path}';g"\
        ${EnKF_ROOT}/perturbation/nml/${NMLFILE} > ${MEMPATH}/${NMLFILE}

    [ ${ENVFILE} ] && ${COPY} ${RUNPATH}/${ENVFILE} ${MEMPATH}/.

    sed "s;^exporter_path=.*$;exporter_path="${exporter_path}";g"\
        ${RUNPATH}/${CONFILE} > ${MEMPATH}/${CONFILE}

    sed -e "s;^MEMNAME=.*$;MEMNAME="${MEMNAME}";g" \
        -e "s;^MEMPATH=.*$;MEMPATH="${MEMPATH}";g" \
        ${RUNPATH}/${DYNFILE} > ${MEMPATH}/${SRCFILE}

#    if [ ${mem} -gt "1" ];then
#        while kill -0 "$XPID"; do
#            echo "Process still running... ${XPID}"
#            sleep 120
#        done
#    fi
## submit nextsim run on docker for the ensemble member #${mem}
#    cd ${MEMPATH}; ID=$( getpid ./${SRUN} ${CONFILE} ${NPROC} ${ENVFILE} )
#    XPID=${ID};echo ${XPID}
done

if [ ${UPDATE} ]; then
# update ensemble
  checkpath ${FILTER}        # check if ensemble filter directory/create
  cd ${FILTER}

# some of the EnKF config files should be modified each analysis cycle
# such as time of the observations/model outputs
# link configuration files
  ${LINK} ${EnKF_CONF}/{enkf,enkf-global,grid,model,obs,obstypes,stats}.prm .
# link makefile, shell script and copy executables
  ${LINK} ${EnKF_CONF}/{Makefile,run_enkf.sh} .
  ${COPY} ${EnKF_EXEC}/enkf_{prep,calc,update} .

# copy/ link the prior.nc of each member into ${FILTER}
  for (( mem=1; mem<=${ESIZE}; mem++ )); do
    touch prior_${ENSEMBLE[${mem}]}.nc
  done

# 'make enkf' will call sequentially enkf_{prep,calc,update} executables
# write a dockerfile to call it
# -CHeCK- calls also clean / we may not want it
fi
