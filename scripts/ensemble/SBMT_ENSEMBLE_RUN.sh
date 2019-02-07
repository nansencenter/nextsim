#!/bin/bash

. neXtSIM_io_paths.env      # paths to relevant directories
. ensemble_run_function.sh  # functions used in this script
. ${DYNFILE}

checkpath ${ENSPATH}        # check if ensemble root directory/create
checkpath ${EnKFDIR}        # check if ensemble filter directory/create

ENSEMBLE+=()                # copy files to ensemble members' paths
for (( mem=1; mem<=${ESIZE}; mem++ )); do
    MEMBER=$(leadingzero 2 ${mem})
        MEMNAME=ENS${MEMBER}
    ENSEMBLE+=([${mem}]=${MEMNAME})
    MEMPATH=${ENSPATH}/${MEMNAME}; checkpath ${MEMPATH}
    SRUN=run_${EXPNAME}_${MEMNAME}.sh

    ${COPY} ${RUNFILE} ${MEMPATH}/${SRUN}
    ${COPY} ${RUNPATH}/${ENVFILE} ${MEMPATH}/.

    sed "s;^exporter_path=.*$;exporter_path="${MEMPATH}";g"\
        ${RUNPATH}/${CONFILE} > ${MEMPATH}/${CONFILE}

    sed -e "s;^MEMNAME=.*$;MEMNAME="${MEMNAME}";g" \
        -e "s;^MEMPATH=.*$;MEMPATH="${MEMPATH}";g" \
        ${RUNPATH}/${DYNFILE} > ${MEMPATH}/${SRCFILE}

    if [ ${mem} -gt "1" ];then
        while kill -0 "$XPID"; do
            echo "Process still running... ${XPID}"
            sleep 120
        done
    fi
# submit nextsim run on docker for the ensemble member #${mem}
    cd ${MEMPATH}; ID=$( getpid ./${SRUN} ${CONFILE} ${NPROC} ${ENVFILE} )
    XPID=${ID};echo ${XPID}
done
