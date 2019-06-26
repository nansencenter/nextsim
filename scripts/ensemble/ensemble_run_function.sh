#!/bin/bash

function leadingzero {
	zeros=$1
	input=$2
	printf "%0${zeros}d" ${input}
	}

# get pid to control work submission
# last argument is the log file; others are the script and its arguments
function getpid() {
  runfile=$1
  confile=$2
  nproc=$3
  envfile=$4
  logfile=sbmt.log
  ${runfile} ${confile} ${nproc} ${envfile} &> ${logfile} &
  echo $!
}

function wait_docker(){

     checkid=$1

     A=($(cat docker.sample.ls | awk '{print $NF}' |\
          sed -n 2,\$p)\
       )
    echo ${A[@]}
    for entry in "${A[@]}"; do
         echo ${entry}
         [[ $checkid == "$entry" ]] && echo "$var present in the array" && break
    done
}
function anywait(){
    for pid in "$@"; do
        while kill -0 "$pid"; do
            echo "Process still running... ${pid}"
            sleep 1
        done
    done
	}

function anywait_w_status2() {
    while true
    do
        alive_pids=()
        for pid in "$@"
        do
            kill -0 "$pid" 2>/dev/null \
                && alive_pids+="$pid "
        done

        if [ ${#alive_pids[@]} -eq 0 ]
        then
            break
        fi

        echo "Process(es) still running... ${alive_pids[@]}"
        sleep 120
    done
    echo 'All processes terminated'
	}

function checkpath { path=$1
	[ -d ${path} ] && echo "!----------- WARNING --------------------!  "\
                       && echo "! Ensemble directory ${path} exists!        "\
                       && echo "!------------- STOP ---------------------!  "\
                       && exit                                               \
                       || echo "${path} is created!                         "\
                       && ${MKDIR} ${path}
	}

# get jobid to control work submission on HPC. Machine-dependent
function jobid {
         output=$($*)
         echo ${output} | head -n1 | cut -d'<' -f2 | cut -d'>' -f1
	}
#ID=$( jobid bsub < ./run_docker_nextsim_bangkok.sh ${NPROC} ); echo ${ID}


