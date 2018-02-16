#! /bin/bash

# Help
if [ "$1" == "-h" -o "$1" == "--help" ]
then
        cat << EOF

usage: CS2_SMOS_pp.sh file_list

where file_list is a list of cs2_smos_ice_thickness files. E.g.

$ ./CS2_SMOS_pp.sh cs2_smos_ice_thickness_2010*
$ ./CS2_SMOS_pp.sh cs2_smos_ice_thickness_*
$ ./CS2_SMOS_pp.sh cs2_smos_ice_thickness_201[12]*
$ ./CS2_SMOS_pp.sh cs2_smos_ice_thickness_????????_????????.nc
$ ./CS2_SMOS_pp.sh cs2smos_ice_thickness_????????_????????_v1.3.nc

for each weekly file, this script creates 7 daily files (by linking)
in the directory where the script is run from

EOF
        exit
fi

# Code starts here
sday=$(( 24*3600 )) # Seconds in the day

# Loop over all the files provided
for file in $*
do
        file0=`basename $file`

        # Check the file name
        echo ${#file0}
        if [ "${file0::22}" != "cs2smos_ice_thickness_" -o ${#file0} -ne 47 ]
        then
                echo "Illegal file name: $file --- skipping"
                continue
        fi

        # Get the dates from the file name
        kernel=`uname -s`
        if [ $kernel == "Darwin" ]
        then
           t0=$( date -j -f %Y%m%d ${file0:22:8} +%s )
           t1=$( date -j -f %Y%m%d ${file0:31:8} +%s )

           # Loop over the dates and link
           # I use the -v option so the script prints out the links as it creates them
           for ((t=$t0;t<=$t1;t=t+$sday))
           do
                   ln -vs $file cs2_smos_ice_thickness_$( date -j -f %s $t +%Y%m%d ).nc
           done
        else
           t0=${file0:22:8}
           for i in `seq 0 6`
           do
              t1=`date --date="$t0 +${i}days" +%Y%m%d`
              ln -vs `readlink -f $file` cs2_smos_ice_thickness_${t1}.nc
           done
        fi
done
