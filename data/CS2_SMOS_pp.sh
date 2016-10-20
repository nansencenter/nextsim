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

EOF
        exit
fi

# Code starts here
sday=$(( 24*3600 )) # Seconds in the day

# Loop over all the files provided
for file in $*
do
        # Check the file name
        if [ "${file::23}" != "cs2_smos_ice_thickness_" -o ${#file} -ne 43 ]
        then
                echo "Illegal file name: $file --- skipping"
                continue
        fi

        # Get the dates from the file name
        t0=$( date -j -f %Y%m%d ${file:23:8} +%s )
        t1=$( date -j -f %Y%m%d ${file:32:8} +%s )

        # Loop over the dates and link
        # I use the -v option so the script prints out the links as it creates them
        for ((t=$t0;t<=$t1;t=t+$sday))
        do
                ln -vs $file cs2_smos_ice_thickness_$( date -j -f %s $t +%Y%m%d ).nc
        done
done

