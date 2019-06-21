#! /bin/bash -x

version=$( git describe --dirty --always --tags )
date=$( date )

cat << EOF > version.hpp
#if !defined (NEXTSIM_VERSION_GIT)
#define NEXTSIM_VERSION_GIT "$version"
#endif

#if !defined (NEXTSIM_BUILD_TIME)
#define NEXTSIM_BUILD_TIME "$date"
#endif
EOF

