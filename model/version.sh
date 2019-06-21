#! /bin/bash -x

version=$( git describe --dirty --always --tags )
date=$( date )

cat << EOF > version.hpp
#ifndef __VERSION_HPP
#define __VERSION_HPP 1

#define NEXTSIM_VERSION_GIT "$version"
#define NEXTSIM_BUILD_TIME "$date"

#endif
EOF

