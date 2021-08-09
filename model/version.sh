#! /bin/bash -x

version=$( git describe --dirty --always --tags )
date=$( date )
commit=$( git rev-parse HEAD )
branch=$( git rev-parse --abbrev-ref HEAD )
cc_path=$( which ${CC} )
cc_version=$( ${CC} -dumpversion )
cxx_path=$( which ${CXX} )
cxx_version=$( ${CXX} -dumpversion )

cat << EOF > version.hpp
#ifndef __VERSION_HPP
#define __VERSION_HPP 1

#define NEXTSIM_VERSION_GIT "$version"
#define NEXTSIM_BRANCH_GIT  "$branch"
#define NEXTSIM_COMMIT_GIT  "$commit"
#define NEXTSIM_BUILD_TIME  "$date"
#define NEXTSIMDIR          "$NEXTSIMDIR"
#define CC_PATH             "$cc_path"
#define CC_VERSION          "$cc_version"
#define CXX_PATH            "$cxx_path"
#define CXX_VERSION         "$cxx_version"

#endif
EOF

