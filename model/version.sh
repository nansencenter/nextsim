#! /bin/bash -x

version=$( git describe --dirty --always --tags )
date=$( date )
commit=$( git rev-parse HEAD )
branch=$( git rev-parse --abbrev-ref HEAD )

cat << EOF > version.hpp
#ifndef __VERSION_HPP
#define __VERSION_HPP 1

#define NEXTSIM_VERSION_GIT "$version"
#define NEXTSIM_BRANCH_GIT "$branch"
#define NEXTSIM_COMMIT_GIT "$commit"
#define NEXTSIM_BUILD_TIME "$date"

#endif
EOF

