#!/bin/bash
################################################################################
##
## Isaac Genome Alignment Software
## Copyright (c) 2010-2014 Illumina, Inc.
## All rights reserved.
##
## This software is provided under the terms and conditions of the
## Illumina Public License 1
##
## You should have received a copy of the Illumina Public License 1
## along with this program. If not, see
## <https://github.com/sequencing/licenses/>.
##
################################################################################
##
## file cwinstallCmakecommon.sh
##
## Installation script for cmake
##
## author Come Raczy
##
################################################################################
REDIST_DIR=$1
INSTALL_DIR=$2
if [[ $# -ge 3 ]] ; then PARALLEL=$3 ; else PARALLEL=1 ; fi

. `dirname "$0"`/common.sh

BUILD_DIR=${INSTALL_DIR}/build
BIN_DIR=${INSTALL_DIR}/bin
LIB_DIR=${INSTALL_DIR}/lib
INCLUDE_DIR=${INSTALL_DIR}/include

CMAKE_MAJOR=2
CMAKE_MINOR=8
CMAKE_PATCH=0
CMAKE_REQUIRED="$CMAKE_MAJOR.$CMAKE_MINOR.$CMAKE_PATCH"
TARBALL_VERSION="2.8.9"
SCRIPT=`basename "$0"`
SOURCE_TARBALL=${REDIST_DIR}/cmake-$TARBALL_VERSION.tar.gz
TARBALL_COMPRESSION=z
SOURCE_DIR=${BUILD_DIR}/cmake-$TARBALL_VERSION
CMAKE_DIR=cmake-$CMAKE_MAJOR.$CMAKE_MINOR

common_options $@

if [[ $CLEAN ]] ; then
    echo removing $SOURCE_DIR >&2
    rm -rf $SOURCE_DIR
    rm -rf ${INSTALL_DIR}/{doc,share}/$CMAKE_DIR
    rm -f ${BIN_DIR}/{ccmake,cmake,cpack,ctest}
    rm -f ${INSTALL_DIR}/man/man1/{ccmake,cmake,cmakecommands,cmakecompat,cmakemodules,cmakeprops,cmakevars,cpack,ctest}.1
    exit 0
fi

AVAILABLE_CMAKE_VERSION=`cmake --version 2> /dev/null`
if [[ "${AVAILABLE_CMAKE_VERSION}" =~ ^cmake\ version\ ([0-9]+)\.([0-9]+)\.([0-9]+) && ! $FORCE ]] ; then
    MAJOR=${BASH_REMATCH[1]}
    MINOR=${BASH_REMATCH[2]}
    PATCH=${BASH_REMATCH[3]}
    if [[ "$MAJOR" -eq "$CMAKE_MAJOR" && ( "$MINOR" -gt "$CMAKE_MINOR" || "$MINOR" -eq "$CMAKE_MINOR" && "$PATCH" -ge "$CMAKE_PATCH"  ) ]] ; then
        echo "${BASH_REMATCH[0]} (>= $CMAKE_REQUIRED) is already installed" >&2
        echo nothing to be done >&2
        exit 1
    fi
fi 

OLD_CMAKE_VERSION=`${BIN_DIR}/cmake --version 2> /dev/null`;
if [[ $OLD_CMAKE_VERSION == "cmake version $TARBALL_VERSION" && ! $FORCE ]] ; then
    echo cmake version \"$TARBALL_VERSION\" is already installed at ${BIN_DIR}/cmake >&2
    echo nothing to be done >&2
    exit 0
elif [[ $OLD_CMAKE_VERSION != "" ]] ; then
    echo unable to install cmake version \"$TARBALL_VERSION\" in ${BIN_DIR} >&2 
    echo cmake version \"$OLD_CMAKE_VERSION\" is in the way. >&2
    echo Please use an empty location to build the product. >&2
    exit 2
fi 


##
## cleanup all existing source directory before proceeding
##
rm -rf $SOURCE_DIR

common_create_source

echo "Extracted cmake version $TARBALL_VERSION source code into $SOURCE_DIR" >&2
echo "Installing cmake using: './bootstrap --prefix=\"${INSTALL_DIR}\" --parallel=\"$PARALLEL\" && make -j \"$PARALLEL\" && make install'" >&2
cd $SOURCE_DIR && ./bootstrap --prefix="${INSTALL_DIR}" --parallel="$PARALLEL" && make -j "$PARALLEL" && make install

if [ $? != 0 ] ; then echo "cmake: build failed: Terminating..." >&2 ; exit 2 ; fi

echo "Cleaning up ${SOURCE_DIR}" >&2
rm -rf ${SOURCE_DIR}

echo CMake installed successfully >&2

exit 0
