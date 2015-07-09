################################################################################
##
## Isaac Genome Alignment Software
## Copyright (c) 2010-2014 Illumina, Inc.
## All rights reserved.
##
## This software is provided under the terms and conditions of the
## GNU GENERAL PUBLIC LICENSE Version 3
##
## You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
## along with this program. If not, see
## <https://github.com/illumina/licenses/>.
##
################################################################################
##
## file common.sh
##
## Definition of functions and variables common to all bootstrap scripts.
##
## author Come Raczy
##
################################################################################

function common_options () {
    TEMP=`getopt -n $SCRIPT -o fcd -- "$@"`
    if [ $? != 0 ] ; then echo $SCRIPT: invalid option  >&2; echo "Terminating..." >&2 ; exit 2 ; fi
    eval set -- "$TEMP"
    FORCE=
    CLEAN=
    DEBUG=
    while true ; do
        case "$1" in
            -f) FORCE=true ; shift ;;
            -c) CLEAN=true ; shift ;;
            -d) DEBUG=true ; shift ;;
            --)              shift ; break ;;
            *) echo "Internal error!" >&2; exit 2 ;;
        esac
    done
}

function common_create_source () {
    if [[ ! -e $SOURCE_TARBALL ]] ; then
        echo $SCRIPT: source tarball $SOURCE_TARBALL not found >&2
        exit 2
    fi  
    echo Decompressing $SOURCE_TARBALL >&2
    mkdir -p ${BUILD_DIR}
    tar -C ${BUILD_DIR} -${TARBALL_COMPRESSION}xf $SOURCE_TARBALL
    
    if [[ ! -d $SOURCE_DIR ]] ; then
        echo $SOURCE_DIR does not exist >&2
        exit 2
    fi
}

