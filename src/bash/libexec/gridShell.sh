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
## file gridShell.sh
##
## make $(SHELL) wrapper that uses GRID_SHELL environment variable to execute
## the loggingShell.sh
##
## author Roman Petrovski
##
################################################################################

#set -x
set -o pipefail
shopt -s compat31 2>/dev/null

LOG_SHELL_ORI=@iSAAC_HOME@@iSAAC_FULL_LIBEXECDIR@/loggingShell.sh

LOG_LEVEL=$1
shift

LOG_DATE_FORMAT=$1
shift

MAKE=$1
shift

FULL_TARGET=$1
shift

ALL_PREREQS=$1
shift

NEWER_PREREQS=$1
shift

SHELL_ORI=$1
shift

[[ "" == "$FULL_TARGET" ]] && echo "ERROR: target cannot be emtpy. Reason: $NEWER_PREREQS\nDependencies: $ALL_PREREQS\nCmd: $@" >&2 && exit 2

# get rid of -c that comes from make
[[ '-c' == $1 ]] && shift

[[ "" == "$QRSH_CMD" ]] && "ERROR: QRSH_CMD environment variable is not set" >&2 && exit 2

removed_multilines=$(echo "$@" | while read -r l; do line_without_newline=$(echo -n "$l"); echo -n "${line_without_newline%\\} "; done)
quoted_doublequotes=${removed_multilines//\"/\"\\\"\"}
escaped_dollars=${quoted_doublequotes//\$/\\\$}
escaped_backticks=${escaped_dollars//\`/\\\`}

$QRSH_CMD $LOG_SHELL_ORI "'$LOG_LEVEL'" "'$LOG_DATE_FORMAT'" "'$MAKE'" "'$FULL_TARGET'" "'$ALL_PREREQS'" "'$NEWER_PREREQS'" "'$SHELL_ORI'" -c "\"$escaped_backticks\""

exit $?
