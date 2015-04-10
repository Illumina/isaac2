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
## file Log.mk
##
## brief Partial makefile providing basic logging support
##
## - Redefines the SHELL
##
## author Roman Petrovski
##
################################################################################

SHELL_LOG_OLD := $(SHELL)

# Important: some targets have a very big number of prerequisites. Using $(wordlist... to
# avoid 'make: execvp: .../loggingshell.sh: Argument list too long'

# in make -n mode use debug sheell to be able to see the information on the prerequisites 
ifeq (n,$(filter n,$(MAKEFLAGS)))
override SHELL = $(if $@,$(warning [$@ ($?)]))$(SHELL_LOG_OLD)
override SHELL_LOG = $(if $@,$(warning [$@ ($?)]))$(SHELL_LOG_OLD)
else

LOG_SHELL:= $(NAME_SHELL) $(LIBEXEC_DIR)/loggingShell.sh
override SHELL = $(if $@,$(LOG_SHELL) '$(iSAAC_LOG_LEVEL)' '$(LOG_DATE_FORMAT)' '$(MAKE)' '$@' '$(wordlist 1, 50, $^)' '$(wordlist 1, 50, $?)' '$(SHELL_LOG_OLD)',$(SHELL_LOG_OLD))
override SHELL_LOG = $(if $@,$(LOG_SHELL) '$(iSAAC_LOG_LEVEL)' '$(LOG_DATE_FORMAT)' '$(MAKE)' '$@' '$(wordlist 1, 50, $^)' '$(wordlist 1, 50, $?)' '$(SHELL_LOG_OLD)',$(SHELL_LOG_OLD))
endif
