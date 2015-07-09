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
## file UnpackReference.mk
##
## brief Defines appropriate rules
##
## author Roman Petrovski
##
################################################################################

# first target needs to be defined in the beginning. Ohterwise includes such as
# Log.mk cause unexpected behavior
.PHONY: firsttarget
firsttarget: all

MAKEFILES_DIR:=@iSAAC_HOME@@iSAAC_FULL_DATADIR@/makefiles

# Import the global configuration
include $(MAKEFILES_DIR)/common/Config.mk

include $(MAKEFILES_DIR)/common/Sentinel.mk

# Import the logging functionalities
include $(MAKEFILES_DIR)/common/Log.mk

# Import the debug functionalities
include $(MAKEFILES_DIR)/common/Debug.mk

include $(MAKEFILES_DIR)/reference/Config.mk

ifeq (,$(INPUT_ARCHIVE))
$(error "INPUT_ARCHIVE is not defined")
endif

UNPACK_REFERENCE_SUBMAKE:=$(MAKEFILES_DIR)/reference/UnpackReferenceSubmake.mk

UNPACKED_SORTED_REFERENCE_XML:=$(TEMP_DIR)/sorted-reference.xml

SORTED_REFERENCE_XML:=sorted-reference.xml

$(UNPACKED_SORTED_REFERENCE_XML): $(INPUT_ARCHIVE) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(TAR) -C $(TEMP_DIR) --touch -xvf $(INPUT_ARCHIVE)

.PHONY: all
all: SHELL:=$(SHELL_LOG_OLD)
all: $(UNPACKED_SORTED_REFERENCE_XML)
	$(MAKE) -f $(UNPACK_REFERENCE_SUBMAKE) $(SORTED_REFERENCE_XML) && \
	$(CMDPREFIX) ($(LOG_INFO) "All done!")


