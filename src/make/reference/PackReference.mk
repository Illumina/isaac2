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
## file PackReference.mk
##
## brief Defines appropriate rules
##
## author Roman Petrovski
##
################################################################################

# first target needs to be defined in the beginning. Ohterwise includes such as
# Log.mk cause unexpected behavior
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

ifeq (,$(REFERENCE_GENOME))
$(error "REFERENCE_GENOME is not defined")
endif

ifeq (,$(OUTPUT_FILE))
$(error "OUTPUT_FILE is not defined")
endif

ORIGINAL_GENOME_FASTA:=$(notdir $(shell $(XSLTPROC) $(GET_ANY_FASTA_PATH_XSL) $(REFERENCE_GENOME)))
GENOME_FASTA:=$(TEMP_DIR)/$(ORIGINAL_GENOME_FASTA)

# in some cases the reference can be spread over multiple .fa files. 
# CONSOLIDATED_REFERENCE_XML ensures there is only one .fa file to pack
CONSOLIDATED_REFERENCE_XML:=$(TEMP_DIR)/sorted-reference.xml

ALL_GENOME_ANNOTATION_PATHS:=$(shell $(XSLTPROC) $(GET_ANNOTATION_FILE_PATHS_XSL) $(REFERENCE_GENOME))

ifeq (original,$(SEED_LENGTHS))
SUPPORTED_SEED_LENGTHS:=$(shell $(XSLTPROC) $(GET_SUPPORTED_SEED_LENGTHS_XSL) $(REFERENCE_GENOME))
else
SUPPORTED_SEED_LENGTHS:=$(SEED_LENGTHS)
endif

ALL_GENOME_NEIGHBORS:=$(foreach s, $(SUPPORTED_SEED_LENGTHS), $(TEMP_DIR)/$(GENOME_NEIGHBORS_PREFIX)$(s)$(GENOME_NEIGHBORS_SUFFIX))

genome_neighbors_seed_length=$(@:$(TEMP_DIR)/$(GENOME_NEIGHBORS_PREFIX)%$(GENOME_NEIGHBORS_SUFFIX)=%)

$(ALL_GENOME_NEIGHBORS): $(REFERENCE_GENOME) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(EXTRACT_NEIGHBORS) --reference-genome $< --seed-length $(genome_neighbors_seed_length) --output-file $(SAFEPIPETARGET)

$(CONSOLIDATED_REFERENCE_XML): $(REFERENCE_GENOME) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(REORDER_REFERENCE) --reference-genome $< --output-directory $(TEMP_DIR) --output-xml $(SAFEPIPETARGET) 

$(OUTPUT_FILE): $(ALL_GENOME_ANNOTATION_PATHS) $(ALL_GENOME_NEIGHBORS) $(CONSOLIDATED_REFERENCE_XML)
	$(CMDPREFIX) $(TAR) -czvO \
		$(foreach ga, $(ALL_GENOME_ANNOTATION_PATHS), -C $(dir $(ga)) $(notdir $(ga))) \
		-C $(CURDIR)/$(TEMP_DIR) \
		$(foreach gn, $(ALL_GENOME_NEIGHBORS), $(notdir $(gn))) \
		$(notdir $(GENOME_FASTA)) \
		$(notdir $(CONSOLIDATED_REFERENCE_XML)) \
		>$(SAFEPIPETARGET) 

all: $(OUTPUT_FILE)

