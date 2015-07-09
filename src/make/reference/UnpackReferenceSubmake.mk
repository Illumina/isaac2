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
## file UnpackReferenceSubmake.mk
##
## brief Defines appropriate rules
##
## author Roman Petrovski
##
################################################################################

MAKEFILES_DIR:=@iSAAC_HOME@@iSAAC_FULL_DATADIR@/makefiles

# Import the global configuration
include $(MAKEFILES_DIR)/common/Config.mk

include $(MAKEFILES_DIR)/common/Sentinel.mk

# Import the logging functionalities
include $(MAKEFILES_DIR)/common/Log.mk

# Import the debug functionalities
include $(MAKEFILES_DIR)/common/Debug.mk

include $(MAKEFILES_DIR)/reference/Config.mk

ifeq (,$(MASK_WIDTH))
$(error "MASK_WIDTH is not defined")
endif


ifeq (yes,$(MOVABLE))
ROOT_PATH:=.
else
ROOT_PATH:=$(CURDIR)
endif


UNPACKED_SORTED_REFERENCE_XML:=$(TEMP_DIR)/sorted-reference.xml
ORIGINAL_GENOME_FASTA:=$(notdir $(shell $(XSLTPROC) $(GET_ANY_FASTA_PATH_XSL) $(UNPACKED_SORTED_REFERENCE_XML)))

UNPACKED_GENOME_FILE:=$(TEMP_DIR)/$(ORIGINAL_GENOME_FASTA)
GENOME_FILE:=$(ROOT_PATH)/$(ORIGINAL_GENOME_FASTA)

MASK_COUNT:=$(shell $(AWK) 'BEGIN{print 2^$(MASK_WIDTH)}')
MASK_LIST:=$(wordlist 1, $(MASK_COUNT), $(shell $(SEQ) --equal-width 0 $(MASK_COUNT)))

MASK_FILE_PREFIX:=kmer-positions-
MASK_FILE_MIDDLE:=-

SORTED_REFERENCE_XML:=sorted-reference.xml

GENOME_ANNOTATION_PATH:=$(shell $(XSLTPROC) $(GET_ANNOTATION_FILE_PATHS_XSL) $(UNPACKED_SORTED_REFERENCE_XML))
ifneq (1,$(words $(GENOME_ANNOTATION_PATH)))
$(error "Expected 1 annotation path, got: $(GENOME_ANNOTATION_PATH)")
endif

ANNOTATION:=$(notdir $(GENOME_ANNOTATION_PATH))
ANNOTATION_XML:=$(TEMP_DIR)/$(ANNOTATION).xml
UNPACKED_ANNOTATION:=$(TEMP_DIR)/$(ANNOTATION)

mask=$(word 2,$(subst $(MASK_FILE_MIDDLE), ,$(@:$(TEMP_DIR)/$(MASK_FILE_PREFIX)%$(MASK_FILE_XML_SUFFIX)=%)))
seed_length=$(word 1,$(subst $(MASK_FILE_MIDDLE), ,$(@:$(TEMP_DIR)/$(MASK_FILE_PREFIX)%$(MASK_FILE_XML_SUFFIX)=%)))

get_format_version=$(shell $(XSLTPROC) $(GET_FORMAT_VERSION_XSL) $(UNPACKED_SORTED_REFERENCE_XML))
ifeq (original,$(SEED_LENGTHS))
get_seed_length_list=$(shell $(XSLTPROC) $(GET_SUPPORTED_SEED_LENGTHS_XSL) $(UNPACKED_SORTED_REFERENCE_XML))
else
get_seed_length_list:=$(SEED_LENGTHS)
endif

ifneq ($(CURRENT_REFERENCE_FORMAT_VERSION),$(get_format_version))
$(error "Unsupported packed reference format $(get_format_version). Version $(CURRENT_REFERENCE_FORMAT_VERSION) is required")
endif

$(GENOME_FILE): $(UNPACKED_GENOME_FILE)
	$(CMDPREFIX) $(CP) $< $(SAFEPIPETARGET)

$(ANNOTATION) : $(UNPACKED_ANNOTATION)
	$(CMDPREFIX) $(CP) $(TEMP_DIR)/$@ $(SAFEPIPETARGET)

$(ANNOTATION_XML): $(ANNOTATION)
	$(CMDPREFIX) echo "\
	<SortedReference>\
		<FormatVersion>$(CURRENT_REFERENCE_FORMAT_VERSION)</FormatVersion>\
		<Annotations>\
			<Annotation Type='k-uniqueness' K='$(patsubst $(TEMP_DIR)/%uniqueness$(ANNOTATION_MASK_FILE_SUFFIX).xml,%,$@)'>\
				<Format>$(ANNOTATION_FORMAT)</Format>\
				<File>$(ROOT_PATH)/$(patsubst $(TEMP_DIR)/%.xml,%,$@)</File>\
			</Annotation>\
		</Annotations>\
	</SortedReference>\
	" >$(SAFEPIPETARGET)

$(CONTIGS_XML): $(GENOME_FILE) $(TEMP_DIR)/.sentinel $(UNPACKED_SORTED_REFERENCE_XML)
	$(CMDPREFIX) $(PRINT_CONTIGS) -g $(GENOME_FILE) --original-metadata $(UNPACKED_SORTED_REFERENCE_XML) >$(SAFEPIPETARGET)

ALL_MASK_XMLS:=$(foreach sl, $(get_seed_length_list), $(foreach m, $(MASK_LIST), $(TEMP_DIR)/$(MASK_FILE_PREFIX)$(sl)$(MASK_FILE_MIDDLE)$(m)$(MASK_FILE_XML_SUFFIX)))
$(ALL_MASK_XMLS): $(CONTIGS_XML) $(ALL_GENOME_NEIGHBORS)
	$(CMDPREFIX) $(SORT_REFERENCE) -r $(CONTIGS_XML) --mask-width $(MASK_WIDTH) --mask $(mask) \
		--output-file $(ROOT_PATH)/$(notdir $(@:%$(MASK_FILE_XML_SUFFIX)=%$(MASK_FILE_SUFFIX))) \
		--repeat-threshold $(REPEAT_THRESHOLD) \
		--seed-length $(seed_length) \
		--genome-neighbors $(TEMP_DIR)/$(GENOME_NEIGHBORS_PREFIX)$(seed_length)$(GENOME_NEIGHBORS_SUFFIX) \
		$(PARALLEL_SORT) >$(SAFEPIPETARGET)

$(SORTED_REFERENCE_XML): $(CONTIGS_XML) $(ALL_MASK_XMLS) $(ANNOTATION_XML)
	$(CMDPREFIX) $(MERGE_REFERENCES) --make-absolute-paths no $(foreach part, $^, -i '$(part)') -o $(SAFEPIPETARGET)

