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
## file SortReference.mk
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

include $(MAKEFILES_DIR)/reference/AnnotateNeighbors.mk

ifeq (,$(MASK_WIDTH))
$(error "MASK_WIDTH is not defined")
endif

ifeq (,$(GENOME_FILE))
$(error "GENOME_FILE is not defined")
endif

ifeq (,$(SEED_LENGTHS))
$(error "SEED_LENGTHS is not defined")
endif

SORTED_REFERENCE_XML:=sorted-reference.xml

MASK_FILE_PREFIX:=kmer-positions-
MASKS:=$(shell $(SEQ) -w 0 $$(( (1<<$(MASK_WIDTH)) - 1)))
ALL_MASKS:=$(foreach s, $(SEED_LENGTHS), $(foreach m, $(MASKS), $(MASK_FILE_PREFIX)$(s)-$(m)$(MASK_FILE_SUFFIX)))
ALL_MASK_XMLS:=$(foreach f, $(ALL_MASKS), $(f:%$(MASK_FILE_SUFFIX)=$(TEMP_DIR)/%$(MASK_FILE_XML_SUFFIX)))

$(ALL_MASK_XMLS):mask=$(lastword $(subst -, ,$(subst $(MASK_FILE_XML_SUFFIX),,$(lastword $(subst $(MASK_FILE_PREFIX), ,$@)))))
$(ALL_MASK_XMLS):seed_length=$(firstword $(subst -, ,$(subst $(MASK_FILE_XML_SUFFIX),,$(lastword $(subst $(MASK_FILE_PREFIX), ,$@)))))
$(ALL_MASK_XMLS):mask_file=$(notdir $(@:%$(MASK_FILE_XML_SUFFIX)=%$(MASK_FILE_SUFFIX)))

$(CONTIGS_XML): $(GENOME_FILE) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(PRINT_CONTIGS) -g $(GENOME_FILE) >$(SAFEPIPETARGET)

ifneq (,$(ANNOTATION_MASKS))
.SECONDEXPANSION:
$(ALL_MASK_XMLS): $(CONTIGS_XML) $(TEMP_DIR)/.sentinel $(TEMP_DIR)/$(GENOME_NEIGHBORS_PREFIX)$$(seed_length)$(GENOME_NEIGHBORS_SUFFIX)
	$(CMDPREFIX) $(SORT_REFERENCE) -r $(CONTIGS_XML) --mask-width $(MASK_WIDTH) --mask $(mask) \
		--seed-length $(seed_length) \
		--output-file $(CURDIR)/$(mask_file) \
		--repeat-threshold $(REPEAT_THRESHOLD) \
		--genome-neighbors $(lastword $^) >$(SAFEPIPETARGET)

ANNOTATION_XML:=$(TEMP_DIR)/annotation.xml

$(ANNOTATION_XML) : $(ANNOTATION_FILE)
	$(CMDPREFIX) echo "\
	<SortedReference>\
		<FormatVersion>$(CURRENT_REFERENCE_FORMAT_VERSION)</FormatVersion>\
		<Annotations>\
			<Annotation Type='k-uniqueness' K='$(NEIGHBORHOOD_DISTANCE)'>\
				<Format>$(ANNOTATION_FORMAT)</Format>\
				<File>$(CURDIR)/$(ANNOTATION_FILE)</File>\
			</Annotation>\
		</Annotations>\
	</SortedReference>\
	" >$(SAFEPIPETARGET)

$(SORTED_REFERENCE_XML): $(CONTIGS_XML) $(ANNOTATION_XML) $(ALL_MASK_XMLS)
	$(CMDPREFIX) $(MERGE_REFERENCES) $(foreach part, $^, -i '$(part)') -o $(SAFEPIPETARGET)

all: $(SORTED_REFERENCE_XML) $(ANNOTATION_WIG)
	$(CMDPREFIX) $(LOG_INFO) "All done!"

else
$(ALL_MASK_XMLS): $(CONTIGS_XML) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(SORT_REFERENCE) -r $(CONTIGS_XML) --mask-width $(MASK_WIDTH) --mask $(mask) \
		--seed-length $(SEED_LENGTH) \
		--output-file $(mask_file) \
		--repeat-threshold $(REPEAT_THRESHOLD) >$(SAFEPIPETARGET)

$(SORTED_REFERENCE_XML): $(CONTIGS_XML) $(ALL_MASK_XMLS)
	$(CMDPREFIX) $(MERGE_REFERENCES) $(foreach part, $^, -i '$(part)') -o $(SAFEPIPETARGET)

all: $(SORTED_REFERENCE_XML)
	$(CMDPREFIX) $(LOG_INFO) "All done!"

endif


