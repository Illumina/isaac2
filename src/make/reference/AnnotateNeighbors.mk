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
## file AnnotateNeighbors.mk
##
## brief File to be included from SortReference.mk
##
## author Roman Petrovski
##
################################################################################

include $(MAKEFILES_DIR)/reference/FindNeighbors.mk

ifeq (,$(SEED_LENGTHS))
$(error "SEED_LENGTHS is not defined")
endif

ifeq (,$(ANNOTATION_FILE))
$(error "ANNOTATION_FILE is not defined")
endif

ifeq (,$(CONTIGS_XML))
$(error "CONTIGS_XML is not defined")
endif

# for setting the neighbors bit in the reference kmers we just need an indicate of whether the number of 1 and 2 distance neighbors is 0 or not
# sum aggregate should work fine
NEIGHBORS_1or2_FILES:=$(foreach s,$(SEED_LENGTHS), $(CURDIR)/$(TEMP_DIR)/neighbors-1or2-$(s)$(NEIGHBORS_FILE_SUFFIX))
.SECONDEXPANSION:
$(NEIGHBORS_1or2_FILES): $(CURDIR)/$(TEMP_DIR)/neighbors-1-$$(target_seed_length)$(NEIGHBORS_FILE_SUFFIX) $(CURDIR)/$(TEMP_DIR)/neighbors-2-$$(target_seed_length)$(NEIGHBORS_FILE_SUFFIX)
	$(CMDPREFIX) $(MERGE_ANNOTATIONS) -r $(CONTIGS_XML) --merged-type neighbor-counts \
	-i $(firstword $^) \
	-i $(lastword $^) --aggregate-function sum -o $(SAFEPIPETARGET)

ALL_GENOME_NEIGHBORS_DAT:=$(foreach s, $(SEED_LENGTHS), $(TEMP_DIR)/$(GENOME_NEIGHBORS_PREFIX)$(s)$(GENOME_NEIGHBORS_SUFFIX))

.SECONDEXPANSION:
$(ALL_GENOME_NEIGHBORS_DAT): $(CURDIR)/$(TEMP_DIR)/neighbors-1or2-$$(target_seed_length)$(NEIGHBORS_FILE_SUFFIX) $(CONTIGS_XML)
	$(CMDPREFIX) $(EXTRACT_NEIGHBORS_FROM_ANNOTATION) -r $(CONTIGS_XML) \
	-i $< -b $(ANNOTATION_BITS) --output-file $(SAFEPIPETARGET)
