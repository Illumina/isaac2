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
## file FindNeighbors.mk
##
## brief Rules for producing KUL annotation file
##
## author Roman Petrovski
##
################################################################################

ifeq (,$(ANNOTATION_MASK_WIDTH))
$(error "ANNOTATION_MASK_WIDTH is not defined")
endif

ifeq (,$(NEIGHBORHOOD_DISTANCE))
$(error "NEIGHBORHOOD_DISTANCE is not defined")
endif

ifeq (,$(CONTIGS_XML))
$(error "CONTIGS_XML is not defined")
endif


ANNOTATION_SEED_LENGTHS:=16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80
ANNOTATION_EDIT_DISTANCES:=0 1 2

ANNOTATION_MASKS:=$(shell $(SEQ) -w 0 $$(( (1<<$(ANNOTATION_MASK_WIDTH)) - 1)))

ANNOTATION_K_MAX:=$(lastword $(ANNOTATION_SEED_LENGTHS))

ANNOTATION_KMER_POSITION_FILE_PREFIX:=neighbor-positions-
ANNOTATION_KMER_POSITION_FILE_SUFFIX:=.dat
ANNOTATION_KMER_POSITION_FILE_XML_SUFFIX:=.xml
ANNOTATION_KMER_POSITION_FILES:=$(foreach s,$(ANNOTATION_SEED_LENGTHS),$(foreach m, $(ANNOTATION_MASKS), $(CURDIR)/$(TEMP_DIR)/$(ANNOTATION_KMER_POSITION_FILE_PREFIX)$(s)-$(m)$(ANNOTATION_KMER_POSITION_FILE_SUFFIX)))
ANNOTATION_KMER_POSITION_XMLS:=$(foreach f, $(ANNOTATION_KMER_POSITION_FILES), $(f:%$(ANNOTATION_KMER_POSITION_FILE_SUFFIX)=%$(ANNOTATION_KMER_POSITION_FILE_XML_SUFFIX)))

$(ANNOTATION_KMER_POSITION_XMLS):mask=$(lastword $(subst -, ,$(subst $(ANNOTATION_KMER_POSITION_FILE_XML_SUFFIX),,$(lastword $(subst $(ANNOTATION_KMER_POSITION_FILE_PREFIX), ,$(notdir $@))))))
$(ANNOTATION_KMER_POSITION_XMLS):seed_length=$(firstword $(subst -, ,$(subst $(ANNOTATION_KMER_POSITION_FILE_XML_SUFFIX),,$(lastword $(subst $(ANNOTATION_KMER_POSITION_FILE_PREFIX), ,$(notdir $@))))))
$(ANNOTATION_KMER_POSITION_XMLS):mask_file=$(@:%$(ANNOTATION_KMER_POSITION_FILE_XML_SUFFIX)=%$(ANNOTATION_KMER_POSITION_FILE_SUFFIX))
$(ANNOTATION_KMER_POSITION_XMLS): $(CONTIGS_XML) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(SORT_REFERENCE) -r $(CONTIGS_XML) --mask-width $(ANNOTATION_MASK_WIDTH) --mask $(mask) \
		--seed-length $(seed_length) \
		--output-file $(mask_file) \
		--repeat-threshold 0 >$(SAFEPIPETARGET)

ANNOTATION_SORTED_REFERENCE_XMLS:=$(foreach s,$(ANNOTATION_SEED_LENGTHS),$(CURDIR)/$(TEMP_DIR)/$(ANNOTATION_KMER_POSITION_FILE_PREFIX)$(s).xml)
.SECONDEXPANSION:
$(ANNOTATION_SORTED_REFERENCE_XMLS): $(CONTIGS_XML) $(foreach m, $(ANNOTATION_MASKS), $(CURDIR)/$(TEMP_DIR)/$(ANNOTATION_KMER_POSITION_FILE_PREFIX)$$(target_seed_length)-$(m)$(ANNOTATION_KMER_POSITION_FILE_XML_SUFFIX))
	$(CMDPREFIX) $(MERGE_REFERENCES)  --merge-annotations no $(foreach part, $^, -i '$(part)') -o $(SAFEPIPETARGET)

edit_distance_permutation_mask_files=\
	$(foreach m,$(ANNOTATION_MASKS),$(CURDIR)/$(TEMP_DIR)/neighbors-$(1)-$(2)-$(m)$(NEIGHBORS_FILE_SUFFIX))
MASK_FILES:=$(foreach seed,$(ANNOTATION_SEED_LENGTHS),\
	$(foreach ed,$(ANNOTATION_EDIT_DISTANCES), \
		$(call edit_distance_permutation_mask_files,$(ed),$(seed))))

target_edit_distance=$(word 2,$(subst -, ,$(notdir $@)))
annotation_seed_length=$(word 3,$(subst -, ,$(word 1,$(subst ., ,$(notdir $(1))))))
target_seed_length=$(call annotation_seed_length,$(notdir $@))

EDIT_DISTANCE_FILES:=$(foreach s,$(ANNOTATION_SEED_LENGTHS),$(foreach d,$(ANNOTATION_EDIT_DISTANCES), $(CURDIR)/$(TEMP_DIR)/neighbors-$(d)-$(s)$(NEIGHBORS_FILE_SUFFIX)))

ifneq (0,$(ANNOTATION_MASK_WIDTH))
target_mask=$(word 4,$(subst -, ,$(subst $(NEIGHBORS_FILE_SUFFIX),,$(notdir $@))))
.SECONDEXPANSION:
$(MASK_FILES): $(CURDIR)/$(TEMP_DIR)/$(ANNOTATION_KMER_POSITION_FILE_PREFIX)$$(target_seed_length).xml $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(FIND_NEIGHBORS) -r $< \
		--seed-length $(target_seed_length) \
		--neighborhood-distance $(target_edit_distance) \
		--mask-width $(ANNOTATION_MASK_WIDTH) \
		--mask $(target_mask) \
		--output-file $(SAFEPIPETARGET)

.SECONDEXPANSION:
$(EDIT_DISTANCE_FILES): $(CONTIGS_XML) $$(call edit_distance_permutation_mask_files,$$(target_edit_distance),$$(target_seed_length))
	$(CMDPREFIX) $(MERGE_ANNOTATIONS) -r $(CONTIGS_XML) --merged-type neighbor-counts \
	-i $(word 2, $^) $(foreach part, $(wordlist 3,$(words $^),$^), -i '$(part)'  --aggregate-function sum) \
		--output-file $(SAFEPIPETARGET)
else
.SECONDEXPANSION:
$(EDIT_DISTANCE_FILES): $(CURDIR)/$(TEMP_DIR)/$(ANNOTATION_KMER_POSITION_FILE_PREFIX)$$(target_seed_length).xml $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(FIND_NEIGHBORS) -r $< \
		--seed-length $(target_seed_length) \
		--neighborhood-distance $(target_edit_distance) \
		--mask-width $(ANNOTATION_MASK_WIDTH) \
		--mask 0 \
		--output-file $(SAFEPIPETARGET)
endif

AGGREGATE_12:=zero-and-leq100
AGGREGATE_012:=zero-and-leq100

NEIGHBORS_12_FILES:=$(foreach s,$(ANNOTATION_SEED_LENGTHS), $(CURDIR)/$(TEMP_DIR)/neighbors-1+2-$(s)$(NEIGHBORS_FILE_SUFFIX))
.SECONDEXPANSION:
$(NEIGHBORS_12_FILES): $(CONTIGS_XML) $(CURDIR)/$(TEMP_DIR)/neighbors-1-$$(target_seed_length)$(NEIGHBORS_FILE_SUFFIX) $(CURDIR)/$(TEMP_DIR)/neighbors-2-$$(target_seed_length)$(NEIGHBORS_FILE_SUFFIX)
	$(CMDPREFIX) $(MERGE_ANNOTATIONS) -r $(CONTIGS_XML) --merged-type neighbor-counts \
	-i $(word 2, $^) \
	-i $(lastword $^) --aggregate-function $(AGGREGATE_12) -o $(SAFEPIPETARGET)

NEIGHBORS_012_FILES:=$(foreach s,$(ANNOTATION_SEED_LENGTHS), $(CURDIR)/$(TEMP_DIR)/neighbors-0+1+2-$(s)$(NEIGHBORS_FILE_SUFFIX))
.SECONDEXPANSION:
$(NEIGHBORS_012_FILES): $(CONTIGS_XML) $(CURDIR)/$(TEMP_DIR)/neighbors-0-$$(target_seed_length)$(NEIGHBORS_FILE_SUFFIX) $(CURDIR)/$(TEMP_DIR)/neighbors-1+2-$$(target_seed_length)$(NEIGHBORS_FILE_SUFFIX)
	$(CMDPREFIX) $(MERGE_ANNOTATIONS) -r $(CONTIGS_XML) --merged-type neighbor-counts \
	-i $(word 2,$^) \
	-i $(lastword $^) --aggregate-function $(AGGREGATE_012) -o $(SAFEPIPETARGET)

ANNOTATION_FILE:=$(NEIGHBORHOOD_DISTANCE)uniqueness$(ANNOTATION_MASK_FILE_SUFFIX)
$(ANNOTATION_FILE): $(CONTIGS_XML) $(NEIGHBORS_012_FILES)
	$(CMDPREFIX) $(MERGE_ANNOTATIONS) -r $(CONTIGS_XML) --merged-type kul \
	-i init-kul-65535 \
	$(foreach nf,$(NEIGHBORS_012_FILES), -i $(nf) --aggregate-function mask-$(call annotation_seed_length,$(nf))) -o $(SAFEPIPETARGET)

