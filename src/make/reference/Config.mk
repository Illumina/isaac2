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
## file Config.mk
##
## brief Common configuration file for all makefiles.
##
## Defines paths, file names, extensions, etc.
##
## author Roman Petrovski
##
################################################################################

MASK_FILE_XML_SUFFIX=.xml
MASK_FILE_SUFFIX=.dat
MASK_TMP_FILE_SUFFIX=.orig

ANNOTATION_BITS:=16
ANNOTATION_MASK_FILE_SUFFIX:=.$(ANNOTATION_BITS)bpb.gz
ANNOTATION_FORMAT:=$(ANNOTATION_BITS)bpb

NEIGHBORS_FILE_SUFFIX:=.16bpb.gz

REPEAT_THRESHOLD:=1000
CURRENT_REFERENCE_FORMAT_VERSION:=6

GENOME_NEIGHBORS_PREFIX:=neighbors-1or2-
GENOME_NEIGHBORS_SUFFIX:=.1bpb

CONTIGS_XML:=$(TEMP_DIR)/contigs.xml
