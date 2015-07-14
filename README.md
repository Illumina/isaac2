Isaac Genome Alignment Software
Copyright (c) 2010-2014 Illumina, Inc.
All rights reserved.

This software is provided under the terms and conditions of the
GNU GENERAL PUBLIC LICENSE Version 3

You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
along with this program. If not, see
<https://github.com/illumina/licenses/>.

---

#Isaac aligner

To configure and install the product, see the [src/INSTALL](src/INSTALL).

The component algorithms for the Isaac aligner are intended for developers to re-use and improve them. 
This version is not commercially supported and provided as is under [GNU GENERAL PUBLIC LICENSE Version 3](https://github.com/illumina/licenses). 
Isaac is included in the [HiSeq Analysis Software described here](http://support.illumina.com/sequencing/sequencing_software/hiseq-analysis-software-v2-0.html) which is commercially supported, and in the Isaac apps on the [BaseSpace platform](https://basespace.illumina.com/home/index). For more information, see our [blog post](http://blog.basespace.illumina.com/2013/06/04/introducing-fast-free-alignment-and-variant-calling-with-the-isaac-human-whole-genome-sequencing-app/).

See [src/Changes](src/Changes) for the release history details.

#Quick start guide

The following example commands will get you through preparing the reference and aligning some example PhiX data. For more examples and details
please see [src/markdown/manual.md](src/markdown/manual.md)

##Prepare reference genome

This step can take long time on bigger genomes. Normally it is done once per reference and the output is reused in further data analyses.

    $ <isaac>/bin/isaac-sort-reference -g <isaac>/data/share/*/data/examples/PhiX/iGenomes/PhiX/NCBI/1993-04-28/Sequence/Chromosomes/phix.fa -o ./PhiX

##Process fastq data

Analyse a pair of fastq files and produce bam output.

    $ <isaac>/bin/isaac-align -r ./PhiX/sorted-reference.xml -b <isaac>/data/share/*/data/examples/PhiX/Fastq -f fastq --use-bases-mask y150,y150 --variable-read-length yes -m10

