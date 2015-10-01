# README #

### PURC: Polyploids Untangling ... ###
PURC is a pipeline for extracting alleles from amplicon-sequencing data (PacBio, Illumina,...etc), and is geared toward analyzing polyploid species complexes.

* Work flow

* Version


### How do I get set up? ###
Usage: ./ppp.py configuration_file > out
Example: ./ppp.py ppp_configuration.txt > summary.txt
For more info, try: ./ppp.py -help

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Citation ###

This script relies heavily on USEARCH, MUSCLE, and BLAST.
If this script assisted with a publication, please cite the following papers
(or updated citations, depending on the versions of USEARCH, etc., used).

PPP: 
-Awesome paper by carl and fay-wei. Awesome journal. Awesome page numbers.

USEARCH/UCLUST: 
-Edgar, R.C. 2010. Search and clustering orders of magnitude faster than BLAST. 
Bioinformatics 26(19), 2460-2461.

UCHIME:
-Edgar, R.C., B.J. Haas, J.C. Clemente, C. Quince, R. Knight. 2011. 
UCHIME improves sensitivity and speed of chimera detection, Bioinformatics 27(16), 2194-2200.

Cutadapt:
-Martin, M. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. 
EMBnet.journal 17:10-12.

MUSCLE:
-Edgar, R.C. 2004. MUSCLE: Multiple sequence alignment with high accuracy and high throughput. 
Nucleic Acids Research 32:1792-1797.

BLAST: 
-Camacho, C., G. Coulouris, V. Avagyan, N. Ma, J. Papadopoulos, et al. 2009. 
BLAST+: Architecture and applications. BMC Bioinformatics 10: 421.

### Who do I talk to? ###

Carl Rothfels
Fay-Wei Li