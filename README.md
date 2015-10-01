# PURC: Polyploids Untangling ... #

### Overview ###
PURC is a pipeline for extracting alleles from amplicon sequencing data (PacBio, Illumina,...etc), and is geared toward analyzing polyploid species complexes. 

Work flow for PacBio amplicon seq:

* Check concatemers and split them if requested (more on concatemers [here](https://github.com/PacificBiosciences/cDNA_primer/wiki/Artificial-concatemers,-PCR-chimeras,-and-fusion-genes))
* Identify barcodes and remove them
* Trim primers and other adapters
* Assign each sequence to specimen based on the barcode and user-specified reference database
* Cluster and remove chimeric sequences, iteratively
* Woop woop

### To setup ###
Purc has several dependencies and we bundled most of them together in the distribution. To setup, cd to the purc directory, and type 
```
#!shell
./install_dependencies.sh
```
If you get "permission denied" error, then 
```
#!shell
chmod +x install_dependencies.sh
```

There are however three other dependencies that you have to install:
* Python - Version 2.7 or later. We have not tested purc on Python 3, and it will probably not work

* Python - Version 2.7 or later. We have not tested purc on Python 3, and it will probably not work.
* [BioPython](http://biopython.org/wiki/Main_Page) - Version 1.6 or later. You need to have [Numpy](http://www.numpy.org) in place before installing BioPython. Please refer to BioPython [manual](http://biopython.org/DIST/docs/install/Installation.htmlall/Installation.html) for installation instruction.
* BLAST+

### To run ###
Usage: ./purc.py configuration_file > out

Example: 
```
#!shell
./purc.py ppp_configuration.txt > summary.txt
```

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