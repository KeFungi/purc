# PURC: Pipeline for Untangling Reticulate Complexes

### Overview ###
PURC is a pipeline for extracting alleles from amplicon sequencing data (PacBio, Illumina,...etc), and is geared toward analyzing polyploid species complexes. 

Workflow for PacBio amplicon seq:

* Check concatemers and split them if requested (more on concatemers [here](https://github.com/PacificBiosciences/cDNA_primer/wiki/Artificial-concatemers,-PCR-chimeras,-and-fusion-genes))
* Identify barcodes and remove them
* Trim primers and other adapters
* Assign each sequence to specimen based on the barcode and user-specified reference database
* Cluster and remove chimeric sequences, iteratively
* Woop woop

Note that PURC works on Mac machines. We have not tested it on Linux and PC. [TODO]

### To setup ###
PURC is consist of purc.py and a number of dependencies. We bundled most of the dependencies (cutadapt, muscle and usearch) together in the distribution. To get the dependencies in place, cd to the purc directory, and type: 
```
#!shell
./install_dependencies.sh
```
If you get "permission denied" error, then try this first:
```
#!shell
chmod +x install_dependencies.sh
```

There are however three other dependencies that you have to install yourself:

* [Python](https://www.python.org) - Version 2.7 or later. We have not tested purc on Python 3, and it will probably not work.
* [BioPython](http://biopython.org/wiki/Main_Page) - Version 1.6 or later. You need to have [Numpy](http://www.numpy.org) in place before installing BioPython. Please refer to BioPython [manual](http://biopython.org/DIST/docs/install/Installation.htmlall/Installation.html) for installation instruction.
* [BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - Version 2.2.30 or later. Place the executables in your PATH. If you are using Mac, the easiest way is to download the .dmg file and follow the installer's instruction. To test if this is installed correctly, open Terminal and type "blastn -h" (without quotations). If you see a bunch of stuff pooped out (i.e. "USAGE: ...blah blah blah"), then you are good to go. However, if you get "command not found" error, then BLAST is not installed correctly.  


### To get files ready ###
You will need to get the following things ready:

* **Configuration file** - This contains all the file names and parameters information that PURC needs. There are several examples come with the PURC distribution (open them in TextWrangler or similar text editor).
* **Barcode sequence file** - fasta format. The sequence name is the barcode name.
```c
#!shell
>BC01
ACTACATATGAGATG
>BC02
TCATGAGTCGACACTA
>BC03
TATCTATCGTATACGC
```

* **Reference sequence file** - fasta format. Each sequence name must specify the locus ('locus=') and group ('group=') that this sequence represents. You can also note where is sequence came from ('ref_taxon='). Each designation is separated by /. For example:
```
#!shell
>locus=ApP/group=A/ref_taxon=G_app_BC2
TGCCACACTGGTGAGTATTGTCTTACTTTTTGTTATCCTTTTTCTTGGTGAGAAAGGGTACTGTGCATGGCATATTCACGTCAGAATCCAAGACCCCCGCTTGGGGCCGAGGGTGACAAGGATGTTTCTGTTGGGTGATACCTGTGATGCCAGTTGGAGCAAGAGTAAAATCAACTTTGTAAACATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAATGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAATATATTGACTACTTCATGCATGGAATTCGAGCAGCTATCACATGTTGATGTTTTTTTTTGAGGCGAGGGGTTCTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATGTGCTATTTTTGCATAAATGCTACGTTACAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTTAGTTGATGATTGCATGCGTACACCTTTATGTTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTTAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTGTATTCCTTTTGGATTAAAGTTTATCTCCCTTGTAATTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT
>locus=ApP/group=A/ref_taxon=G_disj_BC20
TGCCACACTGGTGAGTATTGTCTTACTTTTTGTTATCCTTTTTCTTGGTGAGAAAGGGTACTGTGTATGGCATATTCACGTCATAATCCAAGACCCCCGCTTGGGGCTGGGGGGTGACAAGGATGTTTCTGTTGGGTGATACCTGTGATGCCAGTTGGAGCAAGAGTAAAATCAACTTTGTAAACATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAATGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAATATATTGACTACTTCATGCATGGAATTTGAGCAGCTATCACATGTTGATGTTTTTTTTTCAGGCGAGGGGTTCTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATGTGGGTTTTTTGACATGGCTATTTTTGCATAAATGCTACGTTACAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTTAGTTGATGATTGCATGCGTACACCTTTATGCTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTTAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTTTATTTCTTTTGGATTAAAGTTTATCTCCCTTGTATTTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT
```

* **Map files** - tab delimited text file for linking barcodes and groups to specimens. The first column is the barcode (has to be identical to the barcode seq file), the second column is the group designation, and the third is the specimen name. For example:
```
#!shell
BC01	A	A_jap_8703
BC01	B	C_frag_8111
BC02	A	A_ten_4225
```

### To run ###
PURC can be run by: 
```
#!shell
/Users/fayweili/Programs/purc/purc.py purc_configuration.txt > log.txt
```
This assumes that the purc directory is in /Users/fayweili/Programs. **DO NOT** copy purc.py to your working directory; instead, call purc.py from there. Alternative, you can add the purc directory into your PATH, and in this case, you can run by: 
```
#!shell
purc.py purc_configuration.txt > log.txt
```

### Tutorial 1 - PacBio ###


### Tutorial 2 - Illumina dual barcodes ###


### Citation ###
This script relies heavily on USEARCH, MUSCLE, and BLAST.
If this script assisted with a publication, please cite the following papers
(or updated citations, depending on the versions of USEARCH, etc., used).

PURC: 
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
Fay-Wei Li
Carl Rothfels