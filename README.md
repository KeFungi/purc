# PURC: Pipeline for Untangling Reticulate Complexes

## **Overview** ##
PURC is a pipeline for extracting alleles or homeologs from amplicon sequencing data (PacBio, Illumina, etc), de-multiplexing them (labeling each sequence with its locus and source sample), and cleaning them (clustering sequences, removing chimeras). It is geared toward analyzing polyploid species complexes but is also effective for other applications; the final output of a full run includes an alignment for each locus with each homeolog or allele sequence in the amplicon data labeled with the source sample information and amount of coverage. 

PURC allows users to extract and analyze all the copies of a given locus present in an amplicon pool, so is particularly useful in cases such as the study of polyploid complexes or large gene families, which were historically tractable only via time- and expense-intensive cloning approaches. Within a single run users can analyse as many accessions as they like, limited only by the desired coverage/allele and the number of sequences that can be uniquely linked back to their source samples. PURC can perform this linking using barcode sequences, locus identity, or phylogenetic information. For example, within a single run an individual barcode can be used multiple times if the locus is different, or if the accessions can be distinguished phylogenetically (i.e., of different genera or other clades that can be distinguished using BLAST). While it is most useful in cases where multiple copies of fairly long sequences have been amplified from individual accessions, it also (when used on PacBio sequencing data) provides considerable cost- and time savings for sequencing classic "single-copy" markers, such as those from the plastid or mitochondrion.


Rough outline of PURC's workflow:

* Optional: scan for concatemers and split them into their component sequences (more on concatemers [here](https://github.com/PacificBiosciences/cDNA_primer/wiki/Artificial-concatemers,-PCR-chimeras,-and-fusion-genes))
* Identify and remove barcode sequences
* Trim primers and other adapter sequences
* Assign each read to its source specimen based on the barcode and a user-specified list of reference sequences
* Cluster sequences and remove chimeric sequences, iteratively
* Woop woop

Note that we designed PURC for Mac OSX. We have not tested it on Linux or PCs.

## **Quick Start** ##
### Step 1: setup ###
PURC consists of purc.py (and two other variations-- xx, xx--that we describe below) and a number of dependencies. We bundled most of the dependencies (cutadapt, muscle and usearch) together in the distribution. To get the dependencies in place, cd to the purc directory, and type: 
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

* [Python](https://www.python.org) - Version 2.7 or later (within Python 2 -- we have not tested purc on Python 3, and it will probably not work).
* [BioPython](http://biopython.org/wiki/Main_Page) - Version 1.6 or later. You need to have [Numpy](http://www.numpy.org) in place before installing BioPython. Please refer to BioPython [manual](http://biopython.org/DIST/docs/install/Installation.htmlall/Installation.html) for installation instruction.
* [BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - Version 2.2.30 or later. Place the executables in your PATH. Provided that you're are using a Mac, the easiest way is to download the .dmg file and follow the installer's instruction. To test if BLAST is installed correctly, open Terminal and type "blastn -h" (without quotation marks). If you see a bunch of stuff pooped out (i.e. "USAGE: ...blah blah blah"), then you are good to go. However, if you get "command not found" error, then BLAST is not installed correctly.  

### Step 2: get files ready ###
PURC requires the following files:

* **Configuration file** - This contains the file names and parameters that PURC needs. There are several examples provided with the PURC distribution (open and edit them in TextWrangler or similar text editor).

* **Barcode sequence file** - In fasta format. PURC uses this file to identify barcodes in the sequence reads. The sequence name is the barcode name. E.g.:
         
        >BC01
        ACTACATATGAGATG
        >BC02
        TCATGAGTCGACACTA
        >BC03
        TATCTATCGTATACGC

* **Reference sequence file** - In fasta format. PURC uses these references to assign reads to loci (and optionally to phylogenetic group). Each reference seq name must specify the locus ('locus=') that this sequence represents; you can optionally note where the sequence came from ('ref_taxon='). The ref_taxon information can contain anything (as long as no spaces are special characters are included) and is not used by PURC--it is only useful in allowing users to keep track of the sequences in their reference file more easily. Each designation is separated by /. For example:

        >locus=ApP/ref_taxon=Cystopteris_bulbifera
        TGCCACACTGGTGAGTATTATCTTACTCTTTTTCTTGGTGAGAAAGGGTAGTGTGCATGGCATATTCATGTCAACCCCCGCTTGGGACCGGGGGCGACAAGGATGTTACTGTTGGGTGATACCTGTGATGCCCAATTGGAGCAAGAGTAAAATCAACTTTGTAAATATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTTAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAACGCTTACGATTTGCCTGCAAGGTAAAAGTTGAACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAACATATTGACTACTTTATCAATGGAATTTGAGCAGCTATCACATGTCGATGTTTATTTTGAGGCAAGGGGTTCTTGTTTGCATGCGGTTGTAGAAATGTCTTATCACATTTCTATTTGGGTTTTTTGACATTGCTATTTTTGCATAAATGCTAAGTTACAAATTAGAATTGTTTACTTGTTTGTTTGGAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTCATTTGATGATTGCATGCGTACACCTTTATGTTCATTTAGGCTATGTTTTGTCAGCTCACAAATTTTTGATGTTAAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTTAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTTTTTTCCTTTTGGATTAAAGTTTATCTCCCTTGTAATTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT
        >locus=GAP/ref_taxon=Cystopteris_bulbifera
        AAAGTAATTTGATGCGTTTTTTTGGTGAAGTCAGTGCATTTTATTCTTTGAAATGGTGAGATCTTATGTTTGCAATCCACTGTTTACCAGGTAGTCCACCAGGAGTTTGGAATTGCCGAGGCCCTCATGACAACTGTGCATGCTACCACTGGTACCGTATTGTCTCCAGCCTCTTTGATTCTGATGTGCTGTAGTTTATACAGAGGCTGGAATTGAATACCCTTGTCATGACTAGCTCTTATATGTCGGTATTTGTGCTTCCAGCTACTCAAAAGACTGTGGATGGTCCCTCTGGCAAAGATTGGCGGGGTGGCCGAGGTGCTGGGCAGAATATCATTCCGAGCTCAACCGGTGCTGCCAAGGTACACTGCTAACTTAACCAAAACCTAAAACATCCTACTAATGTGGTGGCTATTTTATTGCCTTGTAAACATGTTGGTCTTTGTGGTAAACTGCAGGCAGTAGGAAAGGTGCTTCCAGAGTTAAATGGAAAGCTCACTGGAATGGCATTCCGAGTGCCTACCCCCAATGTCTCAGTTGTGGACTTGACTTGCCGCCTAGAGAAAGGGGCATCATACGATACCATCAAAGCAGCTGTGAAGTATGAACATTACATGCTCTAGTCATTGAAACAATTTTTTTTGTCTAAATTAATTTATTCTACTGATGCTGTCTAATTTTTAGCCTACCTGAAAGCTACATCTAATGTGATTTATTACACTTTTTGTAGGGCTGCATCTGAAGGTCCCATGAAAGGAATTCTGGGATACACTGAAGATGATGTTGTATCCACGGACTTTGTCGGTGATGAAAGGTAACTTCCACTTTAGAATTGCAGAAGGTAGTTTCAGTCATGCTCGCTCAAATTCCTCTAACGGTGGACTGGGTTGTTTGCAGATCAAGCATATTTGATGCTAAAGCAGGCATTGCCCTCAGTGATCGGTTTGTAAAGCTTGTTTCGT
	

	When sequences from multiple specimens are pooled together and share the same barcode, PURC can demultiplex them using user-provided "group" information. A reference sequence must be provided for each group, and the groups must be sufficiently phylogenetically distinct that BLAST can match each sequences to a reference sequences from the appropriate group. You need to supply the group information in the reference sequences ('group='). In this example, group "A" is the genus Cystopteris, group "B" is the genus Acystopteris, and group "C" is Gymnocarpium:

        >locus=ApP/group=A/ref_taxon=Cystopteris_fragilis_Utah
        TGCCACACTGGTGAGTATTGTCTTACTTTTTGTTATCCTTTTTCTTGGTGAGAAAGGGTACTGTGCATGGCATATTCACGTCAGAATCCAAGACCCCCGCTTGGGGCCGAGGGTGACAAGGATGTTTCTGTTGGGTGATACCTGTGATGCCAGTTGGAGCAAGAGTAAAATCAACTTTGTAAACATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAATGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAATATATTGACTACTTCATGCATGGAATTCGAGCAGCTATCACATGTTGATGTTTTTTTTTGAGGCGAGGGGTTCTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATGTGCTATTTTTGCATAAATGCTACGTTACAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTTAGTTGATGATTGCATGCGTACACCTTTATGTTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTTAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTGTATTCCTTTTGGATTAAAGTTTATCTCCCTTGTAATTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT
        >locus=ApP/group=B/ref_taxon=Acystopteris_japonica_Taiwan
        TGCCACACTGGTGAGTATTGTCTTACTTTTTGTTATCCTTTTTCTTGGTGAGAAAGGGTACTGTGTATGGCATATTCACGTCATAATCCAAGACCCCCGCTTGGGGCTGGGGGGTGACAAGGATGTTTCTGTTGGGTGATACCTGTGATGCCAGTTGGAGCAAGAGTAAAATCAACTTTGTAAACATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAATGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAATATATTGACTACTTCATGCATGGAATTTGAGCAGCTATCACATGTTGATGTTTTTTTTTCAGGCGAGGGGTTCTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATGTGGGTTTTTTGACATGGCTATTTTTGCATAAATGCTACGTTACAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTTAGTTGATGATTGCATGCGTACACCTTTATGCTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTTAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTTTATTTCTTTTGGATTAAAGTTTATCTCCCTTGTATTTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT
        >locus=ApP/group=C/ref_taxon=Gymnocarpium_dryopteris_MA
        TGCCACACTGGTGAGTATTGTCTTACTTTTTGTTATCCTTTTTCTTGGTGAGAAAGGGTACTGTGTATGGCATATTCACGTCATAATCCAAGACCCCCGCTTGGGGCTGGGGGGTGACAAGGATGTTTCTGTTGGGTGATACCTGTGATGCCAGTTGGAGCAAGAGTAAAATCAACTTTGTAAACATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAATGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAATATATTGACTACTTCATGCATGGAATTTGAGCAGCTATCACATGTTGATGTTTTTTTTTCAGGCGAGGGGTTCTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATGTGGGTTTTTTGACATGGCTATTTTTGCATAAATGCTACGTTACAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTTAGTTGATGATTGCATGCGTACACCTTTATGCTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTTAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTTTATTTCTTTTGGATTAAAGTTTATCTCCCTTGTATTTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT

* **Map files** - tab delimited text file, one for each locus. This tells PURC the correspondence between barcode and specimen. The first column is the barcode (has to be identical to the barcode seq file), and the second is the specimen name. For example:

        BC01	Cystopteris_fragilis_Utah
        BC02	Cystopteris_fragilis_Arizona
        BC03	Cystopteris_fragilis_Taiwan
        
	When one barcode contains multiple specimens, you need to add the group designation in the second column. In this case, one Acystopteris and one Cystopteris specimen both were labeled by barcode 1, but they are divergent enough that PURC can pull them apart (based on the reference sequence file):
	        
        BC01	A	Acystopteris_japonica_Taiwan
        BC01	B	Cystopteris_fragilis_Utah
        BC02	A	Acystopteris_japonica_Japan

    If two barcodes are used (one on each primer), then the first two columns are the barcodes and the third the specimen names:

    	BCF1	BCR1	Cystopteris_fragilis_Utah
    	BCF2	BCR1	Cystopteris_fragilis_Arizona
    	BCF1	BCR2	Cystopteris_fragilis_Taiwan
    	BCF2	BCR2	Cystopteris_fragilis_AMagicPlace

### Step 3: run ###
PURC can be run by: 
```
#!shell
/Users/fayweili/Programs/purc/purc.py purc_configuration.txt
```
This assumes that the purc directory is in /Users/fayweili/Programs. **DO NOT** copy purc.py to your working directory; instead, call purc.py from there. Alternative, you can add the purc directory into your PATH, and in this case, you can run by: 
```
#!shell
purc.py purc_configuration.txt
```

### Example 1 - PacBio ###
In this dataset, four loci were amplified from 30 Cystopteris specimens. All the specimens were labeled with different barcodes, and all the PCR reactions were pooled together and sequenced in one PacBio SMRT cell. Because one barcode corresponds to one specimen, in the configuration file we specify: 
```
Multiplex_per_barcode	= 0
```


### Example 2 - PacBio ###
In this dataset, four loci were amplified from xx Cystopteris, xx Acystopteris and xx Gymnocarpium specimens. Because we don't have enough barcoded primers, we multiplex each barcode to cover one Cystopteris, one Acystopteris and one Gymnocarpium. For example, G_dry_7000, A_ten_4225 and C_mil_6761 are all labeled as BC03, but were assigned as different "groups" (A, B, C) in the map and reference files. In the configuration file we need to specify: 
```
Multiplex_per_barcode	= 1
```

And in the map files we need to add the group designation in the second column:
```
BC03	A	G_dry_7000
BC03	B	A_ten_4225
BC03	C	C_mil_6761
...
```
Also for each group of each locus, we need to include the reference sequences, so that purc can tell them apart:
```
>locus=ApP/group=A/ref_taxon=Gymnocarpium
TGCCACACTGGTGAGTATTGTCTTACTTTTTGTTATCCTTTTTCTTGGTGAGAAAGGGTACTGTGCATGGCATATTCACGTCAGAATCCAAGACCCCCGCTTGGGGCCGAGGGTGACAAGGATGTTTCTGTTGGGTGATACCTGTGATGCCAGTTGGAGCAAGAGTAAAATCAACTTTGTAAACATCATCTATTTGAAGGATTAACAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAATGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAATATATTGACTACTTCATGCATGGAATTCGAGCAGCTATCACATGTTGATGTTTTTTTTTGAGGCGAGGGGTTCTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATGTGCTATTTTTGCATAAATGCTACGTTACAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTTAGTTGATGATTGCATGCGTACACCTTTATGTTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTTAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGCTCAAAATTAAAGGAGAGTGTATTCCTTTTGGATTAAAGTTTATCTCCCTTGTAATTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT-------------------------------------------------------------------------------------
>locus=ApP/group=B/ref_taxon=Acystopteris
TGCCACACTGGTGAGTATTGTCTTACTCTTTTTCTTGGTGAGAAAGGGTAAAGTGCATGGCATATTCACGTCAACCCCCGCTTGGGGCCGGGAGTGACAAGGATGTTACTGTTGGGTGATACCTGTGATGCCCAGTTGGAGCAGGAGTAAAATCGACTTTGTAAATATCATTTATTTGAGGGATTAACAGACATGGTATTTAAATTCCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAACGCTTACGATTTGCCTGCAAGGTAAAAGTTGCACAATGCTCAAGGTGGGGCTAGTTCTTGTCACGTAAGCAAGGATCATCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAACATATTGACCACTTCATCCATGGAATTTGAGCAGCTATCACATGTCGATGTTTTTTTTTGAGGCGAGGGGTTCTTGTTTGCATGTGGTTGTAGAAATGTTTTATCACATTTCTATTCGGATTTTTTGACATGGCTATTTTTGCATAAATGCTAAGTTCCAAATTAGAATTGTTTACTTGTTTGTTTGTAGGAAATCTCAAACGACTGTCTTTTGCTCTTGTATGCTCATTTGATGATTGCATGCGTACACCTTTATGTTCATTTCAGGCTATGTTTTGTCAGCTCACAAGTTTTTGATGTTAAACCTAACATGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCATTGTTACCGATCTTGTTTAGAGACTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAAACAGATCAAAATTAAAGAGAGTTTTTTCCTTTTGGATTAAAGTTTATCTCCCTTGTAATTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT------------------------------------------------------------------------------------------
>locus=ApP/group=C/ref_taxon=Cystopteris1
TGCCACACTGGTGAGTATTATCTTACTCTTTTTCTTGGTGAGAAAGGGTAGTGTGCATGGCATATTCATGTCAACCCCCGCTTGGCACCGGGGGTGATAAGGATGTTACTGTTGGGTGATACCTGTGATGCCCAATTGGAGCAAATAGTAAAATCAACTTCGTAAATATCATTTATTTGAAGGATTAATAGACATGGTATTTAAATTCCTCTCACATTCAAAACAGGGTGGTTGCAGAACTGGTATGGCCAAAGTAACGAACGCTTACGATTTGCCTGCAAGGTAAAAGTTGAACAATGCTCAAGGTGGGGCTACTCCTTTTGTCACTTAAGCAAGGATCTTCAAGCATGTAAAATTATTCTCCCTCAACTTTGCTTTACAAAAGAAATTTAACATATTGACTACTTCATCCATGGAATTTGAGCAGCTACACATGTAGATGTTTTTTTTTTGAGGCAAGGGGTTCTTGTTTGCATGTGGTTGTAGAAATGTTTATCACATTTCTATTTGGGTCTTTTGACATGGCTATTTTTGCATAACTGCCAAGTTACAAATTAGAATTGTTCACTTGCTGTTTGTACGAAATCTCAATTGACTGTCTTTTGCTCTTGTATGCTCATTTGATGATTGCATGCGTACACCTTTATGTTCATTTCAGGCTATGTTTTGTCAGCTCACAAGATTTTGATGCTAAACCTAACGTGACAGGAAAGTTATACATACTGTTGGTCCAAGATATGCTGTAAAATATCATACAGCTGCAGAAAATGCTCTAAGTCGTTGCTACCGATCTTGTTTAGAGGCTTTGATTGACTTAGGCCTTCAAAGGTACCAGCTGCTTGTTTAACAGCTCAAAATTAAAGGAGAGTTTTTTCCTTTTGGATTAAAGTTTATCTCCCATGTAATTCTTGCAGCATTGCCCTGGGGTGTATTTACACAGAGTCTAAAGGCTAT---------------------------------------------------------------------------------------
...
```
        
### Example 3 - Illumina ###
PURC can work with Illumina data too. You'd need to merge the pair-end reads first, for example, using [PEAR](https://github.com/xflouris/PEAR). In this dataset, each sequence is labeled by two barcodes, one on the forward and the other on the reverse primer. In the map file, we tell PURC which barcode combination corresponds to which specimen. We invoke dual-barcode mode by:
```
Dual_barcode		= 1
```
And in the map file:
```
BCF1	BCR1	species1
BCF2	BCR1	species2
BCF3	BCR1	species3
...
```

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