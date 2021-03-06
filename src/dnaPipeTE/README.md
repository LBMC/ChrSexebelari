dnaPipeTE - README

v.1.3-23/07/2017

>dnaPipeTE (for De-Novo Assembly & Annotation Pipeline for Transposable Elements), is a pipeline tool designed to find, quantify and get annotation of Transposable Elements in small samples of NGS datasets.
It is very usefull to quantify the proportion of TEs in newly sequenced genomes since it do not requier genome assembly and works on small datasets.

more info at: https://lbbe.univ-lyon1.fr/-dnaPipeTE-?lang=en

******Changelog v1.3********
- fix bug for the fastq samble using in blast1
- remove `bin` folder and replace it with the `init.sh` script so that user can make their own installation.
- Clean git repository of larges files

***********************
******Changelog v1.2********

- Estimation of repeat content is now performed on the ratio of aligned bases (bp) on repeat contig over the total number of base sampled, instead of the number of reads mapping / total of read sampled; this produces a better estimate of the repeat content and reduces potential overestimations. In addition, it allows more accurate estimates if the size of reads used as input is variable.
- If different part of one same read match different repeats contigs (e.g. in case adjacent TEs or TE in TE), all bases are retained instead only the one of the best hit.
- New graph "Bases per component" replaces "reads per component"; is very similar to reads per component graph but represent the total amount of bases aligned over the dnaPipeTE contigs.
- Bug fix: in last version, repbase library was not merged to annotated dnaPipeTE contigs for repeat estimates, now it is (as presented in the pipeline [cartoon](http://gbe.oxfordjournals.org/content/7/4/1192/F1.large.jpg)
- New option: "-Trin_glue" to specify a minimum number of reads supporting the joining of kmer contigs during assembly (Chrysalis step in trinity)
- New option: "-contig_length" to set a minimum size (in bp) to report a contig (default is 200 bp)

***********************



<h1>1 - INSTALLATION</h1>

dnaPipeTE is a pipeline relying on different programs to run. The first step
is to install those dependencies.
Execute the content of the `init.sh` script to do so.

If you encounter some issues during installation, do not hesitate to [ask for help on the forum](http://dnapipete.4rumer.com/) !

##System requirement

To date, dnaPipeTE only runs on Linux x64 environments (tested on ubuntu 14.04 PC and Debian 3.2.57-3 x86_64 cluster).

However, Trinity (used for assembly) uses a lot of RAM ! Here are some examples of RAM usages :

- 100 000 reads ~10 Go RAM (two Trinity iterations)
- 3 000 000 reads ~40 Go RAM (two Trinity iterations)

Thus we recommend to use it on assembly-dedicated servers but it could work (if RAM is sufficient) on a PC.

##Dependencies

To run, dnaPipeTE needs the following programs to be installed:

- **[Python 3](https://www.python.org/download/releases/3.1.1)**, including *argparse*, *configparser*, *os*, *re*, *subprocess*, *time*, *sys*, *random*, *ntpath*
- **[Perl 5](https://www.perl.org/)**
- **[R](http://www.r-project.org/index.html)** version 3.0.2 or later (not tested below) including *ggplot2* package.
- **[TRF](http://tandem.bu.edu/trf/trf.download.html)** (Tandem Repeat Finder, see below for installation)

The following dependancies are provided in the package: (./bin/ folder)

- **[GNU Parallel](http://www.gnu.org/software/parallel)** version 3.
- **[Trinity](http://pbil.univ-lyon1.fr/pub/divers/goubert/trinityrnaseq_r20140413p1.tar.gz)** (RNAseq assembly) vers. 2014-04-13
- **[RepeatMasker](http://repeatmasker.org/RMDownload.html)**, including **[RMblastn](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/LATEST)**
- **[blastn](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)** (from blast+ suite)

###Installation

- Download and unzip the dnaPipeTE package in the folder of your choice

- or clone the github repository
```
git clone https://github.com/clemgoub/dnaPipeTE
```

It will create a new directory called dnaPipeTE with all the necessary files. Do not move or modify any of those file.

- Download the TRF executable for Linux X64 environment in the folder dnaPipeTE/bin/
```bash
cd ./bin/
mv trfXXX.linux64.exe trf
```

- Download the formated Repbase database for RepeatMasker at [GIRI](http://www.girinst.org/server/RepBase/index.php)
unpack and move the "Library" folder into `./bin/RepeatMasker` in the dnaPipeTE folder

Your are ready !!!

###Testing dnaPipeTE installation

- To test the installation, a sample file is available either in .fastq or .fastq.gz ( test_dataset.fastq[.gz]). This file is only provided to test if all the differents coponents of dnaPipeTE works well together.

```
cd ~yourdirectory/dnaPipeTE
python3 ./dnaPipeTE.py -input ./test/test_dataset.fastq -output ~/path/to/the/output_folder -genome_size 1000000 -genome_coverage 0.1 -sample_number 1
```
If the pipeline worked properly you should be able to see the 3 output graphs (piechart, bases_per_component and landscapes) with annotations similar to those provided in the folder `dnaPipeTE/test/test_dataset_example_output`.
Since this is a very small dataset used with at low coverage, it is normal that running this test several times won't produce exactly the same estimate of the total amount of repeats.

<h1>2 - RUNNING dnaPipeTE</h1>

##Input File

The input file must be a **single-end FASTQ or FASTQ.GZ** file of NGS reads}. dnaPipeTE do not handle paired-end (we found chimerism issues in PE assembly of TEs).
Typically, your input file is your cleaned sequencing output. dnaPipeTE will sample it (so you can put deep-sequencing data in input) to produce "low coverage" samples for the run (see in the next section for sample size).
Using .fasta file as input works, however, since dnaPipeTE expects .fastq (4 lines per sequence), only half of the file will be considered (2 lines / 4 wont be read in the sampling step).

>**IMPORTANT: We recommend to remove from your reads mitochondrial and other none nucleic DNA, such as known symbionts or guts' bacterias. We found that if mitochondrial reads are left in the samples, RepeatMasker will annotate the corresponding conting to "Gypsy-12_DVir-I LTR/Gypsy" with however a weak score.**

##Run dnaPipeTE

Move into dnaPipeTE folder and type:

```
cd ~yourdirectory/dnaPipeTE
python3 ./dnaPipeTE.py -input ~/path/to/your/input.fastq[.gz] -output ~/path/to/the/output_folder -cpu N -genome_size N -genome_coverage N -sample_number N [...]
```
 >**/!\ VERY IMPORTANT:** run dnaPipeTE from its installation folder, otherwise the dependant scripts won't run. This advice is important especially if you run it into a computing cluster/server, ask first your job to move to the install folder before executing the  python command /!\

**dnaPipeTE arguments:**

|Argument|Type|Description|
|---|---|---|
|-input | INPUT_FILE_PATH | input fastq or fastq.gz files (single end only). It will be sampled |
|-output | OUTPUT_FOLDER | complete path with name for the outputs |  
|-cpu | INTEGER | maximum number of cpu to use |
|-sample_number | INTEGER | number of trinity iterations |
|-genome_size | INTEGER | size of the genome [use it with -genome_coverage; if used, do not use -sample_size] Ex. 175000000 for 175Mb |
|-genome_coverage | FLOAT | coverage of the genome for each sample [use it with -genome_size; if used, do not use -sample_size] Ex: 0.1 for 0.1X coverage per sample |
|-sample_size | INTEGER | number of reads to sample [use without -genome_size and -genome_coverage] |
|-RM_lib | PATH_TO_FILE.fasta | path to custom repeat library for RepeatMasker. The format myst be a valid .fasta file with for each repeat the following name: >Repeat_name#CLASS/Subclass with CLASS in "DNA, LINE, LTR, SINE, MITE, Helitron, Simple Repeat, Satellite" (if not set, default is to use RepeatMasker database) |
|-species | STRING | RepeatMasker library to use. Must be a valid NCBI for species or clade ex: homo,drosophila, "ciona savignyi". By default "All" is used. Do not used with -RM_lib |
|-RM_t |FLOAT| Annotation threshold: minimal percentage of the query (dnaPipeTE contig) aligned on the repeat to keep the anotation from RepeatMasker. Ex: 0.2 for 20% of query in db |
|-keep_Trinity_output | | Keep Trinity output files at the end of the run. Default files are removed (large and numerous).|
|-Trin_glue | INTEGER | number of overlapping reads to join Inchworm (k-mer) contigs in Trinity (default 1) |
|-contig_length | INTEGER | minimum size of a repeat contig to be retained (default 200bp) |

Continuing a crashed run:
>dnaPipeTE is capable to skip some steps if a run crashes after a checkpoint. For example, if it crashes during the Trinity assembly, the sampling won't be performed again if you launch the run again in the same output folder. The checkpoints are for now 1-sampling of Trinity inputs; 2- Trinity assembly. More to follow...

<h1>3 - dnaPipeTE OUTPUTS</h1>

dnaPipeTE produces a lot of outputs, some of them are very interesting.

The outfolder is divided into the following parts:

- **main folder (output name):**

**important files:**

|File|Description|
|---|---|
|	"Trinity.fasta" | this file contains the dnaPipeTE contigs, this is the last assembly performed with Trinity |
|		"reads\_per\_component\_and\_annotation" | table with the count of reads and bp aligned per dnaPipeTE contigs (from blastn 1), as well as its best RepeatMasker annotation. Col1: counts (#reads); Col2: aligned bases; Col3 dnaPipeTE contig name; col4 RepeatMakser annotation; col5 proportion of the dnaPipeTE contig that received the RM hit |
|		"pieChart.pdf/png" | graph with the relative proportion of the main repeat classes, informs about the estimated proportion of repeats in the genome (from blastn 2 and 3) |
|"Bases\_per\_component.pdf/png" | graph with the number of base-pairs aligned on each dnaPipeTE contig (from blast 1), ordered by genome proportion of the dnaPipeTE contig.|
|		"landscapes.pdf" | TEs landscape graphs (TE age distribution). Plot the blastn divergence distribution between reads and the contigs on which they map. |

**less important files you may like:**

|File|Description|
|---|---|
|"Trinity.fasta.out" | raw RepeatMasker output (not sorted) of Trinity.fasta on repbase libraries.|
|"Counts.txt"| count of bp of the sample aligned for each TE class (used for the pieChart)|
|"Reads\_to\_components\_Rtable.txt"| input file to compute the reads and bp per contig (one line per reads)|
|"reads_landscape"| reads used for the landscape graph, including the blastn divergence from one reads to the contig on which it maps.|

- **"Annotation" folder:**

**important files:**

|File|Description|
|---|---|
|	"one_RM_hit_per_Trinity_contigs"| sorted RepeatMasker output containing the best hit on repbase for each of the dnaPipeTE contigs (Trinity.fasta)|
|	"Best_RM_annot_80_80"| subset of the previous table, including contigs for which at least 80% of the sequence is mapping to at least 80% percent of the target sequence.|
|"Best_RM_annot_partial"| same but for contigs for which at least 80% of the sequence is mapping to less than 80% percent of the target sequence|
|"[repeat-class].fasta"| subsets of the Trinity.fasta file for each repeat type detected by RepeatMasker|
|"unannotated.fasta"| subsets of the Trinity.fasta for contigs that didn't find any match...|

- **"blast_out" folder:**

**important files:**

|File|Description|
|---|---|
|"sorted.reads_vs_Trinity.fasta.blast.out"| best hit per reads from blastn 1|
|"sorted.reads_vs_annotated.blast.out"| best hit per reads from blastn 2|
|"sorted.reads_vs_unannotated.blast.out"| best hit per reads from blastn 3|

**less important files you may like:**

|File|Description|
|---|---|
|		"reads_vs_[anything]"| raw blast out from previous files|

- **Trinity_runX**
	Those files contains the raw Trinity outputs and intermediates files produced during assembly steps. For futher detail see the Trinity documentation (http://trinityrnaseq.sourceforge.net/)
