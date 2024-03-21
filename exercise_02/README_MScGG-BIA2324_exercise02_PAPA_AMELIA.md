---
### BEGIN-Of-YAML-Block ###
#
## ######################################################################################
##
##   README_MScGG-BIA2324_exercise02_SURNAME_NAME.md
##
##   A LaTeX-extended MarkDown template for MScGG-BIA practical exercise submissions.
##
## ######################################################################################
##
##                 CopyLeft 2023/24 (CC:BY-NC-SA) --- Josep F Abril
##
##   This file should be considered under the Creative Commons BY-NC-SA License
##   (Attribution-Noncommercial-ShareAlike). The material is provided "AS IS", 
##   mainly for teaching purposes, and is distributed in the hope that it will
##   be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
##   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
## ######################################################################################
#
# The current execise number
thyexercise: 02
#
# title-meta is title string without any LaTeX format to be used as pdftitle, part of emails subject...
title-meta: MScGG-BIA Exercise 02
#
# title is the big title for the cover page, fully LaTeX formated to fit into a shortstack command...
title: |
  \textsc{MScGG-Advanced Bioinformatics}
  \textsc{Practical Report}
#subtitle:
#
# runtitle is the running header or footer, used i.e. by fancyheadings...
runtitle: |
  MScGG-BIA Exercise 02 Report
#
# author-meta sets the pdfauthor variable...
author-meta: !!str 'Amelia PAPA @ MScGG BIA Advanced Bionformatics'
#
# authors to appear on the title page...
author:
- name: Amelia PAPA
  myemail: apapapap30
  mydomain: alumnes.ub.edu
#
# IMPORTANT: we need to define email as two fields (user/domain)
#            to avoid parsing email by YAML, as in the following example.
#
# - name: Josep F Abril
#   myemail: jabril   #
#   mydomain: ub.edu  # the real complete email address was in this case: jabril@ub.edu
#
# authorshort defines a brief list of authors for headings, i.e.: Abril, JF
authorshort: Papa, AP
#
# template formating variables...
papersize: a4paper
fontsize: 10pt
geometry: margin=1.5cm
toc: true
lof: true
lot: true
colorlinks: true
urlcolor: blue
citecolor: green
linkcolor: red
#
# LaTeX definitions (available for the main document)
further-defs: 
-  \def\BIAVC{\href{https://campusvirtual.ub.edu/course/view.php?id=77505}{BIA Advanced Bioinformatics Virtual Campus at UB}}
-  \def\ScerS{\textit{Saccharomyces cerevisiae} (strain S288C)}
-  \def\scerS{\textit{S.~cerevisiae} (strain S288C)}
-  \def\Scer{\textit{Saccharomyces cerevisiae}}
-  \def\scer{\textit{S.~cerevisiae}}
-  \def\sra{\textsc{SRA}}
-  \def\Sra{\textsc{Short Reads Archive}}
-  \def\SRA{\href{https://www.ncbi.nlm.nih.gov/sra}{\sra}}
#
### End-Of-YAML-Block ###
---

```{=comment}
We do not need the comment LaTeX environment to hide a block,
with pandoc >2.x one can define "comment" blocks, which are even
escaped and not passed to LaTeX and pandoc (achieving the same hidden effect).
\begin{comment}
\begin{comment} 
The \texttt{comment} \LaTeX\ environment allow us to
include any piece of text, code, etc, but it wil not be included into
the final PDF report when compiling the \texttt{MarkDown} file. You
can open/close one of such environments at any time if you need them.
\end{comment}
```


# Introduction

Different high-throughput sequencing experiments have been performed
over genomic DNA samples of \Scer\ and paired-end raw reads were
provided. We were asked to assemble the reads obtained from at least
two of the suggested datasets. Those sequence sets may have
differences in sequencing methodology, but mainly they differ in
whole-genome coverage. We can evaluate if differences in coverage can
have an impact on the final assembled genome as we already have a
reference genome.


## Objectives

* To check qualities and other properties of sequencing reads.
* To run an assembly protocol: cleaning raw reads with `trimmomatic`
  and assembling the contigs to reconstruct a small genome using `SOAPdenovo`.
* To map the reads back to the assembly so we can visualize the
  coverage across the sequences and estimate the insert size.
* To compare our assembly against a reference genome, so we can assess
  the performance of the assembler and the protocol.
* To check completeness of our assembly with `BUSCO`.
* We introduce some \LaTeX\ examples for citing paper references as footnotes.


## Pre-requisites

Here we provide the commands required to install some software
pre-requisites on our machines.

We may need to install `pandoc`, the following command will suffice for today session:

```{.sh}
sudo apt-get install pandoc                    \
                     texlive-latex-recommended \
                     texlive-latex-extra       \
                     texlive-fonts-recommended \
                     texlive-fonts-extra
```

In case you want to play with \LaTeX, I will recommend you to install
the complete set with this command:

```{.sh}
sudo apt-get install texlive-full
```

You can also install optional packages, such a text editor with
programming facilities and extensions, like `emacs` or `geany` (you
can also use `sublime`, `atom`, `gedit`, ...):

```{.sh}
sudo apt-get install emacs geany vim
```

However, if you are experiencing \LaTeX or `pandoc` formatting errors,
try to install the latest `pandoc` version. Just follow the
instructions from: https://pandoc.org/installing.html

```{.sh}
# On a Debian/Ubuntu/Mint box, you will probably visit
# pandoc's releases page, to get the latest version from repository:
#     https://github.com/jgm/pandoc/releases/latest
# Then, you can run those two commands:
wget https://github.com/jgm/pandoc/releases/download/3.1.11.1/pandoc-3.1.11.1-1-amd64.deb

sudo dpkg -i pandoc-3.1.11.1-1-amd64.deb

# On a MacOSX box, you can download the following installer for the last version:
#    https://github.com/jgm/pandoc/releases/download/3.1.11.1/pandoc-3.1.11.1-arm64-macOS.pkg
#                                                          (ARM architecture on new Macs)
#    https://github.com/jgm/pandoc/releases/download/3.1.11.1/pandoc-3.1.11.1-x86_64-macOS.pkg
#                                                          (Intel architecture on old Macs)
```

When using `pandoc` version greater than `2.x` we will be able to
apply further macros and tags on our report file (like embed `LaTeX`
blocks).


## Initialization

First of all we need to download the exercise tarball from the \BIAVC,
unpack such file, modify the files accordingly to the user within the
exercise folder, and set it as the current working directory for the
rest of the exercise...

```{.sh}
# You probably have already done this step.
tar -zxvf MScGG-BIA2324_exercise_02.tgz
cd exercise_02

# Open project vars config file using your text editor of choice
# (for instance vim, emacs, gedit, sublime, atom, ...);
# fix "NAME" and "SURNAME" placeholders on it
# and save those changes before continuing.
geany projectvars.sh &
#50676

# Let's start with some initialization.
source projectvars.sh
echo $WDR $NM

# Rename report file including your "NAME" and "SURNAME"
mv -v README_MScGG-BIA2324_exercise02_SURNAME_NAME.md \
      README_MScGG-BIA2324_exercise02_${NM}.md

# Now you are ready to play with the protocol in the README file
geany README_MScGG-BIA2324_exercise02_yourSurname_yourName.md &
#50703

# You can test if we can compile the MarkDown document.
# You probably must install some dependencies yet...
runpandoc
```

## Installing Bioinformatics software

You can find below examples on how to install a software tool from
different repositories or systems (you should have `emboss` installed
from previous exercises though).

```{.sh}
#################################
# emboss - European molecular biology open software suite

# on a debian/ubuntu/mint linux system (DEBs)
apt-cache search emboss     # to check if there is such a package
sudo apt-get install emboss # to install such a package

# on a redhat/fedora/centos linux system (RPMs)
yum search emboss           # to check if there is such a package
su -c 'yum install emboss'

# on a SUSE/openSuse  linux system
zypper search "emboss"
sudo zypper install emboss

# on a Mac system using anaconda packages (https://conda.io/docs/index.html)
conda search emboss
sudo conda install -c bioconda emboss

# on a Mac system using mac ports (https://guide.macports.org/)
port search emboss
sudo port install emboss

# you can also install the package if available for the CygWin environment
# running on a Windows box (hhtp://www.cygwin.com/)

# add your packaging system here if you have not used any of the above commands...
```

From now on, we assume that you are using a `Debian`-based linux
distribution (like `Ubuntu` or `Mint`), so we will show only the
corresponding set of commands for that distribution. Therefore, you
can find here a series of commands to install the tools to be used on
this exercise.


```{.sh}
#################################
# NCBI SRA Toolkit https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-ubuntu64.tar.gz \
     -O $WDR/bin/sratoolkit.3.0.10-ubuntu64.tar.gz
     
pushd $WDR/bin
tar -zxvf sratoolkit.3.0.10-ubuntu64.tar.gz
cd sratoolkit.3.0.10-ubuntu64
popd

export PATH=$WDR/bin/sratoolkit.3.0.0-ubuntu64/bin:$PATH
# NOTE: you can add the above setting to the projectvars.sh file
```
All necessary software will be installed through the environment.yml file
```{.sh}
# seqtk - sampling, trimming, fastq2fasta, subsequence, reverse complement
brew install seqtk

# jellyfish - count k-mers in DNA sequences
brew install jellyfish

# fastqc - quality control for NGS sequence data
brew install fastqc

# trimmomatic - flexible read trimming tool for Illumina NGS data
brew install trimmomatic

# samtools - processing sequence alignments in SAM and BAM formats
# bamtools - toolkit for manipulating BAM (genome alignment) files
# picard-tools - Command line tools to manipulate SAM and BAM files
brew install samtools 
brew install bamtools
brew install picard-tools

# Downloading the latest version of picard
# see https://github.com/broadinstitute/picard/releases/latest
wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar \
     -O $BIN/picard.jar

# bwa - Burrows-Wheeler Aligner
# bowtie2 - ultrafast memory-efficient short read aligner
brew install bwa bowtie2

# soapdenovo2 - short-read assembly method to build de novo draft assembly
# The github repository for soapdenovo said to install from brew, but it didn't work
# for me. 

conda install -c bioconda soapdenovo2
soapdenovo-63mer # checks installation

# igv - Integrative Genomics Viewer. Already installed.

# ncbi-blast+ - next generation suite of BLAST sequence search tools
brew install blast

# gnuplot-qt - a portable command-line driven interactive data and function plotting utility
#              It is required for mummerplot later on...
brew install qt # gnuplot depends on qt
brew install gnuplot

# graphicsmagick - collection of image processing tools (replacement for imagemagick)
brew install imagemagick 
 
# mummer - Efficient sequence alignment of full genomes
brew reinstall mummer
#
# IMPORTANT: Some fixes are needed prior to run gnuplot within mummerplot
#            (install mummer package first)
#
#
#   + gnuplot fails because some instruction was implemented in later versions
perl -i -pe 's/^\(.*set mouse.*\)$/#\1/
            ' /usr/local/bin/mummerplot
#
#   + mummerplot cannot find gnuplot:
perl -i -pe 's/system (\"gnuplot --version\")/system (\"\/usr\/bin\/gnuplot --version\")/
            ' /usr/local/bin/mummerplot
perl -i -pe 's/my \$cmd = \"gnuplot\";/my \$cmd = \"\/usr\/bin\/gnuplot\";/
            ' /usr/local/bin/mummerplot
#
#   + just in case mummerplot returns error: Inappropriate ioctl for device
sudo ln -vfs /usr/local/bin/gnuplot-qt /etc/alternatives/gnuplot;
#
# Further NOTES:
#   You probably need to run all those fixes if your system has already gnuplot version > 4.
#   MacOS users lacking sed can try with "perl -i -pe" instead of "sed -i".

#
# #       If using conda, you can try:
#             conda create -n busco5 -c conda-forge -c bioconda busco=5.4.3
#       Then, run it after activating the environment:
#             conda activate busco5
#       And remember to set up all project vars again as you enter into a new shell.
#
cd $BIN/
git clone https://gitlab.com/ezlab/busco.git
cd busco/
sudo python3 setup.py install
cd $WDR
```

### Using `conda/mamba` environments:

Another way to install the software required to complete the exercises
is to use `conda` environments. You can install `conda` following the
instructions from [this link](https://conda.io/projects/conda/en/latest/user-guide/install/index.html);
you can also use `mamba` instead, which is a compact and faster
implementation of `conda`, from the instructions at [this
link](https://github.com/conda-forge/miniforge#install). Once you have
one of those environment managers installed, you can follow the
commands in the next code block to create the `BScBI-CG2324_exercises`
environment and activate it. __You probably have the conda environment
created from the previous exercise, then you can jump to the next
block of code.__

```{.sh}
#
# ***Important***: ensure that you run the create command
#                  outside any other environment (even the `base` one),
#                  for a fresh install of the proper dependencies.
#
# If you have conda instead of mamba already installed on your system
# you can just replace 'mamba' by 'conda' on the commands below:
conda env create --file environment.yml

# Now you can run the tools installed on that environment by activating it:
conda activate MScGG-BIA2324_exercises

# Remember that each time you deactivate a conda environment
# all shell variables defined inside will be lost
# (unless they were exported before activating the conda environment).
# Anyway, you can reload project vars with:
source projectvars.sh

# To return to the initial terminal state, you must deactivate the environment:
mamba deactivate
```

You can review the contents of the environment YAML file at the
Appendices (see section \ref{prg:environmentYML} on page
\pageref{prg:environmentYML}),

Let's start with the analyses, and may the shell be with you...


## Datasets

The budding yeast
\href{https://www.ncbi.nlm.nih.gov/genome/?term=txid559292[Organism:noexp]}{\Scer}
is one of the major model organisms for understanding cellular and
molecular processes in eukaryotes. This single-celled organism is also
important in industry, where it is used to make bread, beer, wine,
enzymes, and pharmaceuticals. The \scer\ genome has approximately
\qty{12}{\Mbp}, distributed in 16 chromosomes. Here, we are going to
work on raw reads from one of selected four different sequencing
experiments that have been already submitted to the \Sra\ (\SRA) and
we will compare the resulting assemblies against the reference
\scerS. Table \ref{tbl:srareadsets} summarizes the information about
those sets.

```{=latex}
% if pandoc complains or it is not running LaTeX code as expected,
% this block can also be saved into a docs/fig_kmercountdist.tex file
% and then loaded into main doc using the LaTeX "\input" macro.
\begin{table}[!ht]
\begin{center}
\begin{scriptsize}
\begin{tabular}{l|cllc}
%
Sequence set  &
SRA accession &
Sequencer     &
Run info      &
Run ID        \\\hline
%
Original WT strain &
\href{https://www.ncbi.nlm.nih.gov/sra/SRX3242873[accn]}{SRX3242873} &
Illumina HiSeq 4000 &
2.8M spots, 559.2M bases, 221.2Mb sra &
\href{https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6130428}{SRR6130428} \\
%
S288c &
\href{https://www.ncbi.nlm.nih.gov/sra/SRX1746300[accn]}{SRX1746300} &
Illumina HiSeq 2000 &
5.3M spots, 1.1G bases, 733.6Mb sra &
\href{https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR3481383}{SRR3481383} \\
%
MiSeq PE-sequencing of SAMD00065885 &
\href{https://www.ncbi.nlm.nih.gov/sra/DRX070537[accn]}{DRX070537} &
Illumina MiSeq &
5.6M spots, 3.3G bases, 1.6Gb sra &
\href{https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=DRR076693}{DRR076693} \\
%
S288c genomic DNA library &
\href{https://www.ncbi.nlm.nih.gov/sra/SRX4414623[accn]}{SRX4414623} &
Illumina HiSeq 2500 &
43.1M spots, 8.6G bases, 3.2Gb sra &
\href{https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR7548448}{SRR7548448} \\
%
\end{tabular}
\end{scriptsize}
\parbox{0.75\linewidth}{%
 \caption[Summary of raw reads datasets that we can use from \texttt{SRA} database.]{%
  \label{tbl:srareadsets}\textbf{Summary of raw reads datasets that we can use from \texttt{SRA} database}. Smaller set will be used for the examples, and you can choose another one to complete the exercise.
 }%caption
}%parbox
\end{center}
\end{table}
```

We will get first the reference genome\footnote{Engel et al. "The
Reference Genome Sequence of \textit{Saccharomyces cerevisiae}: Then
and Now". G3 (Bethesda), g3.113.008995v1, 2013
(\href{https://www.ncbi.nlm.nih.gov/pubmed/?term=24374639}{PMID:24374639}).},
\scerS (assembly `R64.4.1`). This genome version has 16 chromosome
sequences, including the mitochondrion genome, totaling
\qty{12157105}{\bp} (see Table \ref{tbl:yeastchrs}).


```{.sh}
mkdir -vp $WDR/seqs

URL="https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases"
wget $URL/S288C_reference_genome_Current_Release.tgz \
     -O $WDR/seqs/S288C_reference_genome_Current_Release.tgz

pushd $WDR/seqs/
tar -zxvf S288C_reference_genome_Current_Release.tgz
popd

# you may consider copying those vars to your projectvars.sh file
export REFDIR="S288C_reference_genome_R64-4-1_20230830"
export REFGEN="S288C_reference_sequence_R64-4-1_20230830"

brew install gawk
gzcat $WDR/seqs/$REFDIR/$REFGEN.fsa.gz | \
  infoseq -only -length -noheading -sequence fasta::stdin  2> /dev/null | \
  gawk '{ s+=$1 } END{ printf "# Yeast genome size (v.R64.4.1): %dbp\n", s }' 
# Yeast genome size (v.R64.4.1): 12157105bp
```

```{.sh}
#
### NOTE ###
#
# You can include a LaTeX table with the chromosome sizes and GC content,
# which can be produced with a command like this:
#
gzcat $WDR/seqs/$REFDIR/$REFGEN.fsa.gz | \
  infoseq -noheading -sequence fasta::stdin  2> /dev/null | \
    gawk 'BEGIN {
            printf "%15s & %10s & %12s & %10s \\\\\n",
                   "Chromosome", "GenBank ID", "Length (bp)", "GC content";
          }
          $0 != /^[ \t]*$/ {
	        L=$0;
	        sub(/^.*\[(chromosome|location)=/,"",L);
	        sub(/\].*$/,"",L);
	        sub(/_/,"\\_",$3);
            printf "%15s & %10s & %12d & %8.2f\\%% \\\\\n",
                   L, $3, $6, $7;
          }' > $WDR/docs/chromosomes_info.tex;
#
###### NOTE ###### -------------------------------------------------------
#
```

```{=latex}
% if pandoc complains or it is not running LaTeX code as expected,
% this block can also be saved into a docs/fig_kmercountdist.tex file
% and then loaded into main doc using the LaTeX "\input" macro.
\begin{table}[!ht]
\begin{center}
\begin{tabular}{rrrr}
%%%
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  YOUR TABLE HERE...  %%%%%%
\csname @@input\endcsname docs/chromosomes_info
%%%
\end{tabular}
\parbox{0.75\linewidth}{%
 \caption[Reference \Scer\ chromosome summary.]{%
  \label{tbl:yeastchrs}\textbf{Reference \Scer\ chromosome summary.} This table will show information about length and GC content for the chromosomes of the \scer\ reference genome.
 }%caption
}%parbox
\end{center}
\end{table}
```

In order to run the exercise, we can just choose the the smaller
dataset, `SRR6130428`, which has a 221.2Mb `sra` file to be
downloaded. You can choose one of the four suggested files, **but be
aware of your hardware limitations**; moreover all the code examples
will be based on the smaller dataset, so that you will need to change
many stuff across the exercise vcode blocks..

```{.sh}
## The example reads dataset: SRR6130428

# 221.2Mb downloads
mkdir -vp $WDR/seqs/SRR6130428
brew install sratoolkit
prefetch -v SRR6130428  -O $WDR/seqs
#
# Ensure that the SRA file is stored as $WDR/seqs/SRR6130428/SRR6130428.sra,
#   not as $WDR/seqs/SRR6130428/SRR6130428/SRR6130428.sra
#   (you can move it to parent folder if needed).
#
# NOTE: Old versions of prefetch do not implement option "-O",
#       so that they can return an error. Just remove that option
#       and its argument, and look for the downloaded files
#       at the following folder: $HOME/ncbi/public/sra
#       You probably will need to move the file with these commands:
#         mkdir -vp $WDR/seqs/SRR6130428/;
#         mv -v $HOME/ncbi/public/sra/SRR6130428.sra $WDR/seqs/SRR6130428/;

fastq-dump -X 2 -Z $WDR/seqs/SRR6130428/SRR6130428.sra
# Read 2 spots for seqs/SRR6130428/SRR6130428.sra
# Written 2 spots for seqs/SRR6130428/SRR6130428.sra
# @SRR6130428.1 1 length=202
# NGACCTTGGCGTTGGTGTAGGTCACCACTCCGATTTTGAGCTAGAATATGAAATTCATCACTGGAATAAGTTTCATAAGGACAAAGATCCAGACGTTAAAGNATGGTAAAAATATATTTTAAGTGCTTTGATTACTTACTACATGCTAATTGACTACATACATAGTGTCTTGAATACTTTCTCAGTCTCAACTATTCATCTT
# +SRR6130428.1 1 length=202
# #AAFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJFJJJJJJJJJJFFJJJJJJJJJJJJJJJJJJJJJJJ#AAA7FAFFJJJFJFF7FJJAJJ<FFJJFAJJJJJJJJFFJFJJ<JJFFJJF-FJFJJJAFJJJJFJJJJF7AJJ-77F7JJ<AFFF-F-J7-FJ-AFJFF
# @SRR6130428.2 2 length=202
# NTGCTACTCTCATGGTCTCAATACTGCCGCCGACATTCTGTCCCACATACTAAATCTCTTCCCGTCATTATCGCCCGCATCCGGTGCCGTAAATGCAAAACNGGACAATTTTTATTATATTGCATCTAATAGCAATAGGATAAGAAAGGTGAAAAAGCAAAAGCAATAGTGCATTGTGATGTGGAGAATAAGGTGCATACGA
# +SRR6130428.2 2 length=202
# #AAAAFFAFFJJJJJJJJJJJJJJJJJJFJJJFJJJJJJJFJJJJAJJJFJ<JJFA<FFAJJFFJJ<JFJJJJJJJFFJ7JJJJFJJJJJFFJAAFJFFJJ#AAFFFJAAFA-FAAFJFJJ-FF-77AJFFJJJJA<<A7-<-<<-<-77FFJJJJJJFJFAJJJ7---<77F--7-7A<FAJJA7A-7JJF--<-----7-

SQSET=SRR6130428
fastq-dump -I --split-files $WDR/seqs/${SQSET}/${SQSET}.sra \
           --gzip --outdir $WDR/seqs/${SQSET}/
#    Read 2768518 spots for seqs/SRR6130428/SRR6130428.sra
# Written 2768518 spots for seqs/SRR6130428/SRR6130428.sra
```

Here we have the commands to download any of the other three reads
datasets (\textbf{but remember that it is not necessary to fulfill
this exercise}).

```{.sh}
## Choose another dataset for an alternative assembly

# 733.6Mb downloads
mkdir -vp $WDR/seqs/SRR3481383
prefetch -v SRR3481383 -O $WDR/seqs

# 1.6Gb downloads
mkdir -vp $WDR/seqs/DRR076693
prefetch -v DRR076693  -O $WDR/seqs

# 3.2Gb downloads
mkdir -vp $WDR/seqs/SRR7548448
prefetch -v SRR7548448 -O $WDR/seqs

for SQSET in SRR3481383 DRR076693 SRR7548448;
  do {
    echo "# Working on $n" 1>&2;
    fastq-dump -I --split-files $WDR/seqs/${SQSET}/${SQSET}.sra \
               --gzip --outdir $WDR/seqs/${SQSET}/;
  }; done
# Working on SRR3481383
Read 5346334 spots for seqs/SRR3481383/SRR3481383.sra
Written 5346334 spots for seqs/SRR3481383/SRR3481383.sra
# Working on DRR076693
Read 5625017 spots for seqs/DRR076693/DRR076693.sra
Written 5625017 spots for seqs/DRR076693/DRR076693.sra
# Working on DRR076693
Read 43091394 spots for seqs/SRR7548448/SRR7548448.sra
Written 43091394 spots for seqs/SRR7548448/SRR7548448.sra
```


# The Assembly Protocol

## Exploratory data analysis of the raw reads

On this initial step, we will run
\href{https://www.bioinformatics.babraham.ac.uk/projects/fastqc/}{\texttt{fastQC}}
over `fastq` files for the forward (R1) and reverse (R2) raw reads
sets. You will need to compare those two sets from three of the
resulting plots: the quality distribution per base position
(boxplots), the base content per position (lineplots), and the read
sequences GC content distribution (lineplot). The program will
generate an `HTML` summary page, as well as a `zip` file containing
all the figures in a folder and some results in tabular format.

```{.sh}
SQSET=SRR6130428

mkdir -vp $WDR/seqs/${SQSET}/${SQSET}.QC

for READSET in 1 2;
  do {
    echo "# Running fastQC on $SQSET R${READSET}..." 1>&2;
    fastqc -t 8 --format fastq                                  \
           --contaminants $BIN/fastqc_conf/contaminant_list.txt \
           --adapters     $BIN/fastqc_conf/adapter_list.txt     \
           --limits       $BIN/fastqc_conf/limits.txt           \
           -o $WDR/seqs/${SQSET}/${SQSET}.QC                    \
           $WDR/seqs/${SQSET}/${SQSET}_${READSET}.fastq.gz      \
           2> $WDR/seqs/${SQSET}/${SQSET}.QC/fastQC_${SQSET}_${READSET}.log 1>&2;
  }; done
```

You can now open the HTML page or the zip file that `fastQC` has
produced for each reads set, and save the requested figures into the
exercise `images` folder. Rename those images in `png` format so they
fit in the table cells that make Figure \ref{fig:fastqcsetA}.

```{=latex}
% if pandoc complains or it is not running LaTeX code as expected,
% this block can also be saved into a docs/fig_kmercountdist.tex file
% and then loaded into main doc using the LaTeX "\input" macro.
\begin{figure}[!ht]
\begin{center}
 \begin{tabular}{c@{}c}
  \includegraphics[width=0.435\linewidth]{{images/fastqc_SRR6130428_quality_R1}.png}     &
  \includegraphics[width=0.435\linewidth]{{images/fastqc_SRR6130428_quality_R2}.png}     \\[-0.75ex]
  \includegraphics[width=0.435\linewidth]{{images/fastqc_SRR6130428_basecontent_R1}.png} &
  \includegraphics[width=0.435\linewidth]{{images/fastqc_SRR6130428_basecontent_R2}.png} \\[-0.75ex]
  \includegraphics[width=0.435\linewidth]{{images/fastqc_SRR6130428_gccontent_R1}.png}   &
  \includegraphics[width=0.435\linewidth]{{images/fastqc_SRR6130428_gccontent_R2}.png}   \\[-0.75ex]
 \end{tabular}
 \parbox{0.8\linewidth}{%
  \caption[Raw reads basic sequence analyses for \texttt{SRR6130428}]{%
   \label{fig:fastqcsetA}\textbf{Raw reads basic sequence analyses for \texttt{SRR6130428}.} Left column shows results for the forward reads (R1), right column for the reverse reads (R2). On top panels we can observe phred scores distribution per base position, mid panels correspond to the base composition per position, while the bottom ones show the GC content distribution across reads.
  }%caption
 }%parbox
\end{center}
\end{figure}
```


## Cleaning and trimming reads with `trimmomatic`

A crucial step before starting the assembly is to remove any
contaminant sequences, like the sequencing adapters, and low quality
segments. For this purpose, we are going to use `trimmomatic`, but
you can take another raw reads cleaner, such as `cutadapt`, to perform
this task and compare the results.

```{.sh}
##
## IMPORTANT ## 
#   If you may run out ouf memory or disk space on following analyses,
#   thus, you may consider subsampling the fasq sequences.
# 
##  You can take the first 100000 sequences with the "head" command:
# 
# gunzip -c $WDR/seqs/${SQSET}/${SQSET}_1.fastq.gz | head -n 400000 | \
#        gzip -c9 - > $WDR/seqs/${SQSET}/${SQSET}_1.subset100k.fastq.gz;
# gunzip -c $WDR/seqs/${SQSET}/${SQSET}_2.fastq.gz | head -n 400000 | \
#        gzip -c9 - > $WDR/seqs/${SQSET}/${SQSET}_2.subset100k.fastq.gz;
# 
##  You can also use "seqtk" program (if installed), as follows:
#
# seqtk sample -s11 $WDR/seqs/${SQSET}/${SQSET}_1.fastq.gz 0.1 | \
#         gzip -c9 - > $WDR/seqs/${SQSET}/${SQSET}_1.subset10pct.fastq.gz;
# seqtk sample -s11 $WDR/seqs/${SQSET}/${SQSET}_2.fastq.gz 0.1 | \
#         gzip -c9 - > $WDR/seqs/${SQSET}/${SQSET}_2.subset10pct.fastq.gz;
#
##  The -s option specifies the seed value for the random number generator,
#   it needs to be the same value for file1 and file2 to keep the paired reads
#   on the sampled output. 0.1 argument sets the percent of seqs to sample (10%).
#
##  Just fix the following commands to work with the proper subset files instead.

# to find where trimmomatic was installed in your system
# find /usr -name trimmomatic.jar
ls -l $(brew --prefix trimmomatic)/bin/trimmomatic
# /usr/local/opt/trimmomatic/bin/trimmomatic

# Take the result of the previous command as the folder to set on TMC:
export TMC=/usr/local/opt/trimmomatic/bin/trimmomatic
# Set this variable to suit your system installation:
# export TMA=/usr/share/trimmomatic/
#
# You can find trimmomatic parameters documentation at:
# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
export TRMPAR="LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30 TOPHRED33";
export TRMPECLP="ILLUMINACLIP:$TMC/TruSeq2-PE.fa:2:30:10";

$TMC  PE                          \
     $WDR/seqs/${SQSET}/${SQSET}_1.fastq.gz \
     $WDR/seqs/${SQSET}/${SQSET}_2.fastq.gz \
     $WDR/seqs/${SQSET}/${SQSET}_1.trimmo_pe.fastq.gz \
     $WDR/seqs/${SQSET}/${SQSET}_1.trimmo_sg.fastq.gz \
     $WDR/seqs/${SQSET}/${SQSET}_2.trimmo_pe.fastq.gz \
     $WDR/seqs/${SQSET}/${SQSET}_2.trimmo_sg.fastq.gz \
     $TRMPECLP $TRMPAR                      \
  2> $WDR/seqs/${SQSET}/${SQSET}.trimmo.log 1>&2 ;

tail -n 2 $WDR/seqs/${SQSET}/${SQSET}.trimmo.log
#  Input Read Pairs:         2768518
#            Both Surviving: 2536084 (91.60%)
#    Forward Only Surviving:  176365 ( 6.37%)
#    Reverse Only Surviving:    4082 ( 0.15%)
#                   Dropped:   51987 ( 1.88%)
# TrimmomaticPE: Completed successfully
# My results below:
# Input Read Pairs: 2768518 
# Both Surviving: 2587195 (93.45%) 
# Forward Only Surviving: 175015 (6.32%) 
# Reverse Only Surviving: 4176 (0.15%) 
# Dropped: 2132 (0.08%)
# TrimmomaticPE: Completed successfully

```

Next we should check the new reads filtered, i.e. by running `fastqc`
over the pair-end reads from `trimmomatic`.

```{.sh}
SQSET=SRR6130428

mkdir -vp $WDR/seqs/${SQSET}/${SQSET}.trimmo.QC

for READSET in 1 2;
  do {
    echo "# Running fastQC on trimmed PE reads $SQSET R${READSET}..." 1>&2;
    fastqc -t 8 --format fastq                                  \
           --contaminants $BIN/fastqc_conf/contaminant_list.txt \
           --adapters     $BIN/fastqc_conf/adapter_list.txt     \
           --limits       $BIN/fastqc_conf/limits.txt           \
           -o $WDR/seqs/${SQSET}/${SQSET}.trimmo.QC             \
	          $WDR/seqs/${SQSET}/${SQSET}_${READSET}.trimmo_pe.fastq.gz \
           2> $WDR/seqs/${SQSET}/${SQSET}.trimmo.QC/fastQC_${SQSET}_${READSET}.log 1>&2;
  }; done
```

```{=latex}
% if pandoc complains or it is not running LaTeX code as expected,
% this block can also be saved into a docs/fig_kmercountdist.tex file
% and then loaded into main doc using the LaTeX "\input" macro.
\begin{figure}[!ht]
\begin{center}
 \begin{tabular}{c@{}c}
  \includegraphics[width=0.435\linewidth]{{images/trimmo_SRR6130428_quality_R1}.png}     &
  \includegraphics[width=0.435\linewidth]{{images/trimmo_SRR6130428_quality_R2}.png}     \\[-0.75ex]
  \includegraphics[width=0.435\linewidth]{{images/trimmo_SRR6130428_basecontent_R1}.png} &
  \includegraphics[width=0.435\linewidth]{{images/trimmo_SRR6130428_basecontent_R2}.png} \\[-0.75ex]
  \includegraphics[width=0.435\linewidth]{{images/trimmo_SRR6130428_gccontent_R1}.png}   &
  \includegraphics[width=0.435\linewidth]{{images/trimmo_SRR6130428_gccontent_R2}.png}   \\[-0.75ex]
 \end{tabular}
 \parbox{0.8\linewidth}{%
  \caption[Cleaned reads basic sequence analyses for \texttt{SRR6130428}]{%
   \label{fig:fastqcsetA}\textbf{Cleaned reads basic sequence analyses for \texttt{SRR6130428}.} Left column shows results for the forward reads (R1), right column for the reverse reads (R2). On top panels we can observe phred scores distribution per base position, mid panels correspond to the base composition per position, while the bottom ones show the GC content distribution across reads.
  }%caption
 }%parbox
\end{center}
\end{figure}
```

We can now estimate the sequencing coverage, as we already have the
size of the reference genome, and the total amount of nucleotides
generated by the sequencing projects.

```{.sh}
# raw PE-reads
gunzip -c $WDR/seqs/${SQSET}/${SQSET}_[12].fastq.gz | \
  infoseq -sequence fastq::stdin -only -length -noheading | \
    gawk '{ s+=$1; n++ }
          END{ printf "# Total %d sequences and %d nucleotides\n", n, s }'
# Total 5537036 sequences and 559240636 nucleotides (my results)

# cleaned PE-reads
gunzip -c $WDR/seqs/${SQSET}/${SQSET}_[12].trimmo_pe.fastq.gz | \
  infoseq -sequence fastq::stdin -only -length -noheading | \
    gawk '{ s+=$1; n++ }
          END{ printf "# Total %d sequences and %d nucleotides\n", n, s }'
# Total 5174390 sequences and 507747193 nucleotides (my results)
```

For raw and cleaned reads we can estimate a sequencing coverage of
$\num{559240636}/\num{12157105}=46.00$X and $\num{559240636}/\num{12157105}=41.09$X
respectively.


## Assembling reads with `SOAPdenovo`

`SOAPdenovo`\footnote{Luo et al. "SOAPdenovo2: an empirically improved
memory-efficient short-read de novo assembler". \textit{GigaScience},
1:18, 2012.} is easy to install and to configure. It has been reported
to perform like other more complex tools, generating in some tests
less chimeric contigs and missasemblies.

```{.r}
# SQSET=SRR6130428
export SPD=$WDR/soapdenovo/$SQSET
mkdir -vp $SPD
  
export GENOMESIZE=12157105;   # estimated or real... ;^D

# Generating the configuration file:
cat > $SPD/${SQSET}_soap.conf <<EOF
max_rd_len=100
[LIB]
avg_ins=450
reverse_seq=0
asm_flags=3
rd_len_cutoff=100
rank=1
pair_num_cutoff=3
map_len=32
q1=$WDR/seqs/${SQSET}/${SQSET}_1.trimmo_pe.fastq.gz
q2=$WDR/seqs/${SQSET}/${SQSET}_2.trimmo_pe.fastq.gz
EOF

cat > $SPD/${SQSET}_soap.conf <<EOF
# Raw reads from SRA experiment ${SQSET}
# Description: Original WT strain, Illumina HiSeq 4000 (2x100bp)
# Insert size was not described on the SRA project, using 450bp
# ---
# maximal read length
max_rd_len=100
#
[LIB]
# average insert size
avg_ins=450
# if sequence needs to be reversed (0 for short insert PE libs)
reverse_seq=0
# in which part(s) the reads are used (3 for both contig and scaffold assembly)
asm_flags=3
# use only first 100 bps of each read
rd_len_cutoff=100
# in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
# minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
# a pair of fastq file, read 1 file should always be followed by read 2 file
# we are going to use pair-end reads after filtering with trimmomatic
q1=$SPD/${SQSET}_1.trimmo_pe.fastq.gz
q2=$SPD/${SQSET}_2.trimmo_pe.fastq.gz
EOF

# Make sure you are in the MScGG-BIA2324_exercises environment, because that is where
# you installed soapdenovo
# Building the k-mers graph
SOAPdenovo-63mer pregraph -s $SPD/${SQSET}_soap.conf -K 63 -R -p 8 \
                           -o $SPD/${SQSET}_k63_graph               \
                           2> $SPD/${SQSET}_k63_pregraph.log 1>&2;
			    
# Contiging stage
SOAPdenovo-63mer contig   -g $SPD/${SQSET}_k63_graph -R -p 8 \
                           2> $SPD/${SQSET}_k63_contig.log 1>&2;

# Mapping reads back over the contigs
SOAPdenovo-63mer map      -s $SPD/${SQSET}_soap.conf -p 8 \
                           -g $SPD/${SQSET}_k63_graph      \
                           2> $SPD/${SQSET}_k63_map.log 1>&2;

# Scaffolding contigs if we have long range reads (i.e. a 2kbp mate-pairs run)
# in our case, as we only have a single pair-ends library, we will get a result
# that will be similar or the same as what we have obtained in the contiging stage
SOAPdenovo-63mer scaff    -F -p 8 -N $GENOMESIZE            \
                           -g $SPD/${SQSET}_k63_graph_prefix \
                           2> $SPD/${SQSET}_k63_scaff.log 1>&2;

# just by looking at the $SPD/${SQSET}_k63_contig.log file,
# we can retrieve info about contigs assembly, such as N50
tail $SPD/${SQSET}_k63_contig.log
# There are 3156 contig(s) longer than 100, sum up 11779878 bp, with average length 3732.
# The longest length is 69946 bp, contig N50 is 14294 bp, contig N90 is 3279 bp.
# 4867 contig(s) longer than 64 output.
```

```{.sh}
# We need to compute some extra stats from the assemblies,
# we can download from github one of the assemblathon scripts for this purpose
# (to facilitate the task, it is already available on the bin folder)
#
export PERL5LIB=$BIN;
#
gzcat $WDR/seqs/$REFDIR/$REFGEN.fsa.gz | \
  $BIN/assemblathon_stats.pl \
            -csv -genome_size 12160000 - \
          > $WDR/stats/assembly_stats_$REFGEN.txt

gzcat $WDR/seqs/$REFDIR/$REFGEN.fsa.gz > $WDR/seqs/$REFDIR/$REFGEN.fsa
$BIN/assemblathon_stats.pl \ 
            -csv -genome_size 12160000 $WDR/seqs/$REFDIR/$REFGEN.fsa \
            > $WDR/stats/assembly_stats_$REFGEN.csv

$BIN/assemblathon_stats.pl \
            -csv -genome_size 12160000 \
            $WDR/soapdenovo/SRR6130428/SRR6130428_k63_graph.contig.fa \
          > $WDR/stats/assembly_stats_SRR6130428_soapdenovo_k63_graph.contig.txt
```

You can check files `*.contig` and `*.scafSeq` in case they were
created, for the fasta sequences assembled contigs and scaffolds
respectively.

__Open questions arise:__

* Can you calculate the contigs lengths and plot that distribution?
* What do you think it will happen when using a larger coverage reads set?
* Do you think knowing the real insert size will improve the assembly?

Some can be answered by comparing with the assembly produced over
another of the suggested raw-reads sets.


## Estimating insert size with `picard`

Let's check whether estimated insert size was an educated guess. We
must first align the PE reads against the reference genome or the
assembly, then we can take the alignments in `bam` format to estimate
insert size with `picard` tool.

```{.r}
SQSET=SRR6130428;
export BWT=$WDR/bowtie;

mkdir -vp $BWT/$SQSET;

# by linking contigs file to one file with a fasta suffix,
# we can facilitate using those sequences on IGV later on...
ln -vs ./${SQSET}_k63_graph.contig \
       $WDR/soapdenovo/${SQSET}/${SQSET}_k63_graph.contig.fa;

# preparing reference sequence databases for bowtie
FSAD=$WDR/seqs/$REFDIR;
bowtie2-build --large-index -o 2 \
              $FSAD/$REFGEN.fsa.gz \
              $BWT/scer_refgenome.bowtiedb \
           2> $BWT/scer_refgenome.bowtiedb.log 1>&2;
	
bowtie2-build --large-index -o 2 \
              $WDR/soapdenovo/${SQSET}/${SQSET}_k63_graph.contig \
              $BWT/${SQSET}/${SQSET}_k63_graph.contig.bowtiedb \
           2> $BWT/${SQSET}/${SQSET}_k63_graph.contig.bowtiedb.log 1>&2;

# mapping pe reads over reference sequences
TMP=$WDR/tmpsort;
mkdir $TMP
PEfileR1=$WDR/seqs/${SQSET}/${SQSET}_1.trimmo_pe.fastq.gz;
PEfileR2=$WDR/seqs/${SQSET}/${SQSET}_2.trimmo_pe.fastq.gz;

BWTBF=$BWT/${SQSET}/${SQSET}-x-scer_refgen.bowtie;
TMPBF=$TMP/${SQSET}-x-scer_refgen.bowtie;

bowtie2 -q --threads 8 -k 5 -L 12                    \
        --local --sensitive-local --no-unal --met 60 \
        --met-file $BWTBF.metrics                    \
        -x         $BWT/scer_refgenome.bowtiedb      \
        -1         $PEfileR1  \
        -2         $PEfileR2  \
        -S         $TMPBF.sam \
        2>         $BWTBF.log 1>&2;
        
( samtools view -Sb -o $TMPBF.bam \
                       $TMPBF.sam;
  samtools sort $TMPBF.bam   \
             -o $TMPBF.sorted.bam;
  cp -v $TMPBF.sorted.bam $BWTBF.sorted.bam
  rm -v $TMPBF.sam $TMPBF.bam
  ) 2> $BWTBF.bowtie2sortedbam.log 1>&2;

samtools index $BWTBF.sorted.bam;

java -jar $BIN/picard.jar CollectInsertSizeMetrics   \
          HISTOGRAM_FILE=$BWTBF.insertsize.hist.pdf  \
                   INPUT=$BWTBF.sorted.bam           \
                  OUTPUT=$BWTBF.insertsize_stats.txt \
           ASSUME_SORTED=true \
              DEVIATIONS=25   \
                      2> $BWTBF.insertsize_stats.log;


# now let's check the assembled contigs

BWTBF=$BWT/${SQSET}/${SQSET}-x-soapk63ctgs.bowtie;
TMPBF=$TMP/${SQSET}-x-soapk63ctgs.bowtie;

bowtie2 -q --threads 8 -k 5 -L 12                    \
        --local --sensitive-local --no-unal --met 60 \
        --met-file $BWTBF.metrics                    \
        -x         $BWT/scer_refgenome.bowtiedb      \
        -1         $PEfileR1  \
        -2         $PEfileR2  \
        -S         $TMPBF.sam \
        2>         $BWTBF.log 1>&2;
	
( samtools view -Sb -o $TMPBF.bam \
                       $TMPBF.sam;
  samtools sort $TMPBF.bam   \
             -o $TMPBF.sorted.bam;
  mv -v $TMPBF.sorted.bam $BWTBF.sorted.bam;
  rm -v $TMPBF.sam $TMPBF.bam
  ) 2> $BWTBF.bowtie2sortedbam.log 1>&2;

samtools index $BWTBF.sorted.bam;

java -jar $BIN/picard.jar CollectInsertSizeMetrics   \
          HISTOGRAM_FILE=$BWTBF.insertsize.hist.pdf  \
                   INPUT=$BWTBF.sorted.bam           \
                  OUTPUT=$BWTBF.insertsize_stats.txt \
           ASSUME_SORTED=true \
              DEVIATIONS=25   \
                      2> $BWTBF.insertsize_stats.log;
```

You can include here the two PDF histograms generated by `picard` from
the reads alignment over the reference and `soapdenovo`
assemblies. You can open the genome set and the aligned reads using a
browser like `igv` or `tablet`.

```{=latex}
\begin{figure}[!ht]
\begin{center}
 \begin{tabular}{c@{}c}
  \includegraphics[width=0.435\linewidth]{{images/hist_picard_ref}.png}     &
  \includegraphics[width=0.435\linewidth]{{images/hist_picard_soap}.png}     \\
 \end{tabular}
 \parbox{0.8\linewidth}{%
  \caption[Histogram of read alignments]{%
   \label{fig:fastqcsetA}\textbf{Histogram of read alignments} Left column shows reads aligned over the reference genome and right column shows reads aligned over soapdenovo assembly. The histograms are generated by picard.
  }%caption
 }%parbox
\end{center}
\end{figure}
```

# Exploring the assemblies

## Filter out contigs mapping to reference chromosomes

On this section we are going to run `dnadiff` from the `MUMmer`
package\footnote{S. Kurtz, A. Phillippy, A.L. Delcher, M. Smoot,
M. Shumway, C. Antonescu, and S.L. Salzberg.\newline\hspace*{1cm}
"Versatile and open software for comparing large genomes."
\textit{Genome Biology}, 5:R12, 2004.} to compare assembled contigs
against the reference or between them. In order to speed up the
example, we will first use `NCBI-BLAST`\footnote{C. Camacho,
G. Coulouris, V. Avagyan, N. Ma, J. Papadopoulos, K. Bealer, and
T.L. Madden.\newline\hspace*{1cm} "BLAST+: architecture and
applications." \textit{BMC Bioinformatics}, 10:421, 2008.} to project
all the assembled contigs into the chosen reference chromosomes, so we
will reduce all the downstream calculations. We need to create a
database for each assembly sequence sets and we will query by the
reference selected chromosomes.

```{.sh}
SQSET=SRR6130428

mkdir -vp $WDR/blast/dbs

makeblastdb -in $WDR/soapdenovo/$SQSET/${SQSET}_k63_graph.contig.fa \
            -dbtype nucl \
	    -title "${SQSET}_SOAPdenovo_k63_contigs" \
	    -out $WDR/blast/dbs/${SQSET}_SOAPdenovo_k63_contigs \
	      2> $WDR/blast/dbs/${SQSET}_SOAPdenovo_k63_contigs.log 1>&2;
```

We will focus on a couple of reference genome chromosomes, `chrI` and
`chrM` (mitochondrion genome). Thus, we need to download them from the
genome repository or to filter them out from the whole genome fasta
file as shown below:

```{.sh}
mkdir -vp $WDR/seqs/chrs;

for SQ in chrI:NC_001133 chrM:NC_001224;
  do {
    SQN=${SQ%%:*}; # get the chr from SQ string
    SQI=${SQ##*:}; # get the refseq id from SQ string
    echo "# Filtering $SQN [$SQI] from whole genome fasta file..." 1>&2;
    samtools faidx \
             $WDR/seqs/$REFDIR/$REFGEN.fsa \
	     'ref|'$SQI'|' | \
      sed 's/^>.*$/>Scer_'$SQN'/;' > $WDR/seqs/chrs/$SQN.fa
  }; done
```

Now we can `BLAST` reference sequences against the newly created
database; we will use `megablast` option as we are comparing sequences
for the same species, and that `BLAST` program has the parameters
optimized for this kind of genomic searches.

```{.sh}
# we define here a custom BLAST tabular output format and export it to become an
# environment variable
BLASTOUTFORMAT='6 qseqid qlen sseqid slen qstart qend sstart send length';
BLASTOUTFORMAT=$BLASTOUTFORMAT' score evalue bitscore pident nident ppos positive';
BLASTOUTFORMAT=$BLASTOUTFORMAT' mismatch gapopen gaps qframe sframe';
export BLASTOUTFORMAT;

for SQ in chrI chrM;
  do {
    echo "# Running MEGABLAST: $SQSET x $SQ ..." 1>&2;
    mkdir -vp $WDR/blast/${SQ}-x-${SQSET};
    blastn -task megablast -num_threads 8                        \
           -db    $WDR/blast/dbs/${SQSET}_SOAPdenovo_k63_contigs \
           -outfmt "$BLASTOUTFORMAT"                             \
           -query $WDR/seqs/chrs/$SQ.fa                          \
           -out   $WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast.out \
             2>   $WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast.log;
  }; done
```

Once we have found the matching contigs, let's filter them out from
the whole assembly fasta file.

```{.sh}
for SQ in chrI chrM;
  do {
    echo "# Running MEGABLAST: $SQSET x $SQ ..." 1>&2;
    OFBN="$WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast";
    # this is a check to ensure we start with an empty contigs fasta file
    if [ -e "$OFBN.fa" ];
      then
        printf '' > "$OFBN.fa"; # rm can be also used here, but dangerous for novice
      fi;
    # get the contig IDs from third column and filter out sequences
    gawk '{ print $3 }' "$OFBN.out" | sed 's/^>//' | sort | uniq | \
      while read SQID;
        do {
           samtools faidx \
                    $WDR/soapdenovo/$SQSET/${SQSET}_k63_graph.contig.fa \
                    "$SQID";
        }; done >> "$OFBN.fa" \
                2> "$OFBN.fa.log";
  }; done
```

Now that we have the filtered out the required contigs, we can start the following sequence comparison procedure based on `dnadiff`:


```{.sh}
# Part 1: DNADIFF Execution and Alignment Plotting
for SQ in chrI chrM; do
    printf "# Running DNADIFF protocol: $SQSET x $SQ ..." 1>&2;
    IFBN="$WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast";
    OFBD="$WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast.dnadiff";
    dnadiff -p $OFBD \
            $WDR/seqs/chrs/$SQ.fa \
            $IFBN.fa \
            2> $OFBD.log;
    printf " DNAdiff..." 1>&2;
    mummerplot --large --layout --fat \
            -t "Alignment Plot: ${SQ}-x-${SQSET}" \
            -p $OFBD \
            --png \
            $OFBD.1delta \
            2> $OFBD.alnplot.log;
    magick convert -verbose $OFBD.ps $OFBD.png;
    printf " ALNplot..." 1>&2;
    printf " DONE\n" 1>&2;
done

# Part 2: Coverage Plotting and Summary
for SQ in chrI chrM; do
    printf "# Running DNADIFF protocol: $SQSET x $SQ ..." 1>&2;
    IFBN="$WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast";
    OFBD="$WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast.dnadiff";
    mummerplot --large --layout --fat --coverage --postscript \
               -t "Coverage Plot: ${SQ}-x-${SQSET}" \
               -p $OFBD.covg \
               $OFBD.1delta \
               2> $OFBD.cvgplot.log;
    magick convert -verbose $OFBD.covg.ps $OFBD.covg.png;
    printf " CVGplot..." 1>&2;
    ( grep "TotalBases"   $OFBD.report;
      grep "AlignedBases" $OFBD.report;
      grep "AvgIdentity"  $OFBD.report
      ) > $OFBD.shortsummary;
    printf " DONE\n" 1>&2;
done

```

```{=latex}
% if pandoc complains or it is not running LaTeX code as expected,
% this block can also be saved into a docs/fig_kmercountdist.tex file
% and then loaded into main doc using the LaTeX "\input" macro.
\begin{figure}[!ht]
\begin{center}
 \begin{tabular}{cc}
  \bf chrI & \bf chrM \\
  \includegraphics[width=0.25\linewidth, angle=-90]{{blast/chrI-x-SRR6130428/chrI-x-SRR6130428_megablast.dnadiff}.png} &
  \includegraphics[width=0.25\linewidth, angle=-90]{{blast/chrM-x-SRR6130428/chrM-x-SRR6130428_megablast.dnadiff}.png} \\
  %% trim={<left> <lower> <right> <upper>}
  \includegraphics[width=0.25\linewidth, trim={0 0 0 21cm}, clip, angle=-90]%
                  {{blast/chrI-x-SRR6130428/chrI-x-SRR6130428_megablast.dnadiff.covg}.png} &
  \includegraphics[width=0.25\linewidth, trim={0 0 0 21cm}, clip, angle=-90]%
                  {{blast/chrM-x-SRR6130428/chrM-x-SRR6130428_megablast.dnadiff.covg}.png} \\
 \end{tabular}
 \parbox{0.75\linewidth}{%
  \caption[\texttt{dnadiff} comparison between two reference chromosomes and \texttt{SRR6130428} contigs]{%
   \label{fig:fastqcsetA}\textbf{\texttt{dnadiff} comparison between two reference chromosomes and \texttt{SRR6130428} contigs.} Top panels show the alignment plots of contigs from \texttt{SRR6130428} assembly mapped over two reference \Scer\ chromosomes, chrI and chrM on left and right panels respectively. Bottom panels show the alignment coverage for the same sequence relations. It is evident from the comparison that contigs aligning to chrM have better contiguity despite overall coverages are quite similar on both reference chromosomes.
  }%caption
 }%parbox
\end{center}
\end{figure}

```

Discuss later on which of the two chromosome assemblies do you think
had a better outcome and try to guess why. We may need to analyze
whether the contigs fully align to the reference chromosomes, __can you
calculate the average coverage of the aligned reads?__


## Assessment of genome completeness with `BUSCO`

`BUSCO`\footnote{M. Manni, M.R. Berkeley, M. Seppey, F.A. Simão, and
E.M. Zdobnov.\newline\hspace*{1cm} "\texttt{BUSCO} Update: Novel and
Streamlined Workflows along with Broader and Deeper Phylogenetic
Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes."
\textit{Molecular Biology and Evolution}, 38(10)4647-–4654, 2021.}
estimates the completeness and redundancy of processed genomic data
based on universal single-copy orthologs. First of all, we need to
check if there is a clade-specific parameters set that fits with our
organism, \Scer\ belong to the _Saccharomycetes_ class, within the
_Ascomycota_ phylum in the Fungi kingdom (see
\href{https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=4932}{NCBI
Taxonomy browser species card}).

```{.sh}
run_BUSCO.py --list-datasets
# 2022-10-24 19:06:49 INFO:	Downloading information on latest versions of BUSCO data...
# 2022-10-24 19:06:52 INFO:	Downloading file 'https://busco-data.ezlab.org/v5/data/information/lineages_list.2021-12-14.txt.tar.gz'
# 2022-10-24 19:06:53 INFO:	Decompressing file '/home/lopep/SANDBOX/BScCG2324/exercise_02/busco_downloads/information/lineages_list.2021-12-14.txt.tar.gz'
# 
# ################################################
# 
# Datasets available to be used with BUSCO v4 and v5:
# 
#  bacteria_odb10
#      - ...
#  archaea_odb10
#      - ...
#  eukaryota_odb10
#      - ...
#      - fungi_odb10
#          - ascomycota_odb10
#              - ...
#              - saccharomycetes_odb10
#              - ...
#          - ...
#      - ...
#  viruses (no root dataset)
#      - ...
# 
```

Luckily for us, there is one class specific set to evaluate the
completeness of our genome assembly `saccharomycetes_odb10`. However,
it is worth that you consider also to run the `BUSCO` commands below
using a more general parameters set, such as `fungi_odb10` or even
`eukaryota_odb10` (or both), and compare the corresponding
results. This can be useful to extrapolate the assessment to other
species genomes for which we will not have as much information as for
this yeast.

```{.sh}
mkdir -vp $WDR/busco/$SQSET

# fix PATH to point bin folders where you installed bbmap and busco programs
export PATH=$BIN:$BIN/bbmap:$BIN/busco/bin:$PATH

# cd $WDR
busco -m genome  --download_base_url  https://busco-data.ezlab.org/v5  \
         -i $WDR/soapdenovo/$SQSET/${SQSET}_k63_graph.contig.fa \
         -o ./busco/$SQSET/${SQSET}_k63_contigs \
         -l saccharomycetes_odb10 -f
## if you get an error "A run with the name xxx already exists"
## you should add command-line option "-f" to force overwriting those files
#    ---------------------------------------------------
#    |Results from dataset saccharomycetes_odb10       |
#    ---------------------------------------------------
#    |  C:98.8%[S:96.7%,D:2.1%],F:0.5%,M:0.7%,n:2137   |
#    |2111  Complete BUSCOs (C)                        |
#    |2067  Complete and single-copy BUSCOs (S)        |
#    |44    Complete and duplicated BUSCOs (D)         |
#    |11    Fragmented BUSCOs (F)                      |
#    |15    Missing BUSCOs (M)                         |
#    |2137  Total BUSCO groups searched                |
#    ---------------------------------------------------
# 2022-10-26 19:55:26 INFO: BUSCO analysis done. Total running time: 440 seconds
```

You can take a step further to plot the resulting completeness values,
using the `generate_plot.py` script that is provided under the scripts
folder of your `busco` installation:

```{.sh}
cd ./busco/$SQSET/
python3 $BIN/busco/scripts/generate_plot.py \
```

**Embed here the corresponding figure**; if you have ran more than one
assembly, you can combine the `busco` results to compare completeness
among them.

```{=latex}
\begin{figure}[htbp]
  \centering
  \includegraphics[width=0.8\textwidth]{busco/SRR6130428/busco_figure.png} 
  \caption{BUSCO Assessment Results}
  \label{fig:example}
\end{figure}
```


# Discussion

__IMPORTANT__ Discuss your results here (around 300 words). And
remember to include in the Appendices section (see page
\pageref{sec:appendices}), any extra script you wrote from this
exercise `bin` folder using the `loadfile` macro. We can take
advantage of the \LaTeX\ referencing capabilities, as described in the
first exercise template.


\clearpage

# Appendices
\label{sec:appendices}

## Software

We have used the following versions:

```{.sh}
uname -a
# Linux aleph 5.15.0-48-generic #54-Ubuntu SMP
# Fri Aug 26 13:26:29 UTC 2022 x86_64 x86_64 x86_64 GNU/Linux

perl -v
# This is perl 5, version 34, subversion 0 (v5.34.0)
# built for x86_64-linux-gnu-thread-multi
# (with 57 registered patches, see perl -V for more detail)

python3 --version
# Python 3.10.6

R --version
# R version 4.1.2 (2021-11-01) -- "Bird Hippie"
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)

wget --version
# GNU Wget 1.21.2 built on linux-gnu.

infoseq -version
# EMBOSS:6.6.0.0

gunzip --version
# gunzip (gzip) 1.10

apt-cache policy fastqc trimmomatic \
                 samtools bamtools picard-tools \
                 bwa bowtie2 soapdenovo2 | \
  gawk '$0 !~ /^ / { printf "#%15s ", $1 }
        $1 == "Installed:" { print $2 }'
#        fastqc: 0.11.9+dfsg-5
#   trimmomatic: 0.39+dfsg-2
#      samtools: 1.13-4
#      bamtools: 2.5.1+dfsg-10build1
#  picard-tools: 2.26.10+dfsg-1
#           bwa: 0.7.17-6
#       bowtie2: 2.4.4-1
#   soapdenovo2: 242+dfsg-2

busco -v
# BUSCO 5.4.3

pandoc --version
# pandoc 2.9.2.1
# Compiled with pandoc-types 1.20, texmath 0.12.0.2, skylighting 0.8.5
```


## Supplementary files
\label{sec:supplfiles}


### `conda` environment dependencies for the exercise

\loadfile{environment.yml}{environment.yml}{prg:environmentYML}


### Project specific scripts

```{=latex}
% now it is commented on the LaTeX document
\loadfile{assemblathon\_stats.pl}{bin/assemblathon_stats.pl}{prg:assemblathonstatsPERL}
```


### Shell global vars and settings for this project

\loadfile{projectvars.sh}{projectvars.sh}{prg:projectvarsBASH}


## About this document

This document was be compiled into a PDF using `pandoc` (see
`projectvars.sh` from previous subsection) and some `LaTeX` packages
installed in this linux system. `synaptic`, `apt-get` or `aptitude`
can be used to retrieve and install those tools from linux
repositories. As the `raw_tex` extension has been provided to the
`markdown_github` and `tex_math_dollars` formats, now this document
supports inline \LaTeX\ and inline formulas!

You can get further information from the following links about the
[Mark Down syntax](http://daringfireball.net/projects/markdown/syntax#link), as
well as from the manual pages (just type `man pandoc` and/or `man
pandoc_markdown`).
