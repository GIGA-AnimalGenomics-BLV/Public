# Proviral Integration Site Calling Pipeline 

## Prerequsits 

* bowtie >= 1.1.2
* bowtie2 >= 2.2.9
* cutadapt >= 1.7.1
* fastx 0.0.13
* samtools >= 0.1.19
* bbmap
* R >= 3.2.2
  * PIC (2.0)
  * dplyr (>= 0.7.6)
  * ggplot2 (>= 2.2.1)
  * tibble (>= 1.4.2)
  * readr (>= 1.1.1)
  * WriteXLS (>= 4.0.0)
  * ShortRead (>= 1.3)
  * stringr (>= 1.3)
  * tidyr (>= 0.8)
  * GenomicRanges (>= 1.32.2)

Versions specified here have been tested succesfully. 

Although the pipeline is relatively light in memory, running it with less than 12 Gb of RAM is not recommended. 

## Introduction

This pipeline is dedicated to the analysis of integration site "clonality" data as described in Rosewick *et al*, 2017. However it can easely repurpose for any other ligation-mediated clonality methods.

Described pipeline has been now cited in:

* Rosewick *et al*, **Nature Communications**, 2017
* Artesi *et al*, **Leukemia**, 2017
* Hahaut *et al*, **Biorxiv**, 2018

A test dataset can be found in: [test data online](https://)

This example is intended to work on HTLV proviruses. Some parameters (index, LTR sequences, ...) have to be adapated for BLV or HIV.


## Pipeline

### Indexes

Before mapping an index need to be created. Three FASTA files have to be provided:

* Complete Viral Sequence
* Host Genome
* LTR sequences 
 * LTRs sequences should comprised the "chromosomes": >LTR3 and >LTR5
 * Provide only the sequences starting from the primer 'start' until the LTR's end.
 
To annotate final results a GTF file downloaded from ENSEMBL is also required [ensembl FTP](https://www.ensembl.org/info/data/ftp/index.html).


#### 1. Bowtie2 viral-host genome

```
viral_genome="path/to/viral/genome.fasta"
host_genome="path/to/host/genome.fasta"

cat $viral_genome $host_genome > viral_host.fa

bowtie2-build viral_host.fa viral_host 
```

#### 2. Bowtie2 LTR sequences

```
LTR_sequences="path/to/LTR.fa"

bowtie2-build $LTR_sequences LTR_sequence 
```

#### 3. Annotation GTF

```
GTF="path/to/ensembl_host.gtf.gz"

zcat $GTF | grep ensembl[[:space:]]gene | awk -F '\t' '{print $9}' | cut -d ';' -f1 | cut -d ' ' -f2 | sed 's/"//g' > geneInfos/geneID.txt
zcat $GTF | grep ensembl[[:space:]]gene | awk -F '\t' '{print $9}' | cut -d ';' -f 3 | awk -F " " '{if($1 == "gene_name"){print $2} else print "NA"}' | sed 's/"//g' > geneInfos/geneName.txt
zcat $GTF | grep ensembl[[:space:]]gene | awk -F '\t' '{print $1"\t"$4"\t"$5"\t"$7}' > geneInfos/genePositions.txt

paste -d '\t' geneInfos/genePositions.txt geneInfos/geneID.txt geneInfos/geneName.txt | sort > geneInfos/geneInfos.sorted.txt
```

### Variable definition

Several information are required to start working on raw data: 

```
# HTLV:

LTR5len=36 # Length of the LTR5 sequence, starting from the LTR to the primer end
LTR3len=45 # Length of the LTR3 sequence, starting from the primer to the LTR end
virus="HTLV_ATK" # Name given to the indexed viral sequence
linkLen=22 # Length of the linker sequence

sampleName="mySample"
r1="path/to/raw/file/r1.fastq"
r2="path/to/raw/file/r2.fastq"
linker="ATGGTGCCAATGGC" # ==> EXAMPLE, should be replaced by the one of the sample

bowtie2_virusHost_index="path/to/bowtie2/virusHost/index/prefix"
bowtie2_LTR_index="path/to/bowtie2/LTR/index/prefix"
geneBedFile="path/to/geneInfos.sorted.txt"

LTR3sequence_firstbp="TTAGTACACA"
LTR5sequence_firstbp="TGACAATGAC"
```

### 1. Select reads harboring linker sequence

In order to speed-up detection of linker sequences, short-read aligner should be prefered over BLAST/BLAT. We create on-the-fly a new linker index for each sample.

```
mkdir linkerBowtie
cd linkerBowtie

echo >linker_$linker >> linker.fa
echo $linker >> linker.fa

bowtie-build -q linker.fa linker

cd ..
```

The linkers is then extracted from each R2 reads based on its hypothetical position (here from position 9 to 9+linkerLen). Sequencing quality is taken into account (-Q33).

```
fastx_trimmer -Q33 -f 9 -l $linkLen -i $r2 -o read_linkers.fastq
```

Sequences are mapped onto the index using bowtie. Read IDs harboring the linker are selected.

```
bowtie -l 7 -p 2 -v 0 linkerBowtie/linker read_linkers.fastq | egrep `echo $linker | sed 's/\;/\|/g'` | awk '{print $1}' > readID_linkerFound.txt
```

Finally only the reads with the right linker are selected

```
bbmap/filterbyname.sh in=$r1 in2=$r2 out=R1_linkerFound.fastq out2=R2_linkerFound.fastq names=readID_linkerFound.txt include=t overwrite=t
```

### 2. Select reads harboring LTR sequence

LTR sequence is located in R1, after complexity sequence (8bp). Length should be provided as argument (LTR3len or LTR5len)

```
fastx_trimmer -Q33 -f 8 -l $LTR3len -i R1_linkerFound.fastq -o R1_virusLTR3_edge.fastq
fastx_trimmer -Q33 -f 8 -l $LTR5len -i R1_linkerFound.fastq -o R1_virusLTR5_edge.fastq
```

Align onto the proviral sequence (LTR edge only) and select IDs from reads mapping onto it with the following arguments:
```
bowtie2 -L 14 -p 2 -N 1 -x $bowtie2_LTR_index -U R1_virusLTR3_edge.fastq | grep "LTR3" | awk '{print $1}' | grep -v "^@"  > readID_LTR3Found.txt
bowtie2 -L 14 -p 2 -N 1 -x $bowtie2_LTR_index -U R1_virusLTR5_edge.fastq | grep "LTR5" | awk '{print $1}' | grep -v "^@"  > readID_LTR5Found.txt
```

-L represent the seed substring (minimum 3 and maximum LTR length), -N the number of mismatch in the initial seed alignment (0 or 1), -x the bowtie2 index, -U to specify that data are unpaired. Reads mapping to the LTR are obtained after running grep "LTR3" | awk '{print $1}' | grep -v "^@".

Intersect the results with R1/R2 reads to only get sequences mapping to the linker and LTR

```
bbmap/filterbyname.sh in=R1_linkerFound.fastq in2=R2_linkerFound.fastq out=R1_LTR3Found_linkerFound.fastq out2=R2_LTR3Found_linkerFound.fastq names=readID_LTR3Found.txt include=t overwrite=t
bbmap/filterbyname.sh in=R1_linkerFound.fastq in2=R2_linkerFound.fastq out=R1_LTR5Found_linkerFound.fastq out2=R2_LTR5Found_linkerFound.fastq names=readID_LTR5Found.txt include=t overwrite=t
```

Go from a 0-base to 1-base system:

```
((LTR5len++))
((LTR3len++))
((linkLen++))
```

### 3. Trim the reads

In order to get a better mapping efficiency alignment onto host genome is performed after trimming of LTR, linker and random tags.

First start by removing the linker sequence from R2:

```
fastx_trimmer -Q33 -f $linkLen -i R2_LTR3Found_linkerFound.fastq -o R2_LTR3Found_linkerFound_Host.fastq
fastx_trimmer -Q33 -f $linkLen -i R2_LTR5Found_linkerFound.fastq -o R2_LTR5Found_linkerFound_Host.fastq
```

Then LTR sequences are removed from R1:

```
fastx_trimmer -Q33 -f $LTR3len -i R1_LTR3Found_linkerFound.fastq -o R1_LTR3Found_linkerFound_Host.fastq
fastx_trimmer -Q33 -f $LTR5len -i R1_LTR5Found_linkerFound.fastq -o R1_LTR5Found_linkerFound_Host.fastq
```

Another important source of noise in clonality libraries is the overlap between R1-R2 reads. This overlap creates softclipping in R1 due to the sequencing of R2-linker and in R2 due to sequencing of R1-LTR. To remove them cutadapt is used. Here HTLV ATK is taken as example.

On R1:

```
# Specify LTR sequences to remove
LTR3forward_firstBp=LTR3sequence_firstbp	# [TTAGTACACA]
LTR3Complement_firstBp=$(echo $LTR3forward_firstBp | tr "[ATGCatgcNn]" "[TACGtacgNn]")	# [AATCATGTG]
LTR5reverse_firstBp=`rev $LTR5sequence_firstbp`	# [TGACAATGAC]
LTR5reverseComplement_firstBp=$(echo $LTR5reverse_firstBp | tr "[ATGCatgcNn]" "[TACGtacgNn]")	# [ACTGTTACT]

cutadapt --quiet -n 2 -m 10 -g $LTR3forward_firstBp R1_LTR3Found_linkerFound_Host.fastq | cutadapt --quiet -n 2 -m 10 -g $LTR3Complement_firstBp - > R1_LTR3Found_linkerFound_Host.trimmed.fastq
cutadapt --quiet -n 2 -m 10 -g $LTR5reverse_firstBp R1_LTR5Found_linkerFound_Host.fastq | cutadapt --quiet -n 2 -m 10 -g $LTR5reverseComplement_firstBp - > R1_LTR5Found_linkerFound_Host.trimmed.fastq
```

On R2:

```
# Specify LTR sequences to remove
LTR3forward_firstBp=LTR3sequence_firstbp	# [TGTGTACTAA]
LTR3reverse_firstbp=`rev $LTR3sequence_firstbp`	# [ACACATGATT]
LTR5reverse_firstbp=`rev $LTR5sequence_firstbp`	# [ACTGTTACTG]
LTR5forward_firstBp=LTR5sequence_firstbp	# [TGACAATGAC]

cutadapt --quiet -m 10 -a $LTR3reverse_firstbp R2_LTR3Found_linkerFound_Host.fastq | cutadapt --quiet -m 10 -a $LTR3forward_firstBp - > R2_LTR3Found_linkerFound_Host.trimmed.fastq
cutadapt --quiet -m 10 -a $LTR5reverse_firstbp R2_LTR5Found_linkerFound_Host.fastq | cutadapt --quiet -m 10 -a $LTR5forward_firstBp - > R2_LTR5Found_linkerFound_Host.trimmed.fastq
```

Finally remaining of linker sequence are trimmed:

```
revLinker=$(echo $link | rev)
revCompLinker=$(echo $revLinker | tr "[ATGCatgcNn]" "[TACGtacgNn]")

cutadapt --quiet -m 10 -a $revLinker R1_LTR3Found_linkerFound_Host.trimmed.fastq | cutadapt -m 10 --quiet -a $revCompLinker - > R1_LTR3Found_linkerFound_Host.trimmed2.fastq
cutadapt --quiet -m 10 -a $revLinker R1_LTR5Found_linkerFound_Host.trimmed.fastq | cutadapt -m 10 --quiet -a $revCompLinker - > R1_LTR5Found_linkerFound_Host.trimmed2.fastq
```

After trimming some reads might have too short R1 or an absence of mate. Pairs are resynchronised: 

```
python resynchronizePairs.pu R1_LTR3Found_linkerFound_Host.trimmed2.fastq R2_LTR3Found_linkerFound_Host.trimmed.fastq
python resynchronizePairs.py R1_LTR5Found_linkerFound_Host.trimmed2.fastq R2_LTR5Found_linkerFound_Host.trimmed.fastq
```

Reassign R1/R2 values to the cleaned FASTQ files:

```
R1_LTR3="R1_LTR3Found_linkerFound_Host.trimmed2.fastq_pairs_R1.fastq"
R2_LTR3="R2_LTR3Found_linkerFound_Host.trimmed.fastq_pairs_R2.fastq"
R1_LTR5="R1_LTR5Found_linkerFound_Host.trimmed2.fastq_pairs_R1.fastq"
R2_LTR5="R2_LTR5Found_linkerFound_Host.trimmed.fastq_pairs_R2.fastq"
```

### 4. Align onto host genome

#### 4.1.1. Extract best quality reads (i.e. proper pairs) and separated first/second strands

Align onto host genome with viral genome as separated chromosome. Allow -N for one mismatch maximum in the initial seed sequence. Although only primary alignments are reported reads mapping to >11 positions are excluded (-k 10).  

```
bowtie2 -p 2 --very-sensitive -N 1 --no-mixed -x $bowtie2_virusHost_index -1 $R1_LTR3 -2 $R2_LTR3 > LTR3_candidateIS_bowtie2_BEST.sam
bowtie2 -p 2 --very-sensitive -N 1 --no-mixed -x $bowtie2_virusHost_index -1 $R1_LTR5 -2 $R2_LTR5 > LTR5_candidateIS_bowtie2_BEST.sam
```

#### 4.1.2. SAM to BAM

```
$samtools view -Sb LTR3_candidateIS_bowtie2_BEST.sam > LTR3_candidateIS_bowtie2_BEST.bam
$samtools sort LTR3_candidateIS_bowtie2_BEST.bam LTR3_candidateIS_bowtie2_BEST.sorted
$samtools index LTR3_candidateIS_bowtie2_BEST.sorted.bam
 
$samtools view -Sb LTR5_candidateIS_bowtie2_BEST.sam > LTR5_candidateIS_bowtie2_BEST.bam
$samtools sort LTR5_candidateIS_bowtie2_BEST.bam LTR5_candidateIS_bowtie2_BEST.sorted
$samtools index LTR5_candidateIS_bowtie2_BEST.sorted.bam
```

#### 4.2.1. Extract alternative position of each reads

```
bowtie2 -p 2 --very-sensitive -k 11 -N 1 --no-mixed -x $bowtie2_virusHost_index -1 $R1_LTR3 -2 $R2_LTR3 > LTR3_candidateIS_ALTERN_bowtie2.sam
bowtie2 -p 2 --very-sensitive -k 11 -N 1 --no-mixed -x $bowtie2_virusHost_index -1 $R1_LTR5 -2 $R2_LTR5 > LTR5_candidateIS_ALTERN_bowtie2.sam
```

#### 4.2.2. SAM to BAM

```
$samtools view -Sb LTR3_candidateIS_bowtie2_ALTERN.sam > LTR3_candidateIS_bowtie2_ALTERN.bam
$samtools sort LTR3_candidateIS_bowtie2_ALTERN.bam LTR3_candidateIS_bowtie2_ALTERN.sorted
$samtools index LTR3_candidateIS_bowtie2_ALTERN.sorted.bam

$samtools view -Sb LTR5_candidateIS_bowtie2_ALTERN.sam > LTR5_candidateIS_bowtie2_ALTERN.bam
$samtools sort LTR5_candidateIS_bowtie2_ALTERN.bam LTR5_candidateIS_bowtie2_ALTERN.sorted
$samtools index LTR5_candidateIS_bowtie2_ALTERN.sorted.bam
```

### 5. Extract reads mapping to the provirus (for quantification purposes):

```
samtools view -F 0x40 -q 30 -S LTR3_candidateIS_bowtie2.sam | awk -v vir=$virus '{if($3 == vir){print $0}}' > LTR3_pureViralReads_bowtie2.sam
samtools view -F 0x40 -q 30 -S LTR5_candidateIS_bowtie2.sam | awk -v vir=$virus '{if($3 == vir){print $0}}' > LTR5_pureViralReads_bowtie2.sam
```

### 6. Extract random tag sequence from the R2 raw data

R2 random tag is located at the first 8bp of R2 reads

```
fastx_trimmer -Q33 -l 8 -i R2.fastq | fastq_to_fasta -Q33 | fasta_formatter -t | awk '{print $1"\t"$3}' > R2_randomTag.txt
```

### 7. Extract integration sites and compute viral abundances:

```
Rscript run_PIC.R "LTR3_candidateIS_bowtie2_BEST.sorted.bam" "LTR5_candidateIS_bowtie2_BEST.sorted.bam" "LTR3_candidateIS_bowtie2_ALTERN.sorted.bam" "LTR5_candidateIS_bowtie2_ALTERN.sorted.bam" "R2_randomTag.txt" $name $geneBedFile "R1.fastq" $virus 30
```

run_PIC.R simply contains the following arguments

```
library(PIC)

args <- commandArgs(trailingOnly = TRUE)

PIC(
LTR3.args = args[1],
LTR5.args = args[2],
LTR3.altern = args[3],
LTR5.altern = args[4],
randomTag.args = args[5],
sampleName.args = args[6],
geneBedFile.args = args[7],
rawFASTQ.args = args[8],
virus.args = args[9],
mapqSTRINGENT = args[10]
)
```

PIC R wrapper function for clonality analysis. To better understand what is achieved behind the hood, please refer yourself to PIC library [PIC R library](https://github.com/GIGA-AnimalGenomics-BLV/PIC)
