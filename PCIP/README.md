# Pooled CRISPR Inverse PCR sequencing 

Pooled CRISPR Inverse PCR sequencing (PCIP-seq) protocol has been designed to sequence the insertion site and its associated provirus using Oxford Nanopore long-reads technology. The protocol can be found in [PCIP-seq](https://www.biorxiv.org/content/10.1101/558130v2).

<p align="center">
  <img src="WORKFLOW/Protocol.jpg">
</p>

## SUMMARY: Analysis

**1.** Raw Data Basecalling (.fast5 to .fastq).
**2.** Mapping to the TARGET-HOST Reference.
**3.** Calling Integration Sites (IS).
**4.** Calling Integration Sites Specific Variants.

**DISCLAIMER:** Given the extreme diversity of proviruses, integration sites and genomic alterations, the following pipeline should not be considered as universal. It can be used to prioritize interesting integration sites (IS) but needs to be adapted to your personnel needs. 

## PREREQUISITES

* albacore (≥ 2.3.1) or guppy (≥ 3.1.5.)
* [porechops (≥ 0.2.4.)](https://github.com/rrwick/Porechop) 
* [samtools (≥ 1.9.)](http://samtools.sourceforge.net/) 
* [minimap2 (≥ 2.10.)](https://github.com/lh3/minimap2) 
* [loFreq (≥ 2.1.2.)](http://csb5.github.io/lofreq/) 
* [R ≥3.5.1](https://www.r-project.org/)
  - [Bioconductor (≥ 3.7.)](https://www.bioconductor.org/install/) 
  - [Rsamtools (≥ 1.32.3.)](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) 
  - [GenomicRanges (≥ 1.32.7.)](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
  - [changepoint (≥ 1.32.3.)](https://cran.r-project.org/web/packages/changepoint/index.html) 
  - [dplyr (≥ 0.7.8.)](https://cran.r-project.org/web/packages/dplyr/index.html)
  - [ggplot2 (≥ 2.3.1.)](https://cran.r-project.org/web/packages/ggplot2/index.html)
  - [magrittr (≥ 1.5.)](https://cran.r-project.org/web/packages/magrittr/index.html)
  - [purrr (≥ 0.2.5.)](https://cran.r-project.org/web/packages/purrr/index.html)
  - [readr (≥ 1.3.1.)](https://cran.r-project.org/web/packages/readr/index.html)
  - [tibble (≥ 2.0.1.)](https://cran.r-project.org/web/packages/tibble/index.html)
  - [tidyr (≥ 0.8.2.)](https://cran.r-project.org/web/packages/tidyr/index.html)

Code has been run on a Linux fedora 7.2 ('Nitrogen'). 

### Genome:

PCIP-Seq pipeline is based on the detection of chimeric TARGET-HOST reads using a chimeric viral-host reference genome. For instance, HTLV-1 - hg38:

```
GENOME="/path/to/TARGETHOST_INDEX/Homo_sapiens.GRCh38.dna.chromosome.fa"
TARGET="/path/to/TARGETHOST_INDEX/HTLV_genome.fa"

mkdir TARGETHOST_INDEX
cd TARGETHOST_INDEX

cat $GENOME $TARGET > HSA_GRCh38_HTLV.fa

TARGETHOST="/path/to/TARGETHOST_INDEX/HSA_GRCh38_HTLV.fa"
```

### Variables

A serie of variables have to be provided as well:

- RAW - nanopore data generated as described in [PCIP-seq](https://www.biorxiv.org/content/10.1101/558130v2) in FAST5 format.
- ANNOTATION - Path to gene annotation in GTF.gz format.
- NAME - Name of the sample
- TARGETNAME - Name of the TARGET chromosome  as the FASTA sequence (i.e., HTLV_ATK, HIV_U1, BLV, ...). 
- OUTDIR - Path to the output directory. 
- THREAD - Number of threads.

```
RAW="/path/to/raw.fast5"
ANNOTATION="/path/to/Homo_sapiens.GRCh38.gtf.gz"
NAME="HTLV_PatientA"
TARGETNAME="HTLV_ATK"
OUTDIR="/path/to/output/HIV_PatientA"
THREAD=4
LENGTHTARGET=9091 #base-pairs
```

### Indexing:

Both Minimap2 and LoFreq require indexing of the TARGETHOST FASTA files:

```
minimap2 -d /path/to/TARGETHOST_INDEX/HSA_GRCh38_HTLV.mni $TARGETHOST
samtools faidx $TARGETHOST

TARGETHOSTINDEX="/path/to/TARGETHOST_INDEX/HSA_GRCh38_HTLV.mni"
```

## PART 1: BASECALLING

Basecalling has been performed using guppy in high-accuracy mode. 

```
# Basecalling
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 --recursive --output_format fastq --input $RAW --save_path /path/to/mySample.fastq -t $THREAD

# Trim Adapters
porechop --discard_middle -i /path/to/mySample.fastq -b /path/to/mySample_trimmed.fastq

FASTQ="/path/to/mySample_trimmed.fastq"
```

## PART 2: MAPPING TO THE REFERENCE

Detection of new integration sites (IS) requires a [Pairwise mApping Format \(PAF\)](https://github.com/lh3/miniasm/blob/master/PAF.md) alignment as input. 

```
minimap2 -cx map-ont -t $THREAD $TARGETHOSTINDEX $FASTQ > "$NAME"_minimap2_TARGETHOST.paf
```

If you want to visualize the alignment results using [Integrated Genome Viewer \(IGV\)](http://software.broadinstitute.org/software/igv/) minimap2 needs to output a [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) file as well.

```
minimap2 -ax map-ont -t $THREAD $TARGETHOSTINDEX $FASTQ > "$NAME"_minimap2_TARGETHOST.sam
samtools view -@ $THREAD -Sb "$NAME"_minimap2_TARGETHOST.sam > "$NAME"_minimap2_TARGETHOST.bam
samtools sort -@ $THREAD "$NAME"_minimap2_TARGETHOST.bam -o "$NAME"_minimap2_TARGETHOST.sorted.bam
samtools index -@ $THREAD "$NAME"_minimap2_TARGETHOST.sorted.bam
```

## PART 3: EXTRACTING INTEGRATION SITES

R functions to extract integration sites and the corresponding clone abundance can be downloaded from this repository (PCIP_1.0.tar.gz). After installing R prerequisites, run the following command in the terminal to install PCIP package:

```
R CMD INSTALL PCIP_1.0.tar.gz
```

Integration site detection is performed using the following R functions. A detailled description of each function can be found at:. The following command can be run using R:

```
# 0. Load functions 
library(PCIP)

# 1. Initialize Variables:
PAF.path = "path/to/"$NAME"_minimap2_TARGETHOST.sorted.paf"
targetName = "HTLV_ATK"
lengthTarget = 9091
distanceLTR = 200

# 2. Read PAF file:
PAF <- readPairwiseAlignmentFile(alignFile = PAF.path)

# 3. Filter data and keep chimeric reads:
PAF.filter <- PCIP_filter(minimap2PAF = PAF, targetName = targetName)

# 4. Get Target-Genome breakpoints:
PAF.breakpoints <- PCIP_getBreakPoints(PAF = PAF.filter, lengthTarget = lengthTarget, targetName = targetName)

# 5. Group viral-host junctions into IS:
PAF.integrationSite <- PCIP_summarise(PCIPbreakpoints = PAF.breakpoints, distanceLTR = distanceLTR, mergeShearSiteDistance = 0)

# 6. Save:
write.table(PAF.integrationSite[[1]], paste0(out.prefix, "-insertionSites.txt"), sep = "\t", row.names = F, quote = F)
write.table(PAF.integrationSite[[2]], paste0(out.prefix, "-splitFASTQ.txt"), sep = "\t", row.names = F, quote = F)
```

## PART 4: VIRAL/HOST ALTERATIONS

R PCIP functions returns two outputs as a list. First, an integration site table (``insertionSites.txt``). Second, the read IDs associated to each IS (``splitFASTQ.txt``). These can be used to extract the reads associated to a specific integration site. 

For each unique IS, extract the target-specific and the host mutations (within +/- 20Kb).

```
awk -F'\t' 'NR>1 {print $2}' path/to/mySample-splitFASTQ.txt | sort | uniq > mySample_uniqueIS.txt

while IFS=$'\t' read -r -a line
do
  ## 1. Create an IS specific folder:
  # ID = chr_start_end
  ID=`echo "${line[0]}" | sed 's/\:/_/g' | sed 's/\-/_/g'`
  mkdir $ID
  cd $ID
  
  ## 2. Extract reads belonging to that particular IS and remap them:
  grep -F "${line[0]}" mySample-splitFASTQ.txt | awk -F'\t' '{print "@"$1}' > "$ID"_readID.fastq
  
  ## 3. Remap them to the viral sequence:
  minimap2 -ax map-ont -t $THREAD $TARGETHOSTINDEX "$ID".fastq > "$ID"_extracted.sam
  samtools view -@ $THREAD -Sb "$ID"_extracted.sam > "$ID"_extracted.bam
  samtools sort -@ $THREAD "$ID"_extracted.bam -o "$ID"_extracted.sorted.bam
  samtools index -@ $THREAD "$ID"_extracted.sorted.bam

  ## 4. Call single-nucleotide polymorphisms:
  
  ### 4.1. Create a .bed file containing the position to call. A 20Kb window is added upstream and downstream:
  win=20000
  IFS=_ read chr start end <<< $ID
  echo -e $chr'\t'`expr $start - $win`'\t'`expr $end + $win` > GENOME.bed
  echo -e $TARGETNAME'\t'1'\t'$LENGTHTARGET > TARGET.bed

  ## 4.2. Call variant on the HOST GENOME
  loFreq call-parallel --pp-threads $THREAD -f $GENOME -l GENOME.bed "$ID"_extracted.sorted.bam > mySample_"$ID"_LoFreq_GENOME.vcf
  
  ## 4.3. Call variant on the TARGET GENOME
  loFreq call-parallel --pp-threads $THREAD -f $TARGET -l TARGET.bed "$ID"_extracted.sorted.bam > mySample_"$ID"_LoFreq_TARGET.vcf

  # Custom scripts or [NGLMR + Sniffles](https://github.com/fritzsedlazeck/Sniffles) can be used to call genomic rearrangements.

done < mySample_uniqueIS.txt

```

----

# CALLING INTEGRATION SITES: DETAILS

