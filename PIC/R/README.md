# R PIC: Get IS Table

This pipeline can be used to process data generated following Rosewick *et al* (2017) or Artesi *et al* (2017). 

## INTRODUCTION

This R package contains the ``PIC()`` wrapper function which detect IS and assign their abundance. The following arguments are required:

* **STRINGENT:** Only the best mapping position of each read (PRIMARY)
* **ALL:** Report up to 11 mapping position per read

|VARIABLE|TYPE|DEFINITION|
|:-:|---|---|
|LTR3.args|**character**|path/to/file.bam reads supporting LTR3 IS (STRINGENT)
|LTR5.args|**character**|path/to/file.bam reads supporting LTR5 IS (STRINGENT)
|LTR3.altern|**character**|path/to/file.bam reads supporting LTR3 IS (ALL)
|LTR5.altern|**character**|path/to/file.bam reads supporting LTR5 IS (ALL)
|randomTag.args|**character**|path/to/file.txt random tag of each reads
|sampleName.args|**character**|sample_name
|geneBedFile.args|**character**|path/to/annotation.bed
|rawFASTQ.args|**character**|path/to/file.fastq R1 fastq ID, for statistics purpose
|winRecall|**numeric**|size of the RECALL window
|virus.args|**character**|viral_chromosome_name
|mapqSTRINGENT|**numeric**|MAPQ of stringent reads

### **NB:** Reads supporting the 5'LTR or 3'LTR are first processed separately

--- 
	
### 1. ``loadClonalityData()``

* **GOALS:** Load the BAM file. 
* **RETURNS:** List of R1 reads, R2 reads and viral reads (mapping only to the virus)

|read_id|flag|chr|pos|mapq|cigar|strand|AS|XS|isPRIMARY|isSTRINGENT|N_ScoreLower|numberAlignment|
|:-:|---|---|---|---|---|---|---|---|---|---|---|---|
|M00991:68:000000000-AU82H:1:1105:21631:18279|147|chr10|2656307|1|12M|-|-20|NA|TRUE|FALSE|0|1|
|M00991:68:000000000-AU82H:1:1101:18425:23433|163|chr10|10614734|1|13M|+|-20|NA|TRUE|FALSE|0|1|
|M00991:68:000000000-AU82H:1:1110:17539:21149|163|chr10|11813034|33|20M|+|-20|NA|TRUE|TRUE|0|1|

In addition of the [SAM fields](https://samtools.github.io/hts-specs/SAMv1.pdf), reads are flagged with:

* **isPRIMARY:** Is it a primary alignment ?
* **isSTRIGNENT:** Is the mapping quality (MAPQ) > mapqSTRINGENT
* **numberAlignment:** How many different alignments are returned by the aligner ?

--- 
### 2. ``getISposition()``

* **GOALS:** Get IS positions and abundances

#### 2.1. ``extractISposition()``:

* For each read, the IS, shear site and random tag is assigned, taking into account remaining soft-clipping (``getCIGARsize()`` & ``splitCIGAR()``).
* Flag reads in ``properPair`` based on their SAM flag (83, 99, 147 or 163)

|read_id|flag|chr|mapq|strand|AS|XS|isPRIMARY|isSTRINGENT|N_ScoreLower|numberAlignment|exactPosition|properPair|shearSite|randomTag|LTR|
|:-:|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|M00991:68:000000000-AU82H:1:1105:21631…|99|chr10|1|+|-33|NA|TRUE|FALSE|0|1|2656307|TRUE|chr10:2656307|ATTTAGCG|LTR3|
|M00991:68:000000000-AU82H:1:1101:18425…|83|chr10|1|-|-33|NA|TRUE|FALSE|0|1|10614746|TRUE|chr10:106147…|AAATGGAC|LTR3|
|M00991:68:000000000-AU82H:1:1110:17539…|83|chr10|33|-|-33|NA|TRUE|FALSE|0|1|11813053|TRUE|chr10:118130…|ACCTGAAA|LTR3|

* ``exactPosition`` corresponds to the viral-host junction.

#### 2.2. Collapse reads into IS:

The IS is called using the STRINGENT reads (in proper pairs and MAPQ > ``mapq.val``). We currently use MAPQ = 33. 

IS abundance is called using the STRINGENT reads. For debugging reasons we call in twice:

* **RAW ABUNDANCE:** All the reads are used.
* **FILTERED ABUNDANCE:** PCR duplicates are removed by collapsing reads displaying the exact same random tag and shear site. 

Although we only use the abundance called with the STRINGENT reads, we define two additional categories for debugging purposes:

1. **PROPER-PAIRS:** Abundance is called with only the reads in proper pairs with a MAPQ > 3
2. **ALL:** Abundance is called with all the reads that support the IS

These two categories are added to the final IS table and can be used to fine-tune the MAPQ parameters.

#### 2.3. Close IS are regrouped:

The same proviral insertion can sometimes be identified by several reads displaying slightly different IS, often due to sequencing or mapping errors. To account for such errors, we regroup IS located within ``mapgap`` up/downstream window (i.e., 75 bp). 

The position supported by the highest abundance is returned. All the reads (raw and filtered) are summed to compute the actual IS abundance. 

#### 2.4. RECALL:

If RECALL=TRUE, the IS abundance is recomputed using all the reads located within a ``winRecall`` up/downstream window, regardless of the MAPQ or SAM flag.
 
--- 
### 3. ``mergeLTRs_V2``

After processing separately the reads supporting the 3'LTR and 5'LTR, the two tables are merged. 

* 3'LTR IS within 25 bp of a 5' LTR are merged. 
	* The position reported is by default the one of the 3'LTR. If not available, the 5'LTR position is reported
	* The final IS abundance is computed using the following function: ``max(LTR5.filtered, LTR3.filtered)``
		* *This process avoid overestimating the IS abundance.
 
--- 
### 4. ``annotateIS()``

Add information about the position of each IS relative to the closest genes or genomic features. 
 
---  
### 5. ``getStatistics()``

Returns the run statistics.
 
--- 
### 6. ``Outputs``

Four outputs are created with ``PIC()``

|OUTPUT|DEFINITION|
|:-:|---|
|``sampleName.args``-clonalityResults.txt|Results without merging LTRs|
|``sampleName.args``-mergedIS.txt|Results after merging LTRs. TXT file|
|``sampleName.args``-mergedIS.xls|Results after merging LTRs. XLS file|
|``sampleName.args``-SIMPLIFIED_mergedIS.txt|Simplified table containing the merged LTRs results. Only the STRINGENT columns are reported|
|``sampleName.args``-statistics.txt|Run statistics|

Description of each field and examples are located in this github at [result field description](https://github.com/GIGA-AnimalGenomics-BLV/Public/tree/master/PIC/R/OutputFields.description.xlsx)
 
 
---

## R: tagContamination() 

Cross-contaminations between libraries as well as non-specific amplification of genomic positions happens in NGS clonality datasets. These need to be filtered out before comparing samples. We define four categories:

1. **UNIQUE:** IS found in only one individual
2. **NON-SPECIFIC:** IS found in > nonSpecific (%) of the sequenced libraries
3. **ENTROPY**: IS found in different individuals for which we can find the actual owner based on its distribution between individuals
	* *ENTROPY_RECURRENCE:* The IS is mainly more often detected in one individual
	* *ENTROPY_ABUNDANCE:* The IS is abundantly detected in one individual 
4. **DUBIOUS:** IS found in different individuals but with not information to retrieve the actual owner
 
DUBIOUS and NON-SPECIFIC IS are excluded from further analysis. To detect such IS we use the following function:

```
tagContamination(IS = NULL, nonSpecific = 6, filt.recurrence = 0.85, filt.abundance = 0.85, minReadMax = 5, mapgap = 5, report = FALSE)
```

### VARIABLES

The function takes an IS table containing at least the following columns:

|seqnames|start|filtered.max|ID|sample|
|:-:|---|---|---|---|
|chr1|4010901|100|sample1|sample1_library1|
|chr2|1103712|2000|sample1|sample1_library1|
|chr1|4010901|1|sample2|sample2_library1|
|chr3|531901|100|sample2|sample2_library1|

Fine-tuning of the entropy parameters is extremely important and will depend on the structure of your dataset (longitudinal samples, sequencing depth, ...). The ``tagContamination()`` is here as an example but might have to be adapted for your needs.

### DETAILS

``tagContamination()`` is divided in four parts.

1. IS are clustered with a small up/down windows (``maxgap``).
	* **WHY?** Exact position of a IS can be slightly shifted due to mapping errors, alternative detection of the 3'LTR or 5'LTR, mapping errors, mismatches, *etc*
2. For each IS, in every individuals, the recurrence (= number times it is detected) and maximal abundance are computed.
3. For each IS, the shannon entropy (log2) of the recurrences and maximal abundances are computed separately
4. Based on the filtering options the IS are separated in 5 CATEGORIES. 
	* We typically use this two sets of parameters:
		* **ENTROPY_RECURRENCE:** entropy < 0.85 (``filt.recurrence``)
		* **ENTROPY_ABUNDANCE:** entropy < 0.85 (``filt.abundance``) & at least one animal needs to have > 5 reads supporting the IS (``minReadMax``)
		* **NON-SPECIFIC:** the IS is detected in > 6% of the samples (``nonSpecific``)
		* **UNIQUE:** IS only detected in one individual
		* **DUBIOUS:** All the other IS
5. The TRUE owner of each IS is assigned.
	* **NON-SPECIFIC & DUBIOUS:** NA
	* **UNIQUE:** Only one
	* **ENTROPY_RECURRENCE:** The owner is the individual with the highest recurrence 
	* **ENTROPY_ABUNDANCE:** The owner is the individual with the maximal abundance
	
	
The function reports either the provided IS table with CATEGORY annotations as a new column (``report = FALSE``) or an intermediate table containing the recurrence, maximal abundances, shannon entropies, *etc* for each IS in each individual. This table can used to determine the right entropy parameters.

|index|ID|recurence|max.abundance|numberSample|numberIndividuals|e.recurrence|e.abundance|CATEGORY|position|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|198833|229|1|1|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390
|198833|233|27|362|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390
|198833|234|1|4|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390
|198833|236|2|2|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390

For is particular case, one IS (``position``) is detected in 31 samples (``numberSample``) coming from 4 individuals (``numberIndividuals`` & ``ID``). This IS is detected 27th times in 233 (``recurrence``) with a maximal abundance of 362 reads (``max.abundance``). Shannon entropies of the abundances (``e.abundance``) and recurrence (``e.recurrence``) show a clear bias for individual 233. 
 
 
---

## R: integrationSiteMotif() 

**GOAL:** Extract nucleotide motifs surrounding a single genomic position.

The function takes an IS table containing at least the following columns:

|seqnames|start|strand|
|:-:|:-:|:-:|
|chr1|4010901|*|
|chr2|1103712|*|

Strand can be +, - or * (for undetermined).

```
integrationSiteMotif(IS = NULL, win = 20, fasta = "path/to/genome.fasta")
```

Nucleotides sequences surrounding each IS are retrieved from the genome fasta file (``fasta``) using an up/downstream window of ``win`` bases.

SeqLogo graphics can be plotted using [ggseqlogo](https://github.com/omarwagih/ggseqlogo)

 