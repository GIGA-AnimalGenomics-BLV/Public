# R PIC: Get IS Table

## INTRODUCTION

This package contains the ``PIC()`` wrapper function which takes the following arguments:

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

Two types of mapping are needed:

* **STRINGENT:** Only the best mapping position of each read (PRIMARY)
* **ALL:** Report up to 11 mapping position per read

Reads supporting the 5'LTR or 3'LTR are processed separately. 
 
### 1. ``loadClonalityData()``

**GOALS:** Load the BAM file. 
**RETURNS:** List of R1 reads, R2 reads and viral reads (mapping only to the virus)

|read_id|flag|chr|pos|mapq|cigar|strand|AS|XS|isPRIMARY|isSTRINGENT|N_ScoreLower|numberAlignment|
|:-:|---|---|---|---|---|---|---|---|---|---|---|---|
|M00991:68:000000000-AU82H:1:1105:21631:18279|147|chr10|2656307|1|12M|-|-20|NA|TRUE|FALSE|0|1|
|M00991:68:000000000-AU82H:1:1101:18425:23433|163|chr10|10614734|1|13M|+|-20|NA|TRUE|FALSE|0|1|
|M00991:68:000000000-AU82H:1:1110:17539:21149|163|chr10|11813034|33|20M|+|-20|NA|TRUE|TRUE|0|1|

In addition of the [SAM fields](https://samtools.github.io/hts-specs/SAMv1.pdf), reads are flagged with:

* isPRIMARY: Is it a primary alignment ?
* isSTRIGNENT: Is the mapping quality (MAPQ) > mapqSTRINGENT

### 2. ``getISposition()``

**GOALS:** Get the IS, shear sites and random tags of each read with ``extractISposition()`` then call IS positions and abundances
**OPTIONS:** [RECALL=TRUE] to compute the RECALL abundances





As ``extractISposition()`` can be slow, with the ``SAVE=TRUE`` an intermediate table summarising IS, shear sites and random tag is saved (under: ``sampleName.args`` "_ISposFixed_QUAL_", ``LTR``, ".txt"). This file will be automatically reused if the function is run twice on the same dataset.


1. 'STRIGENT' reads and 'ALL' reads are loaded using loadClonalityData()
2. IS positions and abundances are called with getISposition()

	* extractISposition(..., STRINGENT_reads): get a table containing readID, IS, shear sites and random.
	* Reads sharing the same IS are regrouped and IS abundance is calculated with either:
		* RAW READS: Abundance based on all the reads including PCR duplicates
		* FILTERED READS: Abundance without PCR duplicates (collapse reads sharing IS/random_tag/shear_site.
	* IS located within a 75bp upstream/downstream window are collapsed as they probably arise from mapping errors.
		* RAW and FILTERED reads are summed
		* Position of the IS with the highest number of FILTERED reads is reported
	* If RECALL=TRUE, the abundance of each IS called with the STRINGENT reads is recalculated using all the reads
		* extractISposition(..., ALL_reads)
		* Cross these IS with the STRINGENT IS
		* For each IS, get the reads within a 600 bp up/down window
		* Compute the RECALL abundance with these new reads (RAW.RECALL and FILTERED.RECALL)

3. Combine LTR5 and LTR3 tables with mergeLTRs_V2()

	* LTR3 IS within 25 bp upstream/downstream window of an LTR5 IS are merged
		* When available LTR3 position is systematically reported
		* To avoid artificial inflation of the IS abundance, the max(LTR3 abundance, LTR5 abundance) is reported

4. Gene annotation is added using annotateIS()

5. Statistics of the run are computed: getStatistics()

6. Final results are outputed.


This function creates several outputs:

1. IS table
2. Statistic table


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
		* *ENTROPY_RECURRENCE:* entropy < 0.85 (``filt.recurrence``)
		* *ENTROPY_ABUNDANCE:* entropy < 0.85 (``filt.abundance``) & at least one animal needs to have > 5 reads supporting the IS (``minReadMax``)
		* *NON-SPECIFIC:* the IS is detected in > 6% of the samples (``nonSpecific``)
		* *UNIQUE:* IS only detected in one individual
		* *DUBIOUS:* All the other IS
5. The TRUE owner of each IS is assigned.
	* *NON-SPECIFIC & DUBIOUS:* NA
	* *UNIQUE:* Only one
	* *ENTROPY_RECURRENCE:* The owner is the individual with the highest recurrence 
	* *ENTROPY_ABUNDANCE:* The owner is the individual with the maximal abundance
	
	
The function reports either the provided IS table with CATEGORY annotations as a new column (``report = FALSE``) or an intermediate table containing the recurrence, maximal abundances, shannon entropies, *etc* for each IS in each individual. This table can used to determine the right entropy parameters.

|index|ID|recurence|max.abundance|numberSample|numberIndividuals|e.recurrence|e.abundance|CATEGORY|position|
|:-:|---|---|---|---|---|---|---|---|---|
|198833|229|1|1|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390
|198833|233|27|362|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390
|198833|234|1|4|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390
|198833|236|2|2|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390

For is particular case, one IS (``position``) is detected in 31 samples (``numberSample``) coming from 4 individuals (``numberIndividuals`` & ``ID``). This IS is detected 27th times in 233 (``recurrence``) with a maximal abundance of 362 reads (``max.abundance``). Shannon entropies of the abundances (``e.abundance``) and recurrence (``e.recurrence``) show a clear bias for individual 233. 

## R: integrationSiteMotif() 

**GOAL:** Extract nucleotide motifs surrounding a single genomic position.

The function takes an IS table containing at least the following columns:

|seqnames|start|strand|
|:-:|---|---|
|chr1|4010901|*|
|chr2|1103712|*|

Strand can be +, - or * (for undetermined).

```
integrationSiteMotif(IS = NULL, win = 20, fasta = "path/to/genome.fasta")
```

Nucleotides sequences surrounding each IS are retrieved from the genome fasta file (``fasta``) using an up/downstream window of ``win`` bases.

SeqLogo graphics can be plotted using [ggseqlogo](https://github.com/omarwagih/ggseqlogo)

```

```


 