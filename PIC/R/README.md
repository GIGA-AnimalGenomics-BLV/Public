# R PIC: Get IS Table

## INTRODUCTION

This R package contains the ``PIC()`` wrapper function which detects integration sites (IS) and computes the relative abundance of the corresponding clone. The following arguments are required:

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

#### **NB:** Reads supporting the 5'LTR or 3'LTR are first processed separately

--- 
	
### 1. ``loadClonalityData()``

* **GOALS:** Load the BAM file. 
* **RETURNS:** List of R1 reads, R2 reads and viral reads (mapping only to the virus)

|read_id|flag|chr|pos|mapq|cigar|strand|AS|XS|isPRIMARY|isSTRINGENT|N_ScoreLower|numberAlignment|
|:-:|---|---|---|---|---|---|---|---|---|---|---|---|
|M00991:68:000000000-AU82H:1:1105:21631:18279|147|chr10|2656307|1|40M|-|-20|NA|TRUE|FALSE|0|1|
|M00991:68:000000000-AU82H:1:1101:18425:23433|163|chr10|10614734|1|35M|+|-20|NA|TRUE|FALSE|0|1|
|M00991:68:000000000-AU82H:1:1110:17539:21149|163|chr10|11813034|33|20M|+|-20|NA|TRUE|TRUE|0|1|

In addition of the [SAM fields](https://samtools.github.io/hts-specs/SAMv1.pdf), reads are flagged with:

* **isPRIMARY:** Is it a primary alignment ?
* **isSTRIGNENT:** Is the mapping quality (MAPQ) > mapqSTRINGENT
* **numberAlignment:** How many different alignments are returned by the aligner ?

--- 
### 2. ``getISposition()``

* **GOALS:** Get IS positions and abundances

#### 2.1. ``extractISposition()``:

* For each read, the IS, shear site and random tag are extracted, taking into account remaining soft-clipping (``getCIGARsize`` & ``splitCIGAR``).
* Flag best quality reads (in **proper pairs**) based on their SAM flag (83, 99, 147 or 163).

|read_id|flag|chr|mapq|strand|AS|XS|isPRIMARY|isSTRINGENT|N_ScoreLower|numberAlignment|exactPosition|properPair|shearSite|randomTag|LTR|
|:-:|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|M00991:68:000000000-AU82H:1:1105:21631…|99|chr10|1|+|-33|NA|TRUE|FALSE|0|1|2656307|TRUE|chr10:2656307|ATTTAGCG|LTR3|
|M00991:68:000000000-AU82H:1:1101:18425…|83|chr10|1|-|-33|NA|TRUE|FALSE|0|1|10614746|TRUE|chr10:106147…|AAATGGAC|LTR3|
|M00991:68:000000000-AU82H:1:1110:17539…|83|chr10|33|-|-33|NA|TRUE|FALSE|0|1|11813053|TRUE|chr10:118130…|ACCTGAAA|LTR3|

* ``exactPosition`` corresponds to the virus-host junction.

#### 2.2. Collapse reads into IS:

Reads sharing the same virus-host junction are grouped as IS. The exact IS is called using the STRINGENT reads (i.e., proper pairs and MAPQ > ``mapq.val``). We currently use MAPQ = 33. 

The  IS abundance is called using the STRINGENT reads:

* **RAW ABUNDANCE:** All the reads are used.
* **FILTERED ABUNDANCE:** PCR duplicates are removed by collapsing reads that show the same random tag and the same shear site. 

Although we only use the abundance called with the STRINGENT reads, we define two additional types of abundance for debugging purposes:

1. **PROPER-PAIRS:** Abundance is called with the reads in proper pairs with a MAPQ > 3.
2. **ALL:** Abundance is called with all the reads that could support the IS.

These two types are added to the final IS table and can be used to fine-tune the MAPQ parameters.

#### 2.3. Nearby IS are grouped:

The same proviral insertion site is sometimes be identified by multiple reads each displaying slightly different virus-host junctions, (i.e., due to sequencing or mapping errors). To account for such errors, we group IS located within ``mapgap`` up/downstream window (i.e., 75 bp). 

The position supported by the highest number of reads is returned as the IS. All the reads are summed to compute the IS abundance. 

#### 2.4. RECALL:

If RECALL=TRUE, the IS abundance is recomputed using all the reads located within a ``winRecall`` up/downstream window, regardless of the MAPQ or SAM flag.
 
--- 
### 3. ``mergeLTRs_V2``

After separately processing  reads supporting the 3'LTR and 5'LTR, the two tables are merged. 

* 3'LTR IS within 25 bp of a 5' LTR are considered as supporting the same IS. 
	* The position reported is by default the position of the 3'LTR. If not available, the 5'LTR position is reported.
	* The final IS abundance is computed using the following function: ``max(LTR5.filtered, LTR3.filtered)``
		* This process avoids overestimating the IS abundance.
 
--- 
### 4. ``annotateIS()``

Adds information about the position of each IS relative to the closest genes or genomic features. 
 
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

Cross-contaminations observed in next-generation clonality sequencing datasets mainly results from carry-over during the sequencing process itself. Recurrent IS due to contamination need to be filtered out. We define three categories of IS:

1. **UNIQUE:** IS found in only one individual
2. **ENTROPY**: IS found in different individuals. The IS distribution is skewed for one individual. The actual 'owner'/sample of these IS can be retrieved. 
	* *ENTROPY_RECURRENCE:* The IS is mainly detected in one individual.
	* *ENTROPY_ABUNDANCE:* The IS is abundantly detected in one individual.
3. **UNCERTAIN:** IS found in several individuals but the information is not sufficient to assign the IS to its 'owner' with confidence.
 
UNCERTAIN IS are excluded from further analysis. To detect such IS we use the following function:

```
tagContamination(IS = NULL, filt.recurrence = 0.85, filt.abundance = 0.85, minReadMax = 5, mapgap = 0, report = FALSE)
```

### VARIABLES

The function takes an IS table containing the following columns:

|seqnames|start|filtered.max|ID|sample|
|:-:|---|---|---|---|
|chr1|4010901|100|sample1|sample1_library1|
|chr2|1103712|2000|sample1|sample1_library1|
|chr1|4010901|1|sample2|sample2_library1|
|chr3|531901|100|sample2|sample2_library1|

Fine-tuning of the entropy parameters is extremely important and will depend on the structure of your dataset (longitudinal samples, sequencing depth, ...). The ``tagContamination()`` is shown here as an example but it may require some adjustments for your needs.

### DETAILS

``tagContamination()`` has four parts.

1. IS are clustered within a small up/down windows (``mapgap``).
	* This option was set to 0 in our current pipeline but could be increased to account for mapping errors.
2. For each IS, in each individual, the recurrence - number of times it is detected - and maximal abundance are computed.
3. For each IS, the shannon entropy (log2) of the recurrences and maximal abundance are computed separately. High recurrence or abundance in one individual can be used to find the IS's owner. 
4. Based on the filtering options the IS are separated in 4 CATEGORIES. 
	* **ENTROPY_RECURRENCE:** entropy < 0.85 (``filt.recurrence``).
	* **ENTROPY_ABUNDANCE:** entropy < 0.85 (``filt.abundance``) and at least one individual with > 5 reads supporting the IS (``minReadMax``).
	* **UNIQUE:** IS only detected in one individual.
	* **UNCERTAIN:** All remaining IS.
5. The TRUE owner of each IS is assigned.
	* **UNCERTAIN:** NA
	* **UNIQUE:** Only one.
	* **ENTROPY_RECURRENCE:** The owner is the individual with the highest recurrence.
	* **ENTROPY_ABUNDANCE:** The owner is the individual with the highest abundance.
	
	
The function reports either the provided IS table with CATEGORY annotations as a new column (``report = FALSE``) or an intermediate table containing the recurrence, maximal abundances, shannon entropies, *etc* for each IS in each individual. This table can be used to determine the appropriate entropy parameters.

|index|ID|recurence|max.abundance|numberSample|numberIndividuals|e.recurrence|e.abundance|CATEGORY|position|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|198833|229|1|1|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390
|198833|233|27|362|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390
|198833|234|1|4|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390
|198833|236|2|2|31|4|0.748|0.162|ENTROPY_RECURRENCE|OAR17:4806385-4806390

Regarding this particular case, one IS (``position``) is detected in 31 samples (``numberSample``) originating from 4 individuals (``numberIndividuals`` & ``ID``). This IS is detected 27 times in individual 233 (``recurrence``) with an abundance of 362 reads (highest) (``max.abundance``). Shannon entropies of the abundance (``e.abundance``) and recurrence (``e.recurrence``) show a clear bias for individual 233. 
 
### FINDING THE RIGHT PARAMETERS:

To find the parameters fitting our dataset we currently use the following procedure:

1. Run the pipeline on the data using default parameters (set report = TRUE).
2. Choose one of the parameters (entropy recurrence or filtered).
3. Order the data by the max(recurrence) or max(filtered-reads) among all the individuals. For instance, max(recurrence):

|IS|Animal_1|Animal_2|Animal_3|MAX|Entropy|
|:-:|:-:|:-:|:-:|:-:|:-:|
|IS-#1|30|1|0|30|0.2771134|
|IS-#2|8|1|1|8|0.6390319|
|IS-#3|4|1|1|4|0.8675632|
|IS-#3|2|1|1|2|1.039721|

4. Go through the list from top-to-bottom, looking at the decision taken by the algorithm. 

**NB:** In this example, there is no doubt that IS-#1 belongs to Animal #1 as it is seen 30 times in it but only once in Animal #2.

5. We currently set an empirical threshold at the first point in the list for which where we are unable to make a  clear decision. In the example above,  IS-#3 probably belongs to Animal #1 but we cannot be fully confident about this.
6. Threshold is set below the value of the last confident call (i.e., 0.85 here). 


 
