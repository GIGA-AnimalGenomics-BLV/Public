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

STRINGENT: Only the best mapping position of each read
ALL: Report up to 11 mapping position per read

## PIC() FUNCTION

LTR3 and LTR5 are first treated separately with the following process:

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


## tagContamination() FUNCTION

## integrationSiteMotif() FUNCTION
