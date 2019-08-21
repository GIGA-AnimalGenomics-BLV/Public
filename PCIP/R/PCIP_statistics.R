#' @title get Statistics for PCIP-seq
#'
#' @description
#' Extract some useful statistics about the run
#'
#' @param pathToResults character. Path/to/results/
#' @param ID character. ID of the sample
#'
#' @return Write a statistics report in Path/to/results/$ID-statistics.txt
#'
#' @keywords PCIP
#'
#' @author Vincent Hahaut
#'
#' @export
PCIP_statistics <- function(pathToResults = NULL, ID = NULL){

	suppressPackageStartupMessages(library(tidyverse))

	raw.path <- list.files(paste0(pathToResults, "/nobackup/RAW/", collapse = ""), full.names = T)
	postFilter.path <- paste0(pathToResults, "/nobackup/2_SITES/", ID, "-readBreakPoints.txt", collapse = "")
	insertionTable.path <- paste0(pathToResults, "/nobackup/2_SITES/", ID, "_210bp-insertionSites.txt", collapse = "")
	variantTarget.path <- paste0(pathToResults,  "/nobackup/4_VARIANTS/loFreq/", ID, "_TARGET_loFreq.vcf", collapse = "")
	variantGenome.path <- paste0(pathToResults, "/nobackup/4_VARIANTS/loFreq/", ID, "_GENOME_loFreq.vcf", collapse = "")
	selectedProvirusForCall.path <- paste0(pathToResults, "/nobackup/3_SPLIT/", ID, "_uniqueIS.txt", collapse = "")
	selectedProvirusForCall_FILTERED.path <- paste0(pathToResults, "/nobackup/3_SPLIT/FASTQ/", collapse = "")
	deletionTable.targetSV <- paste0(pathToResults, "/nobackup/4_VARIANTS/changePoint/", ID, "_TARGET_changePoint.txt", collapse = "")
	deletionTable.sniffles <- paste0(pathToResults, "/nobackup/4_VARIANTS/sniffles/", ID, "_TARGET_Sniffles.vcf", collapse = "")

	# OPEN:
	widths <- as.numeric(system(paste0("zcat ", raw.path, " | awk 'NR%4==2 {print length}'"), intern = TRUE))
	raw <- length(widths)
	meanLength <- round(mean(widths), 2)
	sdLength <- round(sd(widths), 2)

	postFilter <- suppressMessages(read_tsv(postFilter.path))
	insertionTable <- suppressMessages(read_tsv(insertionTable.path))
	variantTarget <- suppressMessages(read_tsv(variantTarget.path, comment = "#", col_names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")))
	variantGenome <- suppressMessages(read_tsv(variantGenome.path, comment = "#", col_names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")))

	if(nrow(variantTarget) > 0){
  	variantTarget <- bind_cols(variantTarget %>%
  		select(-INFO),
  		apply(str_split(variantTarget$INFO, pattern = ";", simplify = T), 1, function(x) str_split(x, "=", simplify = T)[,2]) %>%
  		t() %>%
  		as_data_frame())
  	colnames(variantTarget) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "DP", "AF", "SB", "DP4", "VIRUS", "SAMPLE")
	} else {
	  variantTarget <- data_frame()
	}

	if(nrow(variantGenome) > 0){
	variantGenome <- bind_cols(variantGenome %>%
		select(-INFO),
		apply(str_split(variantGenome$INFO, pattern = ";", simplify = T), 1, function(x) str_split(x, "=", simplify = T)[,2]) %>%
		t() %>%
		as_data_frame())
	colnames(variantGenome) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "DP", "AF", "SB", "DP4", "VIRUS", "SAMPLE")
	} else {
	  variantGenome <- data_frame()
	}

	# SEQUENCING:
	numberReads <- raw
	numberReads_postFilter <- nrow(postFilter)
	percentageReads_useful <- round(numberReads_postFilter/numberReads, 3)

	# INTEGRATION SITES:
	IS <- nrow(insertionTable)
	IS.05 <- nrow(insertionTable %>% filter(count.max >5))
	IS.10 <- nrow(insertionTable %>% filter(count.max >10))
	IS.25 <- nrow(insertionTable %>% filter(count.max >25))
	IS.50 <- nrow(insertionTable %>% filter(count.max >50))

	# VARIANT DETECTION:
	passingFiltersONE <- suppressMessages(nrow(read_tsv(selectedProvirusForCall.path, col_names = F)))
	passingFiltersTWO <- length(list.files(selectedProvirusForCall_FILTERED.path))
	numberProvirusWithVariantsInTarget <- length(unique(variantTarget$VIRUS))
	numberProvirusWithVariantsInGenome <- length(unique(variantGenome$VIRUS))
	target <- nrow(variantTarget)
	target.60 <- nrow(variantTarget %>% filter(AF > 0.6))
	numberProvirusWithVariantsInTarget.60 <- length(unique(variantTarget$VIRUS[variantTarget$AF > 0.6]))

	genome <- nrow(variantGenome)
	genome.60 <- nrow(variantGenome %>% filter(AF > 0.6))

	# DELETIONS
	targetSV <- nrow(read_tsv(deletionTable.targetSV))
	sniffles <- nrow(read_tsv(deletionTable.sniffles))

	# WRITE:
	out <- paste0(pathToResults, "/", ID, "-statistics.txt", collapse = "")
	print(out)

    if( file.exists(out) ){
        file.remove(out )
    }

	write.table("", out, sep = "\t", quote = F, row.names = F, col.names = F, append = T)

  cat(paste0("# Date: ", Sys.time(), "\n"), file = out, append = TRUE)
	cat("\n\n\n##############\n\n\n", file = out, append = TRUE)
	cat("## SEQUENCING\n\n", file = out, append = TRUE)
	cat(paste0("Number reads: ", numberReads, "\n"), file = out, append = TRUE)
	cat(paste0("Number reads chimeric used: ", numberReads_postFilter, "\n"), file = out, append = TRUE)
	cat(paste0("Number reads chimeric used (%): ", 100*percentageReads_useful, "\n"), file = out, append = TRUE)
	cat(paste0("Average read size: ", meanLength, "\n"), file = out, append = TRUE)
	cat(paste0("Standard deviation read size: ", sdLength, "\n"), file = out, append = TRUE)

	cat("\n## INTEGRATION SITES\n\n", file = out, append = TRUE)
	cat(paste0("Number IS detected: ", IS, "\n"), file = out, append = TRUE)
	cat(paste0("Number IS detected (>5): ", IS.05, "\n"), file = out, append = TRUE)
	cat(paste0("Number IS detected (>10): ", IS.10, "\n"), file = out, append = TRUE)
	cat(paste0("Number IS detected (>25): ", IS.25, "\n"), file = out, append = TRUE)
	cat(paste0("Number IS detected (>50): ", IS.50, "\n"), file = out, append = TRUE)

	cat("\n## VARIANT DETECTION\n\n", file = out, append = TRUE)
	cat(paste0("FILTER1, Number provirus with >10 supporting reads (IS): ", passingFiltersONE, "\n"), file = out, append = TRUE)
	cat(paste0("FILTER2, Number provirus with >10 supporting reads after remapping: ", passingFiltersTWO, "\n"), file = out, append = TRUE)
	cat(paste0("Number provirus with variant in TARGET: ", numberProvirusWithVariantsInTarget, "\n"), file = out, append = TRUE)
	cat(paste0("Number variant detected TARGET: ", target, "\n"), file = out, append = TRUE)
	cat(paste0("Number variant detected TARGET (AF > 0.6): ", target.60, "\n"), file = out, append = TRUE)
	cat(paste0("Number provirus with variant in TARGET (AF > 0.6): ", numberProvirusWithVariantsInTarget.60, "\n"), file = out, append = TRUE)

	cat(paste0("Number provirus with variant in GENOME: ", numberProvirusWithVariantsInGenome, "\n"), file = out, append = TRUE)
	cat(paste0("Number variant detected GENOME: ", genome, "\n"), file = out, append = TRUE)
	cat(paste0("Number variant detected GENOME (AF > 0.6): ", genome.60, "\n"), file = out, append = TRUE)

	cat("\n## STRUCTURAL VARIANTS\n\n", file = out, append = TRUE)
	cat(paste0("Number of SV found with TargetSV: ", targetSV, "\n"), file = out, append = TRUE)
	cat(paste0("Number of SV found with Sniffles: ", sniffles, "\n"), file = out, append = TRUE)

}
