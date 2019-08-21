#' @title Call deletions from PCIP-seq Nanopore data
#'
#' @description
#' Primary function of this tool is to detect deletions lying within defined genomic segment.
#' Sudden drops of mean and variance of coverage is monitored through a change point PELT mean-var function.
#' Decreases and increases are use to delimit windows.
#' Windows harboring a ratio deletion / coverage of more than 80% are flagged as DELETION.
#'
#' @param bam.path character. Path/to/myfile.sorted.bam. Index (.bai, Samtools) is expected in the same folder
#' @param viralLength numeric. Size of the provirus
#' @param minseglen numeric. Minimal size of the deletion to be found (shorter than 10 is not recommended)
#' @param out character. Path/to/output/prefix
#' @param minCoverage numeric. Minimum coverage to call deletion
#' @param provirus character. Can be set to the TARGET chromosome name to call deletion only in the TARGET or "all" to avoid any pre-filtering
#' @param results character. If "FILTER", only results passing the FILTER threshold will be returned in an aggregated format.
#' @param ID character. Sample ID
#'
#' @return graphic. Depicts Base-Coverage and Deletion coverage accros the sequence. Significant change points are marked in red and deletions are highlighted in black.
#' @return tibble. Deletion table report.
#' \itemize{
#' \item CHROM: character. Chromosome name
#' \item START: integer. Window start position
#' \item END: integer. Window end position
#' \item WIDTH: integer. Window width
#' \item BASE.SUM: integer. Number of sequenced "ATCG"
#' \item BASE.COVERAGE: double. "ATCG" coverage
#' \item DEL.SUM: integer. Number of sequenced "-"
#' \item DEL.COVERAGE: double. "-" coverage
#' \item COVERAGE: double. Total coverage of the window ("ATCG", "+", "-")
#' \item CATEGORY: character. Currently only "DEL" category exists
#' \item RATIO: double. Ratio between DEL.COVERAGE and BASE.COVERAGE.
#' \item FILTER: character. COVERAGE (coverage is lower than minCoverage), SHORT (the segment is smaller than minseglen), OUT (Base / deletion ratio is bigger than 5), PASS (Base / deletion ratio is smaller than .8), AMBIGUIOUS (other cases to be manually checked).
#' \item ID: character. Sample ID
#' }
#'
#' @keywords PCIP
#'
#' @author Vincent Hahaut
#'
#' @import changepoint
#' @import dplyr
#' @import readr
#' @import purrr
#' @importFrom magrittr %>%
#' @import tibble
#' @import tidyr
#' @import ggplot2
#' @importFrom stats end na.omit sd start
#' @importFrom utils head write.table
#'
#' @export
PCIP_targetSV <- function(bam.path = NULL, viralLength = NULL, minseglen = 50, out = NULL, minCoverage = 5, ID = NULL, provirus = "all", results = "FILTER"){

  #suppressPackageStartupMessages(library(changepoint))
  #suppressPackageStartupMessages(library(tidyverse))

  if(is_empty(bam.path)){ stop('Empty or Absent bam.path argument!') }
  if(is_empty(viralLength)){ stop('Empty or Absent viralLength argument!') }
  if(is_empty(minseglen)){ stop('Empty or Absent minseglen argument!') }
  if(is_empty(out)){ stop('Empty or Absent out argument!') }
  if(is_empty(minCoverage)){ stop('Empty or Absent minCoverage argument!') }
  if(is_empty(provirus)){ stop('Empty or Absent provirus argument!') }
  if(is_empty(results)){ stop('Empty or Absent results argument!') }

  # 1. Load coverage information from BAM file using SAMTOOLS PILEUP
  coverage <- Rsamtools::pileup(bam.path, index = paste0(bam.path, ".bai"), pileupParam = Rsamtools::PileupParam(ignore_query_Ns = FALSE, include_insertions = TRUE))

  if(provirus != "all"){
    coverage <- coverage[coverage$seqnames == provirus,]
  }

  # 2. Tidy Coverage tibble
  # Get base-count DEL (-), INS (+) and BASE (A,T,C,G)
  coverage.tidy <- left_join(
    data_frame(pos = 1:viralLength),
    coverage %>%
      as.data.frame() %>%
      dplyr::filter(nucleotide %in% c("A", "T", "C", "G", "-", "+")) %>%
      mutate(category =
               case_when(
                 nucleotide == "+" ~ "INS",
                 nucleotide == "-" ~ "DEL",
                 TRUE ~ "BASE"
               )
      ) %>%
      group_by(pos, category) %>%
      summarise(
        count = sum(count)
      ),
    by = c("pos")
  ) %>%
    arrange(pos) %>%
    mutate(count = replace_na(count, 0),
           category = replace_na(category, "BASE")) %>%
    spread(category, count) %>%
    mutate(DEL = replace_na(DEL, 0),
           INS = replace_na(INS, 0),
           BASE = replace_na(BASE, 0)) %>%
    rowwise() %>%
    mutate(coverage = DEL + INS + BASE)

  # 3. Compute the change point function
  changePts.coverage.DEL <- cpt.meanvar(coverage.tidy$DEL, method = "PELT", minseglen = minseglen)

  # 4. Summarise results and create the report

  # 4.1. Get windows between changePoints

  # DELETIONS
  windows.DEL <- c(0, changePts.coverage.DEL@cpts) %>%
    as.data.frame() %>%
    dplyr::rename("start" = ".") %>%
    mutate(end = lead(start, default = as.character(viralLength)))

  coverage.windows <- apply(windows.DEL, 1, function(x)
    coverage.tidy[x["start"]:x["end"],] %>%
      mutate(CHROM = provirus,
             START = min(pos),
             END = max(pos),
             WIDTH = END-START,
             BASE.sum = sum(BASE),
             BASE.coverage = BASE.sum / WIDTH,
             DEL.sum = sum(DEL),
             DEL.coverage = DEL.sum / WIDTH,
             INS.sum = sum(INS),
             INS.coverage = INS.sum / WIDTH,
             COVERAGE = (BASE.sum + DEL.sum + INS.sum) / WIDTH,
             RATIO = DEL.coverage / COVERAGE,
             CATEGORY = "<DEL>") %>%
      select(CHROM, START, END, WIDTH, BASE.sum, BASE.coverage, DEL.sum, DEL.coverage, INS.sum, INS.coverage, COVERAGE, RATIO, CATEGORY) %>%
      distinct()
  ) %>%
    bind_rows() %>%
    filter(COVERAGE != "Inf") %>%
    mutate(FILTER =
             case_when(
               WIDTH < minseglen ~ "SHORT",
               COVERAGE < minCoverage ~ "COVERAGE",
               RATIO > 5 ~ "OUT",
               COVERAGE >= minCoverage & RATIO > .8 ~ "PASS",
               TRUE ~ "AMBIGUIOUS"
             ),
           ID = ID
    )

  # INSERTIONS
  coverage.windows.INS <- coverage.tidy %>%
    rowwise() %>%
    mutate(ratio.INS = INS / BASE) %>%
    ungroup() %>%
    filter(ratio.INS > 1.25 & INS > 40)

  coverage.windows.INS <- coverage.windows.INS %>%
    mutate(CHROM = provirus,
           START = pos,
           END = pos,
           WIDTH = 1,
           CATEGORY = "<INS>",
           BASE.sum = BASE,
           BASE.coverage = BASE / WIDTH,
           DEL.sum = DEL,
           DEL.coverage = DEL / WIDTH,
           INS.sum = INS,
           INS.coverage = INS / WIDTH,
           COVERAGE = coverage,
           RATIO = ratio.INS,
           FILTER = "PASS",
           ID = ID) %>%
    select(-pos, -BASE, -DEL, -INS, -coverage, -ratio.INS)


  # 4.4. Group similar windows:
  if(results == "FILTER"){

    # 4.4.1. Remove regions not passing the FILTER
    coverage.windows.filtered <- coverage.windows %>%
      filter(FILTER == "PASS")

    if(nrow(coverage.windows.filtered) > 0){
      # 4.4.2. Aggregate close windows given a 5bp window
      coverage.windows.filtered.gr <- coverage.windows.filtered %>%
        mutate(START = START ,
               END = END) %>%
        GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "CHROM", start.field = "START", end.field = "END", ignore.strand = TRUE)

      windows.gr <- GenomicRanges::reduce(coverage.windows.filtered.gr)
      overlaps <- GenomicRanges::findOverlaps(coverage.windows.filtered.gr, windows.gr, maxgap = 5)

      # 4.4.3. Regroup close windows
      report <- coverage.windows.filtered %>%
        mutate(windows = S4Vectors::subjectHits(overlaps)) %>%
        group_by(CATEGORY, windows) %>%
        summarise(
          CHROM = unique(CHROM),
          START = min(START),
          END = max(END),
          WIDTH = END - START,
          BASE.sum = sum(BASE.sum),
          BASE.coverage = BASE.sum/WIDTH,
          DEL.sum = sum(DEL.sum),
          DEL.coverage = DEL.sum / WIDTH,
          INS.sum = sum(INS.sum),
          INS.coverage = INS.sum / WIDTH,
          COVERAGE = (BASE.sum + DEL.sum + INS.sum) / WIDTH) %>%
        rowwise() %>%
        mutate(
          RATIO = DEL.coverage / COVERAGE,
          FILTER = "PASS",
          ID = ID
        ) %>%
        bind_rows(coverage.windows.INS) %>%
        select(CHROM, START, END, WIDTH, CATEGORY, BASE.sum, BASE.coverage, DEL.sum, DEL.coverage, INS.sum, INS.coverage, COVERAGE, RATIO, FILTER, ID)

    } else {

      report <- data_frame(CHROM = vector(),
                           START = vector(),
                           END = vector(),
                           WIDTH = vector(),
                           BASE.sum = vector(),
                           BASE.coverage = vector(),
                           DEL.sum = vector(),
                           DEL.coverage = vector(),
                           INS.sum = vector(),
                           INS.coverage = vector(),
                           COVERAGE = vector(),
                           CATEGORY = vector(),
                           FILTER = vector(),
                           ID = ID)
    }

    write.table(report, paste0(out, "-DEL_detection.report.txt"), sep = "\t", quote = F, row.names = F)

  } else {

    write.table(report, paste0(out, "-DEL_detection.report.txt"), sep = "\t", quote = F, row.names = F)

  }

  # 5. Plot the results
  coverage.plot <- ggplot(coverage.tidy) +
    geom_point(aes(x = pos, y = coverage.tidy$BASE, color = "Coverage"), alpha = 0.5) +
    scale_x_continuous(breaks = seq(1, viralLength, 1000)) +
    theme_minimal() +
    ylab("Coverage ATCG (# Reads)") +
    xlab("Viral Position (bp)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom",
          legend.title = element_blank()) +
    geom_vline(xintercept = c(filter(report, CATEGORY == "<DEL>")$START, filter(report, CATEGORY == "<DEL>")$END), color = "red") +
    geom_vline(xintercept = c(filter(report, CATEGORY == "<INS>")$START, filter(report, CATEGORY == "<INS>")$END), color = "blue") +
    geom_line(aes(x = pos, y = DEL, color = "Deletions"), alpha = 0.6, size = 1.2) +
    geom_line(aes(x = pos, y = INS, color = "Insertions"), alpha = 0.6, size = 1.2) +
    scale_color_brewer(palette = "Dark2")

  rectangle <- report %>% filter(FILTER == "PASS")
  if(nrow(rectangle) > 0){
    sommet <- max(coverage.tidy$BASE)
    coverage.plot <- coverage.plot +
      geom_rect(data = rectangle %>% filter(CATEGORY == "<DEL>"), aes(xmin = START, xmax = END, ymin = 0.985*sommet, ymax = sommet, fill = CATEGORY)) +
      scale_fill_manual(values = c("#974201"))
  }

  # 6. Save
  ggsave(plot = coverage.plot, filename = paste0(out, "-DEL_detection.coverage.pdf"), device = "pdf", width = 15,height = 6)

}


