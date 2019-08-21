#' @title Read Pairwise Alignment File (PAF) From Minimap2
#'
#' @description Read the first 18th columns from pairwise alignment (PAF) produced by minimap2.
#' @description For more information about the PAF format refers to: \href{https://github.com/lh3/miniasm/blob/master/PAF.md}{Minimap2 PAF format}
#'
#' @param alignFile Path to minimap2 PAF file
#'
#' @keywords PCIP
#'
#' @author Vincent Hahaut
#'
#' @import readr
#' @importFrom stats end na.omit sd start
#' @importFrom utils head write.table
#'
#' @export
readPairwiseAlignmentFile <- function(alignFile = NULL){

  library(tidyverse)

  # 1. LOAD FILE
  print(paste0("Read Pairwise mApping File located at: ", alignFile))

  alignMinimap2 <- suppressWarnings(read_delim(alignFile, delim="\t", col_names = c("qName",
                                                                   "qLength",
                                                                   "qStart",
                                                                   "qEnd",
                                                                   "orientationToRef",
                                                                   "tName",
                                                                   "tLength",
                                                                   "tStart",
                                                                   "tEnd",
                                                                   "matchingBases",
                                                                   "matchingBasesIncludingGaps",
                                                                   "MAPQ",
                                                                   "numberMismatchGaps",
                                                                   "DPscoreMaxScoring",
                                                                   "DPalignmentScore",
                                                                   "numberAmbiguousBases",
                                                                   "alignmentType",
                                                                   "minimizers",
                                                                   "X19", "X20","X21", "X22", "X23"),
                              col_types =
                                cols_only(
                                  qName = col_character(),
                                  qLength = col_integer(),
                                  qStart = col_integer(),
                                  qEnd = col_integer(),
                                  orientationToRef = col_character(),
                                  tName = col_character(),
                                  tLength = col_integer(),
                                  tStart = col_integer(),
                                  tEnd = col_integer(),
                                  matchingBases = col_integer(),
                                  matchingBasesIncludingGaps = col_integer(),
                                  MAPQ = col_integer(),
                                  numberMismatchGaps = col_character(),
                                  DPscoreMaxScoring = col_character(),
                                  DPalignmentScore = col_character(),
                                  numberAmbiguousBases = col_character(),
                                  alignmentType = col_character(),
                                  minimizers = col_character()
                                )
    )
  )

  print(paste0("Number of reads loaded : ", length(unique(alignMinimap2$qName))))

  return(alignMinimap2)

}
