#' @title splitFASTQ
#'
#' @description
#'
#'
#' @param fastq.path character. SAM File path
#' @param mode character. "FULL" (just split the reads by associated IS) or "MAPPED" (split the reads and only keep the mapped parts)
#'
#' @return a tibble containing the .align results
#'
#' @keywords PCIP
#'
#' @author Vincent Hahaut
#'
#' @examples
#' readSAMFile(samFile="my/path/to/samFile.sam")
#'
#' @export
splitFASTQ <- function(breakpoints.path = NULL, readID.cleaned.path = NULL, mode = "breakpoints", PAF.path = NULL, maxgap = 50, out.prefix = NULL){

  suppressPackageStartupMessages(library(ShortRead))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(stringr))

  # 1. Load data:

  print("1. Load Files")

  # 1.1. Cleaned readIDs ( from PCIP_getBreakPoints() )
  readID.cleaned <- read.delim(readID.cleaned.path)

  # 1.2. .align file
  PAF <- readPairwiseAlignmentFile(PAF.path)

  if(mode == "breakpoints"){

    print(paste0("MODE = ", mode))

    # 1.3. Breakpoints
    breakpoints.gr <- read.delim(breakpoints.path) %>%
      rowwise() %>%
      mutate(seqnames = seqnames,
             start = min(edge5.pos, edge3.pos, na.rm = T),
             end = max(edge5.pos, edge3.pos, na.rm = T)) %>%
      select(seqnames, start, end) %>%
      makeGRangesFromDataFrame(keep.extra.columns = T)

  } else if( mode == "windows"){

    print(paste0("MODE = ", mode))

    # 1.3. Breakpoints
    breakpoints.gr <- read.delim(breakpoints.path) %>%
      makeGRangesFromDataFrame(keep.extra.columns = T, ignore.strand = T)

  }

  print("2. Select readIDs from PAF used to call IS")

  # 2. Select readID from .align present in the PCIP_getBreakPoint() function
  PAF.gr <- PAF %>%
    dplyr::filter(qName %in% readID.cleaned$readID) %>%
    select(tName, tStart, tEnd, qName) %>%
    makeGRangesFromDataFrame(
      seqnames.field = "tName",
      start.field = "tStart",
      end.field = "tEnd",
      keep.extra.columns = T,
      ignore.strand = T
    )

  print("3. Associate IS to readIDs from FASTQ")

  # 3. Get readIDs spanning each integration sites
  readIDs <- findOverlaps(PAF.gr, breakpoints.gr, maxgap = maxgap) %>%
    as_data_frame() %>%
    mutate(ID = paste0(breakpoints.gr[subjectHits]),
           readID = PAF.gr$qName[queryHits]) %>%
    select(ID, readID) %>%
    distinct() %>%
    left_join(
      data_frame(ID = paste0(breakpoints.gr)),
      by = c("ID" = "ID")
    ) %>%
    mutate(readID = paste0("@", readID)) %>%
    select(readID, ID)

  if( nchar(out.prefix) > 0 ){
    write.table(readIDs, paste0(out.prefix, "-splitFASTQ_table.txt"), sep = "\t", row.names = F, quote = F)
  }

  return(readIDs)

}




