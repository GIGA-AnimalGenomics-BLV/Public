#' @title Get Read Breakpoints and ShearSites
#'
#' @description
#' Extract TARGET-HOST junction (= Breakpoint) and DNA fragmentation site (= Shearsite) created by sonication .
#'
#' @param PAF tibble. Pairwise Alignment file (PAF) generated with minimap2, red with \code{readPairwiseAlignmentFile} and formated with \code{PCIP_filter}
#' @param targetName character. Name of the TARGET chromosome
#' @param lengthTarget numeric. Length of the TARGET chromosome in base pairs
#' @param gapAlignment numeric. Maximum mean gap tolerated in the alignment between substrings.
#'
#' @return Returns a tibble. Field description:
#' \itemize{
#' \item index (numeric) Line index
#' \item readID (character) ID of the read
#' \item seqnames.genome (character) HOST Chromosome name
#' \item edge5_breakPoint (double) Position of the TARGET 5' edge in the HOST genome
#' \item edge3_breakPoint (double) Position of the TARGET 3' edge in the HOST genome
#' \item shearSite.genome (double) ShearSite position in the HOST genome
#' \item breakPointEdge (character) TARGET edge type detected in this read (edge3, edge5, edge3-edge5 or edge5-edge3)
#' \item distanceEdge3_edge5 (double) Distance between the two TARGET edges if both are detected
#' \item strand.target (character) TARGET substring orientation
#' \item strand.genome (character) GENOME substring orientation
#' \item context (character) TARGET - HOST substrings orientations related to each other
#' \item ligation (character) TARGET - HOST structure of the read as binary code (0 = HOST and 1 = TARGET)
#' \item minDistanceEdge (double) Minimal distance to one TARGET edge
#' }
#'
#' @keywords PCIP
#'
#' @author Vincent Hahaut
#'
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
#' @note
#' minimap2PAF <- readPairwiseAlignmentFile(minimap2PAF="my/path/to/minimap2.align")
#' minimap2PAF.filter <- PCIP_filter(minimap2PAF = minimap2PAF, targetName="HTLV")
#' minimap2PAF.breakpoints <- PCIP_getBreakPoints(PAF = minimap2PAF.filter, lengthTarget = 9091)
#'
#' @export
PCIP_getBreakPoints <- function(PAF = NULL, lengthTarget = NULL, targetName = NULL, gapAlignment = NA, distanceToLTR = NA, returnFILTEREDout = FALSE){

  suppressPackageStartupMessages(library(tidyverse))

  if(is_empty(PAF)){ stop('Empty or Absent PAF argument!') }
  if(is_empty(lengthTarget)){ stop('Empty or Absent lengthTarget argument! Please provide the length of the TARGET chromosome') }
  if(is_empty(targetName)){ stop('Empty or Absent targetName argument! Please provide the TARGET chromosome name') }
  if(isFALSE(all(c("ligation", "mapGap") %in% colnames(PAF)))){ stop('PAF needs to be prepared with PCIP_filter() prior to PCIP_getBreakPoints()!') }

  # 0. OPTION: REMOVE READS WITH MEAN ALIGNMENT GAPS >= gapAlignment
  if(!is.na(gapAlignment)){PAF <- filter(PAF, meanGap < gapAlignment)}

  ###################
  #### 1. TARGET ####
  ###################
  if(nrow(PAF) > 100000){ print('This may take some time, ...') }

  print('1. Summarise TARGET substring information and extract target-genome LTR information')

  targets <- PAF %>%
    select(readID, seqnames, start, end, strand, ligation, readStart_position, readEnd_position) %>%
    # 1.1. Filter out GENOME substrings
    filter(seqnames == targetName) %>%
    # 1.2. Extract the distance to the 5' edge (distanceEdge5) and 3' edge (distanceEdge3)
    # Determine the closest edge (closestEdge) and its distance to the edge of the TARGET (distanceEdgeSub)
    mutate(distanceEdge5 = abs(0 - start),
           distanceEdge3 = abs(as.numeric(lengthTarget) - end)) %>%
    ungroup() %>%
    mutate(minDistanceEdgeSubstring = pmap_dbl(list(distanceEdge5, distanceEdge3), min)) %>%
    # 1.3. Summarise results at the read level
    # LTR5 or LTR3 based on the substring with the closest distance to one LTR
    group_by(readID) %>%
    # 1.3.1.
    mutate(
      closestEdge = ifelse(min(distanceEdge5) < min(distanceEdge3), "edge5", "edge3"),
      strand.target = paste(unique(strand), collapse = ","),
      minDistanceEdge = min(distanceEdge5, distanceEdge3),
      index = which(minDistanceEdgeSubstring == minDistanceEdge)
    ) %>%
    # 1.3.2.
    # To determine where is located the integration in the read after ligation, the position of the LTR-host junction in the read is reported
    # The viral-host junction can be on the left (i.e., 10) or right (i.e., 01) of the read, depending on its position regarding the host genome
    mutate(breakPtsPosition = case_when(
      ligation %in% c("101","1001") ~ ifelse(index == 1, "left", "right"),
      ligation == "01" ~ "right",
      ligation == "10" ~ "left")) %>%
    ungroup() %>%
    # 1.5. Finish
    mutate(strand.target = ifelse(nchar(strand.target) > 1, "*", strand.target))

  # OPTION: Remove reads with large mapping gaps in the target
  # These gaps can come from LTR deletions or mapping issues
  # This option can be used to explore these reads
  if(!is.na(distanceToLTR)){
    print('CAREFUL: distanceToLTR will remove reads with large gaps in the LTRs')
    targets <- filter(targets, minDistanceEdge <= distanceToLTR)

    # OPTION: returnFILTEREDout == TRUE directly returns these reads with mapping gaps
    if(returnFILTEREDout == TRUE){
      return(filter(targets, minDistanceEdge > distanceToLTR))
      stop("returnFILTEREDout == TRUE, return ONLY the reads with large alignment gaps in the provirus")
    }
  }

  # 1001 ligation could be counted twice. Only keep those that have the smallest minDistanceEdge
  if(any(targets$ligation == 1001)){
    index1001 <- targets %>%
    rownames_to_column("rowN") %>%
    filter(ligation == 1001) %>%
    group_by(readID) %>%
    arrange(minDistanceEdgeSubstring) %>%
    dplyr::slice(2)

    targets <- targets[-as.numeric(index1001$rowN),]
  }

  print('2. Summarise GENOME substring information')

  ###################
  #### 2. GENOME ####
  ###################

  genome_concat <- PAF %>%
    select(readID, seqnames, start, end, strand, ligation) %>%
    # 2.1. Filter out TARGET substrings
    filter(seqnames != targetName) %>%
    dplyr::rename("seqnames.genome" = seqnames,
                  "start.genome" = start,
                  "end.genome" = end,
                  "strand.genome" = strand) %>%
    # 2.2. Merge with TARGET informations
    left_join(targets %>%
                group_by(readID) %>%
                select(readID, closestEdge, strand.target, minDistanceEdge, breakPtsPosition) %>%
                distinct(),
              by = c("readID", "readID")
    ) %>%
    # 2.3. Concatenate strand informations
    group_by(readID) %>%
    mutate(strand.genome = paste(unique(strand.genome), collapse = ","),
           strand.genome = if_else(strand.genome %in% c("+", "-"), strand.genome, "*"), # Sometimes minimap2 attribute different strands to the same read: considered as 'unknown'
           context = if_else(strand.target == strand.genome, "concordant", "discordant")) %>%
    # 2.4. Remove ambiguous mapping strands
    filter(strand.genome != "*")

    ####### 2.5. Clean-up
    # 2.5.1. 1001: To avoid inflating counts only keep the 10 or 01 that contains the best integration site
    if(any(genome_concat$ligation == 1001)){
      genome.1001 <- genome_concat %>%
        rownames_to_column("rowN") %>%
        filter(ligation == "1001") %>%
        group_by(readID) %>%
        mutate(index = if_else(breakPtsPosition == "left", rowN[2], rowN[1])) %>%
        ungroup() %>%
        select(index) %>%
        distinct()

      genome_concat <- genome_concat[-as.numeric(genome.1001$index),]
    }

    # 2.5.2. In rare cases substrings have the same distance to LTR. Remove them
    genome_concat <- genome_concat %>%
      group_by(readID) %>%
      filter(n() == 1) %>%
      ungroup()

  # 3. Extract TARGET-HOST ShearSites and Breakpoints
  print('3. Get the Integration site position')

  if(any(genome_concat$closestEdge %in% c("edge3", "edge5"))){

    # Depending on the LTR position in the read and the genome strand, select the right viral-host breakpoints and shear site
    edges <- mutate(genome_concat, category = paste0(strand.genome, ":", breakPtsPosition)) %>% split(.$category)

    if(any(names(edges) %in% "-:left")){
      edges.neg.Left <- mutate(edges[["-:left"]], integrationSite = end.genome, shearSite.genome = start.genome)
    } else{
      edges.neg.Left <- data_frame()
    }

    if(any(names(edges) %in% "+:left")){
      edges.pos.Left <- mutate(edges[["+:left"]], integrationSite = start.genome, shearSite.genome = end.genome)
    } else{
      edges.pos.Left <- data_frame()
    }

    if(any(names(edges) %in% "-:right")){
      edges.neg.right <- mutate(edges[["-:right"]], integrationSite = start.genome, shearSite.genome = end.genome)
    } else{
      edges.neg.right <- data_frame()
    }

    if(any(names(edges) %in% "+:right")){
      edges.pos.right <- mutate(edges[["+:right"]], integrationSite = end.genome, shearSite.genome = start.genome)
    } else{
      edges.pos.right <- data_frame()
    }

  } else {
    stop("No viral-host junctions found!")
  }

  print('4. TIDY')

  read_breakPoints <- bind_rows(edges.neg.Left, edges.pos.Left, edges.neg.right, edges.pos.right) %>%
    select(readID, seqnames.genome, shearSite.genome, closestEdge, integrationSite, strand.target, strand.genome, context, ligation, minDistanceEdge) %>%
    spread(key = closestEdge, value = integrationSite) %>%
    # Can't find the exact shearSites in 10 and 01, put 9999999999 to shearSite to count them all as the same when looking at the abundance. Avoid inflating the counts.
    mutate(shearSite.genome = ifelse(ligation %in% c("10", "01"), 9999999999, shearSite.genome))

  if(any(colnames(read_breakPoints) %in% "edge3")){
    read_breakPoints <- dplyr::rename(read_breakPoints, "edge3_breakPoint" = edge3)
  } else {
    read_breakPoints <- mutate(read_breakPoints, edge3_breakPoint = NA)
  }

  if(any(colnames(read_breakPoints) %in% "edge5")){
    read_breakPoints <- dplyr::rename(read_breakPoints, "edge5_breakPoint" = edge5)
  } else {
    read_breakPoints <- mutate(read_breakPoints, edge5_breakPoint = NA)
  }

  print(paste0('Number of reads with breakpoints detected: ',
               nrow(read_breakPoints),
               " out of ",
               length(unique(PAF$readID)),
               " reads ",
               '(',
               100*round(nrow(read_breakPoints) / length(unique(PAF$readID)), 5),
               '%)'
    )
  )

  return(read_breakPoints)

  print('5. Done! ')
}

# OLD:
#edges <- genome_concat %>%
#  group_by(readID) %>%
#  mutate(integrationSite =
#           case_when(
#             strand.genome == "+" & breakPtsPosition == "left" ~ min(start.genome),
#             strand.genome == "+" & breakPtsPosition == "right" ~ max(end.genome),
#             strand.genome == "-" & breakPtsPosition == "left" ~ max(end.genome),
#             strand.genome == "-" & breakPtsPosition == "right" ~ min(start.genome)
#           )
#  ) %>%
#  mutate(shearSite.genome =
#           case_when(
#             strand.genome == "+" & breakPtsPosition == "left" ~ max(end.genome),
#             strand.genome == "+" & breakPtsPosition == "right" ~ min(start.genome),
#             strand.genome == "-" & breakPtsPosition == "left" ~ min(start.genome),
#             strand.genome == "-" & breakPtsPosition == "right" ~ max(end.genome)
#           )
#  ) %>%
#  ungroup()

