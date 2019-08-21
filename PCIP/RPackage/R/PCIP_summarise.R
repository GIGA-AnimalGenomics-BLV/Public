#' @title Collapse results into an Integration Site table
#'
#' @description
#' Group close viral-host junctions into integration sites.
#' Compute the abundance of each clone and before (raw) or after (filtered) removing PCR duplicates.

#' @param PCIPbreakpoints tibble. TARGET-HOST breakpoints from \code{PCIP_getBreakPoints}
#' @param mergeISdistance numeric. Window (up/down) to regroup close viral-host junctions into an integration site.
#' @param mergeShearSiteDistance numeric. Window (up/down) to collapse reads displaying the same shear site.
#' @param minDistanceLTR numeric. 1001, 01 and 10 reads with large gaps in the LTRs often introduce false positive. Filter those with a distance to the nearest viral edge > minDistanceLTR.
#' @param distanceFLAG numeric. Flag integration sites located too close (within distanceFLAG, up/down).
#'
#' @return list Field description:
#' \enumerate{
#' \item Integration Site Table [[1]].
#' \itemize{
#' \item seqnames (character). Chromosome name.
#' \item strand (character). TARGET orientation related to the HOST genome.
#' \item edge5.pos (double). 5' edge position into HOST genome.
#' \item edge5.raw (double). Number of reads supporting the 5' edge (ALL).
#' \item edge5.count (double). Number of reads supporting the 5' edge (PCR duplicates removed).
#' \item edge3.pos (double). 3' edge position into HOST genome.
#' \item edge3.raw (double). Number of reads supporting the 3' edge (ALL).
#' \item edge3.count (double). Number of reads supporting the 3' edge (PCR duplicates removed).
#' \item ID (character). TARGET location into the HOST genome formated as chr:start-end.
#' \item count.max (double). Maximum number of reads supporting that IS (max(edge5.count, edge3.pos, na.rm = T))
#' \item edge5.meanDistLTR (double). mean distance to the 5' edge in the reads supporting that IS.
#' \item edge5.sdDistLTR (double). sd distance to the 5' edge in the reads supporting that IS.
#' \item edge3.meanDistLTR (double). mean distance to the 3' edge in the reads supporting that IS.
#' \item edge3.sdDistLTR (double). sd distance to the 3' edge in the reads supporting that IS.
#' \item FLAG.DIST (logical). Flag IS with mean distance to one of the two LTR > 100 bp.
#' \item FLAG.READS (logical). Flag IS with less than 10 non-PCR duplicats reads.
#' \item FLAG.OWNER (logical). Flag IS within a distanceFLAG distance of another IS.
#' \item OWNER (numeric). Extremely close IS are often the results of false positive IS due to the read processing. For every IS within a distanceFLAG sliding window, the real IS is the one supported by the most reads.
#'
#' }
#' \item Read description table [[2]].
#' \itemize{
#' \item readID (character). ID of the read.
#' \item ID (character). TARGET location into the HOST genome formated as chr:start-end.
#' \item strand (character). TARGET orientation related to the HOST genome.
#' \item maxIS (character). Position (5' edge or 3' edge) supported by the most reads.
#' }
#' }
#'
#' @keywords PCIP
#'
#' @author Vincent Hahaut
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
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

PCIP_summarise <- function(PCIPbreakpoints = NULL, mergeISdistance = 200, mergeShearSiteDistance = 0, minDistanceLTR = 300, distanceFLAG = 1000){

  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(GenomicRanges))

  if(is_empty(PCIPbreakpoints)){ stop('Empty or Absent PCIPbreakpoints argument!') }
  if(is_empty(mergeISdistance)){ stop('Empty or Absent mergeISdistance argument!') }
  if(isFALSE(all(c("edge5_breakPoint", "edge3_breakPoint", "shearSite.genome") %in% colnames(PCIPbreakpoints)))){ stop('PAF needs to be prepared with PCIP_getBreakPoints() prior to PCIP_summarise()!') }

  breaks <- mutate(PCIPbreakpoints, index = 1:nrow(PCIPbreakpoints)) %>%
    # 01 and 10 are containing large mapping gaps introduce false positive calls. Filter them out.
    filter(ligation %in% c("10", "01", "1001") & minDistanceEdge < minDistanceLTR | ligation %in% c("101"))

  ##############################################
  ## 1. GROUPS BREAKPOINTS/SHEARSITES INTO IS ##
  ##############################################
  print(paste0("1. Group reads lying with a +/- ", mergeISdistance, " bp window."))

  #### ==> IS ####
  # 1.1. Transform the data to GRanges:
  breaks.gr <- breaks %>%
    mutate(start = pmin(.$edge5_breakPoint, .$edge3_breakPoint, na.rm = T),
           end = start) %>%
    select(index, seqnames.genome, start, end) %>%
    makeGRangesFromDataFrame(seqnames.field = "seqnames.genome", ignore.strand = T, keep.extra.columns = T)
  # 1.2. Group close viral-host junctions into integration sites (= groupIS)
  breaks.win <- (breaks.gr + mergeISdistance) %>%
    GenomicRanges::reduce() %>%
    findOverlaps(., breaks.gr) %>%
    as_tibble() %>%
    mutate(index = dplyr::slice(breaks, subjectHits)$index) %>%
    select(index, queryHits) %>%
    dplyr::rename(groupIS = "queryHits")

  #### ==> SHEARSITE ####

  # 2. Regroup close shearsites:

  # 2.1. If a minimal distance between close shear sites is required:
  if(mergeShearSiteDistance > 0){
    print(paste0("mergeShearSiteDistance == ", mergeShearSiteDistance))
    # 2.1.1. Transform the data to GRanges:
    shearSite.gr <- breaks %>%
      mutate(shearSite.genome = ifelse(is.na(shearSite.genome), 0, shearSite.genome), %>%
               start = shearSite.genome,
             end = shearSite.genome
             ) %>%
      select(index, seqnames.genome, start, end) %>%
      makeGRangesFromDataFrame(seqnames.field = "seqnames.genome", ignore.strand = T, keep.extra.columns = T)

    # 2.1.2. Regroup close shear sites given the mergeShearSiteDistance window.
    breaks.shearSite <- (shearSite.gr + mergeShearSiteDistance) %>%
      GenomicRanges::reduce() %>%
      findOverlaps(shearSite.gr) %>%
      as_tibble() %>%
      mutate(index = dplyr::slice(breaks, subjectHits)$index) %>%
      select(index, queryHits) %>%
      dplyr::rename(groupShear = "queryHits")

  } else {
    print("mergeShearSiteDistance == 0")

    # 2.2. If no window is required between shear sites:
    breaks.shearSite <- tibble(
      index = breaks$index,
      groupShear = breaks$shearSite.genome
    ) %>%
      mutate(groupShear = ifelse(is.na(groupShear), 0, groupShear))
  }

  #### ==> REGROUP ####
  # 2.3. Combine informations
  breaks <- breaks %>%
    left_join(breaks.win, by = c("index" = "index")) %>%
    left_join(breaks.shearSite, by = c("index" = "index"))

  ######################################
  ## TIDY AND COMPUTE CLONE ABUNDANCE ##
  ######################################
  print("2. Collapse results into an integration site table.")

  IS.table <- breaks %>%
    group_by(groupIS) %>%
    # 3.1. VIRAL orientation is determined by the orientation of the HOST substring related to the TARGET substring.
    # If enough reads are concordant (+/+), TARGET is oriented "+". If most are discordant (+/- or -/+), TARGET is oriented "-".
    # Rare cases of display reads with no orientation concensus: Flag them as "*" (= Ambigous).
    mutate(
      strand =
        case_when(
          sum(context == "concordant")/n() > 0.8 ~ "+",
          sum(context == "concordant")/n() < 0.2 ~ "-",
          TRUE ~ "*"
        )
    ) %>%
    select(index, seqnames.genome, strand, edge5_breakPoint, edge3_breakPoint, groupShear, groupIS, ligation, minDistanceEdge) %>%
    gather(key = "edge", value = "position", -seqnames.genome, -strand, -groupShear, -groupIS, -index, -ligation, -minDistanceEdge)

  print('3. Compute the clone abundance.')
  # Raw abundance = don't take into account shearsite.
  # count = take the shear site into account.
  abundance <- IS.table %>%
    select(-ligation, -minDistanceEdge) %>%
    group_by(edge, groupIS) %>%
    mutate(count.raw = length(groupIS[!is.na(position)]),
           count = length(unique(groupShear[!is.na(position)]))) %>%
    ungroup() %>%
    select(edge, groupIS, count.raw, count) %>%
    distinct()

  print('4. Identify the integration sites.')
  # Nanopore data are noisy. The integration site is defined as the position supported by the most non PCR-duplicated reads.
  position <- IS.table %>%
    select(-ligation, -minDistanceEdge) %>%
    group_by(edge, groupIS, position) %>%
    # For each IS position of a groupIS, get the number of non-duplicated reads
    mutate(count = length(unique(groupShear))) %>%
    select(edge, groupIS, seqnames.genome, strand, position, count) %>%
    ungroup() %>%
    filter(!is.na(position)) %>%
    group_by(edge, groupIS) %>%
    summarise(
      seqnames = unique(seqnames.genome),
      # In case of equality: select the first.
      position = unique(position[count == max(count, na.rm = T)])[1],
      strand = unique(strand)
    ) %>%
    ungroup()

  # 4. For every IS (5' and 3'), compute the mean and sd of the distance to the closest LTR. Can be used to filter out bad IS.
  meanDistanceLTR <- IS.table %>%
    select(-ligation) %>%
    filter(!is.na(position)) %>%
    group_by(edge, groupIS, position, groupShear) %>%
    summarise(minDistanceLTR = mean(minDistanceEdge)) %>%
    ungroup() %>%
    group_by(edge, groupIS) %>%
    summarise(meanDistanceLTR = mean(minDistanceLTR),
              sdDistanceLTR = sd(minDistanceLTR)) %>%
    ungroup()

  # 5. Regroup results.
  IS.tmp <- left_join(
    left_join(position,
              abundance,
              by = c("edge", "groupIS")),
    meanDistanceLTR,
    by = c("edge", "groupIS"))

  #######################
  # REARRANGE DATAFRAME #
  #######################

  IS.split <- split(IS.tmp, IS.tmp$edge)

  # ==> In case no IS
  if(nrow(IS.table) > 0){

    # ==> In case no 5'LTR
    LTR5.position <- if(is.null(IS.split[["edge3_breakPoint"]])){
      tibble(
        groupIS = IS.split[["edge3_breakPoint"]]$groupIS,
        edge5.pos = rep(NA, nrow(IS.split[["edge3_breakPoint"]])),
        seqnames = IS.split[["edge3_breakPoint"]]$seqnames,
        strand = IS.split[["edge3_breakPoint"]]$strand,
        edge5.raw = rep(NA, nrow(IS.split[["edge3_breakPoint"]])),
        edge5.count = rep(NA, nrow(IS.split[["edge3_breakPoint"]])),
        edge5.meanDistLTR = rep(NA, nrow(IS.split[["edge3_breakPoint"]])),
        edge5.sdDistLTR = rep(NA, nrow(IS.split[["edge3_breakPoint"]]))
      )
    } else {
      IS.split[["edge5_breakPoint"]] %>%
        dplyr::rename("edge5.pos" = position,
                      "edge5.raw" = count.raw,
                      "edge5.count" = count,
                      "edge5.meanDistLTR" = meanDistanceLTR,
                      "edge5.sdDistLTR" = sdDistanceLTR) %>%
        mutate(edge5.meanDistLTR = round(edge5.meanDistLTR, 2),
               edge5.sdDistLTR = round(edge5.sdDistLTR,2)) %>%
        select(-edge)
    }

    # ==> In case no 3'LTR
    LTR3.position <- if(is.null(IS.split[["edge3_breakPoint"]])){
      tibble(
        groupIS = IS.split[["edge5_breakPoint"]]$groupIS,
        edge3.pos = rep(NA, nrow(IS.split[["edge5_breakPoint"]])),
        seqnames = IS.split[["edge5_breakPoint"]]$seqnames,
        strand = IS.split[["edge5_breakPoint"]]$strand,
        edge3.raw = rep(NA, nrow(IS.split[["edge5_breakPoint"]])),
        edge3.count = rep(NA, nrow(IS.split[["edge5_breakPoint"]])),
        edge3.meanDistLTR = rep(NA, nrow(IS.split[["edge3_breakPoint"]])),
        edge3.sdDistLTR = rep(NA, nrow(IS.split[["edge3_breakPoint"]]))
      )
    } else {
      IS.split[["edge3_breakPoint"]] %>%
        dplyr::rename("edge3.pos" = position,
                      "edge3.raw" = count.raw,
                      "edge3.count" = count,
                      "edge3.meanDistLTR" = meanDistanceLTR,
                      "edge3.sdDistLTR" = sdDistanceLTR) %>%
        mutate(edge3.meanDistLTR = round(edge3.meanDistLTR, 2),
               edge3.sdDistLTR = round(edge3.sdDistLTR,2)) %>%
      select(-edge)
    }

  }

  breaks.summerised <- full_join(LTR5.position, LTR3.position, by = c("groupIS", "seqnames", "strand")) %>%
    # CORRECT 1-base system
    mutate(edge3.pos = ifelse(is.na(edge3.pos), NA, edge3.pos+1)) %>%
    rowwise() %>%
    # Extract the maximal count between 5' and 3' LTR
    mutate(count.max = max(edge5.count , edge3.count, na.rm = T),
           ID = paste0(seqnames, ":", min(edge3.pos, edge5.pos, na.rm = T),"-", max(edge3.pos, edge5.pos, na.rm = T))) %>%
    ungroup() %>%
    select(seqnames, strand, edge5.pos, edge5.raw, edge5.count, edge3.pos, edge3.raw, edge3.count, ID, count.max, edge5.meanDistLTR, edge5.sdDistLTR, edge3.meanDistLTR, edge3.sdDistLTR,  groupIS) %>%
    arrange(desc(count.max))

  print('5. Add FLAGS')

  # ADD FLAGS:
  # FLAG.DIST = Is this integration site supported by reads with large mapping gaps ?
  # FLAG.READS = Is this integration site supported by less than 10 non PCR-duplicats reads ?
  breaks.summerised <- mutate(breaks.summerised,
      FLAG.DIST = ifelse(pmin(.$edge5.meanDistLTR, .$edge3.meanDistLTR, na.rm = T)  > 100, F, T),
      FLAG.READS = ifelse(count.max < 10, F, T)
    )

  # FLAG.OWNER:
  # The noise introduce by nanopore sequencing creates dozens exception which cannot all be removed by filtering reads without affecting the final results.
  # False positive integration sites from this pipeline often originate from the inversion between what is considered as the integration site and what is the shear site.
  # This creates IS supported by few reads (i.e., < 5) near large real IS.
  # Removing them without affecting potential close real IS is difficult.
  # We flag them in the final result.
  # Every IS within a distanceFLAG bp (up/down) are regrouped. The IS supported by the most reads is flagged as FLAG.OWNER = TRUE and the others are FLAG.OWNER = FALSE.
  # A OWNER column is created. For every group of close IS, it contains the row number of the IS supported by the highest number of reads.

  # Add a distanceFLAG (up/down) window
  breaks.gr <- GRanges(breaks.summerised$ID) + distanceFLAG
  breaks.gr$count.max <- breaks.summerised$count.max

  # Regroup close IS
  reduced.win <- GenomicRanges::reduce(breaks.gr)

  # Find IS falling in the reduced windows
  overlaps <- findOverlaps(reduced.win, breaks.gr) %>%
    as_tibble() %>%
    group_by(queryHits) %>%
    mutate(maxCount = max(breaks.gr$count.max[subjectHits]),
           FLAG.OWNER = unique(subjectHits[breaks.gr$count.max[subjectHits] == maxCount])[1]) %>% # In case of equality, pick first one
    ungroup() %>%
    select(subjectHits, FLAG.OWNER) %>%
    mutate(subjectHits = as.character(subjectHits)) %>%
    distinct()

  breaks.summerised <- left_join(rownames_to_column(breaks.summerised), overlaps, by = c("rowname" = "subjectHits")) %>%
    mutate(OWNER = FLAG.OWNER,
           FLAG.OWNER = ifelse(FLAG.OWNER == rowname, TRUE, FALSE)) %>%
    select(-rowname)

  #######################
  ## READID - IS TABLE ##
  #######################
  print("6. Extract the read IDs supporting each IS")

  # 5. Return a tibble detailling the read IDs supporting each IS with more than 10 non PCR-duplicats reads.
  fastqID <- select(breaks, readID, groupIS) %>%
    left_join(
      breaks.summerised %>%
        rownames_to_column("index") %>%
        # SELECT ONLY IS WITH 10 non-duplicated reads or for which ownership is certain
        filter(FLAG.READS == TRUE | OWNER == index && FLAG.READS == FALSE) %>%
        select(-index) %>%
        rowwise() %>%
        mutate(maxIS = paste0(seqnames, ":",
                              ifelse(
                                if(is.na(edge5.count)) 0 else edge5.count > if(is.na(edge3.count)) 0 else edge3.count,
                              edge5.pos, edge3.pos))) %>%
        select(groupIS, ID, strand, maxIS),
      by = c("groupIS")) %>%
    dplyr::select(-groupIS) %>%
    filter(!is.na(ID))

  breaks.summerised <- select(breaks.summerised, -groupIS)

  return(list(breaks.summerised, fastqID))

}


