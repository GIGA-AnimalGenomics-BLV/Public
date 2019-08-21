#' @title Collapse results into Integration Site table
#'
#' @description
#' Group and count reads supporting the same breakpoint. Generate an Integration Sites (IS) table.
#'
#' @param PCIPbreakpoints tibble. TARGET-HOST breakpoints from \code{PCIP_getBreakPoints}
#' @param mergeISdistance numerical. Maximal distance (Upstream/Downstream) to group reads supporting the same IS.
#' @param mergeShearSiteDistance numerical. Maximal distance (Upstream/Downstream) to group reads with the same shear site.
#' @param filterIS numerical. Only report readID for IS supported by x reads.
#'
#' @return list Field description:
#' \enumerate{
#' \item Integration Site Table [[1]].
#' \itemize{
#' \item seqnames (character) Chromosome name
#' \item strand (character) TARGET orientation related to the HOST genome
#' \item edge5.pos (double) 5' edge position into HOST genome
#' \item edge5.count (double) Number of reads supporting the 5' edge
#' \item edge3.pos (double) 3' edge position into HOST genome
#' \item edge3.count (double) Number of reads supporting the 3' edge
#' \item ID (character) TARGET location into the HOST genome formated as chr:start-end
#' \item count.max (double) Maximum number of reads supporting that IS (max(edge5.count, edge3.pos, na.rm = T))
#' }
#' \item Read description table [[2]].
#' \itemize{
#' \item readID (character) ID of the read
#' \item ID (character) TARGET location into the HOST genome formated as chr:start-end
#' \item strand (character) TARGET orientation related to the HOST genome
#' \item max.pos (character) Position (5' edge or 3' edge) supported by the most reads and use as anchor to look for variants
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
#' @note
#' minimap2PAF <- readPairwiseAlignmentFile(minimap2PAF="my/path/to/minimap2.align")
#' minimap2PAF.filter <- PCIP_filter(minimap2PAF = minimap2PAF, targetName="HTLV")
#' minimap2PAF.breakpoints <- PCIP_getBreakPoints(PAF = minimap2PAF.filter, lengthTarget = 9091, targetName="HTLV", gapToEdge = 200)
#' myResults <- PCIP_summarise(PCIPbreakpoints = minimap2PAF.breakpoints, mergeISdistance = 200)
#' integrationSite.table <- myResults[[1]]
#' splitFASTQ.table <- myResults[[2]]
#'
#' @export

PCIP_summarise <- function(PCIPbreakpoints = NULL, mergeISdistance = 200, mergeShearSiteDistance = 0, minDistanceLTR = 300, distanceFLAG = 1000){

  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(GenomicRanges))

  if(is_empty(PCIPbreakpoints)){ stop('Empty or Absent PCIPbreakpoints argument!') }
  if(is_empty(mergeISdistance)){ stop('Empty or Absent mergeISdistance argument!') }
  if(isFALSE(all(c("edge5_breakPoint", "edge3_breakPoint", "shearSite.genome") %in% colnames(PCIPbreakpoints)))){ stop('PAF needs to be prepared with PCIP_getBreakPoints() prior to PCIP_summarise()!') }

  breaks <- mutate(PCIPbreakpoints, index = 1:nrow(PCIPbreakpoints)) %>%
    # 01 and 10 are often problematic when then contain large mapping gaps. They create small false positive Integration sites. Remove them
    filter(ligation %in% c("10", "01", "1001") & minDistanceEdge < minDistanceLTR | ligation %in% c("101"))

  ##############################################
  ## 1. GROUPS BREAKPOINTS/SHEARSITES INTO IS ##
  ##############################################
  print(paste0("1. Group reads lying with a +/- ", mergeISdistance, " bp window"))

  #### ==> IS ####
  # 2.1. Prepare Data (toGRANGE)
  breaks.gr <- breaks %>%
    mutate(start = pmin(.$edge5_breakPoint, .$edge3_breakPoint, na.rm = T),
           end = start) %>%
    select(index, seqnames.genome, start, end) %>%
    makeGRangesFromDataFrame(seqnames.field = "seqnames.genome", ignore.strand = T, keep.extra.columns = T)
  # 2.2. Reduce close breakpoints grouping them into Integration Sites
  breaks.win <- (breaks.gr + mergeISdistance) %>%
    GenomicRanges::reduce() %>%
    findOverlaps(., breaks.gr) %>%
    as_tibble() %>%
    mutate(index = dplyr::slice(breaks, subjectHits)$index) %>%
    select(index, queryHits) %>%
    dplyr::rename(groupIS = "queryHits")

  #### ==> SHEARSITE ####

  if(mergeShearSiteDistance > 0){
    print(paste0("mergeShearSiteDistance == ", mergeShearSiteDistance))
    # 2.3. Prepare Data
    shearSite.gr <- breaks %>%
      # REPLACE AMBIGUOUS SHEARSITES BY 0
      # ==> Collapsing PCR duplicates they will be considered as the same and not inflate the IS abundance
      mutate(shearSite.genome = ifelse(is.na(shearSite.genome), 0, shearSite.genome)) %>%
      mutate(
        start = shearSite.genome,
        end = shearSite.genome
      ) %>%
      select(index, seqnames.genome, start, end) %>%
      makeGRangesFromDataFrame(seqnames.field = "seqnames.genome", ignore.strand = T, keep.extra.columns = T)

    # 2.4. Add information about groups of shearsites
    breaks.shearSite <- (shearSite.gr + mergeShearSiteDistance) %>%
      GenomicRanges::reduce() %>%
      findOverlaps(shearSite.gr) %>%
      as_tibble() %>%
      mutate(index = dplyr::slice(breaks, subjectHits)$index) %>%
      select(index, queryHits) %>%
      dplyr::rename(groupShear = "queryHits")

  } else {
    print("mergeShearSiteDistance == 0")
    print("INFO: You can adjust mergeShearSiteDistance > 0")

    breaks.shearSite <- tibble(
      index = breaks$index,
      groupShear = breaks$shearSite.genome
    ) %>%
      mutate(groupShear = ifelse(is.na(groupShear), 0, groupShear))
  }

  #### ==> REGROUP ####
  # 2.5. Combine informations
  breaks <- breaks %>%
    left_join(breaks.win, by = c("index" = "index")) %>%
    left_join(breaks.shearSite, by = c("index" = "index"))

  ######################################
  ## TIDY AND COMPUTE CLONE ABUNDANCE ##
  ######################################
  print("2. Collapse results into Integration Site Table")

  IS.table <- breaks %>%
    group_by(groupIS) %>%
    # 3.1. TARGET orientation is determined by the orientation of the HOST substring related to the TARGET substring.
    # If enough reads are concordant (+/+), TARGET is oriented "+". If most are discordant (+/- or -/+), TARGET is oriented "-".
    # In the rare event of having IS supported with reads harboring multiple orientations, an "Ambiguous" flag is displayed.
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

  # REMOVE PCR DUPLICATES AND GET THE ABUNDANCE
  # Raw abundance = don't take into account shearsite
  # count = take it into account
  abundance <- IS.table %>%
    select(-ligation, -minDistanceEdge) %>%
    group_by(edge, groupIS) %>%
    mutate(count.raw = length(groupIS[!is.na(position)]),
           count = length(unique(groupShear[!is.na(position)]))) %>%
    ungroup() %>%
    select(edge, groupIS, count.raw, count) %>%
    distinct()

  # GET THE BEST DEFINING POSITION FOR EACH IS
  # POSITION = supported by the max(reads-non-duplicated)
  position <- IS.table %>%
    select(-ligation, -minDistanceEdge) %>%
    group_by(edge, groupIS, position) %>%
    mutate(count = length(unique(groupShear))) %>%  # For each IS position of a groupIS, get the number of non-duplicated reads
    select(edge, groupIS, seqnames.genome, strand, position, count) %>%
    ungroup() %>%
    filter(!is.na(position)) %>%
    group_by(edge, groupIS) %>%
    summarise(
      seqnames = unique(seqnames.genome),
      position = unique(position[count == max(count, na.rm = T)])[1],
      strand = unique(strand)
    ) %>%
    ungroup()

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

  # REGROUP
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
    mutate(count.max = max(edge5.count , edge3.count, na.rm = T),
           ID = paste0(seqnames, ":", min(edge3.pos, edge5.pos, na.rm = T),"-", max(edge3.pos, edge5.pos, na.rm = T))) %>%
    ungroup() %>%
    select(seqnames, strand, edge5.pos, edge5.raw, edge5.count, edge3.pos, edge3.raw, edge3.count, ID, count.max, edge5.meanDistLTR, edge5.sdDistLTR, edge3.meanDistLTR, edge3.sdDistLTR,  groupIS) %>%
    arrange(desc(count.max)) %>%
    # ADD FLAGS:
    mutate(
      FLAG.DIST = ifelse(pmin(.$edge5.meanDistLTR, .$edge3.meanDistLTR, na.rm = T)  > 100, F, T),
      FLAG.READS = ifelse(count.max < 10, F, T)
    )

  # FLAG: Wrong attribution of shearSite | integration sites
  # Dozens of exceptions exists. Although they do not usually impact the counts, the can create false positive IS. Often resulting from the exchange between shearSite and IS
  # FLAG Them:
  breaks.gr <- GRanges(breaks.summerised$ID) + distanceFLAG
  breaks.gr$count.max <- breaks.summerised$count.max

  reduced.win <- GenomicRanges::reduce(breaks.gr)

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
  print("3. Generate the final splitFASTQ table")

  # Only selects reads supporting TRUE IS

  # 5. Return a tibble detailling the read IDs supporting IS with more than 10 reads.
  # This tibble is required to extract reads supporting each IS and call variants.
  fastqID <- select(breaks, readID, groupIS) %>%
    left_join(
      breaks.summerised %>%
        rownames_to_column("index") %>%
        # SELECT ONLY IS WITH 10 non-duplicated reads or those with
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
    # Remove reads not attributed to a good IS
    filter(!is.na(ID))

  breaks.summerised <- select(breaks.summerised, -groupIS)

  return(list(breaks.summerised, fastqID))

}


