#' @title PCIP_getTargetHostSequence
#'
#' @description
#' Give a genomic position and create
#'
#' @param fasta.genome Character. Path/to/genome.fa.
#' @param fasta.target Character. Path/to/target.fa.
#' @param position GRanges compatible position. "chrZ:start-end", position where the target is inserted.
#' @param size Numeric. Size of the genomic window to consider.
#' @param orientation "+" or "-". Reverse the target sequence if "-".
#' @param out.path Character. Path/to/, where the new fasta reference will be written.
#'
#' @return Write a FASTA file containing the proviral sequence inserted in the middle of the genomic sequence of interest.
#'
#' @author Vincent Hahaut
#'
#' @examples
#' PCIP_getTargetHostSequence()
#'
#' @export
PCIP_getTargetHostSequence <- function(fasta.genome = NULL, fasta.target = NULL, position = NULL, size = 10000, orientation = "+", out.path = "~/Desktop/"){

  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(Rsamtools))

  gr <- GRanges(position) + size

  # 1. Read fasta files
  genomeSeq <- open(FaFile(fasta.genome))
  targetSeq <- open(FaFile(fasta.target))

  # 2. Sequences
  genome.seq <- getSeq(genomeSeq, gr)
  target.seq <- getSeq(targetSeq)

  if(orientation == "-"){
    target.seq <- reverse(target.seq)
  }

  # 3. Combine
  middle <- floor(width(gr)/2)
  downstream <- substr(genome.seq, 1, middle)
  upstream <- substr(genome.seq, middle, width(gr))

  targetHost <- DNAStringSet(paste(unlist(downstream), unlist(target.seq), unlist(upstream), sep = ""))
  names(targetHost) <- paste0("host_target")

  # 4. Write
  writeXStringSet(targetHost, paste0(out.path, "/", gsub(x=position, pattern = ":|-", replacement = "_"), "_target.fa"))

}
