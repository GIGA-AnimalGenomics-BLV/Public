#' Remove duplicated Integration Sites
#'
#' @param IS tibble. Integration site tibble.
#' @param win numeric. up/down window to add to the position in order to see them in IGV
#' @param mode character. "IS" or "HS" for integration sites or hotspots
#'
#' @return
#'
#' @author Vincent Hahaut
#'
#' @examples
#' VIDO_createGTF_IS()
#'
#' @export
IStoGTF <- function(IS = NULL, win = 1, name = "mySample"){

  if(nrow(IS) > 0){
    # 1. Define colors for the negative and positive strands
    col <- as.character(IS$strand)
    col[col=="*"] <-"#153BC4"
    col[col=="+"] <-"#C31516"
    col[col=="-"] <-"#153BC4"

    # 2. Create the file
    IS$start <- apply(IS, 1, function(x) as.numeric(min(x["edge3.pos"], x["edge5.pos"], na.rm=T)))
    IS$end <- apply(IS, 1, function(x) as.numeric(max(x["edge3.pos"], x["edge5.pos"], na.rm=T)))

    index <- 1:nrow(IS)
    gff <- data.frame(IS$seqnames,".","IS",IS$start, IS$end+win,".", IS$strand,".", paste0(
      "ID=", index,
      ";color=",col,
      ";FilteredReads=",IS$count.max,
      ";SampleName=", name,
      sep=""))

  } else {

    gtf <- data_frame()

  }

  return(gff)

}
