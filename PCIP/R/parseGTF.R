#' Read a GTF file into R
#'
#' @param gtf character. Path to the GTF file.
#' @param type character. Which feature should you extract
#' @param prefix character. Prefix to add to the chromosome names
#'
#' @return
#'
#' @author Vincent Hahaut & Nicolas Rosewick
#'
#' @examples
#' parseGTF()
#'
#' @export
parseGTF <- function(gtf = NULL, type = "exon", prefix = "OAR"){

  # 1. Read table
  gtf <- read.table(gtf, sep="\t", as.is = T)

  # 2. Filter and add prefix
  gtf <- gtf[gtf[,3]==type,]
  gtf[,1] <- paste(prefix,gtf[,1],sep="")

  # 3. Clean file
  g <- do.call(cbind,strsplit(do.call(cbind,strsplit(gtf[,9],";"))[1,]," "))[2,]
  res <- data.frame(do.call(rbind,lapply(split(gtf,g),function(x){
    if(length(grep(x=x[1,9],pattern="gene_name"))>0){
      c(x[1,1],min(x[,4]),max(x[,5]),x[1,7], strsplit(sub("*;","",sub("^\\s+", "",sub('.*gene_id', '', x[1,9])))," ")[[1]][1], strsplit(sub("*;","",sub("^\\s+", "",sub('.*gene_name', '', x[1,9])))," ")[[1]][1],strsplit(sub("*;","",sub("^\\s+", "",sub('.*gene_biotype', '', x[1,9])))," ")[[1]][1])
    }else{
      c(x[1,1],min(x[,4]),max(x[,5]),x[1,7], strsplit(sub("*;","",sub("^\\s+", "",sub('.*gene_id', '', x[1,9])))," ")[[1]][1],NA,strsplit(sub("*;","",sub("^\\s+", "",sub('.*gene_biotype', '', x[1,9])))," ")[[1]][1])
    }
  })),stringAsFactors=F)

  # 4. Rename
  colnames(res)<-c("chr","start","end","orientation","eid","gene","biotype")

  res[,8] <- NULL

  return(res)




}
