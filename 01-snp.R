rm(list=ls())
readSNP <- function(x){
  con <- file(
    x,
    open = 'r')
  xline <- ''
  i <- 1
  SNPs <- list()
  while(length(xline) == 1){
    xline <- readLines(con = con,n = 1)
    if(length(xline) == 0){
      xline <- NULL
    }else{
      if(i == 1){
        header <- strsplit(xline,split = '","')[[1]]
        j <- 1
      }else{
        xout <-  strsplit(xline,split = ',')[[1]]
        sampleID <- xout[1]
        snps <- xout[-1]
        snps[snps == "NA"] <- NA
        out <- list(
          sampleID = sampleID,
          snp = as.numeric(snps)
        )
        SNPs[[j]] <- out
        j <- j + 1
      }
      print(i)
      i <- i + 1
    }
  }
  close(con)
  return(list(header=header,
              snp=SNPs))
}
snp <- readSNP(
  x = './data/transcriptomics_genomics/SNP_Marcotte.csv'
)
matL <- list()
for(i in 1:length(snp$snp)){
  matL[[i]] <- as.matrix(
    snp$snp[[i]]$snp)
}
matDF <- do.call('cbind',matL)
matDF.noNA <- na.omit(matDF)
rm(matDF);rm(matL)
colnames(matDF.noNA) <- sapply(snp$snp,function(x){x$sampleID})
save(snp,file='./data/snp.RData')
mat.snp <- matDF.noNA
save(mat.snp,file='./data/mat.snp.RData')
snp.dist <- dist(t(mat.snp))
pdf('./figures/00-snp.pdf',10,10)
heatmap(as.matrix(snp.dist),main='SNP')
dev.off()

snp.member <- cutree(hclust(snp.dist),k = 2:10)
save(snp.member,
     file = './data/snp.member.RData')