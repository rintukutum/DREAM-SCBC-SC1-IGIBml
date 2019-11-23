rm(list=ls())
cnv <- read.csv(
  './data/transcriptomics_genomics/CNV_Marcotte.csv',
  stringsAsFactors = FALSE,
  row.names = 1
)
cnv.dist <- dist(cnv)
pdf('./figures/00-cnv.pdf',width = 10,height = 10)
heatmap(as.matrix(cnv.dist),main='CNV')
dev.off()
cnv.member <- cutree(hclust(cnv.dist),k = 2:10)
save(cnv.member,
     file = './data/cnv.member.RData')