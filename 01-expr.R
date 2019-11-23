rm(list=ls())
expr <- read.csv(
  './data/transcriptomics_genomics/RNAseq_Marcotte.csv',
  stringsAsFactors = FALSE,
  row.names = 1
)
expr.dist <- dist(as.matrix(expr))
pdf('./figures/00-expr.pdf',10,10)
heatmap(as.matrix(expr.dist),main='EXPR')
dev.off()

expr.member <- cutree(hclust(expr.dist),k = 2:10)
save(expr.member,
     file = './data/expr.member.RData')