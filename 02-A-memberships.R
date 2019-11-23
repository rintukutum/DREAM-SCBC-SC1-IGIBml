rm(list=ls())
load('./data/cnv.member.RData')
load('./data/snp.member.RData')
load('./data/expr.member.RData')
load('./data/phos.conds.member.RData')
dim(phos.conds.member$EGF$t0)
sc1 <- c('AU565','EFM19','HCC2218','LY2','MACLS2','MDAMB436')

phos.cut.2 <- list()
cut <- 2
w <- 1
for(i in 1:length(phos.conds.member)){
  for(j in 1:length(phos.conds.member[[i]])){
    member <- phos.conds.member[[i]][[j]][,cut]
    for(k in 1:length(sc1)){
      CM <- setdiff(names(which(member == member[names(member) == sc1[k]])),
                    sc1[k])
      treatment <- names(phos.conds.member)[i]
      time <- names(phos.conds.member[[i]])[j]
      cCM <- length(CM)
      CM.name <- paste(CM,collapse = ';')
      phos.cut.2[[w]] <- data.frame(
        CL.sc1 = sc1[k],
        treatment = treatment,
        time = as.numeric(gsub('t','',time)),
        cCM = cCM,
        CM.name =CM.name,
        stringsAsFactors = FALSE
      )
      w <- w + 1
    }
  }
}
phos.cut.2.df <- plyr::ldply(phos.cut.2)
dir.create('./data/membership-SC1',showWarnings = FALSE)
save(phos.cut.2.df, file = './data/membership-SC1/phos.cut.2.df.RData')