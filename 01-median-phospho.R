rm(list=ls())
library(ggplot2)
medPHOS <- read.csv(
  './data/median_phospho/median_phospho_data.csv',
  stringsAsFactors = FALSE
)
medPHOSmeta <- data.frame(table(medPHOS[,1:3]))
medPHOSmeta$time <- factor(
  as.numeric(as.character(medPHOSmeta$time)),
  levels = sort(unique(as.numeric(as.character(medPHOSmeta$time))))
)

pdf('./figures/00-medPHOS-cell-line-and-time-and-treatment.pdf',
    width = 14,height = 9)
ggplot(medPHOSmeta,aes(x=time,y=cell_line)) +
  geom_tile(aes(fill=log10(Freq)), width=0.8,height=0.8) +
  facet_grid(.~treatment) +
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
dev.off()


#----------
# Median PHOS dynamics
library(plyr)
timePHOSmed <- dlply(
  medPHOS,
  'time',
  function(x){
    dlply(
      x,
      'treatment'
    )
  }
)
names(timePHOSmed) <- paste0(
  't',
  names(timePHOSmed)
)
timePHOSmed.df <- lapply(
  timePHOSmed,
  function(x){
    lapply(
      x,
      function(x){
        rownames(x) <- as.character(x$cell_line)
        return(x[,-(1:3)])
      })
  }
)
timePHOSmed.cor <- lapply(
  timePHOSmed.df,
  function(x){
    lapply(
      x,
      function(x){
        cor(x = t(x),
            use = "na.or.complete")
      }
    )
  }
)
extractPairDF <- function(mat){
  # mat <- matrix(1:100,nrow = 10)
  # colnames(mat) <- LETTERS[1:10]
  # rownames(mat) <- LETTERS[11:20]
  idx <- lower.tri(mat)
  nCol <- ncol(mat)
  nRow <- nrow(mat)
  colNames <- colnames(mat)
  rowNames <- rownames(mat)
  rowNamesOrder <- rep(rowNames,nRow)[idx]
  
  colNamesOrder <- rep(colNames[-nCol],nRow - (1:(nRow-1)))
  pair.df <- data.frame(
    rowName = rowNamesOrder,
    colName = colNamesOrder,
    value = mat[idx],
    stringsAsFactors = FALSE
  )
  return(pair.df)
}
timePHOSmed.pair <- lapply(
  timePHOSmed.cor,
  function(x){
    lapply(
      x,
      function(x){
        extractPairDF(x)
      }
    )
  }
)
timePHOSmed.pair.cond <- list()
z <- 1
times <- names(timePHOSmed.pair)
for(j in 1:length(times)){
  conds <- names(timePHOSmed.pair[times[j]][[1]])
  for(i in 1:length(conds)){
    df.cond <- timePHOSmed.pair[times[j]][[1]][conds[i]][[1]]
    df.cond$treatment <- conds[i]
    df.cond$time <- times[j]
    timePHOSmed.pair.cond[[z]] <- df.cond
    z <- z + 1
  }
}
idx <- sapply(sapply(timePHOSmed.pair.cond,nrow),is.null)
timePHOSmed.pair.cond.df <- ldply(
  timePHOSmed.pair.cond[!idx]
)

treatTIME.map <- dlply(
  timePHOSmed.pair.cond.df,
  'treatment'
)

for(i in 1:length(treatTIME.map)){
  outname <- paste0(
    './figures/00-medPOS-cor-cell-line-',
    names(treatTIME.map)[i],
    '.pdf'
  )
  pdf(outname,
      width = 22,height = 15)
  mat.df <- treatTIME.map[[i]]
  mat.df$time <- gsub('t','',mat.df$time)
  mat.df$time <- factor(
    as.numeric(as.character(mat.df$time)),
    levels = sort(unique(as.numeric(as.character(mat.df$time))))
  )
  
    
  print(
    ggplot(mat.df,aes(x = rowName,y=colName)) +
  geom_tile(aes(fill=value),width=0.9,height=0.9) +
  facet_wrap(facets = 'time',nrow = 2) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))
  )
  dev.off()
}
#--------------------
timePHOSmed.corL <- list()
z <- 1
times <- names(timePHOSmed.cor)
for(j in 1:length(times)){
  conds <- names(timePHOSmed.cor[times[j]][[1]])
  for(i in 1:length(conds)){
    df.cond <- timePHOSmed.cor[times[j]][[1]][conds[i]][[1]]
    timePHOSmed.corL[[z]] <- list(
      cor = df.cond,
      time = times[j],
      treatment = conds[i]
    )
    z <- z + 1
  }
}
idx.treat <- sapply(
  timePHOSmed.corL,
  function(x){
    x$treatment
  }
)

idx.time <- sapply(
  timePHOSmed.corL,
  function(x){
    x$time
  }
)
treatCOR <- list()
treatments <- unique(idx.treat)
for(i in 1:length(treatments)){
  tCOR <- timePHOSmed.corL[idx.treat == treatments[i]]
  names(tCOR) <- sapply(
    tCOR,
    function(x){
      x$time
    }
  )
  treatCOR[[i]] <- lapply(
    tCOR,
    function(x){
      x$cor
    }
  )
}
names(treatCOR) <- treatments
#-------------
conds <- names(treatCOR)
conds.member <- list()
for(i in 1:length(conds)){
  ts <- names(treatCOR[conds[i]][[1]])
  outname <- paste0(
    './figures/00-median-',
    conds[i],
    '.pdf'
    )
  pdf(outname,width = 10,height = 10)
  member <- list()
  for(j in 1:length(ts)){
    dat_ <- treatCOR[conds[i]][[1]][ts[j]][[1]]
    heatmap(dat_,main = paste0(conds[i],'=',ts[j]))
    member[[j]] <- cutree(hclust(dist(dat_)),k = 2:5)
  }
  names(member) <- ts
  dev.off()
  conds.member[[i]] <- member
}
names(conds.member) <- conds
phos.conds.member <- conds.member
save(phos.conds.member,file = './data/phos.conds.member.RData')
#-------------
library(dendextend)
dendt0 <- treatCOR$EGF$t0 %>% dist %>% hclust(method = "average") %>% as.dendrogram
dendt5.5 <- treatCOR$EGF$t5.5 %>% dist %>% hclust(method = "average") %>% as.dendrogram
dendComb <- dendlist(dendt0,dendt5.5)
tanglegram(dendComb)

###
tangleTREAT <- function(cor0,cor1){
  commonVar <- intersect(colnames(cor1),colnames(cor0))
  dendt0 <- cor0[commonVar,commonVar] %>% dist %>% hclust(method = "average") %>% as.dendrogram
  dendt5.5 <- cor1[commonVar,commonVar] %>% dist %>% hclust(method = "average") %>% as.dendrogram
  dendComb <- dendlist(dendt0,dendt5.5)
  entangleScore <- round(
    entanglement(dendComb), 2
  )
  tanglegram(
    dendComb,
    common_subtrees_color_branches = TRUE
    ) %>% plot(main = paste("entanglement =", entangleScore))
  return(entangleScore)
}
tangleTREAT(
  cor0 = treatCOR$full$t0,
  cor1 = treatCOR$EGF$t0
)
tangleTREAT(
  cor0 = treatCOR$full$t0,
  cor1 = treatCOR$iEGFR$t0
)
tangleTREAT(
  cor0 = treatCOR$full$t0,
  cor1 = treatCOR$iMEK$t0
)
#------------
#
tangleTREAT(
  cor0 = treatCOR$full$t0,
  cor1 = treatCOR$EGF$t0
)
