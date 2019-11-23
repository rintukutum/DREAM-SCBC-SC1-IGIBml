rm(list=ls())
CL.files1 <- list.files(
  path = './data/processed',
  pattern = '.RData',
  full.names = TRUE
)
CL.files2 <- list.files(
  path = './data/SC1/',
  pattern = '.RData',
  full.names = TRUE
)
CL.files <- c(CL.files1,CL.files2)
rmColnames <- c(
  "treatment", "cell_line", 
  "time", "cellID", "fileID",
  # following variables might contain NA
  'p.Akt.Ser473.','p.ERK','p.HER2','p.PLCg2','p.S6'
)
dir.create('./data/pca',showWarnings = FALSE)
library(plyr)
for(i in 1:length(CL.files)){
  load(CL.files[i])
  cat(paste0(i, '. Processing ',CL.files[i],'\n'))
  pca.trTime <- ddply(
    xx,
    'treatment',
    function(x){
      ddply(
        x,
        'time',
        function(x){
          xy <- x[,setdiff(colnames(x),rmColnames)]
          xy.pca <- prcomp(xy)
          prRot <- data.frame(xy.pca$rotation)
          colNames <- colnames(prRot)
          prRot$varName <- rownames(prRot)
          prRot$cell_line <- unique(x$cell_line)
          return(prRot[,c('varName','cell_line',colNames)])
        }
      )
    })
  if(length(grep('processed',CL.files[i])) != 0){
    pca.file <- gsub('processed','pca',CL.files[i])
  }else{
    pca.file <- gsub('SC1','pca',CL.files[i])
  }
  save(pca.trTime,file=pca.file)
}

rm(list=ls())
pca.files <- list.files('./data/pca/',full.names = TRUE)
library(plyr)
pca.tt.out <- list()
for(i in 1:length(pca.files)){
  load(pca.files[i])
  cat(paste0(i, '. Processing ',pca.files[i],'\n'))
  pca.tt <- ddply(pca.trTime,
        'treatment',
        function(x){
          ddply(
            x,
            'time',
            function(x){
              # PC1
              PC <- x[,c('varName','PC1')]
              PC$component <- apply(
                x[,c('cell_line','treatment','time')],
                1,
                function(x){
                  paste(as.character(x),collapse =  '_')
                }
              )
              return(PC)
            })
        })
  pca.tt.out[[i]] <- pca.tt[,-c(1:2)]
}
pca.tt.df <- ldply(
  pca.tt.out
)
save(pca.tt.df,file = './data/comb-PC1/pca.tt.df.RData')
pc1.mat <- reshape2::dcast(
  pca.tt.df,
  component~varName,
  value.var = 'PC1'
)
pc1.df <- pc1.mat[,-1]
rownames(pc1.df) <- pc1.mat[,1]
dir.create('./data/comb-PC1/',showWarnings = FALSE)
save(pc1.df,file = './data/comb-PC1/pc1.df.RData')
rm(list=ls())
load('./data/comb-PC1/pc1.df.RData')
pdf('./figures/03-PC1-dendogram.pdf',width = 6,height = 6)
heatmap(as.matrix(pc1.df))
dev.off()
CLtt.membership <- data.frame(cutree(hclust(dist(pc1.df)),k = 2:20))
load('./data/comb-PC1/pca.tt.df.RData')

oComp <- hclust(dist(pc1.df))$order
oVar <- hclust(dist(t(pc1.df)))$order
library(ggplot2)
pdf('./figures/03-PC1-dendogra-enhanced.pdf',width = 6,height = 200)
ggplot(pca.tt.df,
       aes(x=varName,y=component)) +
  geom_tile(aes(fill=PC1)) +
  scale_fill_gradient2(
    low = '#37c8abff',high = '#d35fbcff'
  ) +
  scale_x_discrete(limit=colnames(pc1.df)[oVar]) +
  scale_y_discrete(limit=rownames(pc1.df)[oComp]) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
dev.off()
  

getCLname <- function(x){
  strsplit(x,split='\\_')[[1]][1]
}
getTTname <- function(x){
  strsplit(x,split='\\_')[[1]][2]
}
getTime <- function(x){
  strsplit(x,split='\\_')[[1]][3]
}

CLtt.membership$cell_line <- sapply(
  rownames(CLtt.membership),getCLname
)
CLtt.membership$treatment <- sapply(
  rownames(CLtt.membership),getTTname
)
CLtt.membership$time <- sapply(
  rownames(CLtt.membership),getTime
)

sc1.files <- list.files(
  './data/SC1/',
  full.names = TRUE)
extractMembership <- function(sc1.data,CLtt.membership){
  combSC1 <- unique(apply(sc1.data[,c('cell_line','treatment','time')],
        1,function(x){
    gsub(' ','',paste(as.character(x),collapse = '_'))
  }))
  combMembers <- list() 
  for(i in 1:length(combSC1)){
    pattern <- combSC1[i]
    out <- strsplit(pattern,split='\\_')[[1]][3]
    if(length(grep("\\.",out)) == 1){
      if(strsplit(out,split = '\\.')[[1]][2] == '0'){
        pattern <- strsplit(pattern,split = '\\.')[[1]][1]
      }
    }
    idx <- grep(pattern,rownames(CLtt.membership))
    mem.vec <- as.numeric(CLtt.membership[idx,1:19])
    idx.max <- which.max(mem.vec)
    mem <- mem.vec[idx.max]
    mem.close <- rownames(CLtt.membership)[CLtt.membership[,idx.max] == mem]
    mem.count <- table(CLtt.membership[,idx.max])[as.character(mem)]
    combMembers[[i]] <- data.frame(
      SC1 = combSC1[i],
      membership = mem,
      mem.count = length(mem.close),
      members = paste(mem.close,collapse = ';'),
      stringsAsFactors = FALSE
    )
  }
  return(plyr::ldply(combMembers))
}
sc1.CLmembers <- list()
for(i in 1:length(sc1.files)){
  load(sc1.files[i])
  cat(paste0(i,'. ',sc1.files[i],'\n'))
  sc1.CLmembers[[i]] <- extractMembership(
    sc1.data = xx,
    CLtt.membership = CLtt.membership
  )
}
dir.create(
  './data/membership-SC1-extended/',
  showWarnings = FALSE)
save(sc1.CLmembers,
     file = './data/membership-SC1-extended/sc1.CLmembers.RData')
rm(list=ls())
load('./data/membership-SC1-extended/sc1.CLmembers.RData')
sc1.CLmembers.df <- plyr::ldply(
  sc1.CLmembers
)
left.comb <- c(
  'HCC2218_EGF_14.0',"HCC2218_EGF_18.0",   "HCC2218_iEGFR_14.0",
  "HCC2218_iEGFR_18.0", "HCC2218_iMEK_14.0",  "HCC2218_iMEK_18.0",
  "HCC2218_iPI3K_14.0", "HCC2218_iPI3K_18.0", "HCC2218_iPKC_14.0",
  "HCC2218_iPKC_18.0" 
)
getCLname <- function(x){
  strsplit(x,split='\\_')[[1]][1]
}
getTTname <- function(x){
  strsplit(x,split='\\_')[[1]][2]
}
getTime <- function(x){
  strsplit(x,split='\\_')[[1]][3]
}
left.cm.data <- list()
for(i in 1:length(left.comb)){
  x <- left.comb[i]
  mm <- sc1.CLmembers.df[grep(x,sc1.CLmembers.df$SC1),'members']
  mms <- strsplit(mm,split='\\;')[[1]]
  pos <- grep(gsub('\\.0','',x),mms)
  subset.mms <- na.omit(mms[setdiff(c((pos-15):(pos+15)),pos)])
  xx <- data.frame(
    cl.name = sapply(subset.mms,getCLname),
    treatment = sapply(subset.mms,getTTname),
    time = sapply(subset.mms,getTime),
    stringsAsFactors = FALSE
  )
  xxDF <- xx[xx$cl.name != 'HCC2218',]
  left.cm.data[[i]] <- xxDF
}
names(left.cm.data) <- left.comb
lapply(left.cm.data,function(x){unique(x$cl.name)})
sapply(left.cm.data,nrow)
dir.create('./data/SC1-left',showWarnings = FALSE)
save(left.cm.data,file = './data/SC1-left/left.cm.data.RData')
