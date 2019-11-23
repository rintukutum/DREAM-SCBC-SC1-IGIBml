#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-C", "--cellLine"), type="character", default=NULL, 
              help="cohort i.e., North/Vadu", metavar="character"),
  make_option(c("-T", "--treatment"), type="character", default=NULL, 
              help="debug mode", metavar="character"),
  make_option(c("-t", "--time"), type="character", default=NULL, 
              help="output file name",
              metavar="character"
  )
)
dir.create('./data/SC1-pred/',showWarnings = FALSE)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
output <- paste0(
  './data/SC1-pred/',
  paste(
    c(opt$cellLine,
      opt$treatment,
      opt$time),
    collapse = '-'
  ),
  '.RData'
)
# opt <- list()
# opt$cellLine <- 'HCC2218'
# opt$treatment <- 'EGF'
# opt$time <- '14'
message('cell-line = ',opt$cellLine)
message('treatment = ',opt$treatment)
message('time      = ',opt$time)
message('OUTPUT    = ',output)
#----
library(plyr)
load('./data/SC1-left/left.cm.data.RData')
source('./funcroom.R')

comb <- paste(c(
  opt$cellLine,
  opt$treatment,
  opt$time),collapse = '_')
message(comb)
idx <- grep(comb,names(left.cm.data))
cm.info <- left.cm.data[idx][[1]]
cl.names <- unique(cm.info$cl.name)
cl.files <- list.files(
  './data/processed',
  full.names = TRUE)
tr.data <- list()
for(i in 1:length(cl.names)){
  r.data.file <- cl.files[grep(cl.names[i],cl.files)]
  load(r.data.file)
  tr.data[[i]] <- xx
  rm(xx)
}
tr.data.df <- ldply(
  tr.data
)
TR.data <- list(
  data = tr.data.df[,setdiff(
    colnames(tr.data.df),
    c('treatment','time','cellID','fileID')
    )],
  metainfo = data.frame(
    cl.name = opt$cellLine,
    treatment = opt$treatment,
    time = opt$time,
    stringsAsFactors = FALSE
  )
)
load('./data/SC1/HCC2218.RData')
tt.data <- xx
idx.treat <- opt$treatment == tt.data$treatment
idx.time <- opt$time == tt.data$time
idx <- idx.treat & idx.time
TT.data <- list(tt.data=tt.data[idx,])
CL.filled <- ML2missing(tr.data=TR.data,tt.data=TT.data)
save(CL.filled, file = output)
