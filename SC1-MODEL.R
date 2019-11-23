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
message('cell-line = ',opt$cellLine)
message('treatment = ',opt$treatment)
message('time      = ',opt$time)
message('OUTPUT      = ',output)
load('./data/CL/CLmetainfo.RData')
load('./data/CL/CLmatrix.RData')
load('./data/membership-SC1/phos.cut.2.df.RData')
#----
library(plyr)
sc1.tr.map <- dlply(
  phos.cut.2.df,
  'CL.sc1',
  function(x){
    dlply(x,
          'treatment',
          function(x){
            dlply(
              x,
              'time'
            )
          })
  }
)
source('./funcroom.R')
CL.tmp <- sc1.tr.map[opt$cellLine][[1]]
cm.data <- CL.tmp[opt$treatment][[1]][opt$time][[1]]
CL.filled <- ml2fill(cm.data = cm.data)
save(CL.filled, file = output)
