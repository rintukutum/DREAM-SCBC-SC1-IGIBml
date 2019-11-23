rm(list=ls())
sc1.files <- list.files(
  path = './data/single_cell_phospo/subchallenge_1/',
  pattern = '.csv',
  full.names = TRUE)
dir.create('./data/SC1',showWarnings = FALSE)

getCLname <- function(x){
  xx <- strsplit(x,split='\\/')[[1]]
  return(gsub('.csv','',xx[length(xx)]))
}
for(i in 1:length(sc1.files)){
  print(i)
  xx <- read.csv(
    sc1.files[i],
    stringsAsFactors = FALSE
  )
  CLout <- paste0(
    './data/SC1/',
    getCLname(sc1.files[i]),
    '.RData'
  )
  save(xx, file = CLout)
}

