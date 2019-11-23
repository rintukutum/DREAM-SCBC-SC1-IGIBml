# install.packages("https://cran.r-project.org/src/contrib/Rcpp_1.0.2.tar.gz",repos=NULL,type="source")
extractTRclData <- function(
  cm.data
){
  #CM.name <- "184A1;BT549;HCC1954;HCC2157;MDAMB468"
  cm.data <- data.frame(cm.data)
  CM.name <- cm.data$CM.name
  cl.id <- strsplit(CM.name,split = ';')[[1]]
  cl.id <- intersect(names(CLmatrix),cl.id)
  cl.tmp.data <- CLmatrix[cl.id]
  cl.tmp.meta <- CLmetainfo[cl.id]
  message("______________________________")
  message("______EXTRACT TRAIN DATA______")
  message(paste0("______TREATMENT ::: ",cm.data$treatment))
  message(paste0("______TIME ::: ",cm.data$time))
  xx.data <- list()
  j <- 1
  idx.any <- c()
  for( i in 1:length(cl.tmp.meta)){
    xx.tmp <- cl.tmp.meta[[i]]
    idx.treat <- xx.tmp$treatment == cm.data$treatment # "iPKC"
    idx.time <- xx.tmp$time == as.numeric(cm.data$time) # 60
    idx.data <- idx.time & idx.treat
    idx.any[i] <- any(idx.data)
    if(any(idx.data) != 0){
      xx.data[[j]] <- cl.tmp.data[[i]][idx.data,]
      j <- j + 1
    }
  }
  names(xx.data) <- cl.id[idx.any]
  message("______JOIN TRAIN DATA")
  xx.dataDF <- ldply(xx.data)
  colnames(xx.dataDF)[1] <- 'CL'
  cm.data$New.cCM <- length(xx.data)
  cm.data$New.CM.name <- paste(cl.id[idx.any],collapse = ';')
  out <- list(
    data = xx.dataDF,
    metainfo = cm.data
  )
  message("______TRAIN DATA READY")
  
  return(out)
}

extractTTclData <- function(cm.data){
  sc1.files <- list.files(
    path = './data/SC1/',
    full.names = TRUE
  )
  sc1.CLname <- gsub(
    '.RData','',
    sapply(
      sc1.files,
      function(x){
        xx <- strsplit(x,split = '\\/')[[1]]
        return(xx[length(xx)])
      })
  )
  load(names(which(sc1.CLname == cm.data$CL.sc1)))
  idx.treat <- xx$treatment == cm.data$treatment
  idx.time <- xx$time == cm.data$time
  idx <- idx.treat & idx.time
  out <- list(
    tt.data = xx[idx,],
    metainfo = cm.data
  )
  return(out)
}
retainExtremeSample <- function(mini.y){
  names(mini.y) <- 1:length(mini.y)
  y.sort <- sort(mini.y)
  # 6% top and bottom
  n12TB <- round((length(y.sort)/100)*6)
  extTB <- c(head(names(y.sort),n12TB),
    tail(names(y.sort),n12TB))
  bias.loop <- setdiff(
    names(y.sort),
    extTB
  )
  set.seed(767)
  fold5 <- caret::createFolds(
    y=mini.y[bias.loop],
    k=5
  )
  fold5.yName <- foreach(f = 1:length(fold5))%do%{
    tmp.y.name <- c(
      names(mini.y[unlist(fold5[-f])]),
      extTB
    )
  }
  return(fold5.yName)
}
ML2missing <- function(
  tr.data,
  tt.data,
  model.path='./data/SC1-MODS/',
  plan='A'){
  TT.data <- tt.data$tt.data
  unqCODE <- apply(TT.data[,1:5],1,function(x){
    gsub(' ','',paste(as.character(x),collapse = '-'))
  })
  rownames(TT.data) <- unqCODE
  sc1.Y <- c('p.Akt.Ser473.','p.ERK','p.HER2','p.PLCg2','p.S6')
  if(plan == 'A'){
    TR.data <- tr.data$data
    sc1.X <- setdiff(colnames(TR.data)[-1],sc1.Y)
  }
  library('doMC')
  library('caret')
  library('glmnet')
  MODs.TT <- foreach(i =1:length(sc1.Y))%do%{
    if(plan != 'A'){
      notPlanA <- dlply(
        tr.data$planB,
        'Y')
      x.cls <- unique(notPlanA[sc1.Y[i]][[1]]$cl.name)
      TR.data <- ldply(tr.data$data[x.cls])[,-1]
      sc1.X <- setdiff(colnames(TR.data)[-c(1:5)],sc1.Y)
      idx.na.y <- is.na(TR.data[,sc1.Y[i]])
      mini.y <- TR.data[!idx.na.y,sc1.Y[i]]
    }else{
      mini.y <- TR.data[,sc1.Y[i]]
    }
    
    message(paste0("__________MODEL :: ",i))
    message(paste0("__________",sc1.Y[i],"__________"))
    names(mini.y) <- 1:length(mini.y)
    bias.fold5 <- retainExtremeSample(mini.y)
    if(plan != 'A'){
      mini.x <- as.matrix(TR.data[!idx.na.y,sc1.X])
    }else{
      mini.x <- as.matrix(TR.data[,sc1.X])
    }
    
    rownames(mini.x) <- 1:length(mini.y)
    registerDoMC(5)
    MODs <- foreach(j = 1:length(bias.fold5))%dopar%{
      #---------------FOLD SPACE
      trset <- bias.fold5[[j]]
      yfold <- mini.y[trset]
      idx.na.y <- is.na(yfold)
      if(any(idx.na.y)){
        trset <- setdiff(trset,trset[idx.na.y])
        yfold <- mini.y[trset]
      }
      xfold <- mini.x[trset,]
      #------------
      ttset <- setdiff(names(mini.y), trset)
      tt.y <- mini.y[ttset]
      idx.na.y <- is.na(tt.y)
      if(any(idx.na.y)){
        ttset <- setdiff(ttset,ttset[idx.na.y])
        tt.y <- mini.y[ttset]
      }
      tt.x <- mini.x[ttset,]
      #---------------MODELS
      #--- LASSO
      fit.lasso <- cv.glmnet(
        xfold, yfold,
        family="gaussian",
        alpha=1)
      pred.lasso <- predict(
        fit.lasso,
        s=fit.lasso$lambda.1se,
        new=tt.x)
      TT.lasso <- predict(
        fit.lasso,
        s=fit.lasso$lambda.1se,
        new=as.matrix(TT.data[,colnames(mini.x)])
      )
      #---- RIDGE
      fit.ridge <- cv.glmnet(
        xfold, yfold, family="gaussian", 
        alpha=0)
      pred.ridge <- predict(
        fit.ridge,
        s=fit.ridge$lambda.1se,
        new=tt.x)
      TT.ridge <- predict(
        fit.ridge,
        s=fit.ridge$lambda.1se,
        new=as.matrix(TT.data[,colnames(mini.x)])
      )
      #---- ELASTIC NET
      fit.elnet <- cv.glmnet(
        xfold, yfold,
        family="gaussian",
        alpha=0.5
      )
      pred.elnet <- predict(
        fit.elnet,
        s=fit.elnet$lambda.1se,
        new=tt.x)
      TT.elnet <- predict(
        fit.elnet,
        s=fit.elnet$lambda.1se,
        new=as.matrix(TT.data[,colnames(mini.x)])
      )
      perf <- data.frame(
        lasso = postResample(pred=pred.lasso,obs = tt.y),
        ridge = postResample(pred=pred.ridge,obs = tt.y),
        elnet = postResample(pred=pred.ridge,obs = tt.y)
      )
      TT.pred <- data.frame(
        lasso = TT.lasso[,1],
        ridge = TT.ridge[,1],
        elnet = TT.elnet[,1]
      )
      message('DONE')
      list(
        perf = perf,
        TT.pred = TT.pred,
        model = list(
          lasso = fit.lasso,
          ridge = fit.ridge,
          elnet = fit.elnet
        )
      )
    }
    registerDoMC(1)
    pred.TT <- lapply(MODs,function(x){
      idx.min <- names(which.min(x$perf['MAE',]))
      x$TT.pred[,idx.min]
    })
    pred.TT.BIND <- do.call('cbind',pred.TT)
    pred.TT.FINAL <- apply(pred.TT.BIND,1,mean)
    
    message(paste0("__________DONE"))
    list(model = MODs,
         pred.TT = pred.TT.FINAL,
         sc1.y = sc1.Y[i]
    )
  }
  modOUT <- apply(
    tr.data$metainfo[,1:3],
    1,
    function(x){paste(x,collapse='-')}
  )
  dir.create(model.path,showWarnings = FALSE)
  modOUT <- paste0(
    model.path,
    modOUT,
    '.RData')
  message("MODEL output: ",modOUT)
  PRED.TT <- lapply(MODs.TT,function(x){x$pred.TT})
  names(PRED.TT) <- sapply(MODs.TT,function(x){x$sc1.y})
  copyTT <- TT.data
  for(i in 1:length(sc1.Y)){
    copyTT[,sc1.Y[i]] <- PRED.TT[sc1.Y[i]][[1]]
  }
  ttColNames <- c(
    "treatment", "cell_line", "time", "cellID", "fileID",
    sc1.Y)
  copyTT <- copyTT[,ttColNames]
  modOUTPUT <- list(
    models = MODs.TT,
    predTT = copyTT
  )
  save(modOUTPUT,
       file = modOUT)
  return(copyTT)
}
ml2fill <- function(cm.data){
  tr.data <- extractTRclData(
    cm.data = cm.data
  )
  tt.data <- extractTTclData(
    cm.data = cm.data
  )
  TT.fill <- ML2missing(tr.data,tt.data)
  return(TT.fill)
}
