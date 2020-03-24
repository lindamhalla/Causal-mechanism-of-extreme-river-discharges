#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The Danube dataset as well as the declustering functions can be downloaded from the supplementary material of
# Asadi, P., Davison, A. C. and Engelke, S. (2015) Extremes on river networks. The Annals of Applied Statistics, 9, 2023–2050.
setwd("AOAS863_CodesAndData")
source("Codes/Functions.R")

### Loading Data ###
FileList <- list.files(path="Data")
for (i in 1:length(FileList)){
  load(paste("Data",FileList[i],sep="/"))
}

StsChos       <- c(1:31)
NoSt          <- length(StsChos)
CatchEucDis   <- as.matrix(dist(CatCntr))
CatchEucDisWt <- as.matrix(dist((CatCntrWt)))

### Declustering and obtaining events ###
Years          <- unique(as.numeric(substr(ComTSs[,1],1,4)))
U              <- 0.1
Lag            <- 4 #leading to a 9-day period
Plotting       <- 0
StNames        <- as.character(StsChos)
YearsWithEvent <- numeric()
AllEvMat       <- matrix(0,length(StNames),1)
AllRes         <- list()
for (i in 1:length(Years)){
  Res       <- ObtainMultVarEvents(TSs=ComTSs,U=U,Lag=Lag,Year=Years[i],StNames=StNames,
                                   Plotting=Plotting,mfrow=c(NoSt,1))  ## For Censored Likelihood
  Events    <- Res$AllEvMat
  Locations <- Res$Location
  if (length(Events)>0){
    YearsWithEvent <- c(YearsWithEvent,rep(Years[i],dim(Events)[2]))
    AllEvMat       <- cbind(AllEvMat,Events)
  }
  AllRes[[i]] <- Res 
}

DataEvents           <- t(AllEvMat[,-1])
rownames(DataEvents) <- YearsWithEvent 

### Final dataset of declustered data, i.e., the matrix of 428 independent events (rows) and 31 sites (columns)
TSNew           <- DataEvents
rownames(TSNew) <- YearsWithEvent 