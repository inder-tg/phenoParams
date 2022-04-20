LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}

# --- this function takes a time series -x- and returns lengthPeriod x length(x)/lengthPeriod matrix
get_pixel_matrix <- function(x,lenPeriod=23){
  output <- matrix(nrow=length(x)/lenPeriod, ncol=lenPeriod)
  
  for(i in seq_len(nrow(output))){
    output[i,] <- x[((i-1) * lenPeriod + 1):(i * lenPeriod)]
  }
output
}
# ---

# --- this function fill gaps in first 3 years of the series stored in m
# m is a matrix nrow=nYears, ncol=23
fill_gap <- function(m){
apply(m[-1,1:3],MARGIN=2,FUN=median)
}
# ---

# --- this function applies haRmonics to each row of m
get_smoothing <- function(m, method="harmR", numFreq=4, delta=0.01){
  smooth_m <- matrix(nrow = nrow(m), ncol = ncol(m))
  for(i in 1:nrow(m)){
    smooth_m[i,] <- haRmonics(y=m[i,], method=method, 
                              numFreq=numFreq, delta=delta)$fit
  }
  
smooth_m
}

# --- Funcion que utiliza funcion tsclust para hacer clusters...
##### NO SE SI VALE LA PENA LA FUNCION ASI porque solo se dejan
##### parametros por default (solo reduce parametros)

tsclusters_centroid <- function(m, type = 'h', seed, distance = "L2"){
  clust <- tsclust(m, type = type, seed = seed, distance = distance, 
                   centroid = shape_extraction, trace = T, 
                   control = hierarchical_control(method = "average")) 
clust
}



# --- this function returns a harmonic fit at the time t as a function of amplitude,
# --- phase angle and number of observations per period

# NOTE: When amp and pha have been obtained from haRmonics(), then pha[1] is always 
# equal to zero. Hence, amp[1]*cos(0 - pha[1])=amp[1]. 

fitHarmonic <- function(amp,pha,L,t){
  sum(sapply( 0:(length(amp)-1),
             function(k) amp[k+1] * cos( 2 * pi * k * t / L - pha[k+1]/(180/pi) )  
             )
      )
}

harmonicFit <- function(amp,pha,L,t){
  sapply(1:length(t), function(s) fitHarmonic(amp=amp, pha=pha, L=L, t=s) )
}

# --- this function computes a continuous harmonic fit for each row in m

# get_smoothing <- function(m, method="harmR", numFreq=4, delta=0.01){
#   smooth_m <- matrix(nrow = nrow(m), ncol = ncol(m))
#   for(i in 1:nrow(m)){
#     smooth_m[i,] <- haRmonics(y=m[i,], method=method, 
#                               numFreq=numFreq, delta=delta)$fit
#   }
#   
#   smooth_m
# }


get_harmonicFit <- function(samples=100, m, method="harmR", numFreq=4, delta=0.1){
  
  smooth_m_aug <- matrix(nrow=samples, ncol=nrow(m))
  for(i in 1:ncol(smooth_m_aug)){
    TEMP <- haRmonics(y=m[i,], method=method, numFreq=numFreq, delta = delta)
    
    smooth_m_aug[,i] <- sapply(seq(0, (ncol(m)-1), length=samples),
                               function(s) fitHarmonic(amp=TEMP$amplitude,
                                                       pha=TEMP$phase,
                                                       t=s, L=ncol(m)))
  }
  
smooth_m_aug
}


# ---

getFPCA_byPixels_v5 <- function(seriesByCluster, mAug){
  
  flag <- F # se usaron menos de 17 series 
  seriesInCluster <- lapply(1:length(unique(seriesByCluster)), 
                            function(s) which(seriesByCluster == s))
  
  if(length(seriesInCluster[[1]]) >= 10) series <- seriesInCluster[[1]]
  if(length(seriesInCluster[[2]]) >= 10) series <- seriesInCluster[[2]]
  
  AUX <- ifelse(any(table(seriesByCluster) >= 10), T, F)
  if(AUX == F) {
    series <- 1:22
    flag <- T
  }
  
  output <- FPCA::FPCAEst(DATA=mAug[,series], k=3)
  
list(fpca=output$f, usedTotal=flag)
}

# --- update on FEB 15, 2022
# --- We got to take into account that we are only using 3 harmonics

firstDerivative <- function(amp,pha,L){
  function(t) amp[2] * (2*pi/L) * (-1) *
    sin( 2 * pi * t / L - pha[2]/(180/pi) ) +
    amp[3] * (2*pi*(2)/L) * (-1) *
    sin( 2 * pi * t * (2) / L - pha[3]/(180/pi) ) +
    amp[4] * (2*pi*(3)/L) * (-1) *
    sin( 2 * pi * t * (3) / L - pha[4]/(180/pi) )
  # +
  #   amp[5] * (2*pi*(4)/L) * (-1) *
  #   sin( 2 * pi * t * (4) / L - pha[5]/(180/pi) )
}

secondDerivative <- function(amp,pha,L,t){
  
  function(t) amp[2] * (2*pi*(2-1)/L)^2 * (-1) * 
    cos( 2 * pi * (2-1) * t / L - pha[2]/(180/pi) ) +
    amp[3] * (2*pi*(3-1)/L)^2 * (-1) * 
    cos( 2 * pi * (3-1) * t / L - pha[3]/(180/pi) ) +
    amp[4] * (2*pi*(4-1)/L)^2 * (-1) * 
    cos( 2 * pi * (4-1) * t / L - pha[4]/(180/pi) )
  # +
  #   amp[5] * (2*pi*(5-1)/L)^2 * (-1) * 
  #   cos( 2 * pi * (5-1) * t / L - pha[5]/(180/pi) )
}

thirdDerivative <- function(amp,pha,L){
  function(t) amp[2] * (2*pi*(2-1)/L)^3 *  
    sin( 2 * pi * (2-1) * t / L - pha[2]/(180/pi) ) +
    amp[3] * (2*pi*(3-1)/L)^3 *  
    sin( 2 * pi * (3-1) * t / L - pha[3]/(180/pi) ) +
    amp[4] * (2*pi*(4-1)/L)^3 *  
    sin( 2 * pi * (4-1) * t / L - pha[4]/(180/pi) )
  # +
  #   amp[5] * (2*pi*(5-1)/L)^3 *  
  #   sin( 2 * pi * (5-1) * t / L - pha[5]/(180/pi) )
}

# --- Added on March 30, 2022
# --- following functions consider 4 harmonics

firstDerivative_4 <- function(amp,pha,L){
  function(t) amp[2] * (2*pi/L) * (-1) *
    sin( 2 * pi * t / L - pha[2]/(180/pi) ) +
    amp[3] * (2*pi*(2)/L) * (-1) *
    sin( 2 * pi * t * (2) / L - pha[3]/(180/pi) ) +
    amp[4] * (2*pi*(3)/L) * (-1) *
    sin( 2 * pi * t * (3) / L - pha[4]/(180/pi) ) +
    amp[5] * (2*pi*(4)/L) * (-1) *
    sin( 2 * pi * t * (4) / L - pha[5]/(180/pi) )
}

secondDerivative_4 <- function(amp,pha,L){
  
  function(t) amp[2] * (2*pi*(2-1)/L)^2 * (-1) * 
    cos( 2 * pi * (2-1) * t / L - pha[2]/(180/pi) ) +
    amp[3] * (2*pi*(3-1)/L)^2 * (-1) * 
    cos( 2 * pi * (3-1) * t / L - pha[3]/(180/pi) ) +
    amp[4] * (2*pi*(4-1)/L)^2 * (-1) * 
    cos( 2 * pi * (4-1) * t / L - pha[4]/(180/pi) ) +
    amp[5] * (2*pi*(5-1)/L)^2 * (-1) *
    cos( 2 * pi * (5-1) * t / L - pha[5]/(180/pi) )
}

thirdDerivative_4 <- function(amp,pha,L){
  function(t) amp[2] * (2*pi*(2-1)/L)^3 *  
    sin( 2 * pi * (2-1) * t / L - pha[2]/(180/pi) ) +
    amp[3] * (2*pi*(3-1)/L)^3 *  
    sin( 2 * pi * (3-1) * t / L - pha[3]/(180/pi) ) +
    amp[4] * (2*pi*(4-1)/L)^3 *  
    sin( 2 * pi * (4-1) * t / L - pha[4]/(180/pi) ) +
    amp[5] * (2*pi*(5-1)/L)^3 *
    sin( 2 * pi * (5-1) * t / L - pha[5]/(180/pi) )
}

# ---

getGlobalMaxMin <- function(f, crtVal, interval){
  # this function assumes an integer interval
  
  MAX <- crtVal[which.max(f(crtVal))]
  MIN <- crtVal[which.min(f(crtVal))]
  
  truncMAX <- ifelse(MAX%%1 < 0.5, floor(MAX), ceiling(MAX))
  truncMIN <- ifelse(MIN%%1 < 0.5, floor(MIN), ceiling(MIN))
  
  list(maxDate=truncMAX, minDate=truncMIN, origMax=MAX, origMIN=MIN)  
}

# --- auxiliary functions to apply getPhenoParams

get_fpcaFiles <- function(typeFPCA=c("dtw", "l2")){
  if(typeFPCA=="dtw"){
    fpcaFILES <- list.files(paste0(getwd(),"/RData_sTBDF/fpca_dtw_4harm"),
                            pattern = ".RData",full.names = TRUE)
    
  } else {
    fpcaFILES <- list.files(paste0(getwd(),"/RData_sTBDF/fpca_l2_4harm"),
                            pattern = ".RData",full.names = TRUE)
  }
  fpcaFILES
}

get_fpca_object <- function(typeFPCA=c("dtw", "l2"), file){
  fpca_temp <- LoadToEnvironment(file)
  if(typeFPCA == "dtw"){
    fpca_object <- fpca_temp$fpca_dtw
  } else {
    fpca_object <- fpca_temp$fpca_l2
  }
  fpca_object
}

# ---

getPhenoParams <- function(x, numFreq=4, fitErrorTol=0.01, delta=0.1,
                           L=365, interval=seq(0,364,length=365)){
  
  fit <- haRmonics(y=x,numFreq=numFreq,fitErrorTol=fitErrorTol,
                   delta=delta)
  
  fprime <- firstDerivative(amp=fit$amplitude,
                            pha=fit$phase,L=L)
  
  
  fbiprime <- secondDerivative(amp=fit$amplitude,
                               pha=fit$phase,L=L)
  
  fthrprime <- thirdDerivative(amp=fit$amplitude,
                               pha=fit$phase,L=L)
  
  
  crtFrstD <- uniroot.all(fprime, c(interval[1],interval[length(interval)]) )
  
  crtSndD <- uniroot.all(fbiprime, c(interval[1],interval[length(interval)]))
  
  crtThrD <- uniroot.all(fthrprime, c(interval[1],interval[length(interval)]))
  
  
  crtFrstD_test <- uniroot.all(fprime_test, c(interval[1],interval[length(interval)]) )
  
  crtSndD_test <- uniroot.all(fbiprime_test, c(interval[1],interval[length(interval)]))
  
  crtThrD_test <- uniroot.all(fthrprime_test, c(interval[1],interval[length(interval)]))
  
  
  
  maxsSndD <- vector("numeric",2)
  minsSndD <- vector("numeric",2)
  minmaxFrstD <- vector("numeric",2)
  
  globalMaxMinSndD <- getGlobalMaxMin(f=fbiprime, crtVal=crtThrD,
                                      interval=interval)
  maxsSndD[1] <- globalMaxMinSndD$maxDate
  minsSndD[1] <- globalMaxMinSndD$minDate
  
  crtThrDAux <- crtThrD[-c(which.max(fbiprime(crtThrD)),which.min(fbiprime(crtThrD)))]
  localMaxMinSndD <- getGlobalMaxMin(f=fbiprime, crtVal=crtThrDAux,
                                     interval=interval)
  maxsSndD[2] <- localMaxMinSndD$maxDate
  minsSndD[2] <- localMaxMinSndD$minDate
  
  globalMaxMinFrstD <- getGlobalMaxMin(f=fprime, crtVal = crtSndD,
                                       interval=interval)
  minmaxFrstD <- c(globalMaxMinFrstD$minDate, globalMaxMinFrstD$maxDate)
  
  L <- list(GU=maxsSndD[1], SoS=minmaxFrstD[2], Mat=minsSndD[1], 
            Sen=minsSndD[2], EoS=minmaxFrstD[1], Dorm=maxsSndD[2])
  
  phenoDates <- L
  
  cortes <- unlist(phenoDates)
  ordenados <- sort(cortes)
  
  if(which(names(ordenados) == "GU") == 1)
  {names(ordenados) <-  c('GU','SoS','Mat','Sen','EoS','Dorm')}
  if(which(names(ordenados) == "GU") == 2)
  {names(ordenados) <-  c('Dorm','GU','SoS','Mat','Sen','EoS')}
  if(which(names(ordenados) == "GU") == 3)
  {names(ordenados) <-  c('EoS','Dorm','GU','SoS','Mat','Sen')}
  if(which(names(ordenados) == "GU") == 4)
  {names(ordenados) <-  c('Sen','EoS','Dorm','GU','SoS','Mat')}
  if(which(names(ordenados) == "GU") == 5)
  {names(ordenados) <-  c('Mat','Sen','EoS','Dorm','GU','SoS')}
  if(which(names(ordenados) == "GU") == 6)
  {names(ordenados) <-  c('SoS','Mat','Sen','EoS','Dorm','GU')}
  
  phenoDates <- ordenados
  phenoDates
}

getPhenoParamsNew <- function(x,numFreq=4,fitErrorTol=0.01,delta=0.1,
                           L=365,interval=seq(0,364,length=365)){
  
  fit <- haRmonics(y=x,numFreq=numFreq,fitErrorTol=fitErrorTol,
                   delta=delta)
  
  fprime <- firstDerivative(amp=fit$amplitude,
                            pha=fit$phase,L=L)
  
  fbiprime <- secondDerivative(amp=fit$amplitude,
                               pha=fit$phase,L=L)
  
  fthrprime <- thirdDerivative(amp=fit$amplitude,
                               pha=fit$phase,L=L)
  
  crtFrstD <- uniroot.all(fprime, c(interval[1],interval[length(interval)]) )
  
  crtSndD <- uniroot.all(fbiprime, c(interval[1],interval[length(interval)]))
  
  crtThrD <- uniroot.all(fthrprime, c(interval[1],interval[length(interval)]))
  
  maxsSndD <- vector("numeric",2) # both local and global maxs
  minsSndD <- vector("numeric",2) # both local and global mins
  minmaxFrstD <- vector("numeric",2)
  
  globalMaxMinSndD <- getGlobalMaxMin(f=fbiprime, crtVal=crtThrD,
                                      interval=interval)
  maxsSndD[1] <- globalMaxMinSndD$maxDate # GU
  minsSndD[1] <- globalMaxMinSndD$minDate # Mat
  
  crtThrDAux <- crtThrD[-c(which.max(fbiprime(crtThrD)),which.min(fbiprime(crtThrD)))]
  localMaxMinSndD <- getGlobalMaxMin(f=fbiprime, crtVal=crtThrDAux,
                                     interval=interval)
  maxsSndD[2] <- localMaxMinSndD$maxDate # Dorm
  minsSndD[2] <- localMaxMinSndD$minDate # Sen
  
  globalMaxMinFrstD <- getGlobalMaxMin(f=fprime, crtVal = crtSndD,
                                       interval=interval)
  # c(EoS,SoS) 
  minmaxFrstD <- c(globalMaxMinFrstD$minDate, globalMaxMinFrstD$maxDate) # min and max global
  
  L <- list(GU=maxsSndD[1], SoS=minmaxFrstD[2], Mat=minsSndD[1], 
            Sen=minsSndD[2], EoS=minmaxFrstD[1], Dorm=maxsSndD[2])
  
  unlist(L)
}



# ---

addid <- function(name,id){
  zeros <- ifelse(nchar(id)==1, "0000", 
                  ifelse(nchar(id)==2, "000", 
                         ifelse(nchar(id)==3, "00",
                                ifelse(nchar(id)==4, "0", ""))))
  paste0(name, "_", zeros, id)
}

# --- to get SpiralPlot

library(spiralize)
library(vcd)


# --- Colors used in spiral and phenoParam plots

cgu <- rgb(173/255,221/255,142/255)
csos <- rgb(120/255,198/255,121/255)
cmat <- rgb(49/255, 163/255,84/255)
csen <- rgb(217/255, 95/255, 14/255)
ceos <- rgb(254/255, 153/255, 41/255)
cdor <- rgb(208/255, 209/255, 230/255)

colores <- c(cgu,csos,cmat,csen,ceos,cdor)

# ---

meses <- c('Ene','Feb','Mar','Abr','May','Jun','Jul','Ago','Sep','Oct','Nov',
           'Dic')


#' Returns a vector with parameter estimate phenoParam
#' 
#' LIST is a list
#' phenoParam a character indicating the parameter estimate to get
#' 
#' This function lacks of an input checking
#' 
#' Value
#' 
#' Returns a vector 
#' 
getDist_phenoParam <- function(LIST, phenoParam=c("GU", "SoS", "Mat","Sen",
                                                  "EoS","Dorm")){
  # phenoParam <- match.arg(phenoParam)
  unlist(sapply(LIST, 
                function(x) as.numeric(as.character(x[[1]][which(x[[1]][[2]] == phenoParam),1]))))
  
}

# source("auxiliaryObjects.R")

#' Get spiral plot of phenological parameters
#' 
#' LIST is a list 
#' 
getSpiralPlot <- function(LIST, height=0.2, LABELS, ...){
  
  gu <- getDist_phenoParam(LIST=LIST,phenoParam="GU")
  sos <- getDist_phenoParam(LIST=LIST,phenoParam="SoS")
  mat <- getDist_phenoParam(LIST=LIST,phenoParam="Mat")
  sen <- getDist_phenoParam(LIST=LIST,phenoParam="Sen")
  eos <- getDist_phenoParam(LIST=LIST,phenoParam="EoS")
  dorm <- getDist_phenoParam(LIST=LIST,phenoParam="Dorm")
  
  season <- cbind(gu, sos, mat, sen, eos, dorm)
  
  spiral_initialize(xlim=c(0, 360*nrow(season)), 
                    start=360+90, end=360*(nrow(season)+1)+90, 
                    reverse=TRUE, ...)
  
  spiral_track(height=height)
  spiral_axis(major_at=seq(0, 360*8, by=30)[1:12], #curved_labels=TRUE,
              labels=LABELS, labels_gp=gpar(cex=1.25, fontface=2)) #
  
  X <- c()
  for(i in 1:nrow(season)){
    X <- c(X, 360 * (season[i,]-1)/364 + (i-1) * 360)
  }
  
  spiral_points(x=X, y=0.5, pch=18, 
                gp=gpar(col=rep(colores,8), cex=1.5))
  
  # grid_legend(1.15,0.25, pch=18, col=colores,
  #             frame=FALSE,
  #             labels = c("GU","SoS","Mat","Sen","EoS","Dorm"),
  #             title="Params")
}



