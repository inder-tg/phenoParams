
# --- to get paths to datasets

path_polygons <- paste0(getwd(), "/data/polygons")

path_shape <- paste0(getwd(), "/data/shapes")

MASTER <- raster( paste0(getwd(), "/data/master.tif") )

# ---

listFiles_shapes <- list.files(path=path_shape,
                               pattern=".RData",
                               full.names=TRUE)


listFiles_polygons <- list.files(path=path_polygons,
                                 pattern=".RData",
                                 full.names=TRUE)


# --- to get spiral plot

cgu <- rgb(173/255,221/255,142/255)
csos <- rgb(120/255,198/255,121/255)
cmat <- rgb(49/255, 163/255,84/255)
csen <- rgb(217/255, 95/255, 14/255)
ceos <- rgb(254/255, 153/255, 41/255)
cdor <- rgb(208/255, 209/255, 230/255)

# colores <- c(cgu,csos,cmat,csen,ceos,cdor)

# --- to get x-axis for phenoParams plot 

DAYS <- c(31,28,31,30,31,30,31,31,30,31,30,31)

DAYS_AT_YEAR <- c(1,cumsum(DAYS))

middleMonth <- sapply(1:length(DAYS_AT_YEAR), 
                      function(s) median(c(DAYS_AT_YEAR[s],DAYS_AT_YEAR[s+1])) )

middleMonth <- ceiling(middleMonth[!is.na(middleMonth)])

MONTHS <- c("Jan", "Feb", "March", "Abr", "May", "Jun", "Jul",
            "Aug", "Sep", "Oct", "Nov", "Dec")

# ---

rastermat <- function(rast,shape){
  rCrop <- crop(rast,shape)
  rMask <- mask(rCrop,shape)
  shpR <- rasterize(shape,rMask)
  
  rasterToPoints(shpR)
}

getCRS <- function(lat, long){
  sampleData <- data.frame(num = 1, lat = lat, long = long)
  
  coordinates(sampleData) <- ~ long + lat
  proj4string(sampleData) <- "+init=epsg:4326"
  
  sampleData <- spTransform(sampleData, 
                            CRS = "+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" )
  
  sampleData
}
# 

textCRS <- function(lat,long){
  coor <- extent(getCRS(lat=lat,long=long))
  
  v <- c(coor[1], coor[3])
  
  v
}


get_timeSeries_byClicking <- function(toPlot, df){
  nRow <- length(unlist(toPlot)) / 2
  
  mat_toPlot <- matrix(as.numeric(unlist(toPlot)), nrow = nRow)
  
  dX <- matrix(nrow = nrow(df))
  
  dY <- matrix(nrow = nrow(df))
  
  aproxX <- numeric(nRow)
  
  aproxY <- numeric(nRow)
  
  dX <- sapply(1:nRow, function(s) abs(df[,1] - mat_toPlot[s,1]))
  
  aproxX <- sapply(1:nRow, function(s) df[which.min(dX[,s]),1] )
  
  dY <- sapply(1:nRow, function(s) abs(df[,2] - mat_toPlot[s,2]))
  
  aproxY <- sapply(1:nRow, function(s) df[which.min(dY[,s]),2] )
  
  toExtract <- matrix(nrow = nRow, ncol = 2)
  
  toExtract[,1] <- aproxX
  toExtract[,2] <- aproxY
  # 
  IND <- 1:length(df)
  xTemp <- which(df[,1] == toExtract[1,1])
  yTemp <- which(df[xTemp,2] == toExtract[1,2])
  # 
  xyRow <- xTemp[yTemp] # df[xTemp[yTemp],1:2]
  
  list(coord = xyRow)
  # xyRow
}


# --- to randomly trim dataset of large polygon

getDataSetLargePolygon <- function(percent, seed=101){
  
  if(percent > 1 | percent <0){
    stop("percent must be in the interval (0,1)")
  }
  
  numberPOLYGON <- 2
  POLYGON_SHP <- LoadToEnvironment(listFiles_shapes[numberPOLYGON])$shp
  
  rasterPOLYGON <- rastermat(rast=MASTER, shape=POLYGON_SHP)
  
  set.seed(seed = seed)
  samplePOLYGON <- sample(1:nrow(rasterPOLYGON),
                          ceiling( nrow(rasterPOLYGON) * percent ))
  set.seed(NULL)
  
  rasterPOLYGON <- rasterPOLYGON[samplePOLYGON,]
  
  paramsPOLYGON <- LoadToEnvironment(listFiles_params_dtw[numberPOLYGON])$param_list
  paramsPOLYGON_new <- LoadToEnvironment(listFiles_params_dtw_new[numberPOLYGON])$param_list
  
  paramsPOLYGON <- paramsPOLYGON[samplePOLYGON]
  paramsPOLYGON_new <- paramsPOLYGON_new[samplePOLYGON]
  
  
  # POLYGON_dataPolygon <- LoadToEnvironment(listFiles_polygons[POLYGON$NUMBER])$poly[[1]] * 1e-4
  # POLYGON$dataPolygon <- POLYGON_dataPolygon[POLYGON$SAMPLE,]
  # 
  # POLYGON_data_fpca_dtw <- LoadToEnvironment(listFiles_fpca_dtw[POLYGON$NUMBER])$fpca_dtw
  # POLYGON$data_fpca_dtw <- POLYGON_data_fpca_dtw[POLYGON$SAMPLE]
  
  
  dataPOLYGON <- LoadToEnvironment(listFiles_polygons[numberPOLYGON])$poly[[1]] * 1e-4
  dataPOLYGON <- dataPOLYGON[samplePOLYGON,]
  
  data_fpca_dtw_POLYGON <- LoadToEnvironment(listFiles_fpca_dtw[numberPOLYGON])$fpca_dtw
  data_fpca_dtw_POLYGON <- data_fpca_dtw_POLYGON[samplePOLYGON]
  
  list(POLYGON_NUMBER=numberPOLYGON, POLYGON_SAMPLE=samplePOLYGON, 
       POLYGON_RASTER=rasterPOLYGON, 
       POLYGON_PARAMS=paramsPOLYGON, POLYGON_PARAMS_new=paramsPOLYGON_new,
       POLYGON_dataPolygon=dataPOLYGON, POLYGON_data_fpca_dtw=data_fpca_dtw_POLYGON)
}

# --- Added on March 30, 2022
# --- This part might be better included in auxFunctions.R

getPhenoParamsTemp <- function(x, numFreq=3, fitErrorTol=0.01, delta=0.1,
                               L=365, interval=seq(1,365,length=365)){
  
  
  if(numFreq==3){
    fit <- haRmonics(y=x, numFreq=numFreq, fitErrorTol=fitErrorTol,
                     delta=delta)
    
    fprime <- firstDerivative(amp=fit$amplitude,
                              pha=fit$phase,L=L)
    
    fbiprime <- secondDerivative(amp=fit$amplitude,
                                 pha=fit$phase,L=L)
    
    fthrprime <- thirdDerivative(amp=fit$amplitude,
                                 pha=fit$phase,L=L)
  }
  
  if(numFreq==4){
    fit <- haRmonics(y=x, numFreq=numFreq, fitErrorTol=fitErrorTol,
                     delta=delta)
    
    fprime <- firstDerivative_4(amp=fit$amplitude,
                                pha=fit$phase,L=L)
    
    fbiprime <- secondDerivative_4(amp=fit$amplitude,
                                   pha=fit$phase,L=L)
    
    fthrprime <- thirdDerivative_4(amp=fit$amplitude,
                                   pha=fit$phase,L=L)
  }
  
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



