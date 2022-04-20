
POLYGON <- reactiveValues()

observeEvent(input$selectPolygon,{
  
  isolate({
    if(input$selectPolygon=="Santiago Amoltepec"){
      POLYGON$SHP <- LoadToEnvironment(listFiles_shapes[3])$shp
      POLYGON$NUMBER <- 3
      POLYGON$RASTER <- rastermat(rast=MASTER, shape=POLYGON$SHP)
    }
    
    if(input$selectPolygon=="Boquerón de Tonalá"){
      POLYGON$SHP <- LoadToEnvironment(listFiles_shapes[1])$shp
      POLYGON$NUMBER <- 1
      POLYGON$RASTER <- rastermat(rast=MASTER, shape=POLYGON$SHP)
    }
    
    if(input$selectPolygon=="Tehuacán-Culcatlán"){
      POLYGON$SHP <- LoadToEnvironment(listFiles_shapes[2])$shp
      POLYGON$NUMBER <- 2
    }
    
    POLYGON$RASTER <- rastermat(rast=MASTER, shape=POLYGON$SHP)
  })
  
})


output$polygon <- renderLeaflet({
  crs <- CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  xyTEST <- SpatialPoints(cbind(POLYGON$RASTER[,1],
                                POLYGON$RASTER[,2]), proj4string=crs)
  longlatTEST <- spTransform(xyTEST, CRS("+proj=longlat +datum=WGS84"))
  
  POINTS <- data.frame(x=longlatTEST@coords[,1], y=longlatTEST@coords[,2],
                       row = 1:nrow(POLYGON$RASTER))
  
  mp_points <- leaflet(POINTS) %>%
    addTiles(group="OSM") %>%
    addProviderTiles(providers$Esri.WorldImagery, group="Esri") %>%
    addProviderTiles(providers$OpenTopoMap, group="Topo") %>%
    addProviderTiles(providers$CartoDB.Positron, group="Carto") %>%
    addLayersControl(baseGroups = c("OSM", "Esri", "Topo", "Carto")) %>%
    addCircles(~x, ~y, radius= ~1, #color = ~palCol(pendiente),
               fillOpacity = 1, 
               popup = paste0("<b>lon: </b>", round(POINTS$x, digits=5), "<br>",
                              "<b>lat: </b>", round(POINTS$y, digits=5), "<br>",
                              "<b>row: </b>", POINTS$row) )
})

output$spiral <- renderPlot({
  POLYGON$PARAMS <- LoadToEnvironment(listFiles_params_dtw[POLYGON$NUMBER])$param_list
  POLYGON$PARAMS_old <- LoadToEnvironment(listFiles_params_dtw_old[POLYGON$NUMBER])$param_list
  # POLYGON$PARAMS_new <- LoadToEnvironment(listFiles_params_dtw_new[POLYGON$NUMBER])$param_list
  
  getSpiralPlot(LIST=POLYGON$PARAMS, LABELS=MONTHS,
                vp_param=list(width=0.5, height=0.7))
  
  grid_legend(1.215, 0.125, pch=18, col=colores,
              frame=FALSE,
              labels=c("GU","SoS","Mad","Sen","EoS","Dorm"),
              title="Params")
})

pixel <- reactiveValues()

observeEvent(input$polygon_click,{
  if(length(input$polygon_click$lng)==0){
    return()
  } else {
    
    isolate({
      LONG <- input$polygon_click$lng #input$longitude_x LONG is y-axis
      LAT <- input$polygon_click$lat #input$latitude_y LAT is x-axis
      
      AUX <- textCRS(lat=LAT, long=LONG)
      ROW <- get_timeSeries_byClicking(toPlot=AUX, df=POLYGON$RASTER)$coord
      
      pixel$ROW <- ROW
      
      POLYGON$dataPolygon <- LoadToEnvironment(listFiles_polygons[POLYGON$NUMBER])$poly[[1]] * 1e-4
      
      ndviTemp <- as.numeric(POLYGON$dataPolygon[ROW,])
      
      ndviMat <- get_pixel_matrix(x=ndviTemp)
      
      ndviMat[1,1:3] <- fill_gap(m=ndviMat)
      
      pixel$ndviMat <- ndviMat
      
      ndviSmoothedSampled <- get_harmonicFit(samples=100, m=ndviMat,
                                             method="harmR", numFreq=4, delta=0.1)
      
      pixel$ndviSmoothedSampled <- ndviSmoothedSampled
      
      POLYGON$data_fpca_dtw <- LoadToEnvironment(listFiles_fpca_dtw[POLYGON$NUMBER])$fpca_dtw
      
      ndviFPCA <- POLYGON$data_fpca_dtw[[ROW]]
      
      fit <- haRmonics(y=ndviFPCA$fpca, numFreq=4, delta=0.1)
      
      pixel$fit <- fit
      
      ndviParams <- POLYGON$PARAMS[ROW]
      
      # ndviParams_new <- POLYGON$PARAMS_new[ROW]
      
      ndviParams_old <- POLYGON$PARAMS_old[ROW]
      
      pixel$ndviParams <- ndviParams
      
      # pixel$ndviParams_new <- ndviParams_new
      
      pixel$ndviParams_old <- ndviParams_old
      
      ndvi365 <- harmonicFit(amp=fit$amplitude, pha=fit$phase,
                             L=365, t=seq(0,364,length=365))
      
      pixel$ndvi365 <- ndvi365
    })
  }
})


output$the_three_plots <- renderUI({
  div(
    style="position:relative",
    tabBox(id="the_three_plots",
           width=NULL,
           height=300,
           tabPanel("Raw data",
                    fluidRow(
                      column(12,
                             wellPanel(
                               plotOutput("raw_data", width="100%", height=400)
                             )
                      )
                    )
           ),
           tabPanel("Smoothed & FPCA",
                    fluidRow(
                      column(12,
                             wellPanel(
                               plotOutput("smoothed_fpca", width="100%", height=400)
                             )
                      )
                    )
           ),
           tabPanel("Phenological Params",
                    fluidRow(
                      column(12,
                             wellPanel(
                               plotOutput("phenoParams", width="100%", height=800)
                             )
                      )
                    )
           ),
           tabPanel("Scatterplot (all-in)",
                    fluidRow(
                      column(12,
                             wellPanel(
                               plotOutput("scatterplot", width="100%", height=600)
                             )
                      )
                    )
           )
    )
  )
})

output$raw_data <- renderPlot({
  
  if(length(input$polygon_click$lng)==0){
    return() # ensures that plotting area is empty when app is loaded
  } else {
    
    isolate({
      par(mar=c(2,4,2,1))
      yRan <- range(pixel$ndviMat)
      plot(pixel$ndviMat[1,], col=1, type="l", ylab="NDVI", ylim=yRan)
      legend("topleft", legend=c(paste0("row: ", pixel$ROW)), bty="n",
             cex=2)
      for(i in 2:nrow(pixel$ndviMat)){
        lines(1:23, pixel$ndviMat[i,], col=i)
      }
    })
    
  }
})

output$smoothed_fpca <- renderPlot({
  par(mar=c(2,4,2,1))
  yRan <- range(pixel$ndviSmoothedSampled)
  plot(pixel$ndviSmoothedSampled[,1], col=1, type="l", ylab="NDVI", ylim=yRan)
  for(i in 2:22){
    lines(1:100, pixel$ndviSmoothedSampled[,i ], col=i)
  }
  par(new=TRUE)
  plot(pixel$fit$fitted, col="red", type="l", lwd=4, ylab="NDVI", ylim=yRan)
})

output$phenoParams <- renderPlot({
  
  fechas <- c()
  tipos <- c()
  if(is.factor(pixel$ndviParams_old[[1]]$params$Fecha)){
    fechas <- as.vector(pixel$ndviParams_old[[1]]$params$Fecha)
    tipos <- as.vector(pixel$ndviParams_old[[1]]$params$Tipo)
  } else {
    fechas <- pixel$ndviParams_old[[1]]$params$Fecha
    tipos <- pixel$ndviParams_old[[1]]$params$Tipo
  }
  
  # fechas_new <- c()
  # tipos_new <- c()
  # if(is.factor(pixel$ndviParams_new[[1]]$params$Fecha)){
  #   fechas_new <- as.vector(pixel$ndviParams_new[[1]]$params$Fecha)
  #   tipos_new <- as.vector(pixel$ndviParams_new[[1]]$params$Tipo)
  # } else {
  #   fechas_new <- pixel$ndviParams_new[[1]]$params$Fecha
  #   tipos_new <- pixel$ndviParams_new[[1]]$params$Tipo
  # }
  
  fechas_4harm <- c()
  tipos_4harm <- c()
  if(is.factor(pixel$ndviParams[[1]]$params$Fecha)){
    fechas_4harm <- as.vector(pixel$ndviParams[[1]]$params$Fecha)
    tipos_4harm <- as.vector(pixel$ndviParams[[1]]$params$Tipo)
  } else {
    fechas_4harm <- pixel$ndviParams[[1]]$params$Fecha
    tipos_4harm <- pixel$ndviParams[[1]]$params$Tipo
  }
  
  
  # COLORES <- c()
  #
  # if(which( tipos == "GU") == 1)
  # {COLORES <-  c(cgu, csos, cmat, csen, ceos, cdor)}
  # if(which( tipos == "GU") == 2)
  # {COLORES <-  c(cdor, cgu, csos, cmat, csen, ceos)}
  # if(which( tipos == "GU") == 3)
  # {COLORES <-  c(ceos, cdor, cgu, csos, cmat, csen)}
  # if(which( tipos == "GU") == 4)
  # {COLORES <-  c(csen, ceos, cdor, cgu, csos, cmat)}
  # if(which( tipos == "GU") == 5)
  # {COLORES <-  c(cmat, csen, ceos, cdor, cgu, csos)}
  # if(which( tipos == "GU") == 6)
  # {COLORES <-  c(csos, cmat, csen, ceos, cdor, cgu)}
  
  # fechas <- pixel$ndviParams[[1]]$params$Fecha
  # COLORES <- colores
  
  # par(mar=c(2,4,2,1),mfrow=c(1,1))
  # yRan <- range(pixel$ndviSmoothedSampled)
  # plot(pixel$ndvi365, col="red", type="l", lwd=4, ylab="NDVI", ylim=yRan,
  #      xaxt="n")
  # axis(side=1, at=DAYS_AT_YEAR, labels=c(MONTHS,""))
  # abline(v=fechas, lwd=3, col=COLORES)
  
  fprime <- firstDerivative(amp=pixel$fit$amplitude,
                            pha=pixel$fit$phase, L=365)
  
  fbiprime <- secondDerivative(amp=pixel$fit$amplitude,
                               pha=pixel$fit$phase, L=365)
  
  
  par(mar=c(2,4,2,1),mfrow=c(3,2))
  yRan <- range(pixel$ndviSmoothedSampled)
  plot(pixel$ndvi365, col="red", type="l", lwd=4, ylab="NDVI", ylim=yRan,
       xaxt="n", cex.lab=1.25, font.lab=4)
  axis(side=1, at=DAYS_AT_YEAR, labels=c(MONTHS,""))
  abline(v=fechas, lwd=3, col=colores)
  plot(pixel$ndvi365, col="red", type="l", lwd=4, ylab="NDVI", ylim=yRan,
       xaxt="n", cex.lab=1.25, font.lab=4)
  axis(side=1, at=DAYS_AT_YEAR, labels=c(MONTHS,""))
  abline(v=fechas_4harm, lwd=3, col=colores)
  
  yRan <- range(fprime(seq(0,364,length=365)))
  yRan[1] <- yRan[1]-0.0025
  yRan[2] <- yRan[2]+0.0025
  plot(fprime(seq(0,364,length=365)), ylim=yRan, ylab="NDVI'", xaxt="n",
       lwd=4, col="blue", cex.lab=1.25, font.lab=4)
  axis(side=1, at=DAYS_AT_YEAR, labels=c(MONTHS,""))
  abline(v=fechas, lwd=3, col=colores)
  plot(fprime(seq(0,364,length=365)), ylim=yRan, ylab="NDVI'", xaxt="n",
       lwd=4, col="blue", cex.lab=1.25, font.lab=4)
  axis(side=1, at=DAYS_AT_YEAR, labels=c(MONTHS,""))
  abline(v=fechas_4harm, lwd=3, col=colores)
  
  
  yRan <- range(fbiprime(seq(0,364,length=365)))
  yRan[1] <- yRan[1]-0.00005
  yRan[2] <- yRan[2]+0.00005
  plot(fbiprime(seq(0,364,length=365)), ylim=yRan, ylab="NDVI''", xaxt="n",
       lwd=4, col="magenta", cex.lab=1.25, font.lab=4)
  axis(side=1, at=DAYS_AT_YEAR, labels=c(MONTHS,""))
  abline(v=fechas, lwd=3, col=colores)
  plot(fbiprime(seq(0,364,length=365)), ylim=yRan, ylab="NDVI''", xaxt="n",
       lwd=4, col="magenta", cex.lab=1.25, font.lab=4)
  axis(side=1, at=DAYS_AT_YEAR, labels=c(MONTHS,""))
  abline(v=fechas_4harm, lwd=3, col=colores)
  
})

output$scatterplot <- renderPlot({
  
  getGU <- getDist_phenoParam(LIST=POLYGON$PARAMS, phenoParam="GU")
  getSoS <- getDist_phenoParam(LIST=POLYGON$PARAMS, phenoParam="SoS")
  getMat <- getDist_phenoParam(LIST=POLYGON$PARAMS, phenoParam="Mat")
  getSen <- getDist_phenoParam(LIST=POLYGON$PARAMS, phenoParam="Sen")
  getEoS <- getDist_phenoParam(LIST=POLYGON$PARAMS, phenoParam="EoS")
  getDorm <- getDist_phenoParam(LIST=POLYGON$PARAMS, phenoParam="Dorm")
  
  params <- data.frame(GU=getGU, SoS=getSoS, Mat=getMat,
                       Sen=getSen, EoS=getEoS, Dorm=getDorm)
  
  pairs(~GU+SoS+Mat+Sen+EoS+Dorm, params)
  
})
