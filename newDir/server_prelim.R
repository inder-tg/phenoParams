
POLIGONO <- reactiveValues()

observeEvent(input$polygonSelect,{
  
    isolate({
      if(input$polygonSelect=="Boqueron de Tonala"){
        POLIGONO$NUMBER <- 1
      }

      if(input$polygonSelect=="Tehuacan-Culcatlan"){
        POLIGONO$NUMBER <- 2
      }

      if(input$polygonSelect=="Santiago Amoltepec"){
        POLIGONO$NUMBER <- 3
      }
    })

  POLIGONO$SHP <- LoadToEnvironment(listFiles_shapes[POLIGONO$NUMBER])$shp
  POLIGONO$RASTER <- rastermat(rast=MASTER, shape=POLIGONO$SHP)

})

output$poligono <- renderLeaflet({
  
  crs <- CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  xyTEST <- SpatialPoints(cbind(POLIGONO$RASTER[,1],
                                POLIGONO$RASTER[,2]), proj4string=crs)
  longlatTEST <- spTransform(xyTEST, CRS("+proj=longlat +datum=WGS84"))
  
  POINTS <- data.frame(x=longlatTEST@coords[,1], y=longlatTEST@coords[,2],
                       row = 1:nrow(POLIGONO$RASTER))
  
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

pixelio <- reactiveValues()

observeEvent(input$poligono_click,{
  if(length(input$poligono_click$lng)==0){
    return()
  } else {

    isolate({
      LONG <- input$poligono_click$lng
      LAT <- input$poligono_click$lat #input$latitude_y LAT is x-axis

      AUX <- textCRS(lat=LAT, long=LONG)
      ROW <- get_timeSeries_byClicking(toPlot=AUX, df=POLIGONO$RASTER)$coord

      pixelio$ROW <- ROW

      POLIGONO$dataPolygon <- LoadToEnvironment(listFiles_polygons[POLIGONO$NUMBER])$poly[[1]] * 1e-4

      ndviTemp <- as.numeric(POLIGONO$dataPolygon[ROW,])

      ndviMat <- get_pixel_matrix(x=ndviTemp)

      ndviMat[1,1:3] <- fill_gap(m=ndviMat)

      pixelio$ndviMat <- ndviMat
    })

  }
})

output$test <- renderUI({
  div(
    style="position:relative",
    tabBox(id="test",
           width=NULL,
           height=300,
           tabPanel("Raw data",
                    fluidRow(
                      column(12,
                             wellPanel(
                               plotOutput("raw_data_prelim", width="100%", height=400)
                             )
                      )
                    )
           ),
           tabPanel("Smoothed & FPCA",
                    fluidRow(
                      column(12,
                        wellPanel(
                          uiOutput("parameters_fpca", width="100%", height=150)
                        )
                      )
                    ),
                    fluidRow(
                      column(6,
                             wellPanel(
                               plotOutput("smoothed_fpca_prelim", width="100%", height=400)
                             )
                      ),
                      column(6,
                        wellPanel(
                          plotOutput("clustering", width="100%", height=400)
                        )
                      )
                    )
           ),
           # tabPanel("Clustering & MS",
           #   fluidRow(
           #     column(12,
           #       wellPanel(
           #         plotOutput("clustering", width="100%", height=400)
           #       )
           #     )
           #   )
           # ),
           tabPanel("Phenological Params",
                    fluidRow(
                      column(8,
                        wellPanel(
                          uiOutput("parameters_phenoParams", width="100%", height=150)
                        )
                      )
                    ),
                    fluidRow(
                      column(6,
                             wellPanel(
                               plotOutput("phenoParams_prelim", width="100%", height=600)
                             )
                      ),
                      column(6,
                        wellPanel(
                          plotOutput("phasePlot", width="100%", height=600)
                        )
                      )
                    )
           )
    )
  )

})

output$raw_data_prelim <- renderPlot({

  if(length(input$poligono_click$lng)==0){
    return() # ensures that plotting area is empty when app is loaded
  } else {

    isolate({
      par(mar=c(2,4,2,1))
      yRan <- range(pixelio$ndviMat)
      plot(pixelio$ndviMat[1,], col=1, type="l", ylab="NDVI", ylim=yRan)
      legend("topleft", legend=c(paste0("row: ", pixelio$ROW)), bty="n",
             cex=2)
      for(i in 2:nrow(pixelio$ndviMat)){
        lines(1:23, pixelio$ndviMat[i,], col=i)
      }
    })

  }
})

output$parameters_fpca <- renderUI({
  tabPanel("Smoothed & FPCA",
           fluidRow(
             column(1,
                    numericInput("numFreq", "Frequencies:", value=3,
                                 min=1, max=floor(nrow(pixelio$ndviMat)-1)/2, width=75)#,
                    # checkboxInput("getFPCA", label="FPCA fit", value=FALSE)
             ),
             column(1,
                    checkboxInput("getFPCA", label="FPCA fit", value=FALSE)
             )
           )
  )
})

output$smoothed_fpca_prelim <- renderPlot({

  if(length(input$numFreq)==0){
    return()
  } else {
    mat_smooth <- get_smoothing(m=pixelio$ndviMat, numFreq = input$numFreq)

    # comment made on March 1, 2022, should I use m=mat_pixel?
    mat_augSmooth <- get_harmonicFit(samples=100, m=mat_smooth,
                                     method="harmR", numFreq=input$numFreq,
                                     delta=0.1)

    # par(mar=c(2,4,2,1))
    # yRan <- range(mat_augSmooth)
    # plot(mat_augSmooth[,1], col=1, type="l", ylab="NDVI", ylim=yRan)
    # for(i in 2:22){
    #   lines(1:100, mat_augSmooth[,i], col=i)
    # }

    par(mar=c(2,4,2,1))
    yRan <- range(mat_augSmooth)
    plot(mat_augSmooth[,1], col=1, type="l", ylab="NDVI", ylim=yRan)
    for(i in 2:22){
      lines(1:100, mat_augSmooth[,i], col=i)
    }

    pixelio$matSmooth <- mat_smooth

  }

})

output$clustering <- renderPlot({

  if(length(input$numFreq)==0){
    return() # ensures that plotting area is empty when app is loaded
  } else {
    clust_DTWxstack <- tsclusters_centroid(pixelio$matSmooth, seed = 11,
                                           distance = "dtw2")

    cmd_distmat_DTW <- cmdscale(clust_DTWxstack@distmat)
    cmd_distmat_DTW <- as.data.frame(cmd_distmat_DTW)

    # clusterDTW <- as.factor(clust_DTWxstack@cluster)

    cmd_distmat_DTW$cluster <- as.factor(clust_DTWxstack@cluster)
    cmd_distmat_DTW$year <- c("'00","'01","'02","'03","'04","'05","'06",
                              "'07","'08","'09","'10","'11","'12","'13",
                              "'14","'15","'16", "'17", "'18", "'19", "'20",
                              "'21")

    names(cmd_distmat_DTW) <- c('V1','V2','cluster','year')
    esc_multidim_DTW <- NULL
    esc_multidim_DTW <- rbind(esc_multidim_DTW, cmd_distmat_DTW)

    # GRAFICA DE ESCALADO MULTIDIMENSIONAL
    ggplot(esc_multidim_DTW,
           aes(x=V1,y=V2,
               label=year,color=cluster)) +
      geom_text(aes(colour = cluster), size=5, alpha=0.7)+
      theme(legend.position = "none") +
      scale_color_manual(values = c("blue", "red")) +
      theme(plot.title = element_text(hjust = 0.5))
  }
})


output$parameters_phenoParams <- renderUI({
  tabPanel("Phenological Params",
           fluidRow(
             column(2,
                    numericInput("freqNum", "Frequencies:", value=3,
                                 min=3, max=4, width=100)
             )
           )
  )
})


output$phenoParams_prelim <- renderPlot({

  if(length(input$freqNum)==0 | length(input$poligono_click$lng)==0){
    return()
  } else {

    print(input$freqNum)

    mat_smooth <- get_smoothing(m=pixelio$ndviMat, numFreq = input$freqNum)

    mat_augSmooth <- get_harmonicFit(samples=100, m=mat_smooth,
                                     method="harmR", numFreq=input$freqNum,
                                     delta=0.1)

    # clust_L2xstack <- tsclusters_centroid(mat_smooth, seed=11, distance="L2")
    clust_DTWxstack <- tsclusters_centroid(mat_smooth, seed=11, distance="dtw_basic")

    # stack_fpca_l2 <- getFPCA_byPixels_v5(seriesByCluster=clust_L2xstack@cluster,
    #                                           mAug=mat_augSmooth)

    ndvi_fpca_dtw <- getFPCA_byPixels_v5(seriesByCluster=clust_DTWxstack@cluster,
                                               mAug=mat_augSmooth)

    fit <- haRmonics(y=ndvi_fpca_dtw$fpca, numFreq=input$freqNum, delta=0.1)

    pixelio$fit <- fit

    ndvi365 <- harmonicFit(amp=fit$amplitude, pha=fit$phase,
                           L=365, t=seq(1,365,length=365))

    pixelio$ndvi365 <- ndvi365

    # print(pixelio$fit$amplitude)

    if(input$freqNum==3){

      fprime <- firstDerivative(amp=pixelio$fit$amplitude,
                                pha=pixelio$fit$phase, L=365)

      fbiprime <- secondDerivative(amp=pixelio$fit$amplitude,
                                   pha=pixelio$fit$phase, L=365)

      fthrprime <- thirdDerivative(amp=pixelio$fit$amplitude,
                                   pha=pixelio$fit$phase, L=365)
    }

    if(input$freqNum==4){
      fprime <- firstDerivative_4(amp=pixelio$fit$amplitude,
                                pha=pixelio$fit$phase, L=365)

      fbiprime <- secondDerivative_4(amp=pixelio$fit$amplitude,
                                   pha=pixelio$fit$phase, L=365)

      fthrprime <- thirdDerivative_4(amp=pixelio$fit$amplitude,
                                   pha=pixelio$fit$phase, L=365)
    }

    pixelio$fprime <- fprime

    pixelio$fbiprime <- fbiprime

    pixelio$fthrprime <- fthrprime

    # ---
    # print(ndvi_fpca_dtw$fpca)

    phenoDates <- getPhenoParamsTemp(x=ndvi_fpca_dtw$fpca, numFreq=input$freqNum)

    # print(names(phenoDates))
    fec_tipo <- cbind(as.vector(phenoDates), names(phenoDates))%>%as.data.frame()
    names(fec_tipo) <- c("Fecha", "Tipo")
    param_list <- list(params=fec_tipo, number_params=phenoDates)

    print(param_list)
    #---

    fechas <- c()
    tipos <- c()
    if(is.factor(param_list$params$Fecha)){
      fechas <- as.vector(param_list$params$Fecha)
      tipos <- as.vector(param_list$params$Tipo)
    } else {
      fechas <- param_list$params$Fecha
      tipos <- param_list$params$Tipo
    }

    # ---

    pixelio$ndviSmoothedSampled <- mat_augSmooth

    par(mar=c(2,4,2,1), mfrow=c(3,1))
    yRan <- range(pixelio$ndviSmoothedSampled)

    plot(pixelio$ndvi365, col="red", type="l", lwd=4, ylab="NDVI", ylim=yRan,
         xaxt="n", cex.lab=1.25, font.lab=4)
    axis(side=1, at=DAYS_AT_YEAR, labels=c(MONTHS,""))
    abline(v=fechas, lwd=3, col=colores)

    yRan <- range(fprime(seq(1,365,length=365)))
    yRan[1] <- yRan[1]-0.0025
    yRan[2] <- yRan[2]+0.0025
    plot(fprime(seq(1, 365, length=365)), ylim=yRan, ylab="NDVI'", xaxt="n",
         lwd=4, col="blue", cex.lab=1.25, font.lab=4)
    axis(side=1, at=DAYS_AT_YEAR, labels=c(MONTHS,""))
    abline(v=fechas, lwd=3, col=colores)
    abline(h=0, lty=3, col="blue")


    yRan <- range(fbiprime(seq(1,365,length=365)))
    yRan[1] <- yRan[1]-0.00005
    yRan[2] <- yRan[2]+0.00005
    plot(fbiprime(seq(1, 365, length=365)), ylim=yRan, ylab="NDVI''", xaxt="n",
         lwd=4, col="magenta", cex.lab=1.25, font.lab=4)
    axis(side=1, at=DAYS_AT_YEAR, labels=c(MONTHS,""))
    abline(v=fechas, lwd=3, col=colores)
    abline(h=0, lty=3, col="blue")

  }

})

output$phasePlot <- renderPlot({
  if(length(input$poligono_click$lng)==0 | length(input$freqNum)==0){
    return()
  } else {

    par(mfrow=c(2,2))

    plot(pixelio$ndvi365,
         pixelio$fprime(t=seq(1, 365, length=365)),
         type="p", xlab="NDVI", ylab="NDVI'")
    abline(v=0, lty=3, col="blue")
    abline(h=0, lty=3, col="blue")

    plot(pixelio$fprime(t=seq(1, 365, length=365)),
         pixelio$fbiprime(t=seq(1, 365, length=365)),
         type="p", xlab="NDVI'", ylab="NDVI''")
    abline(v=0, lty=3, col="blue")
    abline(h=0, lty=3, col="blue")

    plot(pixelio$fprime(t=seq(1, 365, length=365)),
         pixelio$fthrprime(t=seq(1, 365, length=365)),
         type="p", xlab="NDVI'", ylab="NDVI'''")
    abline(v=0, lty=3, col="blue")
    abline(h=0, lty=3, col="blue")

    plot(pixelio$fbiprime(t=seq(1, 365, length=365)),
         pixelio$fthrprime(t=seq(1, 365, length=365)),
         type="p", xlab="NDVI''", ylab="NDVI'''")
    abline(v=0, lty=3, col="blue")
    abline(h=0, lty=3, col="blue")
    
  }

})

