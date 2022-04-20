
library(shiny)
library(shinydashboard)
library(leaflet)
library(leafpop)
library(mapview)
library(spiralize)
library(vcd)
library(geoTS)
library(raster)

library(pryr)
# library(geoTS)
library(dtwclust)
library(TSclust)
library(FPCA)
library(TSA)
library(ggplot2)
library(foreach)
library(doParallel)
library(gtools)
library(eBsc)

library(ggplot2)

library(plyr)
library(dplyr)
library(rlang)
library(stringr)
library(rootSolve)

# source("C:/Users/inder/OneDrive/Desktop/proyectoFanny/Rscripts_sTBDF/auxFunctions.R")
# source("C:/Users/inder/OneDrive/Desktop/shinyRStudio/phenoParams/auxiliaryObjects.R")
source(paste0(getwd(), "/auxFunctions.R"))
source(paste0(getwd(), "/auxObjects.R"))

ui <- dashboardPage(
  dashboardHeader(title = 'phenoParams'),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Pre-Analysis", tabName = "prelim"),
      # menuItem("Visualization", tabName = "visualization"),
      menuItem("Contact us", tabName = "contactUs")
    )
  ),
  
  dashboardBody(
    tabItems(
      
      tabItem("prelim",
        
        fluidRow(
          column(2,
            wellPanel(
              selectInput("polygonSelect",
                "Polygon:", width=175,
                choices=c("Santiago Amoltepec", "Boqueron de Tonala", 
                          "Tehuacan-Culcatlan")
              )
            )
          )
        ),
        
        fluidRow(
          column(12,
            wellPanel(
              leafletOutput("poligono", width="100%", height=500)
            )
          )
        ),
        
        fluidRow(
          column(12,
            # wellPanel(
              uiOutput("test")
            # )
          )
        )
      ),
      
      # tabItem("visualization",
      #         
      #         fluidRow(
      #           column(3,
      #             # wellPanel(
      #               selectInput("selectPolygon",
      #                 "Polygon:", width=175, 
      #                 choices=c("Santiago Amoltepec", "Boquerón de Tonalá", "Tehuacán-Culcatlán")
      #               )
      #             # )
      #           )
      #         ),
      #         
      #         fluidRow(
      #           column(7,
      #                  wellPanel(
      #                    leafletOutput("polygon", width="100%", height=500)
      #                  )
      #                  
      #           ),
      #           column(5,
      #                  wellPanel(
      #                    plotOutput("spiral", width="100%", height=500)
      #                  )
      #                  
      #           )
      #         ),
      #         
      #         fluidRow(
      #           column(12,
      #                  # wellPanel(
      #                    uiOutput("the_three_plots", inline=TRUE)
      #                  # )
      #           )
      #         )
      # ),
      
      tabItem("contactUs",
              fluidRow(
                column(8,
                       p("For comments on this shiny app please contact",
                         a(href="https://irt466.wixsite.com/inder", "Inder Tecuapetla"),
                         "(itecuapetla@conabio.gob.mx)"
                       )
                )
              )
      )
    )
  )
)
