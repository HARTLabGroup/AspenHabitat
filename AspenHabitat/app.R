library(tidyverse)
library(shiny)

library(sf)
library(leaflet)

# Create the UI, where the visual appearance of the app will be stored.
ui <- fluidPage(
  titlePanel("Aspen Habitat"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId ="Time",
                  label = "Choose time period",
                  choices = c("1981-2010", "2011-2040", "2041-2070", "2071-2011"),
                  selected = NULL
      ),
      selectInput(inputId ="SSP",
                  label = "Choose Emissions Scenario",
                  choices = c("historical", "SSP1-2.6", "SSP3-7.0", "SSP5-8.5"),
                  selected = NULL
      )  
    ),
    mainPanel(
      leafletOutput(outputId = 'map')
    )
  )
)
# Create the server, which will contain the code used.
server = function(input, output){

}

# Run the application itself
shinyApp(ui, server)