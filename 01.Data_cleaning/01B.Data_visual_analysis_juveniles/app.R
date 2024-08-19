#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#--- import libraries ---#
library(readr)
library(tidyverse)
library(shinythemes)
library(shiny)
library(forcats)
library(bslib)
#--- import data ---# 
df_jresp <- read_csv("resp_results_juveniles.csv")

df_jresp$Population <-  fct_collapse(df_jresp$Population, 
                                      `Vlassof cay`= c("Vlassof reef", "Vlassof", "Vlassof Cay", "Vlassof cay"), 
                                      `Arlington reef` = c("Arlington reef","Arlginton reef")) 

df_jresp$Female <-  fct_collapse(df_jresp$Female, 
                                  `CARL359`= c("CARL359", "CARL59")) 



#--- data manipulation ---# 
df_jresp2 <-  df_jresp |> 
  unite("F0", c("Male","Female"), sep="_", remove=FALSE) |>
  mutate(across(1:7, factor), 
         Temperature = factor(Temperature), 
         True_resting = factor(True_resting)) 

df_jresp2_rest <- df_jresp2 |> 
  filter(True_resting == "Y")


# Define UI for application that draws a histogram
ui <- fluidPage( 
  theme = shinytheme("spacelab"),
  headerPanel("Data visualisation"), 
  tabsetPanel( 
    tabPanel('Resting oxygen uptake', 
             sidebarLayout( 
               sidebarPanel(
                 div(style="display:inline-block", 
                     selectInput( 
                       inputId = "F0", 
                       label = "Unique parentage:", 
                       choices = levels(df_jresp2$F0), 
                       width= 200 
                       
                     ))),
                 
               
             mainPanel(title="Resting oxygen uptake", 
                       plotOutput("plot_clutch"),
                       plotOutput("plot_mass"))
               
               )
             ), 
    tabPanel('Maximum oxygen uptake', 
             sidebarLayout( 
               sidebarPanel(
                 div(style="display:inline-block", 
                     selectInput( 
                       inputId = "F0", 
                       label = "Unique parentage:", 
                       choices = levels(df_jresp2$F0), 
                       width= 200 
                       
                     ))),
               
               
               mainPanel(title="Maximum oxygen uptake", 
                         plotOutput("plot_clutch_mmr"),
                         plotOutput("plot_mass_mmr"))
               
             )
    ), 
    tabPanel('Absolute aerobic scope', 
             sidebarLayout( 
               sidebarPanel(
                 div(style="display:inline-block", 
                     selectInput( 
                       inputId = "F0", 
                       label = "Unique parentage:", 
                       choices = levels(df_jresp2$F0), 
                       width= 200 
                       
                     ))),
               
               
               mainPanel(title="Absolute aerobic scope", 
                         plotOutput("plot_clutch_aas"),
                         plotOutput("plot_mass_aas"))
               
             )
    )
    )
  )



# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  Datafinder_rest <- reactive({ 
    req(input$F0)
    filter(df_jresp2_rest, 
           F0 %in% input$F0
           )
    }) 
  
  Datafinder2 <- reactive({ 
    req(input$F0)
    filter(df_jresp2, 
           F0 %in% input$F0
    )
  })

    output$plot_clutch <- renderPlot({
        input$F0 
        isolate({ 
          ggplot(data= Datafinder_rest(), aes(x=Clutch, y=Resting, color=Clutch)) + 
            geom_point(data=df_jresp2, color="grey66", alpha=0.2) +
            geom_point(size=3) + 
            theme_classic()
        })
          }) 
        
    output$plot_mass <- renderPlot({ 
      input$F0 
      isolate({ 
        ggplot(data=Datafinder_rest(), aes(x=Dry_mass, y=Resting, color=Clutch)) + 
          geom_point(data=df_jresp2, color="grey66", alpha=0.2) + 
          geom_point(size=3) + 
          scale_x_continuous(limits = c(0, 0.001), breaks = seq(0, 0.001, 0.0001)) +
          theme_classic()
        })
      }) 
    
    output$plot_clutch_mmr <- renderPlot({
      input$F0 
      isolate({ 
        ggplot(data= Datafinder2(), aes(x=Clutch, y=Max, color=Clutch)) + 
          geom_point(data=df_jresp2, color="grey66", alpha=0.2) +
          geom_point(size=3) + 
          theme_classic()
      })
    }) 
    
    output$plot_mass_mmr <- renderPlot({ 
      input$F0 
      isolate({ 
        ggplot(data=Datafinder2(), aes(x=Dry_mass, y=Max, color=Clutch)) + 
          geom_point(data=df_jresp2, color="grey66", alpha=0.2) + 
          geom_point(size=3) + 
          theme_classic()
      })
    })

    output$plot_clutch_aas <- renderPlot({
      input$F0 
      isolate({ 
        ggplot(data= Datafinder2(), aes(x=Clutch, y=AAS, color=Clutch)) + 
          geom_point(data=df_jresp2, color="grey66", alpha=0.2) +
          geom_point(size=3) + 
          theme_classic()
      })
    }) 
    
    output$plot_mass_aas <- renderPlot({ 
      input$F0 
      isolate({ 
        ggplot(data=Datafinder2(), aes(x=Dry_mass, y=AAS, color=Clutch)) + 
          geom_point(data=df_jresp2, color="grey66", alpha=0.2) + 
          geom_point(size=3) + 
          theme_classic()
      })
    })
        
    }

# Run the application 
shinyApp(ui = ui, server = server)
