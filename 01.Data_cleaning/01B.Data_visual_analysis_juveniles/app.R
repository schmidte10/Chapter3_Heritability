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
  mutate(across(1:11, factor), 
         Temperature = factor(Temperature), 
         True_resting = factor(True_resting))

df_jresp2_rest <- df_jresp2 |> 
  filter(True_resting == "Y")

#--- test code for running figure ---# 

df_jresp2 |> 
  filter(True_resting == "Y") |>
  ggplot(aes(x=Clutch, y=Resting, color=Resting)) +
  geom_point() +
  scale_colour_gradient(low="black", high="red")+
  theme_classic()





# Define UI for application that draws a histogram
ui <- page_navbar(
  title="Respiration (juveniles) Data Quailty Checks", 
  bg="#2D89C8", 
  inverse=TRUE, 
  
  nav_panel(title="Resting oxygen uptake", 
              plotOutput("plot_clutch"),
              plotOutput("plot_mass")),
  nav_panel(title="Maximum oxygen uptake", 
            plotOutput("plot_clutch_mmr"),
            plotOutput("plot_mass_mmr")),
  nav_panel(title="Absolute aerobic scope", 
            plotOutput("plot_clutch_aas"),
            plotOutput("plot_mass_aas")),
  
   

    
    sidebarLayout(
        sidebarPanel(
          div(style="display:inline-block", 
              
              selectInput( 
                inputId = "F0", 
                label = "Unique parentage:", 
                choices = levels(df_jresp2$F0)
                ))
          
        ),

       
        navset_card_underline(title="Data quality check")
           
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
            geom_point() + 
            theme_classic()
        })
          }) 
        
    output$plot_mass <- renderPlot({ 
      input$F0 
      isolate({ 
        ggplot(data=Datafinder_rest(), aes(x=Dry_mass, y=Resting, color=Clutch)) + 
          geom_point(data=df_jresp2, color="grey66", alpha=0.2) + 
          geom_point() + 
          theme_classic()
        })
      }) 
    
    output$plot_clutch_mmr <- renderPlot({
      input$F0 
      isolate({ 
        ggplot(data= Datafinder2(), aes(x=Clutch, y=Max, color=Clutch)) + 
          geom_point(data=df_jresp2, color="grey66", alpha=0.2) +
          geom_point() + 
          theme_classic()
      })
    }) 
    
    output$plot_mass_mmr <- renderPlot({ 
      input$F0 
      isolate({ 
        ggplot(data=Datafinder2(), aes(x=Dry_mass, y=Max, color=Clutch)) + 
          geom_point(data=df_jresp2, color="grey66", alpha=0.2) + 
          geom_point() + 
          theme_classic()
      })
    })

    output$plot_clutch_aas <- renderPlot({
      input$F0 
      isolate({ 
        ggplot(data= Datafinder2(), aes(x=Clutch, y=AAS, color=Clutch)) + 
          geom_point(data=df_jresp2, color="grey66", alpha=0.2) +
          geom_point() + 
          theme_classic()
      })
    }) 
    
    output$plot_mass_aas <- renderPlot({ 
      input$F0 
      isolate({ 
        ggplot(data=Datafinder2(), aes(x=Dry_mass, y=AAS, color=Clutch)) + 
          geom_point(data=df_jresp2, color="grey66", alpha=0.2) + 
          geom_point() + 
          theme_classic()
      })
    })
        
    }

# Run the application 
shinyApp(ui = ui, server = server)
