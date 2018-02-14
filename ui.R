library(leaflet)
library(shinyjs)
library(shinyBS)
library(shinyCustom)
library(shinydashboard)
library(shinyWidgets)
library(tidyverse)
library(rvest)

# last commit date
dt <- read_html('https://github.com/SCCWRP/SGRRMP/commits/master') %>% 
  html_nodes(".commit-group-title") %>% 
  html_text %>% 
  .[1] %>% 
  gsub('^.*Commits on (.*)\\n.*$', '\\1', .)

# watersheds to select from data folder
shds <- list.files('data') %>% 
  gsub('\\.RData$|^scrs_|^spat_', '', .) %>% 
  unique %>% 
  tolower

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  theme = 'styles.css',
  useShinyjs(),
  useShinyCustom(slider_delay = '1500'),
  
  # Application title
  h1('Stream classification in the Santa Margarita Watershed', 
     style = "font-family: 'Volkhov';
        font-weight: 500; line-height: 1.1"),
  
  fluidRow(
  
    column(width = 3, img(src = "sccwrp_logo.jpg", width = '150px'), align = 'center', style = "margin-top: 0px;"),
    
    column(width = 9, 
           h5('This application can be used to explore stream and site classifications for the Santa Margarita Watershed.  Classications are based on the relationship of field CSCI scores at a site to biological expectations for the stream reach.  Expectations are based on user-defined parameters for CSCI thresholds and confidence in the biological expectation. Site classifications for CSCI scores are defined as over-performing, expected, and under-performing.  Stream reach expectations are defined as likely constrained, undetermined, or likely unconstrained. Last updated:', dt)
    ),

    column(width = 12, 
      h4('Created by Marcus W. Beck,', a('marcusb@sccwrp.org', href = 'mailto:marcusb@sccwrp.org'), ", Raphael D. Mazor,", a('raphaelm@sccwrp.org', href = 'mailto:raphaelm@sccwrp.org'))
      ), 
    
    column(width = 4, 
           selectInput('shd',
                       label = h6("Select watershed:"), 
                       choices = shds, 
                       selected = shds[2]
           )
      )

  ),
    
  tabsetPanel(
    
    tabPanel('Maps',

      h5('These two maps show stream reach classifications by COMID and CSCI scores at monitoring stations.  The', strong('left map'), 'shows the predicted CSCI score for a COMID and measured CSCI score at a station from from field data.  The', strong('right map'), 'shows the CSCI score expectation for a COMID and the site classification (or performance) of a monitoring station.'),   

      column(width = 12, 
             
        h5('These sliders control the aesthetics in both maps. Use them to change the point/line sizes and apply a jitter for repeat visits as the same station.  The average CSCI across all samples at a size is shown if jitter is set to zero.'), 
        
        # select point radius
        column(width = 4,
              customSliderInput("pt_sz", 
                          label = h6("Point size:"), 
                          min = 0, 
                          max = 15,
                          value = 4, 
                          step = 1, 
                          width = '400px', 
                          ticks = FALSE
              )
        ),
        
        # select line size
        column(width = 4,
              customSliderInput("ln_sz", 
                          label = h6("Line size:"), 
                          min = 0, 
                          max = 5,
                          value = 1, 
                          step = 0.1, 
                          width = '400px', 
                          ticks = FALSE
              )
        ),          
        
        column(width = 4,
               
               customSliderInput('jitr', 
                           label = h6("Jitter overlaps:"), 
                           min = 0, 
                           max = 500,
                           value = 0, 
                           step = 25, 
                           width = '400px', 
                           ticks = FALSE
               )
               
        )
        
      ),
         
      column(width = 6, 
        
        h5('The toggle switch for the', strong('left map'), 'controls how the CSCI scores at the stations (points) are displayed.  The observed scores from field samples are shown when the switch is off and the differences between the observed scores and the stream reach predictions are shown when the switch is on.'),      

        # show csci differences   
        materialSwitch('difr', 
          label = h6('CSCI observed - predicted:'), 
          status = 'primary',
          right = F
        )
        
      ), 
      
      column(width = 6,
          
        h5('These controls change the attributes in the',  strong('right map'), '. The first slider controls the CSCI threshold and the second slider controls the certainty range of the predicted CSCI scores at each stream reach. Overlap of the certainty range with the CSCI threshold determines the expectation of a reach and performance of the CSCI score at a sampling station. See the plot tab for more.'),      
             
        # select CSCI threshold, master
        sliderTextInput(
          inputId = "thrsh",
          label = h6("CSCI reference threshold:"),
          grid = FALSE,
          force_edges = TRUE,
          selected = '10% (0.79)',
          choices = c('1% (0.63)', '10% (0.79)', '30% (0.89)'),
          width = '600px'
        ),

        # selected tails, master
        sliderTextInput(
          inputId = "tails", 
          label = h6("Confidence range (+/-):"),  
          grid = FALSE, 
          force_edges = TRUE,
          choices = c('More certain (0.45)', '0.40', '0.35', '0.30', '0.25', '0.20', '0.15', '0.10', 'Less certain (0.05)'), 
          width = '600px'
        )

      ),
    
      # map output
      column(width = 6,
        
        leafletOutput('map', width = '100%', height = 550), 
        h3()
             
      ),
     
      # map output
      column(width = 6,
             
        leafletOutput('map_exp', width = '100%', height = 550), 
        h3()

      ) 

    ),
    
    tabPanel('Plot',
      
      h5('This plot shows the CSCI score expectations for every stream reach with CSCI sampling stations.  The CSCI threshold and confidence range define the reach expectation and the CSCI performance for the sampling stations.  Toggle the sliders to see how these change on the plot, including the maps and table in the other tabs.'),
             
      column(width = 4,
             
        # select CSCI threshold
        sliderTextInput(
          inputId = "thrsh2", 
          label = h6("CSCI reference threshold:"),  
          grid = FALSE, 
          force_edges = TRUE,
          selected = '10% (0.79)',
          choices = c('1% (0.63)', '10% (0.79)', '30% (0.89)'), 
          width = '600px'
        )
        
      ),   
      
      column(width = 4, 
             
        # selected tails
        sliderTextInput(
          inputId = "tails2", 
          label = h6("Confidence range (+/-):"),  
          grid = FALSE, 
          force_edges = TRUE,
          choices = c('More certain (0.45)', '0.40', '0.35', '0.30', '0.25', '0.20', '0.15', '0.10', 'Less certain (0.05)'), 
          width = '400px'
        )
        
      ),
        
      column(width = 2,
             
        # order by site
        materialSwitch('bysta', 
                       label = h6('Order by site:'), 
                       status = 'primary',
                       right = F
        )

      ),
    
      column(width = 2,
             
             # order by site
             materialSwitch('nocon', 
                            label = h6('No context:'), 
                            status = 'primary',
                            right = F
             )
             
      ),
      
      # plot output
      column(width = 12,
            
        plotOutput('plo_exp', width = '90%', height = 850)
            
      ) 
             
    ), 
        
    tabPanel('Table', 

      h5('This table summarizes the sampling station performance for CSCI scores shown in the maps and plot in the other tabs. The "types" are finer divisions that further categorize sites relative to the performance.'),
             
      column(width = 6,
            
            # select CSCI threshold
            sliderTextInput(
             inputId = "thrsh3", 
             label = h6("CSCI reference threshold:"),  
             grid = FALSE, 
             force_edges = TRUE,
             selected = '10% (0.79)',
             choices = c('1% (0.63)', '10% (0.79)', '30% (0.89)'), 
             width = '600px'
            )

      ),   
      
      column(width = 6, 
          
            # selected tails
            sliderTextInput(
              inputId = "tails3", 
              label = h6("Confidence range (+/-):"),  
              grid = FALSE, 
              force_edges = TRUE,
              choices = c('More certain (0.45)', '0.40', '0.35', '0.30', '0.25', '0.20', '0.15', '0.10', 'Less certain (0.05)'), 
              width = '600px'
            )
            
      ),

      # table output
      column(width = 12, 
      
        DT::dataTableOutput('tab_sum'), 
        HTML('<p></p>')
                   
      )
    
    )
                  
  )
      
))
            

