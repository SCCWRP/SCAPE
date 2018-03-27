library(leaflet)
library(shinyjs)
library(shinyBS)
library(shinyCustom)
library(shinydashboard)
library(shinyWidgets)
library(ShinyDash)
library(tidyverse)
library(mapview)
library(rvest)

# last commit date
dt <- read_html('https://github.com/SCCWRP/Biological_Constraints/commits/master') %>% 
  html_nodes(".commit-group-title") %>% 
  html_text %>% 
  .[1] %>% 
  gsub('^.*Commits on (.*)\\n.*$', '\\1', .)

# watersheds to select from data folder
shds <- list.files('data') %>% 
  gsub('\\.RData$|^scrs_|^spat_', '', .) %>% 
  unique

# column padding global
pad <- 'padding:0px;'
  
# Define UI for application
shinyUI(fluidPage(
  
  theme = 'styles.css',
  useShinyjs(),
  useShinyCustom(slider_delay = '1500'),
  
  # Application title
  h1(HTML('<h1>Landscape constraints on stream biological integrity</h1>'), 
     style = "font-family: 'Volkhov';
     font-weight: 500; line-height: 1.1"),
  
  fluidRow(
    
    column(width = 12, 
           h4('Created by Marcus W. Beck,', a('marcusb@sccwrp.org', href = 'mailto:marcusb@sccwrp.org'), ", Raphael D. Mazor,", a('raphaelm@sccwrp.org', href = 'mailto:raphaelm@sccwrp.org'))
    ),
    
    column(width = 3, img(src = "sccwrp_logo.jpg", width = '150px'), align = 'center', style = "margin-top: 0px;"),
    
    column(width = 9, 
           h5(HTML('This application can be used to explore landscape constraints on biological integrity of streams.  The application provides context for evaluating stream health by estimating an expectation of biological condition at a given stream reach relative to landscape drivers. Biological condition from field data using the California Stream Condition Index can be compared to the reach expectation.  The process begins by identifying stream classifications and expectations from user-defined parameters for CSCI thresholds and confidence in the biological expectation. Stream classifications of expected biological constraints are defined as likely unconstrained, possibly unconstrained, possibly constrained, and likely constrained.  Observed CSCI scores at a site are then characterized relative to the reach expectations as over-scoring, expected, or under-scoring.  Relative site scores given the expectation can be used to recommend priorities for management actions. View source <a href="https://github.com/SCCWRP/SGRRMP">here</a>. Last updated:'), dt)
    ), 
    
    column(width = 4, 
           selectInput('shd',
                       label = h6("Select watershed:"), 
                       choices = shds, 
                       selected = shds[2]
           )
    )
    
  ),
  
  # master widgets    
  column(width = 12, 
         
         h4('Controlling the app:'),
         
         h5('These controls determine how stream expectations and relative site scores are evaluated. The first slider controls the CSCI threshold and the second slider controls the certainty range of the expected CSCI scores at each stream reach. Overlap of the certainty range with the CSCI threshold determines the expectation of a reach and the relative CSCI score at a station (see step 2). The third switch determines if results are averaged for each site across all visits, or if results for all visits are shown.  Turning the switch to the right will jitter repeat visits on the maps and all tabular summaries (step 4 and 5) will consider each visit as a unique event.'),
         
         # select CSCI threshold, master
         column(width = 4,    
                sliderTextInput(
                  inputId = "thrsh",
                  label = h6("CSCI reference threshold:"),
                  grid = FALSE,
                  force_edges = TRUE,
                  selected = '10% (0.79)',
                  choices = c('1% (0.63)', '10% (0.79)', '30% (0.89)'),
                  width = '400px'
                )
         ),
         
         # selected tails, master
         column(width = 4,
                sliderTextInput(
                  inputId = "tails", 
                  label = h6("Confidence range (+/-):"),  
                  grid = FALSE, 
                  force_edges = TRUE,
                  choices = c('More certain (0.45)', '0.40', '0.35', '0.30', '0.25', '0.20', '0.15', '0.10', 'Less certain (0.05)'), 
                  width = '400px'
                )
         ), 
         
         # apply jitter 
         column(width = 4, 
                # jitr switch   
                materialSwitch('jitr', 
                               label = h6(HTML('Show individual samples at each site:<br/><br/></br>')), 
                               status = 'primary',
                               right = F, 
                               width = '400px'
                )
                
         )
         
  ),
  
  tabsetPanel(id = 'alltabs', 
              
              tabPanel('(1) View maps',
                       
                       h5('These maps show stream reach classifications and CSCI scores at monitoring stations.  The', strong('left map'), 'shows the predicted median CSCI score for a reach and observed CSCI score at a station from field data.  The', strong('right map'), 'shows the CSCI score expectation for a reach and the relative CSCI score at a station for the expectation. See the plot tab (step 2) for more details on how expectations and relative site scores are determined. The toggle switch controls how the CSCI scores at the stations (points) on the left map are displayed.  The observed scores from field samples are shown when the switch is off and the differences between the observed scores and the stream reach median expectations are shown when the switch is on.'),
                       
                       # show csci differences   
                       materialSwitch('difr', 
                                      label = h6('CSCI observed - predicted:'), 
                                      status = 'primary',
                                      right = F
                       ),
                       
                       # map output
                       column(width = 6,
                              
                              leafletOutput('map_med', width = '100%', height = 550), 
                              h3()
                              
                       ),
                       
                       # map output
                       column(width = 6,
                              
                              leafletOutput('map_exp', width = '100%', height = 550), 
                              h3()
                              
                       ) 
                       
              ),
              
              tabPanel('(2) View reach summary',
                       
                       h5('This plot shows the range of CSCI score expectations for every stream reach with CSCI sampling stations.  The CSCI threshold and confidence range define the reach expectation and the relative CSCI score for the sampling stations.  The median for the expected range of CSCI scores at a reach is shown as a white tick. Toggle the sliders to see how these change on the plot, including the maps and table in the other tabs.'),
                       
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
              
              tabPanel('(3) Tabulate reach summary', 
                       
                       h5('This table summarizes the sampling stations for the relative CSCI scores shown in the maps and plot in steps 1 and 2. The "types" are finer divisions that further categorize sites the relative score and CSCI threshold.  The types are based on relative score and location relative to the selected CSCI threshold. The types can be used to recommend priorities for management actions in step 4.'),
                       
                       # table output
                       column(width = 12, 
                              
                              DT::dataTableOutput('tab_sum'), 
                              HTML('<p></p>')
                              
                       )
                       
              ),
              
              tabPanel('(4) Set reach priorities',
                       
                       h5('This plot can be used to identify potential management actions for each site type.  The plot on the right shows a graphical depiction of the types in the able in step 3.  Priorities for each type can be selected from none to many using the drop-down menus for each type on the left.  These priorities can then be viewed with the maps in step 5.'),
                       
                       h5('The default priorities were based on recommendations from a stakeholder group with familiarity of the watershed.  The priorities are generalized into three categories to recommend actions in addition to baseline monitoring and maintenance:'),
                       
                       HTML('<ul>
                            <li><strong>Investigate</strong>: Additional monitoring or review of supplementary data (e.g., aerial imagery)</li>
                            <li><strong>Protect</strong>: Additional scrutiny of proposed development and/or projects</li>
                            <li><strong>Restore</strong>: Targeted action for causal assessment and/or restoration funds</li>
                            </ul>'),
                       
                       # plot output legend
                       column(width = 12,
                              
                              plotOutput('plo_leg', width = '100%', height = 100)
                              
                       ),
                       
                       # site priority selectors
                       column(width = 2,
                              
                              div(style = 'padding:11px;'),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 1", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = c('Investigate', 'Protect'),
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 2", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = NULL,
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 3", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = c('Investigate'),
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 4", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = c('Investigate', 'Restore'),
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 5", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = c('Investigate', 'Protect'),
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 6", label = NULL, choices = c('Protect', 'Investigate', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = NULL,
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 7", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = c('Investigate'),
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 8", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = c('Investigate', 'Restore'),
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 9", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = c('Protect'),
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 10", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = NULL,
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 11", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = NULL,
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 12", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = NULL,
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 13", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = c('Protect'),
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 14", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = NULL,
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 15", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = NULL,
                                              multiple = TRUE
                                  )
                              ),
                              
                              div(style = pad,
                                  pickerInput(inputId = "Site 16", label = NULL, choices = c('Investigate', 'Protect', 'Restore'),
                                              options = list(`actions-box` = TRUE, size = 20), selected = NULL,
                                              multiple = TRUE
                                  )
                              )
                              
                       ),
                       # plot output
                       column(width = 10,
                              
                              plotOutput('plo_exp2', width = '100%', height = 900)
                              
                       )
                       
                       ),
              
              tabPanel('(5) View priorities',
                       
                       h5("These maps show the location of recommended priority actions defind for each site type in step 4. Each site can have more than one priority."),
                       
                       # investigate map
                       column(width = 4,
                              
                              htmlWidgetOutput(outputId = 'cnts_inv', HTML('<h3>Investigate: <b><span id="Investigate"></span></b></h3>')),
                              leafletOutput('bs_inv', width = '100%', height = 550),
                              h3()
                              
                       ),
                       
                       # protect map
                       column(width = 4,
                              
                              htmlWidgetOutput(outputId = 'cnts_pro', HTML('<h3>Protect: <b><span id="Protect"></span></b></h3>')),
                              leafletOutput('bs_pro', width = '100%', height = 550),
                              h3()
                              
                       ),
                       
                       # restore map
                       column(width = 4,
                              
                              htmlWidgetOutput(outputId = 'cnts_res', HTML('<h3>Restore: <b><span id="Restore"></span></b></h3>')),
                              leafletOutput('bs_res', width = '100%', height = 550),
                              h3()
                              
                       )
                       
              )
              
              )
  
  ))



