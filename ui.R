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
dt <- read_html('https://github.com/SCCWRP/SCAPE/commits/master') %>% 
  html_nodes(".commit-group-title") %>% 
  html_text %>% 
  .[1] %>% 
  gsub('^.*Commits on (.*)\\n.*$', '\\1', .)
# # watersheds to select from data folder
# shds <- list.files('data') %>% 
#   gsub('\\.RData$|^scrs_|^spat_', '', .) %>% 
#   unique

# column padding global
pad <- 'padding:0px;'
  
# Define UI for application
shinyUI(fluidPage(
  
  theme = 'styles.css',
  useShinyjs(),
  useShinyCustom(slider_delay = '1500'),
  
  # Application title
  h1(HTML('<h1><b>SCAPE</b>: <b>S</b>tream <b>C</b>lassification <b>A</b>nd <b>P</b>riority <b>E</b>xplorer</h1>'), 
     style = "font-family: 'Volkhov';
     font-weight: 500; line-height: 1.1"),
  
  fluidRow(
    
    column(width = 3, img(src = "logo2.jpg", width = '450px'), align = 'center', style = "margin-top: 0px;")
    
    # column(width = 4, 
    #        selectInput('shd',
    #                    label = h6("Select watershed:"), 
    #                    choices = shds, 
    #                    selected = shds[2]
    #        )
    # )
    
  ),
  
  HTML('<br></br>'),
  
  # master widgets    
  column(width = 12, 

         
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
  
  tabsetPanel(id = 'alltabs', type = 'pills', 
              
              tabPanel('Instructions',
                       
                       tabsetPanel(#type = 'pills',
                       
                                   tabPanel('Overview',
                       
                                     column(width = 12, 
                                            
                                            h5(HTML(' View source <a href="https://github.com/SCCWRP/SGRRMP">here</a>. Last updated:'), dt),
                                            
                                            h5('The SCAPE application can be used to explore landscape constraints on biological integrity of streams.  The application provides context for evaluating stream health by estimating an expectation of biological condition at a given stream reach relative to landscape drivers. Biological condition from field data using the California Stream Condition Index can be compared to the reach expectation.  The process begins by identifying stream classifications and expectations from user-defined parameters for CSCI thresholds and confidence in the biological expectation. Stream classifications of expected biological constraints are defined as likely unconstrained, possibly unconstrained, possibly constrained, and likely constrained.  Observed CSCI scores at a site are then characterized relative to the reach expectations as over-scoring, expected, or under-scoring.  Relative site scores given the expectation can be used to recommend priorities for management actions.')
                                            
                                      )
                                     
                                     ),
                                   
                                   tabPanel('Controlling the app',
                                     
                                     column(width = 12, 
                                            
                                            h5(HTML('Changing the master controls at the top will propogate changes to each map, figure, and table. These controls determine how stream expectations and relative site scores are evaluated by SCAPE. The first slider controls the <b>CSCI reference threshold</b> and the second slider controls the <b>certainty range</b> of the expected CSCI scores at each stream reach. Overlap of the certainty range with the CSCI threshold determines the expectation of a reach and the relative CSCI score at a station (see step 2). The third switch determines if results are averaged for each site across all visits, or if <b>individual samples</b> for repeat visits to the same site are shown.  Turning the switch to the right will jitter repeat visits on the maps and all tabular summaries (step 4 and 5) will consider each visit as a unique event.')),
                                            
                                            h5('Follow these steps to use SCAPE:'),
                                            
                                            h5(HTML('<b>(1) View maps</b> that show stream reach classifications and CSCI scores at monitoring stations.')),
                                            
                                            h5(HTML('<b>(2) View reach summary</b> of CSCI score expectations for every stream reach with CSCI sampling stations.')), 
                                            
                                            h5(HTML('<b>(3) Tabulate reach summary</b> of CSCI score expectations for every stream reach with CSCI sampling stations.')), 
                                            
                                            h5(HTML('<b>(4) Set reach priority</b> to identify potential management actions for each site type. ')), 
                                            
                                            h5(HTML('<b>(5) View priorities</b> by location of recommended actions defined for each site type.'))
                                            
                                            )
                                   ), 
                                   
                                   tabPanel('Definitions',
                       
                                     column(width = 12, 
                                            
                                            h5(HTML('<b>Classification</b>: Determining which streams are likely to be constrained and predicting the range of likely scores within those constraints.')),
                                            
                                            h5(HTML('<b>Constrained</b>: The lack of potential for a stream to support good biological condition. This is based on a defined expectation and does not mean a stream cannot be fully restored or further degraded.')),
                                            
                                            h5(HTML('<b>CSCI</b>: The California Stream Condition Index as a measure of stream biological integrity.  CSCI scores of 0.63, 0.79, and 0.83 are nominal descriptions of potential biological objectives that can be used to classify and prioritize.')),
                                            
                                            h5(HTML('<b>Expected score</b>: A CSCI score from SCAPE, usually reported as a range of low to high expectations or as the median or most likely score.')),
                                            
                                            h5(HTML('<b>Investigate</b>: Additional monitoring or review of supplementary data (e.g., aerial imagery).')),
                                            
                                            h5(HTML('<b>Management priority</b>: A priority recommended by the user to investigate, protect, and/or restore a reach type.  These are actions in addition to baseline maintanence and monitoring that occurs at all sites.')),
                                            
                                            h5(HTML('<b>Observed score</b>: A CSCI score using field observations.')),
                                            
                                            h5(HTML('<b>Prioritization</b>: Identifying where certain types of management should occur.')),
                                            
                                            h5(HTML('<b>Protect</b>: Additional scrutiny of proposed development and/or projects.')),
                                            
                                            h5(HTML('<b>Relative score</b>: Difference of the observed CSCI score from the expected.  Used to define a site as over scoring, expected, or under scoring.')),
                                               
                                            h5(HTML('<b>Restore</b>: Targeted action for causal assessment and/or restoration funds.')), 
                                            
                                            h5(HTML('<b>Type</b>: An additional site description based on the relative CSCI score, the stream, and location to the selected CSCI reference threshold.'))                                            
                                            
                                            )
                                   )
                       )
              ),
              
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
                       
                       h5("These maps show the location of recommended priority actions defined for each site type in step 4. Each site can have more than one priority."),
                       
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



