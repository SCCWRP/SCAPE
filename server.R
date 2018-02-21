library(shiny)
library(sf)
library(tidyverse)
library(leaflet)
library(mapview)
library(stringr)
library(scales)
library(leaflet.minicharts)
library(manipulateWidget)
library(RColorBrewer)
source('R/funcs.R')

prj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

# color palette for stream expectations
pal_exp <- colorFactor(
  palette = RColorBrewer::brewer.pal(9, 'Paired')[c(2, 1, 5, 6)],
  na.color = 'yellow',
  levels = c('likely unconstrained', 'possibly unconstrained', 'possibly constrained', 'likely constrained'))

# color palette for CSCI scoring performance
pal_prf <- colorFactor(
  palette = c(
    colorRampPalette(brewer.pal(9, 'Blues'))(100)[c(90, 74, 58)],
    colorRampPalette(brewer.pal(9, 'Blues'))(100)[c(42, 26, 10)],
    colorRampPalette(brewer.pal(9, 'Reds'))(100)[c(42, 26, 10)],
    colorRampPalette(brewer.pal(9, 'Reds'))(100)[c(90, 74, 58)]
    ),
  na.color = 'yellow',
  levels = c(
    'over scoring (lu)', 'expected (lu)', 'under scoring (lu)',
    'over scoring (pu)', 'expected (pu)', 'under scoring (pu)',
    'over scoring (pc)', 'expected (pc)', 'under scoring (pc)',
    'over scoring (lc)', 'expected (lc)', 'under scoring (lc)')
  )

# server logic
server <- function(input, output, session) {
  
  # spatial polylines from watershed selection
  spat <- reactive({
    
    shd <- input$shd
    spat <- paste0('spat_', shd)
    load(file = paste0('data/', spat, '.RData'))
    get(spat)
    
  })
  
  # csci scores from watershed selection
  scrs <- reactive({
    
    shd <- input$shd
    scrs <- paste0('scrs_', shd)
    load(file = paste0('data/', scrs, '.RData'))
    get(scrs)
    
  })
  
  # csci scores from watershed selection, as mapview object
  scrs_mv <- reactive({
    
    scrs() %>% 
      st_as_sf(coords = c('long', 'lat')) %>% 
      st_set_crs(prj) %>% 
      mapview(layer.name = 'reset') %>% 
      .@map %>% 
      syncWith('maps')
    
  })
  
  # color palette for csci scores
  pal <- reactive({
  
    # color domain, csci scores and expectations
    dmn <- spat() %>% 
      select(matches('full')) %>% 
      data.frame %>% 
      select(-geometry) %>% 
      gather('var', 'val') %>% 
      .$val %>% 
      c(., scrs()$csci) %>% 
      range(na.rm = T)  
    
    colorNumeric(
      palette = c('#d7191c', '#abd9e9', '#2c7bb6'),
      na.color = 'yellow',
      domain = dmn
      )
    
  })
  
  # color palette for csci score differences
  pal_difr <- reactive({
  
    dmn_difr <- spat() %>% 
      select(matches('full')) %>% 
      data.frame %>% 
      select(-geometry) %>% 
      gather('var', 'val') %>% 
      .$val %>% 
      range(na.rm = T)
    dmn_difr <- c(min(scrs()$csci) - dmn_difr[2], max(scrs()$csci) - dmn_difr[1]) %>% 
      abs %>% 
      max
    dmn_difr <- c(-1 * dmn_difr, dmn_difr)
    
    colorNumeric(
      palette = c('black', 'purple', 'white', 'darkgreen', 'black'),
      na.color = 'yellow',
      domain = dmn_difr)
    
  })
  
  # data to plot, polylines with score expections
  dat <- reactive({
 
    # change percentile column name
    out <- spat() %>% 
      rename(lns = `full0.50`)
    
    # output
    out
    
  })
  
  # tails input as reactive, passed to multiple
  tlinp <- reactive({
    
    tails <- input$tails %>% 
      gsub('More certain|Less certain|\\(|\\)|\\s+', '', .) %>% 
      as.numeric
    tails <- 0.5 - tails
    return(tails)
    
  })
  
  # data to plot, polylines with condition expectations
  dat_exp <- reactive({
    
    thrsh <- input$thrsh
    
    # get biological condition expectations
    cls <- getcls2(spat(), thrsh = thrsh, tails = tlinp(), modls = 'full')

    # join with spatial data
    out <- spat() %>% 
      left_join(cls, by = 'COMID')

    out
    
  })
  
  # CSCI scores, take difference from expectation if difr is true
  csci <- reactive({

    jitr <- input$jitr

    # get csci difference
    out <- dat() %>% 
      select(COMID, lns) %>% 
      mutate(COMID = as.character(COMID)) %>% 
      left_join(scrs(), ., by = 'COMID') %>% 
      mutate(csci_difr = csci - lns)
  
    # jitter scores with overlapping lat/lon
    if(jitr != 0){
    
      out <- out %>% 
        mutate(
          lat = ifelse(duplicated(lat), jitter(lat, factor = jitr), lat),
          long = ifelse(duplicated(long), jitter(long, factor = jitr), long)
        )

    # take average csci is jitter is zero
    } else {
      
      out <- out %>% 
        group_by(COMID, StationCode, lat, long) %>% 
        summarise(
          csci = mean(csci, na.rm = TRUE), 
          csci_difr = mean(csci_difr, na.rm = TRUE)
          ) %>% 
        ungroup
        
    }
    
    return(out)
      
  })
  
  # CSCI scores and stream condition expectations, maps only 
  scr_exp_map <- reactive({

    thrsh <- input$thrsh
    
    # process
    incl <- site_exp(spat(), csci(), thrsh = thrsh, tails = tlinp(), modls = 'full') %>% 
      select(-lat, -long) %>% 
      group_by(StationCode) %>% 
      nest
    
    # assign csci station locations for jittr
    out <- csci() %>% 
      select(StationCode, lat, long) %>% 
      group_by(StationCode) %>% 
      nest %>% 
      mutate(StationCode = factor(StationCode, levels = levels(incl$StationCode))) %>% 
      left_join(incl, ., by = 'StationCode') %>% 
      unnest
    
    # add additional perf column for multicolor by strcls (pal_prf)
    out <- get_perf_mlt(out)
    
    return(out)
    
  })

  # CSCI scores and stream condition expectations, not on maps
  # these data are never averaged by station averaged for CSCI
  scr_exp <- reactive({
    
    thrsh <- input$thrsh
    
    # process
    incl <- site_exp(spat(), scrs(), thrsh = thrsh, tails = tlinp(), modls = 'full')
    
    # add additional perf column for multicolor by strcls (pal_prf)
    out <- get_perf_mlt(incl)
    
    return(out)
    
  })
  
  ##
  # these update inputs that are duplicated across tabs
  trs <- ''
  tls <- ''
  
  # thrsh
  observe({
    if (trs != input$thrsh){
      updateSliderInput(session, "thrsh2", value = input$thrsh)
      updateSliderInput(session, "thrsh3", value = input$thrsh)
      trs <<- input$thrsh
    }
  })
  observe({
    if (trs != input$thrsh2){
      updateSliderInput(session, "thrsh", value = input$thrsh2)
      updateSliderInput(session, "thrsh3", value = input$thrsh2)
      trs <<- input$thrsh2
    }
  })
  observe({
    if (trs != input$thrsh3){
      updateSliderInput(session, "thrsh", value = input$thrsh3)
      updateSliderInput(session, "thrsh2", value = input$thrsh3)
      trs <<- input$thrsh3
    }
  })
  
  # tails
  observe({
    if (trs != input$tails){
      updateSliderTextInput(session, "tails2", selected = input$tails)
      updateSliderTextInput(session, "tails3", selected = input$tails)
      trs <<- input$tails
    }
  })
  observe({
    if (trs != input$tails2){
      updateSliderTextInput(session, "tails", selected = input$tails2)
      updateSliderTextInput(session, "tails3", selected = input$tails2)
      trs <<- input$tails2
    }
  })
  observe({
    if (trs != input$tails3){
      updateSliderTextInput(session, "tails", selected = input$tails3)
      updateSliderTextInput(session, "tails2", selected = input$tails3)
      trs <<- input$tails3
    }
  })
  
  # non-reactive base map
  output$map <- renderLeaflet(scrs_mv())
  
  # non-reactive base map, condition expectations
  output$map_exp <- renderLeaflet(scrs_mv())
  
  ##
  # reactive maps
  observe({

    # other inputs
    ptsz <- input$pt_sz
    lnsz <- input$ln_sz
    typs <- input$typs
    difr <- input$difr
    
    # reactives
    dat <- dat()
    dat_exp <- dat_exp()
    scr_exp_map <- scr_exp_map()
    
    # score expectations
    exp <- leafletProxy("map", data = dat) %>%
      clearMarkers() %>%
      clearShapes() %>% 
      clearControls() %>% 
      addPolylines(opacity = 1, weight = lnsz, color = ~pal()(lns), 
                   label = ~paste0(COMID, ', Likely score:', as.character(round(lns, 2)))
      ) %>% 
      addLegend("topright", pal = pal(), values = ~lns,
                title = "Reach prediction (lines)",
                opacity = 1
      )
    
    # csci scores if false, otherwise differences
    if(difr){

      exp <- exp %>% 
        addCircleMarkers(data = csci(), lng = ~long, lat = ~lat, radius = ptsz, weight = 0.9, fillOpacity = 0.8, 
                         label = ~paste0(StationCode, ', CSCI: ', as.character(round(csci_difr, 2))),
                         fillColor = ~pal_difr()(csci_difr), color = 'black'
        ) %>% 
        addLegend("topright", pal = pal_difr(), values = csci()$csci_difr,
                  title = "CSCI difference",
                  opacity = 1
        )
      
    } else {
      
      exp <- exp %>% 
        addCircleMarkers(data = csci(), lng = ~long, lat = ~lat, radius = ptsz, weight = 0.9, fillOpacity = 0.8, 
                         label = ~paste0(StationCode, ', CSCI: ', as.character(round(csci, 2))),
                         fillColor = ~pal()(csci), color = 'black'
        ) %>% 
        addLegend("topright", pal = pal(), values = csci()$csci,
                  title = "CSCI observed (points)",
                  opacity = 1
        )
      
    }
      
    # condition expectations
    exp_bs <- leafletProxy("map_exp", data = dat_exp) %>%
      clearMarkers() %>%
      clearShapes() %>% 
      clearControls()%>% 
      addLegend("topright", pal = pal_exp, values = ~strcls,
                title = "Expected classification (lines)",
                opacity = 1
      ) %>% 
      addPolylines(opacity = 1, weight = lnsz, color = ~pal_exp(strcls), 
                   label = ~paste0(COMID, ', Stream class:', strcls)
      ) %>% 
      addCircleMarkers(data = scr_exp_map, lng = ~long, lat = ~lat, radius = ptsz, weight = 0.9, fillOpacity = 0.9, 
                       label = ~paste0(StationCode, ', CSCI: ', as.character(round(csci, 2)), ', ', perf_mlt),
                       fillColor = ~pal_prf(perf_mlt), color = 'black'
      ) %>% 
      addLegend("topright", pal = pal_prf, values = scr_exp_map$perf_mlt,
                title = "CSCI scoring (points)",
                opacity = 1
      )
  
  # sync the maps
  combineWidgets(exp, exp_bs)
  
  })
  
  # plot of csci scores and expectations by station code
  output$plo_exp <- renderPlot({

    bysta <- input$bysta
    nocon <- input$nocon
    thrsh <- input$thrsh
    
    # CSCI scores and expectations
    toplo1 <- scr_exp() %>% 
      select(COMID, StationCode, datcut, strcls, csci, perf, typelv, perf_mlt) %>% 
      unnest %>% 
      mutate(strcls = factor(strcls, levels = rev(levels(strcls)))) %>% 
      rename(
        `Stream Class` = strcls,
        `Relative\nscoring` = perf_mlt,
        Type = typelv
        )

    # total expected range
    toplo2 <- scr_exp() %>% 
      select(COMID, StationCode, data, strcls) %>% 
      unnest %>% 
      mutate(strcls = factor(strcls, levels = rev(levels(strcls)))) %>% 
      rename(`Stream Class` = strcls)

    # arrange by station if true
    if(bysta){
      
      toplo1 <- toplo1 %>% 
        mutate(StationCode = as.character(StationCode)) %>% 
        arrange(StationCode)%>% 
        mutate(
          StationCode = factor(StationCode),
          StationCode = factor(StationCode, levels = rev(levels(StationCode)))
        )
      
      toplo2 <- toplo2 %>% 
        mutate(StationCode = as.character(StationCode)) %>% 
        arrange(StationCode) %>% 
        mutate(
          StationCode = factor(StationCode),
          StationCode = factor(StationCode, levels = rev(levels(StationCode)))
        )
      
    }
    
    # bare bones plot if true
    if(nocon){
      
      # plot
      p <- ggplot(toplo1, aes(y = StationCode, x = val)) + 
        theme_bw(base_family = 'serif', base_size = 18) +
        theme(
          axis.text.y = element_text(size = 10)
        ) +
        scale_x_continuous('CSCI') +
        scale_y_discrete('Site') +
        geom_point(aes(x = csci), fill = 'white', shape = 21, size = 6, alpha = 0.8) +
        geom_vline(xintercept = thrsh, linetype = 'dashed', size = 1)
           
    # otherwise full
    } else {
      
      # plot
      p <- ggplot(toplo1, aes(y = StationCode, x = val)) + 
        geom_line(data = toplo2, aes(x = val, colour = `Stream Class`), alpha = 0.1, size = 2) +
        geom_line(aes(colour = `Stream Class`), alpha = 0.6, size = 2) + 
        theme_bw(base_family = 'serif', base_size = 18) +
        theme(
          axis.text.y = element_text(size = 10)
        ) +
        scale_x_continuous('CSCI') +
        scale_y_discrete('Site') +
        scale_colour_manual(values = pal_exp(levels(toplo1$`Stream Class`))) +
        geom_point(aes(x = csci, fill = `Relative\nscoring`), shape = 21, size = 6, alpha = 0.8) +
        geom_vline(xintercept = thrsh, linetype = 'dashed', size = 1) +
        scale_fill_manual(values = pal_prf(levels(toplo1$`Relative\nscoring`)), na.value = 'yellow')
      
    }
    
    print(p)
    
  })
  
  # summary tables
  output$tab_sum <- DT::renderDataTable({

    thrsh <- input$thrsh

    # summary table by csci type
    totab <- get_tab(scr_exp(), thrsh = thrsh, tails = tlinp())

    return(totab)

  }, rownames = F, options = list(dom = 't', pageLength = 16))
   
}
