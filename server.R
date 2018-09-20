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
library(gridExtra)
source('R/funcs.R')

prj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
lnsz <- 1.25
ptsz <- 4.5

# color palette for stream expectations
pal_exp <- colorFactor(
  palette = RColorBrewer::brewer.pal(9, 'Paired')[c(2, 1, 5, 6)],
  na.color = 'yellow',
  levels = c('likely unconstrained', 'possibly unconstrained', 'possibly constrained', 'likely constrained'))

# icons for map, created externally
crsz <- 11
trsz <- 15
mapicons <- iconList(
  `over scoring (lu)` = makeIcon("www/overlu.png", trsz, trsz),
  `expected (lu)` = makeIcon("www/explu.png", crsz, crsz),
  `under scoring (lu)`= makeIcon("www/underlu.png", trsz, trsz),
  `over scoring (pu)` = makeIcon("www/overpu.png", trsz, trsz),
  `expected (pu)` = makeIcon("www/exppu.png", crsz, crsz),
  `under scoring (pu)`= makeIcon("www/underpu.png", trsz, trsz),
  `over scoring (pc)` = makeIcon("www/overpc.png", trsz, trsz),
  `expected (pc)` = makeIcon("www/exppc.png", crsz, crsz),
  `under scoring (pc)`= makeIcon("www/underpc.png", trsz, trsz),
  `over scoring (lc)` = makeIcon("www/overlc.png", trsz, trsz),
  `expected (lc)` = makeIcon("www/explc.png", crsz, crsz),
  `under scoring (lc)`= makeIcon("www/underlc.png", trsz, trsz)
)

sts <- paste('Site', seq(1:16))

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

# color palette for stream expectations
pal_pri <- colorFactor(
  palette = RColorBrewer::brewer.pal(9, 'Greys')[c(8, 5, 2)],
  na.color = 'yellow',
  levels = c('Protect', 'Investigate', 'Restore'))

# color palette for CSCI type
pal_typ <- colorFactor(
  palette = RColorBrewer::brewer.pal(11, 'Spectral'),#hue_pal()(100), 
  na.color = 'yellow',
  domain = paste0('Type', sprintf('%02d', seq(1:16)))
)

# example data, csci scores
scrs_ex <- data.frame(
  Site = factor(sts, levels = sts),
  csci = c(1.25, 0.98, 0.81, 0.68, 1.14, 0.87, 0.73, 0.58, 0.98, 0.83, 0.68, 0.5, 0.88, 0.76, 0.57, 0.39)
)

# example data, stream predictions
exps_ex <- data.frame(
  Site = factor(sts, levels = sts),
  minv = rep(c(0.84, 0.68, 0.58, 0.43), each = 4),
  maxv = rep(c(1.14, 0.98, 0.88, 0.73), each = 4),
  stringsAsFactors = F
) 

# # have to reassign sgr data for reactives
# data(scrs)
# scrs_tmp <- scrs
# data(spat)
# spat_tmp <- spat

# server logic
server <- function(input, output, session) {
  
  # spatial polylines from watershed selection
  spat <- reactive({

    # na.omit(spat_tmp)
    shd <- input$shd
    spat <- paste0('spat_', shd)
    load(file = paste0('data/', spat, '.RData'))
    get(spat)
    
  })
  
  # csci scores from watershed selection
  scrs <- reactive({

    # na.omit(scrs_tmp)
    shd <- input$shd
    scrs <- paste0('scrs_', shd)
    load(file = paste0('data/', scrs, '.RData'))
    get(scrs)
    
  })
  
  
  # base mapview
  scrs_mv <- reactive({
    
    scrs() %>% 
      st_as_sf(coords = c('long', 'lat')) %>% 
      st_set_crs(prj) %>% 
      mapview(layer.name = 'reset') %>% 
      .@map
    
  })
  
  # base mapview sync, first tab
  scrs_mv1 <- reactive({
    scrs_mv() %>% 
      syncWith('maps1')
  })
  
  # base mapview sync, final tab
  scrs_mv2 <- reactive({
    scrs_mv() %>% 
      syncWith('maps2')
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
    
    out <- spat()
    # change percentile column name
    names(out)[names(out) %in% 'full0.50'] <- 'lns'
    
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
    
    # get biological condition expectations
    cls <- getcls2(spat(), thrsh = thrsh(), tails = tlinp(), modls = 'full')
    
    # join with spatial data
    out <- spat() %>% 
      left_join(cls, by = 'COMID')
    
    out
    
  })
  
  # example data to plot on priorities tab
  plot_ex <- reactive({
    
    # output
    out <- proc_all(exps_ex, scrs_ex, thrsh = 0.79, tails = 0.05) %>% 
      mutate(perf = factor(perf, levels = rev(levels(perf))))
    
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
    if(jitr){

      out <- out %>% 
        mutate(
          lat = ifelse(duplicated(lat), jitter(lat, factor = 300), lat),
          long = ifelse(duplicated(long), jitter(long, factor = 300), long)
        )
      
    # take average csci if jitter is zero
    } else {

      out <- out %>% 
        group_by(COMID, StationCode, lat, long) %>% 
        summarise(
          csci = mean(csci, na.rm = TRUE), 
          csci_difr = mean(csci_difr, na.rm = TRUE), 
          SampleDate = ifelse(length(SampleDate) > 1, 'multiple dates', as.character(SampleDate))
        ) %>% 
        ungroup
      
    }
    
    return(out)
    
  })
  
  # CSCI scores and stream condition expectations, maps only 
  scr_exp_map <- reactive({
    
    # process
    incl <- site_exp(spat(), csci(), thrsh = thrsh(), tails = tlinp(), modls = 'full') %>% 
      select(-lat, -long) %>% 
      group_by(StationCode) %>% 
      nest
    
    # assign csci station locations for jittr
    out <- csci() %>% 
      select(StationCode, lat, long, SampleDate) %>% 
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
    
    # process
    incl <- site_exp(spat(), scrs(), thrsh = thrsh(), tails = tlinp(), modls = 'full')
    
    # add additional perf column for multicolor by strcls (pal_prf)
    out <- get_perf_mlt(incl)
    
    return(out)
    
  })
  
  # CSCI thresold reactive input
  thrsh <- reactive({
    
    input$thrsh %>%
      gsub('^.*\\(|\\)$', '', .) %>% 
      as.numeric
    
  })
  
  # site priorities from user selections
  scr_pri <- reactive({
    
    out <- get_pri_inp(input, plot_ex(), scr_exp_map()) 
    return(out)
    
  })
  
  # all priority, type counts
  allcnts <- reactive({
    
    # to join to get all priority categories (for those with zero)
    allpri <- data.frame(Priority = c('Protect', 'Investigate', 'Restore'), stringsAsFactors = F)
    
    # to join to get all types (for those with zero)
    alltyp <- data.frame(Type = paste0('Type', sprintf('%02d', seq(1, 16))), stringsAsFactors = F)
    
    # get priority counts, join with allpri for all priority categories
    out <- scr_pri() %>% 
      unnest %>% 
      rename(Type = typelv) %>% 
      group_by(Priority, Type) %>% 
      nest %>% 
      mutate(n = map(data, nrow)) %>%
      select(-data) %>% 
      left_join(allpri, ., by = 'Priority') %>% 
      group_by(Priority) %>% 
      nest %>% 
      mutate(data = map(data, ~ left_join(alltyp, .x, by = 'Type'))) %>% 
      unnest %>% 
      mutate(
        n = map(n, ~ ifelse(is.null(.x), 0, .x))
      )
    
    return(out)
    
  })
  
  # cnts of prioritites across types
  cnts <- reactive({
    out <- allcnts() %>%
      select(Priority, n) %>%
      unnest %>%
      group_by(Priority) %>%
      summarise(n = sum(n)) %>%
      mutate(n = as.character(n)) %>%
      deframe %>%
      as.list
    
    return(out)
    
  })
  output$cnts_inv <- reactive({cnts()['Investigate']})
  output$cnts_pro <- reactive({cnts()['Protect']})
  output$cnts_res <- reactive({cnts()['Restore']})
  
  # plot of csci scores and expectations by station code
  output$plo_exp <- renderPlot({
    
    bysta <- input$bysta
    nocon <- input$nocon
    
    # CSCI scores and expectations
    toplo1 <- scr_exp_map() %>%
      select(COMID, StationCode, datcut, strcls, csci, perf, typelv, perf_mlt) %>%
      unnest %>%
      mutate(
        strcls = factor(strcls, levels = rev(levels(strcls))),
        perf = factor(perf, levels = rev(levels(perf)))
        ) %>%
      rename(
        `Stream Class` = strcls,
        `Relative\nscore` = perf_mlt,
        Type = typelv
      )

    # total expected range
    toplo2 <- scr_exp_map() %>%
      select(COMID, StationCode, data, strcls) %>%
      unnest %>%
      mutate(strcls = factor(strcls, levels = rev(levels(strcls)))) %>%
      rename(`Stream Class` = strcls)
    
    # median expectation
    toplo3 <- scr_exp_map() %>%
      select(COMID, StationCode, datcut) %>%
      unnest %>%
      filter(grepl('0\\.50$', var))
    
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
        geom_point(aes(x = csci), fill = 'white', shape = 21, size = 4, alpha = 0.8) +
        geom_vline(xintercept = thrsh(), linetype = 'dashed', size = 1)
      
    # otherwise full
    } else {

      # plot
      p <- ggplot(toplo1, aes(y = StationCode, x = val)) +
        geom_line(data = toplo2, aes(x = val, colour = `Stream Class`), alpha = 0.1, size = 2) +
        geom_line(aes(colour = `Stream Class`), alpha = 0.6, size = 2) +
        geom_point(data = toplo3, colour = 'white', size = 1, alpha = 1, shape = 15) +
        theme_bw(base_family = 'serif', base_size = 18) +
        theme(
          axis.text.y = element_text(size = 10)
        ) +
        scale_x_continuous('CSCI') +
        scale_y_discrete('Site') +
        scale_colour_manual(values = pal_exp(levels(toplo1$`Stream Class`))) +
        geom_point(aes(x = csci, fill = `Stream Class`, shape = perf), size = 4, alpha = 0.8) +
        geom_vline(xintercept = thrsh(), linetype = 'dashed', size = 1) +
        scale_fill_manual(values = pal_exp(levels(toplo1$`Stream Class`)), na.value = 'yellow', guide = F) + 
        scale_shape_manual('Relative site score', values = c(24, 21, 25))
      
    }
    
    print(p)
    
  })
  
  # summary tables
  output$tab_sum <- DT::renderDataTable({
    
    # summary table by csci type
    totab <- get_tab(scr_exp_map(), thrsh = thrsh(), tails = tlinp())
    
    return(totab)
    
  }, rownames = F, options = list(dom = 't', pageLength = 16))
  
  # the selection priority plot
  siteplo <- reactive({
    
    mythm <- theme_minimal(base_family = 'serif', base_size = 18) +
      theme(
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.text.x = element_blank(),
        legend.position = 'top'
      )
    
    p <- ggplot(plot_ex(), aes(x = typelv)) +
      geom_errorbar(aes(ymin = minv, ymax = maxv, colour = `Stream class`), width = 0, size = 2, alpha = 0.2) +
      geom_errorbar(aes(ymin = minv_qt, ymax = maxv_qt, colour = `Stream class`), width = 0, size = 2, alpha = 0.7) +
      geom_point(aes(y  = `CSCI score`, fill = `Stream class`, shape = perf), size = 7, alpha = 0.8) +
      geom_hline(yintercept = 0.79, linetype = 'dashed') +
      scale_colour_manual(values = pal_exp(levels(plot_ex()$`Stream class`)),
                          guide = guide_legend(direction = 'vertical', title.position = 'left')) +
      scale_fill_manual(values = pal_exp(levels(plot_ex()$`Stream class`)),
                        guide = F) +
      scale_shape_manual('Relative site score', values = c(24, 21, 25), 
                         guide = guide_legend(direction = 'veritical', title.position = 'left')) +
      scale_x_discrete(limits = rev(levels(plot_ex()$typelv))) +
      mythm +
      coord_flip()
    
    # get legend
    pleg <- g_legend(p)
    p <- p +
      mythm %+replace%
      theme(legend.position = 'none')
    
    # output as list
    plo_ls <- list(pleg, p)
    return(plo_ls)
    
  })
  
  # expectation plot legend
  output$plo_leg <- renderPlot({
    grid.arrange(siteplo()[[1]])
  })
  
  # expectation plot
  output$plo_exp2 <- renderPlot({
    siteplo()[[2]]
  })
  
  ######
  # maps
  
  # non-reactive base map for csci observed
  output$map_med <- renderLeaflet(scrs_mv1())
  
  # non-reactive base map, condition expectations
  output$map_exp <- renderLeaflet(scrs_mv1())
  
  ##
  # reactive maps, all steps
  observe({
    
    input$alltabs
    
    # other inputs
    difr <- input$difr
    
    # reactives
    dat <- dat()
    dat_exp <- dat_exp()
    scr_exp_map <- scr_exp_map()

    # score expectations
    exp_med <- leafletProxy("map_med", data = dat) %>%
      clearMarkers() %>%
      clearShapes() %>%
      clearControls() %>%
      addPolylines(opacity = 1, weight = lnsz, color = ~pal()(lns),
                   label = ~paste0(COMID, ', Likely score:', as.character(round(lns, 2)))
      ) %>%
      addLegend("topright", pal = pal(), values = ~lns,
                title = "Reach CSCI prediction",
                opacity = 1, na.label = "not in StreamCat"
      )
    
    # csci scores if false, otherwise differences
    if(difr){
      
      exp_med <- exp_med %>%
        addCircleMarkers(data = csci(), lng = ~long, lat = ~lat, radius = ptsz, weight = 1, fillOpacity = 1,
                         label = ~paste0(StationCode, ', CSCI: ', as.character(round(csci_difr, 2)), ' (', SampleDate, ')'),
                         fillColor = ~pal_difr()(csci_difr), color = 'black'
        ) %>%
        addLegend("topright", pal = pal_difr(), values = csci()$csci_difr,
                  title = "CSCI difference",
                  opacity = 1
        )
      
    } else {
      
      exp_med <- exp_med %>%
        addCircleMarkers(data = csci(), lng = ~long, lat = ~lat, radius = ptsz, weight = 1, fillOpacity = 1,
                         label = ~paste0(StationCode, ', CSCI: ', as.character(round(csci, 2)), ' (', SampleDate, ')'),
                         fillColor = ~pal()(csci), color = 'black'
        ) %>%
        addLegend("topright", pal = pal(), values = csci()$csci,
                  title = "Observed CSCI (points)",
                  opacity = 1
        )
      
    }
    
    # condition expectations
    exp_cls <- leafletProxy("map_exp", data = dat_exp) %>%
      clearMarkers() %>%
      clearShapes() %>%
      clearControls()%>%
      addLegend("topright", pal = pal_exp, values = ~strcls,
                title = "Reach classification",
                opacity = 1, na.label = "not in StreamCat"
      ) %>%
      addPolylines(opacity = 1, weight = lnsz, color = ~pal_exp(strcls),
                   label = ~paste0(COMID, ', Stream class:', strcls)
      ) %>%
      addMarkers(data = scr_exp_map, lng = ~long, lat = ~lat,
                 label = ~paste0(StationCode, ', CSCI: ', as.character(round(csci, 2)), ', ', perf_mlt),
                 icon = ~mapicons[perf_mlt]
                 
      )
    
  })
  
  # non-reactive base map, protect priority
  output$bs_pro <- renderLeaflet(scrs_mv2())
  
  # non-reactive base map, investigate priority
  output$bs_inv <- renderLeaflet(scrs_mv2())
  
  # non-reactive base map, restore priority
  output$bs_res <- renderLeaflet(scrs_mv2())
  
  ##
  # reactive maps, all steps
  observe({
    
    input$alltabs
    
    # reactives
    dat <- dat()
    dat_exp <- dat_exp()
    scr_exp_map <- scr_exp_map()
    scr_pri <- scr_pri()
    
    # get seperate priorities from scr_pri_map
    dat_pro <- filter(scr_pri, Priority %in% 'Protect')$value
    dat_inv <- filter(scr_pri, Priority %in% 'Investigate')$value
    dat_res <- filter(scr_pri, Priority %in% 'Restore')$value
    
    # base protect map
    pri_pro <- leafletProxy("bs_pro", data = dat_exp) %>%
      clearMarkers() %>%
      clearShapes() %>%
      clearControls() %>% 
      addPolylines(opacity = 1, weight = lnsz, color = ~pal_exp(strcls), 
                   label = ~paste0(COMID, ', Stream class:', strcls)
      )
    
    # base investigate map
    pri_inv <- leafletProxy("bs_inv", data = dat_exp) %>%
      clearMarkers() %>%
      clearShapes() %>%
      clearControls() %>% 
      addPolylines(opacity = 1, weight = lnsz, color = ~pal_exp(strcls), 
                   label = ~paste0(COMID, ', Stream class:', strcls)
      ) 
    
    # base restore map
    pri_res <- leafletProxy("bs_res", data = dat_exp) %>%
      clearMarkers() %>%
      clearShapes() %>%
      clearControls() %>% 
      addPolylines(opacity = 1, weight = lnsz, color = ~pal_exp(strcls), 
                   label = ~paste0(COMID, ', Stream class:', strcls)
      )
    
    # add protect points if not empty
    if(length(dat_pro) > 0){
      
      pri_pro <- pri_pro %>% 
        addMarkers(data = dat_pro[[1]], lng = ~long, lat = ~lat,
                   label = ~paste0(StationCode, ', CSCI: ', as.character(round(csci, 2)), ', ', perf_mlt, ', ', typelv),
                   icon = ~mapicons[perf_mlt]
                   
        )
      
    }
    
    # add investigate points if not empty
    if(length(dat_inv) > 0){
      
      pri_inv <- pri_inv %>% 
        addMarkers(data = dat_inv[[1]], lng = ~long, lat = ~lat,
                   label = ~paste0(StationCode, ', CSCI: ', as.character(round(csci, 2)), ', ', perf_mlt, ', ', typelv),
                   icon = ~mapicons[perf_mlt]
                   
        )
        
    }
    
    # add restore points if not empty
    if(length(dat_res) > 0){
      
      pri_res <- pri_res %>% 
        addMarkers(data = dat_res[[1]], lng = ~long, lat = ~lat,
                   label = ~paste0(StationCode, ', CSCI: ', as.character(round(csci, 2)), ', ', perf_mlt, ', ', typelv),
                   icon = ~mapicons[perf_mlt]
                   
        ) 
      
    }
    
  })
  
}