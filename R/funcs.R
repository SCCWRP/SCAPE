#' Get biological expectation with changing threshold and likelihoods
#' 
#' @param datin sf object with stream COMIDS and quantile expectations
#' @param thrsh numeric for CSCI scoring thresholds
#' @param tails numeric for tails to truncate expectations for overlap with thrsh
#' @param modls chr string for selecting core (simple) or full model for expectations
#' @param lbs chr string labels for interval classifications
#' 
#' @return a nested data frame sorted by increaesing median value of expected score of COMID and nested columns as original data, cut data by tails, and sream classification (strcls).  The strcls column indicates if the ranges in datcut are within, above, or below those defind by thrsh.
getcls <- function(datin, thrsh = 0.79, tails = 0.05,  modls = c('core', 'full'), lbs = list('likely unconstrained' = 0, 'undetermined' = 1, 'likely constrained' = 2)){
  
  # sanity check
  if(tails >= 0.5)
    stop('tails must be less 0.5')
  
  # models argument
  modls <- match.arg(modls)
  
  dat <- datin
  st_geometry(dat) <- NULL
  dat <- dat %>% 
    select(matches(paste0('^COMID$|^', modls, '0'))) %>% 
    gather('var', 'val', -COMID) %>% 
    arrange(COMID, var) %>% 
    group_by(COMID) %>% 
    nest %>% 
    mutate(
      datcut = map(data, function(x){
        
        # get quantile labels to filter
        lovl <- 100 * tails
        hivl <- 100 * (1 -  tails) 
        vls <- c(lovl, hivl) %>% 
          str_pad(2, pad = '0') %>% 
          paste0(modls, '0.', .)
        
        # filter by quantile labels and median
        x <- x %>% 
          filter(var %in% vls)
        
        return(x)
        
        
      }),
      
      medv = map(data, ~ filter(.x, var %in% paste0(modls, '0.50')) %>% .$val), 
      strcls_int = map(datcut, function(x){
        
        # return NA if any zero values in predictions
        if(any(is.na(x$val))){
          
          cls <- NA
          return(cls)
          
        } 
        
        # get threshold range
        rngs <- x$val %>% 
          range
        
        # find if above/below or covering thrsh
        cls <- findInterval(thrsh, rngs)
        
        return(cls)
        
      })
      
    ) %>% 
    unnest(medv) %>% 
    arrange(medv) %>% 
    mutate(COMID = factor(COMID, levels = COMID)) %>% 
    unnest(strcls_int) 
  
  # subset lbs by those in interval
  lbs <- unique(dat$strcls_int) %>% 
    na.omit %>% 
    as.numeric %>% 
    match(unlist(lbs)) %>% 
    lbs[.]
  
  # strcls as correct factor levels
  dat <- dat %>%
    mutate(strcls = factor(strcls_int, levels = unlist(lbs), labels = names(lbs)))
  
  return(dat)
  
}

#' Get biological expectation with changing threshold and likelihoods, same as getcls but only returns class designation
#' 
#' @param datin sf object with stream COMIDS and quantile expectations
#' @param thrsh numeric for CSCI scoring thresholds
#' @param tails numeric for tails to truncate expectations for overlap with thrsh
#' @param modls chr string for selecting core (simple) or full model for expectations
#' @param lbs chr string labels for interval classifications
#' 
#' @return a nested data frame sorted by increaesing median value of expected score of COMID and nested columns as original data, cut data by tails, and sream classification (strcls).  The strcls column indicates if the ranges in datcut are within, above, or below those defind by thrsh.
getcls2 <- function(datin, thrsh = 0.79, tails = 0.05, modls = c('core', 'full'), lbs = list('likely unconstrained' = 0, 'undetermined' = 1, 'likely constrained' = 2)){
  
  # sanity check
  if(tails >= 0.5)
    stop('tails must be less than 0.5')
  
  # models argument
  modls <- match.arg(modls)
  
  dat <- datin
  st_geometry(dat) <- NULL
  dat <- dat %>% 
    select(matches(paste0('^COMID$|^', modls, '0'))) %>% 
    gather('var', 'val', -COMID) %>% 
    arrange(COMID, var) %>% 
    group_by(COMID) %>% 
    nest %>% 
    mutate(
      strcls_int = map(data, function(x){
        
        # return NA if any zero values in predictions
        if(any(is.na(x$val))){
          
          cls <- NA
          return(cls)
          
        } 
        
        # get quantile labels to filter
        lovl <- 100 * tails
        hivl <- 100 * (1 -  tails) 
        vls <- c(lovl, hivl) %>% 
          str_pad(2, pad = '0') %>% 
          paste0(modls, '0.', .)
        
        # filter by quantile labels and median
        rngs <- x %>% 
          filter(var %in% vls) %>% 
          .$val %>% 
          range
        
        cls <- findInterval(thrsh, rngs)
        
        return(cls)
        
      })
      
    ) %>% 
    select(-data) %>% 
    unnest 
  
  # subset lbs by those in interval
  lbs <- unique(dat$strcls_int) %>% 
    na.omit %>% 
    as.numeric %>% 
    match(unlist(lbs)) %>% 
    lbs[.] %>% 
    .[names(sort(unlist(.)))]
  
  # strcls as correct factor levels
  dat <- dat %>%
    mutate(strcls = factor(strcls_int, levels = unlist(lbs), labels = names(lbs)))
  
  return(dat)
  
}

#' Get CSCI StationCode expectations and performance classication
#'
#' @param datin sf object with stream COMIDS and quantile expectations
#' @param scrs CSCI scores by COMID and StationCode
#' @param thrsh numeric for CSCI scoring thresholds
#' @param tails numeric for tails to truncate expectations for overlap with thrsh
#' @param lbs chr string labels for site performance as over, expected, or under performing
#' @param ... additional arguments passed to getcls
#' 
site_exp <- function(datin, scrs, thrsh = 0.79, tails = 0.05, lbs = list('over performing' = 2, 'expected' = 1, 'under performing' = 0),
                     ...
){
  
  # site csci scores
  scrs <- scrs %>% 
    mutate(COMID = as.character(COMID)) %>% 
    select(COMID, StationCode, csci, lat, long)
  
  # filter comids with csci scores, classify, join with scores
  incl <- datin %>% 
    filter(COMID %in% scrs$COMID) %>% 
    getcls(thrsh = thrsh, tails = tails, ...) %>% 
    mutate(COMID = as.character(COMID)) %>% 
    left_join(scrs, by = 'COMID') %>% 
    arrange(medv) %>% 
    mutate(StationCode = factor(StationCode, levels = unique(StationCode))) %>% 
    select(-medv)
  
  # get CSCI performance (over/under)
  incl <- incl %>% 
    mutate(
      perf = pmap(list(datcut, csci), function(datcut, csci){
        
        # return NA if any zero values in predictions
        if(any(is.na(datcut$val))){
          
          prf <- NA
          return(prf)
          
        } 
        
        # within datcut interval
        prf <- findInterval(csci, datcut$val)
        
        return(prf)
        
      })
    ) %>% 
    unnest(perf) %>%
    mutate(bythrsh = ifelse(csci < thrsh, 0, 1)) %>%
    unite('typeprf', strcls_int, perf, bythrsh, remove = FALSE) %>% 
    mutate(
      typelv = typ_lbs(typeprf, thrsh = thrsh, tails = tails),
      typeoc = typ_lbs(typeprf, thrsh = thrsh, tails = tails, obs_sc = T)
    ) %>% 
    dplyr::select(-typeprf, -bythrsh)
  
  # subset lbs by those in interval
  lbs <- unique(incl$perf) %>% 
    na.omit %>% 
    as.numeric %>% 
    match(unlist(lbs)) %>% 
    lbs[.] %>% 
    .[names(sort(unlist(.)))]
  
  # perf as correct factor levels
  incl <- incl %>%
    mutate(perf = factor(perf, levels = unlist(lbs), labels = names(lbs)))
  
  return(incl)
  
}

#' Get type labels from three level code 
#'
#' @param vec chr string vector of codes to typify 
#' @param thrsh numeric for CSCI scoring thresholds
#' @param tails numeric for tails to truncate expectations for overlap with thrsh
#' @param obs_sc logical if observed score text qualifiers are returned
#' @param get_cds logical indicating if three level codes as list is returned
#' @details  The three level codes for type are (0, 1, 2), (0, 1, 2), and (0, 1).  The first level is likely unconstrained (0), undetermined (1), and likely constrained (2); the second level is under-performing (0), as expected (1), and over-performing (2); the third level is below threshold (0) and above threshold (1)
typ_lbs <- function(vec = NULL, thrsh = 0.79, tails = 0.05, obs_sc = FALSE, get_cds = FALSE){
  
  # type labels from codes in vec
  lbs <- list(
    Type01 = '0_2_1', 
    Type02 = '0_1_1',
    Type03 = '0_0_1',
    Type04 = '0_0_0',
    Type05 = '1_2_1',
    Type06 = '1_1_1',
    Type07 = '1_1_0',  
    Type08 = '1_0_0',
    Type09 = '2_2_1',
    Type10 = '2_2_0',
    Type11 = '2_1_0',
    Type12 = '2_0_0'
  )
  
  if(get_cds) return(lbs)
  
  # kill NA entries
  vec <- ifelse(grepl('NA', vec), NA, vec)
  
  # get tails in chr format
  lovl <- 100 * tails
  hivl <- 100 * (1 -  tails) 
  vls <- c(lovl, hivl) %>% 
    round(., 0) %>% 
    paste0(., 'th')
  
  # subset lbs by those in vec
  lbs <- unique(vec) %>% 
    na.omit %>%
    match(unlist(lbs)) %>% 
    lbs[.] %>% 
    .[sort(names(.))]
  
  # assign factor levels to vec
  vec <- factor(vec, levels = unlist(lbs), labels = names(lbs))
  
  if(!obs_sc) return(vec)
  
  # observed score chr types
  obs_sc <- list(
    Type01 = paste('>=', vls[2]), 
    Type02 = paste(vls[1], 'to', vls[2]),
    Type03 = paste(thrsh, 'to', vls[1]),
    Type04 = paste('<', thrsh),
    Type05 = paste('>=', vls[2]),
    Type06 = paste(thrsh, 'to', vls[2]),
    Type07 = paste(vls[1], 'to', thrsh), 
    Type08 = paste('<', vls[1]),
    Type09 = paste('>=', thrsh),
    Type10 = paste(vls[2], 'to', thrsh),
    Type11 = paste(vls[1], 'to', vls[2]),
    Type12 = paste('<', vls[1])
  ) %>% 
    enframe('vec', 'obs_sc') %>% 
    unnest
  
  # assign factor levels to vec
  vec <- vec %>% 
    data.frame(vec = .) %>% 
    mutate(vec = as.character(vec)) %>% 
    left_join(., obs_sc, by = 'vec') %>% 
    .$obs_sc
  
  return(vec)
  
}

#' Summarize data from site_exp by counts and all types, used for table in app
#'
#' @param datin output data.frame from site_exp
#' @param thrsh numeric for CSCI scoring thresholds
#' @param tails numeric for tails to truncate expectations for overlap with thrsh
#' @param obs_sc logical if observed score text qualifiers are returned
#' @param lbs_str chr string labels for stream comid expectation as likely constrained, undetermind, and likely unconstrained
#' @param lbs_sta chr string labels for site/station performance as over, expected, or under performing
get_tab <- function(datin, thrsh = 0.79, tails = 0.05, lbs_str = list('likely unconstrained' = 0, 'undetermined' = 1, 'likely constrained' = 2), lbs_sta = list('over performing' = 2, 'expected' = 1, 'under performing' = 0)){ 
  
  # typeoc labels
  typeoc <- typ_lbs(get_cds = T) %>% 
    unlist %>% 
    typ_lbs(., thrsh = thrsh, tails = tails, obs_sc = T)
  
  # type labels from codes
  lbs <- typ_lbs(get_cds = T) %>%
    enframe('typelv', 'codes') %>%
    unnest %>%
    separate(codes, c('strcls', 'perf', 'thrsh'), sep = '_', remove = FALSE) %>% 
    mutate(
      typelv = factor(typelv, levels = typelv),
      strcls = factor(strcls, levels = unlist(lbs_str), labels = names(lbs_str)),
      perf = factor(perf, levels = unlist(lbs_sta), labels = names(lbs_sta))
    ) %>% 
    mutate(typeoc = typeoc) %>% 
    dplyr::select(-codes, -thrsh)
  
  # get type levels and convert all to character
  lvs <- levels(lbs$typelv)
  lbs <- lbs %>% 
    mutate_all(as.character)
  
  # types from observed data, join with complete table
  totab <- datin %>%
    dplyr::select(strcls, perf, typeoc, typelv) %>%
    group_by(strcls, perf, typeoc, typelv) %>%
    summarise(Sites = n()) %>%
    na.omit %>%
    ungroup %>% 
    mutate_all(as.character) %>% 
    left_join(lbs, ., by = c('strcls', 'perf', 'typeoc', 'typelv')) %>% 
    mutate(
      Sites = ifelse(is.na(Sites), 0, Sites), 
      Sites = as.numeric(Sites),
      typelv = factor(typelv, levels = lvs)
    ) %>% 
    arrange(`typelv`) %>%
    rename(
      `Reach expectation` = strcls,
      `Site performance` = perf,
      `Observed score` = typeoc,
      Type = typelv
    ) %>%
    mutate(
      Type = gsub('^Type|^Type0', '', Type),
      Type = as.numeric(Type)
      ) %>% 
    dplyr::select(`Reach expectation`, `Site performance`, `Observed score`, Type, Sites)
  
  return(totab)
  
}

#' Get performance multi classifications
#' Input is scr_exp() reactive from app
get_perf_mlt <- function(scr_exp, lbs = c('over performing (lu)', 'expected (lu)', 'under performing (lu)',
                                          'over performing (u)', 'expected (u)', 'under performing (u)',
                                          'over performing (lc)', 'expected (lc)', 'under performing (lc)')
){
  
  # format perf_mlt as column combos of perf and strcls
  out <- scr_exp %>% 
    mutate(
      torm = map(strcls, function(chr){
        
        if(!is.na(chr)){
          as.character(chr) %>% 
            strsplit(., ' ') %>% 
            lapply(function(x) substr(x, 1, 1) %>% paste(collapse = '') %>% paste0('(', ., ')')) %>% 
            .[[1]]
        } else {
          NA
        }
        
      }),
      torm = unlist(torm)
    ) %>% 
    unite('perf_mlt', perf, torm, sep = ' ', remove = F) %>% 
    mutate(perf_mlt = ifelse(grepl('^NA|NA$', perf_mlt), NA, perf_mlt))
  
  # subset lbs by those in data
  lbs2 <- unique(out$perf_mlt) %>% 
    na.omit %>% 
    match(lbs) %>% 
    lbs[.] 
  
  # get correct order
  lbs <- lbs[lbs %in% lbs2]
  lbs2 <- lbs2[match(lbs, lbs2)]
  
  # reassign factor levels
  out <- out %>% 
    mutate(perf_mlt = factor(perf_mlt, levels = lbs2))
  
  return(out)
  
}
