library(tidyverse)
library(sf)
library(lubridate)
library(rgdal)
library(readxl)

prstr <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

######
# processing data files for use with classify app

## 
# whole state data

# flw_pth <- 'ignore/SMR_NHDPlus.shp'
# shd_pth <- 'ignore/calwater_SWAMP3Code.shp'
flw_pth <- 'S:/Spatial_Data/SMCBasefiles/Frames/CompleteFrame_Flow_LU.shp'
shd_pth <- 'S:/Spatial_Data/SMCBasefiles/Boundaries/SMCSheds/SMCSheds2009.shp'
scr_pth <- 'ignore/csci_061917.csv'
exp_pth <- 'ignore/comid_statewide.Rdata'

# extra SGR sites from KW
sgr_ex <- read_excel('ignore/SGRRMP Missing SCAPE Sites CSCI Scores 05.22.2019.xlsx')

# stream expectations, named comid
load(file = exp_pth)

# csci scores
scrs <- read.csv(scr_pth, header = T, stringsAsFactors = F) %>% 
  dplyr::select(SampleID, StationCode, New_Lat, New_Long, COMID, SampleDate, CSCI) %>% 
  rename(
    csci = CSCI, 
    lat = New_Lat, 
    long = New_Long
  ) %>% 
  mutate(
    SampleDate = dmy(SampleDate), 
    COMID = as.character(COMID)
  )

# add extra SGR sites to scrs
sgr_ex <- sgr_ex %>% 
  rename(
    lat = Latitude, 
    long = Longitude, 
    csci = Result
  ) %>% 
  mutate(
    SampleDate = ymd(SampleDate), 
    COMID = as.character(COMID), 
    rep = 1
  ) %>% 
  unite('SampleID', StationCode, SampleDate, rep, sep = '_', remove = F) %>% 
  dplyr::select(SampleID, StationCode, lat, long, COMID, SampleDate, csci)
scrs <- scrs %>% 
  bind_rows(sgr_ex)

# watersheds
shed <- readOGR(shd_pth) %>% 
  spTransform(CRS(prstr)) %>% 
  st_as_sf

# flowlines
spat <- readOGR(flw_pth) %>% 
  spTransform(CRS(prstr)) %>% 
  st_as_sf %>% 
  st_intersection(shed) %>% 
  select(COMID)  


##
# process separate spatial and score files for each watershed

# sheds to process, appended to file names
shds <- shed$SMC_Name

# process and save files for each shed
for(shd in shds){
  
  # counter
  cat(shd, which(shds %in% shd), 'of', length(shds), '\n')
  
  # filter watershed for intersect
  shd_tmp <- shed %>% 
    filter(SMC_Name %in% shd)
  
  # create spatial polyines from shed intersect, left_join with csci scrs
  sel <- st_covered_by(spat, shd_tmp, sparse = F) 
  spat_tmp <- spat %>% 
    filter(sel[, 1]) %>% 
    left_join(comid, by = 'COMID') %>% 
    select(COMID, matches('^core0'))
  
  # csci scores
  scrs_tmp <- scrs %>% 
    filter(COMID %in% spat_tmp$COMID) %>% 
    unique
  
  # assign unique names to scrs and spat
  scrs_shd <- paste0('scrs_', shd)
  spat_shd <- paste0('spat_', shd)
  assign(scrs_shd, scrs_tmp)
  assign(spat_shd, spat_tmp)
  
  # save unique scrs, spat
  save(list = scrs_shd, file = paste0('data/', scrs_shd, '.RData'))
  save(list = spat_shd, file = paste0('data/', spat_shd, '.RData'))
  
}

######
# spat with wgs projection

data(spat)

prj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
spat <- st_transform(spat, crs = prj)

save(spat, file = 'data/spat.RData')