## getting temp data for Loblolly Marsh Wetland
library(tidyverse)
library(RCurl)
library(ncdf4)


## download NCEP North American Regional Reanalysis: NARR 8x daily air surface temperature data
url = "ftp://ftp.cdc.noaa.gov/Datasets/NARR/monolevel/"

filenames = getURL(url, dirlistonly = TRUE) %>%
  strsplit(., "\\n") %>%
  unlist(.)

download_files <- c()
i = 1979
while(i < 2021) {
  download_files <- append(download_files, paste("air.sfc.", i, ".nc", sep = ""))
  i = i + 1
}

for (filename in download_files) {
  download.file(paste(url, filename, sep = ""), paste("/Volumes/TimeMachine/air.sfc/", filename,
                                                      sep = ""))
}

# find closest coordinates to location 
# 40.554861, -85.030452
filename <- "/Volumes/TimeMachine/air.sfc/air.sfc.1979.nc"

ncfile <- nc_open(filename)
lat <- ncvar_get(ncfile, "latitude")
long <- ncvar_get(ncfile, "longitude")
nc_close(ncfile)

# extract temperature data for location from each file 
loc_cumulative <- data.frame(matrix(ncol=2))
colnames(loc_cumulative) <- c("date", "temp_value")

i = 1979
while(i < 2021) {
  filename <- paste("/Volumes/TimeMachine/air.sfc/air.sfc.", i, ".nc", sep = "")
  
  ncfile <- nc_open(filename)
  
  date <- ncvar_get(ncfile, "date_number")
  arr.anom <-ncvar_get(ncfile, "temperature")
  
  ## close the file
  nc_close(ncfile)
  
  
  
}


