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
lat <- ncvar_get(ncfile, "lat")
lon <- ncvar_get(ncfile, "lon")
nc_close(ncfile)

# this file has a weird coordinate system that isn't linear
# first find the absolute difference between our lat and lon and every lat and lon in the coordinate system 
lat <- abs(lat - (40.554861))
lon <- abs(lon - (-85.030452))

# then add the differences together
dists <- lat + lon

# then bingo! the smallest number will be the closest point
bingo <-  which(dists == min(dists), arr.ind = TRUE)

# extract temperature data for location from each file 
loc_cumulative <- data.frame(matrix(ncol=2))
colnames(loc_cumulative) <- c("date", "temp_value")

i = 1979
while(i < 2021) {
  filename <- paste("/Volumes/TimeMachine/air.sfc/air.sfc.", i, ".nc", sep = "")
  
  ncfile <- nc_open(filename)
  date <- ncvar_get(ncfile, "time") ## hours since 1800-1-1 00:00
  air <- ncvar_get(ncfile, "air")
  nc_close(ncfile)
  
  temp_value <- air[bingo[1], bingo[2],]
  
  loc <- data.frame(date, temp_value)
  loc_cumulative <- rbind(loc_cumulative, loc)

  i = i + 1
}


# convert from degrees K to degrees C
loc_cumulative <- loc_cumulative %>%
  mutate(temp_value = temp_value - 273.15) 



# save:
write.csv(loc_cumulative, "./data-processed/time-series_loblolly-marsh.csv", row.names = FALSE)