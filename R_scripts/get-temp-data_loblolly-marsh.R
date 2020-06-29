## getting temp data for Loblolly Marsh Wetland
library(tidyverse)
library(RCurl)

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


