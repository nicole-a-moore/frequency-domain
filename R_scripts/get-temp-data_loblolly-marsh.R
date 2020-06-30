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


loc_cumulative <- read.csv("./data-processed/time-series_loblolly-marsh.csv") %>%
  .[1:119800,] ## chop to even year


# convert date so plotting is easier 
years <- rep(1979:2019, each = 365*8)
# insert extra day (8 elements) on leap years
leap_years = unique(years[which(years %% 4 == 0)])
z = 1
while (z < length(leap_years) + 1) {
  years = append(years, rep(leap_years[z], 8), after = first(which(years == leap_years[z])))
  z = z+1
}
# add decimal representing fraction of year
fraction_normalyear <- c((1:(365*8))/(365*8))
fraction_leapyear <- c((1:(366*8))/(366*8))
i = 1979
new_date <- c()
while (i < 2020) {
  if (i %% 4 == 0) {
    new_date <- append(new_date, fraction_leapyear)
  }
  else {
    new_date <- append(new_date, fraction_normalyear)
  }
  i = i + 1
}
new_date <- new_date + years
new_date = format(new_date, nsmall = 8)

loc_cumulative$date <- as.numeric(new_date)


# plot the time series
ts <- ggplot(loc_cumulative, aes(x = date, y = temp_value)) +
  geom_line() +
  scale_x_continuous(expand = c(0,0)) + 
  labs(x = "Year", y = "Temperature (°C)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(ts, filename = "/time-series_LoblollyMarsh.png", path = "./figures" , dpi = 300, device = "png", width = 10.71, height = 6.47)


# make power spectrum:
temp <- loc_cumulative$temp_value
L <- length(temp)

# Fourier transform the series: 
dft <- fft(temp)/L

amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2)	# omit first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]				# remove second half of amplitudes (negative half)

frequency <- 1:(L/2)/L ## sampling frequency (1 minute, 2 minutes, 3 minutes etc..) / length of time series 

plotdata <- data.frame(frequency, amp)

# plot magnitudes against frequencies
gg = ggplot(plotdata, aes(x = frequency, y = log10(amp))) + 
  geom_line() +  
  labs(x = "Frequency (1/minutes)", y = "log amplitude") + 
  scale_x_continuous(trans = 'log10', breaks =  c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1), 
                     labels = c("0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1")) +
  scale_y_continuous(expand = c(0,0), limits = c(-6, 1), breaks = c(-6:1)) +
  annotate("text", label = "1 year", x = 1/(365.25*8), y = -5.3, size = 3) +
  annotate("segment", x = 1/(365.25*8), xend = 1/(365.25*8), y = -6, yend = -5.5) +
  annotate("text", label = "1 day", x = 1/(8.5), y = -5.3, size = 3) +
  annotate("segment", x = 1/(8), xend = 1/(8), y = -6, yend = -5.5) +
  annotate("text", label = "1 week", x = 1/(8*7), y = -5.3, size = 3) +
  annotate("segment", x = 1/(8*7), xend = 1/(8*7), y = -6, yend = -5.5) +
  annotate("text", label = "1 month", x = 1/(8*30), y = -5.3, size = 3) +
  annotate("segment", x = 1/(8*30), xend = 1/(8*30), y = -6, yend = -5.5) +
  annotate("text", label = "1/2 day", x = 1/(3.5), y = -5.3, size = 3) +
  annotate("segment", x = 1/(4), xend = 1/(4), y = -6, yend = -5.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(gg, filename = "power-spectrum_LoblollyMarsh.png", path = "./figures", dpi = 300, device = "png", width = 10.71, height = 6.47)




