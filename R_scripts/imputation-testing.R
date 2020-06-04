## exploring imputations
library(tidyverse)
library(imputeTS)
library(spectral)
library(lubridate)

series1 <- read.csv("./data-raw/BarkleyCanyon_BarkleyCanyonMid-East_ConductivityTemperatureDepth_Temperature_20090915T195933Z_20110728T053801Z-NaN_clean.csv", skip = 1, header = T, comment.char = "#")
series2 <- read.csv("./data-raw/BarkleyCanyon_BarkleyCanyonMid-East_ConductivityTemperatureDepth_Temperature_20110728T053901Z_20130703T163146Z-NaN_clean.csv", skip = 1, header = T, comment.char = "#")
series3 <- read.csv("./data-raw/BarkleyCanyon_BarkleyCanyonMid-East_ConductivityTemperatureDepth_Temperature_20130703T163246Z_20141231T235936Z-NaN_clean.csv", skip = 1, header = T, comment.char = "#")
series4 <- read.csv("./data-raw/BarkleyCanyon_BarkleyCanyonMid-East_ConductivityTemperatureDepth_Temperature_20150101T000036Z_20150112T233836Z-NaN_clean.csv", skip = 1, header = T, comment.char = "#")

colnames(series1) <- c("time", ## yyyy-mm-ddThh:mm:ss.sssZ
                       "temp_val",
                       "qc") ## TT.TTTT
colnames(series2) <- colnames(series1)
colnames(series3) <- colnames(series1)
colnames(series4) <- colnames(series1)

series <- rbind(series1, series2, series3, series4)
series <- series %>%
  select(-3)

rm(series1, series2, series3, series4)



## subset to just 2013-2014:
start <- which(as.Date(series$time) == "2013-01-01")
stop <- length(as.Date(series$time) == "2014-01-01")
series_short <- series[start[1]:stop,]

## visualize data gap:
plotNA.gapsize(series_short$temp_val)
plotNA.distribution(series_short$temp_val)
plotNA.distributionBar(series_short$temp_val, breaks = 30)

## plot the series:
date.num <- decimal_date(as.Date(series$time)) 
split <- str_split_fixed(series$time, pattern = "T", n = 2)
split2 <- str_split_fixed(split[,2], pattern = "Z", n = 2)
time <- split2[,1]
time <- str_split_fixed(time, pattern = ":", n = 3) ## split into hours, minutes, seconds

hour <- as.numeric(time[,1])/24
minute <- (as.numeric(time[,2]) + hour*60) / (60*24)
second <- (as.numeric(time[,3]) + minute*60) / (60*24*60)

time.num <-  hour + minute + second ## in fraction of a day where 1 day = 1.0
time.num <- time.num / 365.25 ## convert to fraction of a year

time <- as.numeric(format(date.num + time.num, nsmall = 8))

ts_plotdata <- data.frame(time, series$temp_val)

## visualize the gaps:
ts <- ggplot(ts_plotdata, aes(x = time, y = series.temp_val)) +
  geom_line() +
  scale_x_continuous(limits = c(2013, 2014), breaks = c(2013:2014), 
                     labels = c("2013", "2014")) +
  labs(x = "Year", y = "Temperature (°C)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))



### ------- IMPUTATION 1: KALMAN ------------------

imp_kal_struct <- na_kalman(series$temp_val, smooth = TRUE, model = "StructTS")

imput <- data.frame(series$time, imp_kal_struct)

ts_imput <- ggplot(imput, aes(x = time, y = imp_kal_struct)) +
  geom_line() +
  scale_x_continuous(limits = c(2013, 2014), 
                     breaks = c(2013:2014), 
                     labels = c("2013", "2014")) +
  labs(x = "Year", y = "Temperature (°C)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))


imp_seadec <- na_seadec(series$temp_val, find_frequency = TRUE, algorithm = "interpolation")

imput <- data.frame(series$time, imp_seadec)

ts_imput <- ggplot(imput, aes(x = time, y = imp_seadec)) +
  geom_line() +
  scale_x_continuous(limits = c(2013, 2014), 
                     breaks = c(2013:2014), 
                     labels = c("2013", "2014")) +
  labs(x = "Year", y = "Temperature (°C)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))


imp_kal_auto <- na_kalman(series$temp_val, smooth = TRUE, model = "auto.arima")

imput <- data.frame(series$time, imp_kal_auto)

ts_imput <- ggplot(imput, aes(x = time, y = imp_kal_auto)) +
  geom_line() +
  scale_x_continuous(limits = c(2013, 2014), 
                     breaks = c(2013:2014), 
                     labels = c("2013", "2014")) +
  labs(x = "Year", y = "Temperature (°C)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))

