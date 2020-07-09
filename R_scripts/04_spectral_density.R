#### Version of colour-of-noise script "04_spectral_density.R" using Loblolly Marsh data
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())


## temperature data
place1 <- read.csv("./data-processed/time-series_loblolly-marsh.csv") %>%
  .[1:119800,] %>%
  rename(time = date, sst = temp_value)

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

place1$time <- as.numeric(new_date)


### plot the time series
place1 %>% 
  ggplot(aes(x = time, y = sst)) + geom_line() +
  ylab("Temperature (°C)") + xlab("Date") +
  theme(text = element_text(size=14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

##ggsave("./figures/temp-time-series_LoblollyMarsh.png", width = 6, height = 3)


### plot the correlogram
sst_acf <- acf(place1$sst, lag.max = 365*8) ## sampling frequency was once every 3 hours = 8 times/day

correlation <- data.frame(correlation = sst_acf$acf, lag = sst_acf$lag)

correlation %>% 
  ggplot(aes(x = lag, y = correlation)) + geom_line() +
  geom_hline(yintercept = 0) + ylab("Autocorrelation") + xlab("Lag (days)") +
  scale_x_continuous(breaks = c(0,800,1600,2400), labels = c("0", "100", "200", "300"))

##ggsave("./figures/sst-correlogram_LoblollyMarsh.png", width = 6, height = 3)


eighth <- place1[seq(1, nrow(place1), 8),]
sst_acf <- acf(eighth$sst, lag.max = 365) ## look at only every eighth temp to give daily sampling frequency

correlation <- data.frame(correlation = sst_acf$acf, lag = sst_acf$lag)

correlation %>% 
  ggplot(aes(x = lag, y = correlation)) + geom_line() +
  geom_hline(yintercept = 0) + ylab("Autocorrelation") + xlab("Lag (days)") 

ggsave("./figures/sst-correlogram-eighth_LoblollyMarsh.png", width = 6, height = 3)


### estimate the spectral density
mspect <- spectrum(place1$sst, spans=c(2,2), plot=FALSE)
mspect$freq = mspect$freq/180 ## correct frequency to sampling frequency by dividing by 180mins

delta <- 1
specx <- mspect$freq/delta
specy <- 2*mspect$spec

spectral <- data.frame(specx = mspect$freq/delta, specy = specy)

#### plot
spectral %>% 
  ggplot(aes(x = specx, y = specy)) + geom_line() +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  scale_x_log10(breaks =  c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1), 
                labels = c("0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1")) +
  geom_vline(xintercept = 1/(365*24*60), color = "cadetblue", size = 1) + ## 1 year
  geom_vline(xintercept = 1/(182.5*24*60), color = "green", size = 1) + ## 1/2 year
  geom_vline(xintercept = 1/(30*24*60), color = "turquoise", size = 1) + ## 1 month
  geom_vline(xintercept = 1/(365*10*24*60), color = "purple", size = 1) + ## 10 years
  geom_vline(xintercept = 1/(7*24*60), color = "pink", size = 1) + ## 1 week
  geom_vline(xintercept = 1/(24*60), color = "orange", size = 1) + ## 1 day
  geom_vline(xintercept = 1/(12*60), color = "yellow", size = 1) + ## 1/2 day
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  ylab("Spectral density") +
  xlab("Frequency (1/minutes)")

##ggsave("./figures/sst-spectral-slope-normalred-noise_LoblollyMarsh.png", width = 6, height = 3)


### mess around and plot some of the component sine waves
x <- seq(0,365,length.out=1000)
y <- 9263*sin(x*(2*pi/(365)))
df <- data.frame(x = x, y = y)

x2 <- seq(0,365,length.out=1000)
y2 <- 2000*sin(x2*(2*pi/20))
df2 <- data.frame(x = x2, y = y2)

x3 <- seq(0,365,length.out=1000)
y3 <- 100*sin(x3*(2*pi/7))
df3 <- data.frame(x = x3, y = y3)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + geom_line(aes(x = x, y = y), data = df, color = "cadetblue",size = 1) +
  geom_line(aes(x = x, y = y), data = df2, color = "turquoise", size = 1) +
  geom_line(aes(x = x, y = y), data = df3, color = "pink", size = 1) +
  xlim(0, 365) + ylab("Temperature (°C)") + xlab("Time (days)") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  theme(text = element_text(size=14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_blank())

##ggsave("./figures/sine-waves_LoblollyMarsh.png", width = 6, height = 3)
