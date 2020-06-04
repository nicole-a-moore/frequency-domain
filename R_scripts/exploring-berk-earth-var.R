## visualizing variation in Berkeley Earth TMAX time series'
library(tidyverse)

## read in all of the TMAX temp data for populations in intratherm:
temps <- read_csv("~/Documents/arr/data-processed/arr_temp-data.csv")


## visualize a couple time series to choose a diverse set to examine:
plotdata1 <- data.frame(time = temps$date, temp_val = temps[,600])
colnames(plotdata1) <- c("time", "temp_val")
plotdata2 <- data.frame(time = temps$date, temp_val = temps[,3])
colnames(plotdata2) <- c("time", "temp_val")

ts_combined <- ggplot(plotdata2, aes(x = time, y = temp_val)) +
  geom_line(colour = "red") +
  geom_line(data = plotdata1, aes(x = time, y = temp_val)) +
  scale_x_continuous(limits = c(1960, 1970)) +
  labs(x = "Year", y = "Temperature (°C)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(ts_combined, filename = "time-series_BerkEarth_600and3.png", path = "./figures", dpi = 300, device = "png", width = 8.802532, height = 6.47)

ts_single <- ggplot(plotdata2, aes(x = time, y = temp_val)) +
  geom_line(colour = "red") +
  scale_x_continuous(limits = c(1960, 1970)) +
  labs(x = "Year", y = "Temperature (°C)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(ts2, filename = "time-series_BerkEarth_3.png", path = "./figures", dpi = 300, device = "png", width = 8.802532, height = 6.47)

## 3:  30 to -10
## 70: 30 to 10

## 450: 40 to 0
## 340: 40 to 10

## 616: 25 to -20
## 600: 35 to 20


## plot some power spectra:
## capture length 
L <- length(plotdata1$temp_val)

# Fourier transform the series: 
dft1 <- fft(plotdata1$temp_val)/L
dft2 <- fft(plotdata2$temp_val)/L
amp1 <- sqrt(Im(dft1[-1])^2 + Re(dft1[-1])^2)	
amp1 <- amp1[1:(L/2)]				
amp2 <- sqrt(Im(dft2[-1])^2 + Re(dft2[-1])^2)	
amp2 <- amp2[1:(L/2)]				

frequency <- 1:(L/2)/L ## sampling frequency of 1 day

plotdata_freq1 <- data.frame(frequency, amp1)
plotdata_freq2 <- data.frame(frequency, amp2)

# plot magnitudes against frequencies
g <- ggplot(plotdata_freq1, aes(x=frequency, y=log(amp1))) + 
  geom_line() +
  geom_line(data = plotdata_freq2, aes(x=frequency, y=log(amp2)), colour = "red") +
  labs(x = "Frequency (1/days)", y = "log amplitude") + 
  scale_x_continuous(trans = 'log10', breaks =  c(0.0001, 0.001, 0.01, 0.1), 
                     labels = c("0.0001", "0.001", "0.01", "0.1")) +
  scale_y_continuous(expand = c(0,0), limits = c(-12, 2), breaks = c(-8, -4, 0)) +
  annotate("text", label = "2 days", x = 1/2, y = -11, size = 3) + 
  annotate("segment", x = 1/2, xend = 1/2, y = -11.25, yend = -12) +
  annotate("text", label = "1 year", x = 1/365.25, y = -11, size = 3) + 
  annotate("segment", x = 1/1/365.25, xend = 1/1/365.25, y = -11.25, yend = -12) +
  annotate("text", label = "20 years", x = 1/(365.25*20), y = -11, size = 3) + 
  annotate("segment", x = 1/(365.25*20), xend = 1/(365.25*20), y = -11.25, yend = -12) +
  annotate("text", label = "1 week", x = 1/7, y = -11, size = 3) + 
  annotate("segment", x = 1/7, xend = 1/7, y = -11.25, yend = -12) +
  annotate("text", label = "1 month", x = 1/(365.25/12), y = -11, size = 3) +
  annotate("segment", x = 1/(365.25/12), xend = 1/(365.25/12), y = -11.25, yend = -12) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(g, filename = "power-spectrum_BerkEarth_600and3.png", path = "./figures", dpi = 300, device = "png", width = 8.802532, height = 6.47)



## Parseval’s theorem: the variance of the time series equals the sum of the squared absolute values of the output of the DFT
get.total.SS = function(series){
  L = length(series)
  # Fourier transform series and capture phases and amplitudes of series
  dft<-  fft(series)/L
  amp<-sqrt(Im(dft[-1])^2+Re(dft[-1])^2)	# omit this first term
  amp<-amp[1:(L/2)]				# snag just first half of amplitudes
  # squares of each
  ampt = (2*amp)^2/2
  X = sum(ampt)
  return(X)
}

## find out what portion of variance is seasonal, annual and interannual:
get.seasonal.SS = function(series){
  L = length(series)
  
  dft<-  fft(series)/L
  amp<-sqrt(Im(dft[-1])^2+Re(dft[-1])^2)	
  amp<-amp[1:(L/2)]			
  
  amp <- amp[141:length(amp)] ## subset to amplitudes between 1/365 and year and 1 day
  # squares of each
  ampt = (2*amp)^2/2
  X = sum(ampt)
  return(X)
}

get.annual.SS = function(series){
  L = length(series)
  
  dft<-  fft(series)/L
  amp<-sqrt(Im(dft[-1])^2+Re(dft[-1])^2)	
  amp<-amp[1:(L/2)]			
  
  amp <- amp[140] ## subset to amplitudes between 1/365 and year and 1 day
  # squares of each
  ampt = (2*amp)^2/2
  X = sum(ampt)
  return(X)
}

get.interannual.SS = function(series){
  L = length(series)
  
  dft<-  fft(series)/L
  amp<-sqrt(Im(dft[-1])^2+Re(dft[-1])^2)	
  amp<-amp[1:(L/2)]			
  
  amp <- amp[1:139] ## subset to amplitudes between 1/365 and year and 1 day
  # squares of each
  ampt = (2*amp)^2/2
  X = sum(ampt)
  return(X)
}

sumvar1 = get.total.SS(plotdata1$temp_val) ## 600
sumvar2 = get.total.SS(plotdata2$temp_val) ## 3

seasonalvar1 = get.seasonal.SS(plotdata1$temp_val)
seasonalvar2 = get.seasonal.SS(plotdata2$temp_val)

annualvar1 = get.annual.SS(plotdata1$temp_val)
annualvar2 = get.annual.SS(plotdata2$temp_val)

interannualvar1 = get.interannual.SS(plotdata1$temp_val)
interannualvar2 = get.interannual.SS(plotdata2$temp_val)

ratio1s = seasonalvar1/sumvar1
ratio2s = seasonalvar2/sumvar2

ratio1a = annualvar1/sumvar1
ratio2a = annualvar2/sumvar2

ratio1ia = interannualvar1/sumvar1
ratio2ia = interannualvar2/sumvar2

sd1 <- sd(plotdata1$temp_val)
sd2 <- sd(plotdata2$temp_val)

plotdata_freq1$amp1[140]
plotdata_freq2$amp2[140]



