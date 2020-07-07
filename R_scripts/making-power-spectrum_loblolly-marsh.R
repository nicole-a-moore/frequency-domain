## make power spectrum for loblolly marsh
library(tidyverse)
library(gridExtra)

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

frequency <- 1:(L/2)/(L*3*60) ## in 1/minutes, sampling frequency: 1/3 hours

plotdata <- data.frame(frequency, amp)

# plot magnitudes against frequencies
gg = ggplot(plotdata, aes(x = frequency, y = log10(amp))) + 
  geom_line() +
  labs(x = "Frequency (1/minutes)", y = "log amplitude") + 
  scale_x_continuous(trans = 'log10', breaks =  c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1), 
                     labels = c("0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1")) +
  scale_y_continuous(expand = c(0,0), limits = c(-6, 1), breaks = c(-6:1)) +
  annotate("text", label = "1 year", x = 1/(365.25*24*60), y = -5.3, size = 3) +
  annotate("segment", x = 1/(365.25*24*60), xend = 1/(365.25*24*60), y = -6, yend = -5.5) +
  annotate("text", label = "1 day", x = 1/(24*60) - 0.00005, y = -5.3, size = 3) +
  annotate("segment", x = 1/(24*60), xend = 1/(24*60), y = -6, yend = -5.5) +
  annotate("text", label = "1 week", x = 1/(24*60*7), y = -5.3, size = 3) +
  annotate("segment", x = 1/(24*60*7), xend = 1/(24*60*7), y = -6, yend = -5.5) +
  annotate("text", label = "1 month", x = 1/(30*24*60), y = -5.3, size = 3) +
  annotate("segment", x = 1/(30*24*60), xend = 1/(30*24*60), y = -6, yend = -5.5) +
  annotate("text", label = "1/2 day", x = 1/(12*60), y = -5.3, size = 3) +
  annotate("segment", x = 1/(12*60), xend = 1/(12*60), y = -6, yend = -5.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(gg, filename = "power-spectrum_LoblollyMarsh.png", path = "./figures", dpi = 300, device = "png", width = 10.71, height = 6.47)


## read in lifespans to add
lifespans <- read.csv("./data-raw/species-list_LoblollyMarsh_completed.csv") %>%
  mutate(lifespan = as.numeric(ifelse(lifespan == "unk", NA, as.character(lifespan)))) %>%
  mutate(lifespan = round(lifespan, digits = 0)) %>%
  mutate(lifespan_minutes = lifespan*24*60) %>%
  mutate(lifespan_frequency = 1/lifespan_minutes)
  

gg_withlife <- ggplot(plotdata, aes(x = frequency, y = log10(amp))) + 
  labs(x = "Frequency (1/minutes)", y = "log amplitude") + 
  scale_x_continuous(trans = 'log10', breaks =  c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1), 
                     labels = c("0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1")) +
  scale_y_continuous(expand = c(0,0), limits = c(-6, 1), breaks = c(-6:1)) +
  annotate("text", label = "1 year", x = 1/(365.25*24*60), y = -5.3, size = 3) +
  annotate("segment", x = 1/(365.25*24*60), xend = 1/(365.25*24*60), y = -6, yend = -5.5) +
  annotate("text", label = "1 day", x = 1/(24*60) - 0.00005, y = -5.3, size = 3) +
  annotate("segment", x = 1/(24*60), xend = 1/(24*60), y = -6, yend = -5.5) +
  annotate("text", label = "1 week", x = 1/(24*60*7), y = -5.3, size = 3) +
  annotate("segment", x = 1/(24*60*7), xend = 1/(24*60*7), y = -6, yend = -5.5) +
  annotate("text", label = "1 month", x = 1/(30*24*60), y = -5.3, size = 3) +
  annotate("segment", x = 1/(30*24*60), xend = 1/(30*24*60), y = -6, yend = -5.5) +
  annotate("text", label = "1/2 day", x = 1/(12*60), y = -5.3, size = 3) +
  annotate("segment", x = 1/(12*60), xend = 1/(12*60), y = -6, yend = -5.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) +
  annotate("segment", x = lifespans$lifespan_frequency, 
           xend = lifespans$lifespan_frequency,
           y = -6, 
           yend = 1,
           colour = "red",
           size = 0.15) +
  geom_line() 

## manually create legend 
gg_legend <- ggplot() + 
  geom_blank() +
  theme_void() +
  theme(panel.spacing.x=unit(1, "lines")) +
  scale_y_continuous(limits = c(75, 75)) +
  scale_x_continuous(limits = c(50,90)) +
  annotate("segment", x = 60, xend = 62, y = 75, yend = 75, colour = "red") +
  annotate("text", label = "Lifespan of species in community", x = 72, y = 75)

lay <- rbind(c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(2,2,2,2,2))

gg_complete <- grid.arrange(gg_withlife, gg_legend, layout_matrix = lay)

ggsave(gg_complete, filename = "power-spectrum-with-lifespan_Loblolly-Marsh.png", path = "./figures", dpi = 300, device = "png", height = 6.47, width = 8.369152)

