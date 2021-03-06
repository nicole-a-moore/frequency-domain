## make power spectrum for loblolly marsh
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(png)
library(grid)

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
  mutate(lifespan_frequency = 1/lifespan_minutes) %>%
  group_by(group) %>%
  mutate(group_mean = mean(lifespan, na.rm=TRUE)) %>%
  arrange(desc(group_mean)) %>% 
  ungroup() %>%
  filter(!is.na(lifespan)) %>%
  mutate(group = ifelse(group == "Damselflies", "Damselflies & Dragonflies", as.character(group))) %>%
  mutate(group = ifelse(group == "Dragonflies", "Damselflies & Dragonflies", as.character(group))) %>%
  mutate(group = ifelse(group == "Insects", "Other insects", as.character(group))) %>%
  mutate(group = ifelse(group == "Bees", "Other insects", as.character(group))) %>%
  mutate(group = ifelse(group == "Stoneflies", "Other insects", as.character(group))) 

cols <- data.frame(group = unique(lifespans$group)) %>%
  mutate(colour = brewer.pal(n = 5, name = "Spectral"))

lifespans <- left_join(lifespans, cols)

gg_withlife <- ggplot(plotdata, aes(x = frequency, y = log10(amp))) + 
  annotate("segment", x = lifespans$lifespan_frequency, 
           xend = lifespans$lifespan_frequency,
           y = -6, 
           yend = 1,
           colour = lifespans$colour,
           size = 0.5) +
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
  geom_line() 

## manually create legend 
x <- c(55:62)

gg_legend <- ggplot() + 
  geom_blank() +
  theme_void() +
  theme(panel.spacing.x=unit(1, "lines")) +
  scale_y_continuous(limits = c(75, 75)) +
  scale_x_continuous(limits = c(50,90)) +
  geom_line(aes(x=x, y=75, colour = x)) +
  scale_color_gradientn(colours = unique(lifespans$colour)) +
  annotate("text", label = "Lifespan of species in community", x = 72, y = 75) +
  theme(legend.position = "none")


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

ggsave(gg_complete, filename = "power-spectrum-with-lifespan_Loblolly-Marsh_changed-groups.png", path = "./figures", dpi = 300, device = "png", height = 7.47, width = 10.369152)

## read in and properly colour images
frog <- as.raster(png::readPNG("./data-raw/frog.png")) 
frog[frog != "#FFFFFF"] <- lifespans$colour[first(which(lifespans$group=="Amphibia"))]
frog[frog == "#FFFFFF"] <- "#FFFFFF00"
gg_frog <- rasterGrob(frog, interpolate=TRUE, x = .1, y = .9, height = 0.15, width = .1)

dragonfly <- as.raster(png::readPNG("./data-raw/dragonfly.png"))
dragonfly[dragonfly != "#00000000"] <- lifespans$colour[first(which(lifespans$group=="Damselflies & Dragonflies"))]
gg_dragonfly <- rasterGrob(dragonfly, x = .6, y = .9, height = 0.15, width = .15)

butterfly <- as.raster(png::readPNG("./data-raw/butterfly.png"))
butterfly[butterfly == "#020202FF"] <- lifespans$colour[first(which(lifespans$group=="Butterflies"))]
butterfly[butterfly == "#020202F0"] <- lifespans$colour[first(which(lifespans$group=="Butterflies"))]
butterfly[butterfly == "#47704C00"] <- "#FFFFFF00"
butterfly[butterfly == "#010101A7"] <- lifespans$colour[first(which(lifespans$group=="Butterflies"))]
butterfly[butterfly == "#010101D9"] <- lifespans$colour[first(which(lifespans$group=="Butterflies"))]
butterfly[butterfly == "#010101C2"] <- lifespans$colour[first(which(lifespans$group=="Butterflies"))]
butterfly[butterfly == "#0000002B"] <- "#FFFFFF00"
butterfly[butterfly == "#00000086"] <- "#FFFFFF00"
butterfly[butterfly == "#00000069"] <- "#FFFFFF00"
butterfly[butterfly == "#0000004E"] <- "#FFFFFF00"
gg_butterfly <- rasterGrob(butterfly, interpolate=TRUE, x = .3, y = .9, height = 0.15, width = .1)

reptile <- as.raster(png::readPNG("./data-raw/reptile.png")) 
reptile[reptile != "#00000000"] <- lifespans$colour[first(which(lifespans$group=="Reptilia"))]
reptile[reptile == "#00000000"] <- "#FFFFFF00"
gg_reptile <- rasterGrob(reptile, interpolate=TRUE, x = .2, y = .85, height = 0.2, width = .09)

insects <- as.raster(png::readPNG("./data-raw/insect.png")) 
insects[insects != "#FFFFFF"] <- lifespans$colour[first(which(lifespans$group=="Other insects"))]
insects[insects == "#FFFFFF"] <- "#FFFFFF00"
gg_insects <- rasterGrob(insects, interpolate=TRUE, x = .45, y = .9, height = 0.125, width = .125)


gg_withpics <- ggplot(plotdata, aes(x = frequency, y = log10(amp))) + 
  annotate("segment", x = lifespans$lifespan_frequency, 
           xend = lifespans$lifespan_frequency,
           y = -6, 
           yend = 2,
           colour = lifespans$colour,
           size = 0.5) +
  labs(x = "Frequency (1/minutes)", y = "log amplitude") + 
  scale_x_continuous(trans = 'log10', breaks =  c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1), 
                     labels = c("0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1")) +
  scale_y_continuous(expand = c(0,0), limits = c(-6, 2), breaks = c(-6:2)) +
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
  geom_line() +
  annotation_custom(gg_frog) +
  annotation_custom(gg_dragonfly) +
  annotation_custom(gg_butterfly) +
  annotation_custom(gg_reptile) +
  annotation_custom(gg_insects)

## manually create legend 
x <- c(55:62)

gg_legend <- ggplot() + 
  geom_blank() +
  theme_void() +
  theme(panel.spacing.x=unit(1, "lines")) +
  scale_y_continuous(limits = c(75, 75)) +
  scale_x_continuous(limits = c(50,90)) +
  geom_line(aes(x=x, y=75, colour = x)) +
  scale_color_gradientn(colours = unique(lifespans$colour)) +
  annotate("text", label = "Lifespan of species in community", x = 72, y = 75) +
  theme(legend.position = "none")


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

gg_complete <- grid.arrange(gg_withpics, gg_legend, layout_matrix = lay)

ggsave(gg_complete, filename = "power-spectrum-with-lifespan_Loblolly-Marsh_with-pics.png", path = "./figures", dpi = 300, device = "png", height = 7.47, width = 10.369152)
