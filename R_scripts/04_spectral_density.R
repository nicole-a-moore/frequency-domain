#### Version of colour-of-noise script "04_spectral_density.R" using Loblolly Marsh data
library(tidyverse)
library(cowplot)
library(gridExtra)
library(grid)
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

##ggsave("./figures/sst-correlogram-eighth_LoblollyMarsh.png", width = 6, height = 3)


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




#####################################################################################
## making new version of spectral plot with real lifespan data, labels and pictures
#####################################################################################

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

## choose colours
cols <- c("#D7191C", "#FDAE61", "#66C2A5", "#E6F598", "#3288BD")

cols <- data.frame(group = unique(lifespans$group)) %>%
  mutate(colour = cols)
  ##mutate(colour = c("orange", "darkolivegreen2", "plum", "cadetblue3", "pink2"))
  

lifespans <- left_join(lifespans, cols)

## read in and properly colour images
frog <- as.raster(png::readPNG("./data-raw/frog.png")) 
frog[frog != "#FFFFFF"] <- lifespans$colour[first(which(lifespans$group=="Amphibia"))]
frog[frog == "#FFFFFF"] <- "#FFFFFF00"
gg_frog <- rasterGrob(frog, interpolate=TRUE, x = .225, y = .9, height = 0.15, width = .1)

dragonfly <- as.raster(png::readPNG("./data-raw/dragonfly.png"))
dragonfly[dragonfly != "#00000000"] <- lifespans$colour[first(which(lifespans$group=="Damselflies & Dragonflies"))]
gg_dragonfly <- rasterGrob(dragonfly, x = .65, y = .9, height = 0.15, width = .15)

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
gg_butterfly <- rasterGrob(butterfly, interpolate=TRUE, x = .35, y = .9, height = 0.15, width = .1)

reptile <- as.raster(png::readPNG("./data-raw/reptile.png")) 
reptile[reptile != "#00000000"] <- lifespans$colour[first(which(lifespans$group=="Reptilia"))]
reptile[reptile == "#00000000"] <- "#FFFFFF00"
gg_reptile <- rasterGrob(reptile, interpolate=TRUE, x = .1, y = .875, height = 0.2, width = .09)

insects <- as.raster(png::readPNG("./data-raw/insect.png")) 
insects[insects != "#FFFFFF"] <- lifespans$colour[first(which(lifespans$group=="Other insects"))]
insects[insects == "#FFFFFF"] <- "#FFFFFF00"
gg_insects <- rasterGrob(insects, interpolate=TRUE, x = .5, y = .9, height = 0.125, width = .125)


## plot 
gg_spectral <- ggplot(spectral, aes(x = specx, y = specy)) + 
  annotate("segment", x = lifespans$lifespan_frequency, 
           xend = lifespans$lifespan_frequency,
           y = 0, 
           yend = 100000000,
           colour = lifespans$colour,
           size = 0.5) + 
  geom_line() +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), limits = c(0.01,100000000), expand = c(0,0)) +
  scale_x_log10(breaks =  c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1), 
                labels = c("0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1")) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  ylab("Spectral density") +
  xlab("Frequency (1/minutes)") +
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

gg_complete <- grid.arrange(gg_spectral, gg_legend, layout_matrix = lay)

##ggsave(gg_complete, path = "./figures", filename = "sst-spectral-slope-normalred-noise_LoblollyMarsh_fancy.png", width = 6, height = 4, device = "png")



## another version of plot with some lifespans labelled
## plot 
gg_insects <- rasterGrob(insects, interpolate=TRUE, x = .35, y = .825, height = 0.1, width = .1)
gg_butterfly <- rasterGrob(butterfly, interpolate=TRUE, x = .5, y = .825, height = 0.125, width = .08)
gg_frog <- rasterGrob(frog, interpolate=TRUE, x = .225, y = .825, height = 0.125, width = .08)
gg_dragonfly <- rasterGrob(dragonfly, x = .65, y = .825, height = 0.125, width = .125)
gg_reptile <- rasterGrob(reptile, interpolate=TRUE, x = .1, y = .8, height = 0.175, width = .075)

## 365: 1 year - 17
## 6: 6 days - 54
## 30: 1 month - 52
## 3650: 10 years - 12
## 14600: 40 years - 2
## 170: 170 days - 27
ticks <- lifespans[c(54, 52, 27, 17, 12, 2),] %>%
  mutate(label = c("6 days", "1 month", "     ~1/2 year", "1 year   ", "10 years", "40 years"))

gg_spectral <- ggplot(spectral, aes(x = specx, y = specy)) + 
  annotate("segment", x = lifespans$lifespan_frequency, 
           xend = lifespans$lifespan_frequency,
           y = 0, 
           yend = 100000000,
           colour = lifespans$colour,
           size = 0.5) + 
  annotate("segment", x = ticks$lifespan_frequency, 
           xend = ticks$lifespan_frequency,
           y = 0, 
           yend = 400000000,
           colour = ticks$colour,
           size = 0.5) +
  annotate("text", label = ticks$label, 
           x = ticks$lifespan_frequency,
           y = 800000000,
           colour = ticks$colour,
           size = 3) +
  geom_line() +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), limits = c(0.01,1000000000), expand = c(0,0)) +
  scale_x_log10(breaks =  c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1), 
                labels = c("0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1")) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  ylab("Spectral density") +
  xlab("Frequency (1/minutes)") +
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

gg_complete <- grid.arrange(gg_spectral, gg_legend, layout_matrix = lay)

##ggsave(gg_complete, path = "./figures", filename = "sst-spectral-slope-normalred-noise_LoblollyMarsh_fancy-with-ticks.png", width = 6, height = 4, device = "png")
