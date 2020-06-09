## trying to make a power spectrum graph - woohoo!!!
library(tidyverse)
library(imputeTS)
library(spectral)
library(lubridate)
library(gridExtra)
## time series data: from Ocean Networks Canada


## ---------------------------------------------------------------
## Barkley Canyon Mid-East benthic marine community 
## expect little yearly and seasonal temperature fluctuations since it is so deep and cold
## from 2009-09-15 to 2014-12-31
## presence-absence data collected in Barkley Canyon by Matabos et al. 2014
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


## plot the series:
ggplot(series, aes(x=as.Date(time), y=temp_val)) +
  geom_line() + 
  xlab("") 
## uh oh, gaps in the series!!!!! 

## ______
## inspect them to decide what to do
plotNA.gapsize(series$temp_val)
## seems like gaps are uncommon but when they occur they span a relatively large amount of time
plotNA.distribution(series$temp_val)
plotNA.distributionBar(series$temp_val, breaks = 30)
## looks like they're mostly in the first part of the data
## ______


## chop that part off evenly:
start <- which(as.Date(series$time) == "2011-01-01")
stop <- which(as.Date(series$time) == "2015-01-01")
series <- series[start[1]:stop[1],]

## ______
## check out the NAs now:
plotNA.distribution(series$temp_val)
plotNA.gapsize(series$temp_val)
plotNA.distributionBar(series$temp_val, breaks = 30)
## ______

## impute the missing values following methods from https://arxiv.org/pdf/1510.03924.pdf
## see how autocorrelated data are so we can choose the correct imputation method. if very highly autocorrelated, future temperature values rely a lot on the past values 
acf(series$temp_val, na.action = na.pass)
## a strong positive autocorrelation

## impute missing values in the time series to allow transformation  
BC = na_kalman(series$temp_val)
BC <- readRDS("./data-processed/BC.rds") ## saved RDS


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

ts_plotdata_BC <- data.frame(time, BC)

ts_BC <- ggplot(ts_plotdata_BC, aes(x = time, y = BC)) +
  geom_line() +
  scale_x_continuous(breaks = c(2011:2015), 
                     labels = c("2011", "2012", "2013", "2014", "2015")) +
  labs(x = "Year", y = "Temperature (°C)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))


ggsave(ts_BC, filename = "/time-series_BarkleyCanyon.png", path = "./figures" , dpi = 300, device = "png", width = 10.71, height = 6.47)



## capture length 
L <- length(BC)

# Fourier transform the series: 
dft_BC <- fft(BC)/L

amp_BC <- sqrt(Im(dft_BC[-1])^2 + Re(dft_BC[-1])^2)	# omit first term (represents DC component - y axis shift)
amp_BC <- amp_BC[1:(L/2)]				# remove second half of amplitudes (negative half)

frequency_BC <- 1:(L/2)/L ## sampling frequency (1 minute, 2 minutes, 3 minutes etc..) / length of time series 

plotdata_BC <- data.frame(frequency_BC, amp_BC)

# plot magnitudes against frequencies
g_BC = ggplot(plotdata_BC, aes(x = frequency_BC, y = log10(amp_BC))) + 
  geom_line() +  
  labs(x = "Frequency (1/minutes)", y = "log amplitude") + 
  scale_x_continuous(trans = 'log10', breaks =  c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1), 
                     labels = c("0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1")) +
  scale_y_continuous(expand = c(0,0), limits = c(-9, -1), breaks = c(-8:-2)) + 
  annotate("text", label = "1 day", x = 1/(60*24), y = -8, size = 3) +
  annotate("segment", x = 1/(60*24), xend = 1/(60*24), y = -9, yend = -8.25) +
  annotate("text", label = "1 year", x = 1/(60*24*365), y = -8, size = 3) + 
  annotate("segment", x = 1/(60*24*365), xend = 1/(60*24*365), y = -9, yend = -8.25) +
  annotate("text", label = "6 hours", x = 1/(60*6), y = -8, size = 3) + 
  annotate("segment", x = 1/(60*6), xend = 1/(60*6), y = -9, yend = -8.25) +
  annotate("text", label = "1 week", x = 1/(60*24*7), y = -8, size = 3) + 
  annotate("segment", x = 1/(60*24*7), xend = 1/(60*24*7), y = -9, yend = -8.25) +
  annotate("text", label = "1 month", x = 1/(60*24*365/12), y = -8, size = 3) +
  annotate("segment", x = 1/(60*24*365/12), xend = 1/(60*24*365/12), y = -9, yend = -8.25) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))


ggsave(g_BC, filename = "power-spectrum_BarkleyCanyon.png", path = "./figures", dpi = 300, device = "png", width = 10.71, height = 6.47)







## ---------------------------------------------------------------
## Folger Pinnacle reef community:
## expect more yearly and seasonal temperature fluctuations since it is less deep and more exposed to tides
## from 2012-10-19 to 2014-12-31
filenames <-c("20121019T225611Z_20121025T220554Z-NaN_clean.csv",
              "20121025T220555Z_20121106T115245Z-NaN_clean.csv",
              "20121106T115246Z_20121118T013946Z-NaN_clean.csv",
              "20121118T013947Z_20121129T152635Z-NaN_clean.csv",
              "20121129T152636Z_20121211T051328Z-NaN_clean.csv",
              "20121211T051329Z_20121222T190018Z-NaN_clean.csv",
              "20121222T190019Z_20130103T084657Z-NaN_clean.csv",
              "20130103T084658Z_20130114T223439Z-NaN_clean.csv",
              "20130114T223440Z_20130126T122338Z-NaN_clean.csv",
              "20130126T122339Z_20130207T021237Z-NaN_clean.csv",
              "20130207T021238Z_20130218T160136Z-NaN_clean.csv",
              "20130218T160137Z_20130302T055034Z-NaN_clean.csv",
              "20130302T055035Z_20130313T193909Z-NaN_clean.csv",
              "20130313T193910Z_20130325T092808Z-NaN_clean.csv",
              "20130325T092809Z_20130405T231651Z-NaN_clean.csv",
              "20130405T231652Z_20130417T130550Z-NaN_clean.csv",
              "20130417T130551Z_20130429T025449Z-NaN_clean.csv",
              "20130429T025450Z_20130510T164348Z-NaN_clean.csv",
              "20130510T164349Z_20130522T063246Z-NaN_clean.csv",
              "20130522T063247Z_20130602T202109Z-NaN_clean.csv",
              "20130602T202110Z_20130614T101008Z-NaN_clean.csv",
              "20130614T101009Z_20130625T235906Z-NaN_clean.csv",
              "20130625T235907Z_20130707T134805Z-NaN_clean.csv",
              "20130707T134806Z_20130719T033643Z-NaN_clean.csv",
              "20130719T033644Z_20130730T172336Z-NaN_clean.csv",
              "20130730T172337Z_20130811T071032Z-NaN_clean.csv",
              "20130811T071033Z_20130822T205736Z-NaN_clean.csv",
              "20130822T205737Z_20130903T104440Z-NaN_clean.csv",
              "20130903T104441Z_20130915T003159Z-NaN_clean.csv",
              "20130915T003200Z_20130926T141928Z-NaN_clean.csv",
              "20130926T141929Z_20131008T040702Z-NaN_clean.csv",
              "20131008T040703Z_20131019T175430Z-NaN_clean.csv",
              "20131019T175431Z_20131031T074217Z-NaN_clean.csv",
              "20131031T074218Z_20131111T213116Z-NaN_clean.csv",
              "20131111T213117Z_20131123T112015Z-NaN_clean.csv",
              "20131123T112016Z_20131205T010913Z-NaN_clean.csv",
              "20131205T010914Z_20131216T145613Z-NaN_clean.csv",
              "20131216T145614Z_20131228T044259Z-NaN_clean.csv",
              "20131228T044300Z_20140108T182945Z-NaN_clean.csv",
              "20140108T182946Z_20140120T081632Z-NaN_clean.csv",
              "20140120T081633Z_20140131T220319Z-NaN_clean.csv",
              "20140131T220320Z_20140212T115005Z-NaN_clean.csv",
              "20140212T115006Z_20140224T013653Z-NaN_clean.csv",
              "20140224T013654Z_20140307T152339Z-NaN_clean.csv",
              "20140307T152340Z_20140319T051025Z-NaN_clean.csv",
              "20140319T051026Z_20140330T185715Z-NaN_clean.csv",
              "20140330T185716Z_20140411T084409Z-NaN_clean.csv",
              "20140411T084410Z_20140422T223052Z-NaN_clean.csv",
              "20140422T223053Z_20140504T121738Z-NaN_clean.csv",
              "20140504T121739Z_20140516T020424Z-NaN_clean.csv",
              "20140516T020425Z_20140527T155113Z-NaN_clean.csv",
              "20140527T155114Z_20140608T053759Z-NaN_clean.csv",
              "20140608T053800Z_20140619T192448Z-NaN_clean.csv",
              "20140619T192449Z_20140701T091018Z-NaN_clean.csv",
              "20140701T091019Z_20140712T225705Z-NaN_clean.csv",
              "20140712T225706Z_20140804T105740Z-NaN_clean.csv",
              "20140804T105755Z_20141231T235955Z-NaN_clean.csv"
)

## combine all files into a single data frame
megafile <- data.frame(matrix(ncol = 3))
colnames(megafile) <- c("time", "temp_val", "qc") 

file  = 1
while(file < length(filenames) + 1) {
  current <- read.csv(paste("/Volumes/TimeMachine/search13905844/FolgerPassage_FolgerPinnacle_ConductivityTemperatureDepth_Temperature_", filenames[file], sep = ""), skip = 1, header = T, 
                      comment.char = "#")
  
  colnames(current) <- c("time", "temp_val", "qc")
  megafile <- rbind(megafile, current)
  file = file + 1
}

## save
megafile <- readRDS("./data-processed/megafile.rds")

series <- megafile[,1:2] 
series <- series[-1:-195,] ## remove elements necessary to get to an even minute mark

rm(current)
rm(megafile)

time_minutes <- substr(series$time, 1 , nchar(series$time)-8)## get rid of decimal places representing seconds in the time vector so all minutes have same time chr

## group by minute and find mean of each group of temp_vals 
series_new <- series %>%
  mutate(time = time_minutes) %>%
  group_by(time) %>%
  mutate(temp_val = mean(temp_val)) %>%
  ungroup() %>%
  filter(!duplicated(time))

series = readRDS("./data-processed/series_PR")

start <- which(as.Date(series$time) == "2012-12-31")
stop <- which(as.Date(series$time) == "2014-12-31")
series <- series[start[1]:stop[1],]


## impute missing values in the time series to allow transformation  
PR = na_seadec(series$temp_val, find_frequency = TRUE, algorithm = "interpolation")
PR <- readRDS("./data-processed/PR.rds") ## saved RDS



## plot the series:
date.num <- decimal_date(as.Date(series$time)) 
split <- str_split_fixed(series$time, pattern = "T", n = 2)
time <- split[,2]
time <- str_split_fixed(time, pattern = ":", n = 2) ## split into hours + minutes

hour <- as.numeric(time[,1])/24
minute <- (as.numeric(time[,2]) + hour*60) / (60*24)

time.num <-  hour + minute  ## in fraction of a day where 1 day = 1.0
time.num <- time.num / 365.25 ## convert to fraction of a year

options(digits=12)
time <- as.numeric(format(date.num + time.num, nsmall = 8))

ts_plotdata_PR <- data.frame(time, PR)

ts_PR <- ggplot(ts_plotdata_PR, aes(x = time, y = PR)) +
  geom_line() +
  scale_x_continuous(limits = c(2013, 2015),
                     breaks = c(2013:2015), 
                     labels = c("2013", "2014", "2015")) +
  labs(x = "Year", y = "Temperature (°C)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(ts_PR, filename = "/time-series_FolgerPinnacle.png", path = "./figures" , dpi = 300, device = "png", width = 10.71, height = 6.47)


L = length(PR)

# Fourier transform the series: 
dft_PR <- fft(PR)/L
dft_PR <- readRDS("./data-processed/dft_PR.rds") ## saved RDS

L <- length(dft_PR)

amp_PR <- sqrt(Im(dft_PR[-1])^2 + Re(dft_PR[-1])^2)	# omit first term (represents DC component - y axis shift)
amp_PR <- amp_PR[1:(L/2)]				# remove second half of amplitudes (negative half)

frequency_PR <- 1:(L/2)/L ## sampling frequency (1 minute, 2 minutes, 3 minutes etc..) / length of time series 

plotdata_PR <- data.frame(frequency_PR, amp_PR)

# plot magnitudes against frequencies
g_PR = ggplot(plotdata_PR, aes(x = frequency_PR, y = log(amp_PR))) + 
  geom_line() +  
  labs(x = "Frequency (1/minutes)", y = "log amplitude") +
  scale_x_continuous(trans = 'log10', limits = c(0.0000007, 0.05), 
                     breaks =  c(0.000001, 0.00001, 0.0001, 0.001, 0.01), 
                     labels = c("0.000001", "0.00001", "0.0001", "0.001",
                                "0.01")) +
  scale_y_continuous(expand = c(0,0), limits = c(-22.5, 5), breaks = seq(-20, 5, 5)) +
  annotate("text", label = "1 day", x = 1/(60*24), y = -19.25, size = 3) +
  annotate("segment", x = 1/(60*24), xend = 1/(60*24), y = -22.5, yend = -20) +
  annotate("text", label = "1 year", x = 1/(60*24*365), y = -19.25, size = 3) + 
  annotate("segment", x = 1/(60*24*365), xend = 1/(60*24*365), y = -22.5, yend = -20) +
  annotate("text", label = "6 hours", x = 1/(60*6), y = -19.25, size = 3) + 
  annotate("segment", x = 1/(60*6), xend = 1/(60*6), y = -22.5, yend = -20) +
  annotate("text", label = "1 week", x = 1/(60*24*7), y = -19.25, size = 3) + 
  annotate("segment", x = 1/(60*24*7), xend = 1/(60*24*7), y = -22.5, yend = -20) +
  annotate("text", label = "1 month", x = 1/(60*24*365/12), y = -19.25, size = 3) +
  annotate("segment", x = 1/(60*24*365/12), xend = 1/(60*24*365/12), y = -22.5, yend = -20) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))



ggsave(g_PR, filename = "power-spectrum_FolgerPinnacle.png", path = "./figures", dpi = 300, device = "png", width = 10.71, height = 6.47)






### creating graph with lifespans:
BC_species <- read.csv("./data-raw/BarkleyCanyon_species-list.csv")

species <- BC_species %>%
  filter(lifespan_days != "unk") %>%
  mutate(lifespan_minutes = as.numeric(as.character(lifespan_days))*24*60) %>% 
  mutate(lifespan_as_frequency = 1/lifespan_minutes) %>%
  mutate(species = as.character(species), 
         genus = as.character(genus), 
         higher_tax_extracted = as.character(higher_tax_extracted)) %>%
  mutate(species = ifelse(is.na(species), species_extracted, species)) %>%
  mutate(genus = ifelse(is.na(genus), genus_extracted, genus)) %>%
  mutate(species = ifelse(is.na(species), "", species)) %>%
  mutate(genus = ifelse(is.na(genus), higher_tax_extracted, genus)) 

## visualize
ggplot(species, aes(x = reorder(lifespan_days_source, lifespan_as_frequency), 
                    y = lifespan_as_frequency)) + 
  geom_col() + 
  theme_bw() +
  labs(x = "Taxonomic classification", y = "Lifespan (1/minutes)") +
  theme(axis.text.x = element_text(angle = 90, size = 7)) + 
  scale_x_discrete(labels = c(paste(species$genus, species$species)))
## uh oh, many have smaller frequecies (larger lifepsans) than the time series data we have 

## plot with all lifespans displayed 
g_BC_lifespan <- g_BC + 
  annotate("segment", x = species$lifespan_as_frequency, 
           xend = species$lifespan_as_frequency,
           y = -9, 
           yend = -8.25,
           colour = "red") 

fitonplot <- species$lifespan_as_frequency[which(species$lifespan_as_frequency > min(plotdata_BC$frequency_BC))]

## plot with only lifespans that fit
g_BC_lifespan_actual <- g_BC + 
  annotate("segment", x = fitonplot, 
           xend = fitonplot,
           y = -9, 
           yend = -8.25,
           colour = "red") 

## manually create legend 
gg <- ggplot() + 
  geom_blank() +
  theme_void() +
  theme(panel.spacing.x=unit(1, "lines")) +
  scale_y_continuous(limits = c(75, 75)) +
  scale_x_continuous(limits = c(50,90)) +
  annotate("segment", x = 62, xend = 64, y = 75, yend = 75, colour = "red") +
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

g1 <- grid.arrange(g_BC_lifespan_actual, gg, layout_matrix = lay)
g2 <- grid.arrange(g_BC_lifespan, gg, layout_matrix = lay)

## save 
ggsave(g1, filename = "power-spectrum-with-lifespan_BarkleyCanyon.png", path = "./figures", dpi = 300, device = "png", height = 6.47, width = 8.369152)
ggsave(g2, filename = "power-spectrum-with-lifespan_all_BarkleyCanyon.png", path = "./figures", dpi = 300, device = "png", height = 6.47, width = 8.369152)
