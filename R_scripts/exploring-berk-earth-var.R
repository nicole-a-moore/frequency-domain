## visualizing variation in Berkeley Earth TMAX time series'
library(tidyverse)
library(imputeTS)
library(gridExtra)


## read in all of the TMAX temp data for populations in intratherm:
temps <- read_csv("~/Documents/arr/data-processed/arr_temp-data.csv")
temps <- as.data.frame(temps)


## visualize a couple time series to choose a diverse set to examine:
plotdata1 <- data.frame(time = temps$date, temp_val = temps[,600])
plotdata1<- plotdata1[1:50770,]
colnames(plotdata1) <- c("time", "temp_val")
plotdata2 <- data.frame(time = temps$date, temp_val = temps[,3])
plotdata2<- plotdata2[1:50770,]
colnames(plotdata2) <- c("time", "temp_val")

bothplots <- rbind(plotdata1, plotdata2) %>%
  mutate(location = "Location 3") 

bothplots$location[1:(nrow(plotdata1))] = "Location 600"

ts_combined <- ggplot(bothplots, aes(x = time, y = temp_val, colour = location)) +
  geom_line() +
  scale_colour_manual(values = c("red3", "black")) +
  scale_x_continuous(limits = c(1960, 1970)) +
  labs(x = "Year", y = "Temperature (°C)", colour = "Collection location:") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"),  
        legend.position = "right")

ggsave(ts_combined, filename = "time-series_BerkEarth_600and3.png", path = "./figures", dpi = 300, device = "png", width = 8.802532, height = 6.47)

ts_single <- ggplot(plotdata2, aes(x = time, y = temp_val)) +
  geom_line() +
  scale_x_continuous(limits = c(1960, 1970)) +
  labs(x = "Year", y = "Temperature (°C)") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(ts_single, filename = "time-series_BerkEarth_3.png", path = "./figures", dpi = 300, device = "png", width = 8.802532, height = 6.47)

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
colnames(plotdata_freq1) <- c("frequency", "amp")
colnames(plotdata_freq2) <- colnames(plotdata_freq1) 

bothplots_freq <- rbind(plotdata_freq1, plotdata_freq2) %>%
  mutate(location = "Location 3") 

bothplots_freq$location[1:(nrow(plotdata_freq1))] = "Location 600"


# plot magnitudes against frequencies
g <- ggplot(bothplots_freq, aes(x=frequency, y=log(amp), colour = location)) + 
  geom_line() +
  scale_colour_manual(values = c("red3", "black")) +
  labs(x = "Frequency (1/days)", y = "log amplitude", colour = "Collection location:") + 
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

ggsave(g, filename = "power-spectrum_BerkEarth_600and3.png", path = "./figures", dpi = 300, device = "png", width = 10.802532, height = 6.47)



## Parseval’s theorem: the variance of the time series equals the sum of the squared absolute values of the output of the DFT
get.total.SS = function(series){
  if(length(which(is.na(series))) >= 1) {
    series <- na_seadec(series, find_frequency = TRUE, algorithm = "interpolation")
  }
  
    L = length(series)
    # Fourier transform series and capture phases and amplitudes of series
    dft <-  fft(series) / L
    amp <- sqrt(Im(dft[-1]) ^ 2 + Re(dft[-1]) ^ 2)	# omit this first term
    amp <- amp[1:(L / 2)]				# snag just first half of amplitudes
    # squares of each
    ampt = (2 * amp) ^ 2 / 2
    X = sum(ampt, na.rm = TRUE)
    return(X)
  
}

## find out what portion of variance is seasonal, annual and interannual:
get.seasonal.SS = function(series){
  if(length(which(is.na(series))) >= 1) {
    series <- na_seadec(series, find_frequency = TRUE, algorithm = "interpolation")
  }
 
    L = length(series)
    
    dft<-  fft(series)/L
    amp<-sqrt(Im(dft[-1])^2+Re(dft[-1])^2)	
    amp<-amp[1:(L/2)]			
    
    annual <- which(amp == max(amp))
    
    amp <- amp[annual+1:length(amp)] ## subset to amplitudes between 1/365 and year and 1 day
    # squares of each
    ampt = (2*amp)^2/2
    X = sum(ampt, na.rm = TRUE)
    return(X)
  
}

get.annual.SS = function(series){
  if(length(which(is.na(series))) >= 1) {
    series <- na_seadec(series, find_frequency = TRUE, algorithm = "interpolation")
  }
  
    L = length(series)
    dft<-  fft(series)/L
    amp<-sqrt(Im(dft[-1])^2+Re(dft[-1])^2)	
    amp<-amp[1:(L/2)]			
    
    annual <- which(amp == max(amp))
    
    amp <- amp[annual] ## subset to amplitudes between 1/365 and year and 1 day
    # squares of each
    ampt = (2*amp)^2/2
    X = sum(ampt, na.rm = TRUE)
    return(X)
  
}

get.interannual.SS = function(series){
  if(length(which(is.na(series))) >= 1) {
    series <- na_seadec(series, find_frequency = TRUE, algorithm = "interpolation")
  }
  
    L = length(series)
    
    dft<-  fft(series)/L
    amp<-sqrt(Im(dft[-1])^2+Re(dft[-1])^2)	
    amp<-amp[1:(L/2)]			
    
    annual <- which(amp == max(amp))
    
    amp <- amp[1:(annual-1)] ## subset to amplitudes between 1/365 and year and 1 day
    # squares of each
    ampt = (2*amp)^2/2
    X = sum(ampt, na.rm=TRUE)
    return(X)
  
}

get.halfyearly.SS = function(series){
  if(length(which(is.na(series))) >= 1) {
    series <- na_seadec(series, find_frequency = TRUE, algorithm = "interpolation")
  }

    L = length(series)
    
    dft<-  fft(series)/L
    amp<-sqrt(Im(dft[-1])^2+Re(dft[-1])^2)	
    amp<-amp[1:(L/2)]			
    
    annual <- which(amp == max(amp))
    
    amp <- amp[annual*2] ## subset to amplitudes between 1/365 and year and 1 day
    # squares of each
    ampt = (2*amp)^2/2
    X = sum(ampt, na.rm=TRUE)
    return(X)
}

var1 = get.total.SS(plotdata1$temp_val) ## 600
var2 = get.total.SS(plotdata2$temp_val) ## 3

seasonalvar1 = get.seasonal.SS(plotdata1$temp_val)
seasonalvar2 = get.seasonal.SS(plotdata2$temp_val)

annualvar1 = get.annual.SS(plotdata1$temp_val)
annualvar2 = get.annual.SS(plotdata2$temp_val)

interannualvar1 = get.interannual.SS(plotdata1$temp_val)
interannualvar2 = get.interannual.SS(plotdata2$temp_val)

ratio1s = seasonalvar1/var1
ratio2s = seasonalvar2/var2

ratio1a = annualvar1/var1
ratio2a = annualvar2/var2

ratio1ia = interannualvar1/var1
ratio2ia = interannualvar2/var2

sd1 <- sd(plotdata1$temp_val)
sd2 <- sd(plotdata2$temp_val)



## visualize:
compare <- data.frame(location = c("Location 600","Location 600","Location 600", 
                                   "Location 3","Location 3","Location 3"),
                      var_type = c("Seasonal", "Annual", "Interannual"),
                      proportion = c(ratio1s, ratio1a, ratio1ia, ratio2s, ratio2a, ratio2ia),
                      sd = c(sd1, sd1, sd1, sd2, sd2, sd2),
                      var = c(var1,var1,var1,var2,var2,var2))

c <- ggplot(compare, aes(x = location, y = proportion, fill = var_type)) +
  geom_col(width = .60) + 
  scale_fill_manual(values = hcl.colors(n=3)) +
  labs(x = "Collection location", y = "Proportion of total variance of temperature", 
       fill = "Type of variation:") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))

c2 <- ggplot(compare, aes(x = reorder(var_type, -proportion), y = proportion, fill = location)) +
  geom_col(width = .60, position = "dodge") + 
  scale_fill_manual(values = c("red3", "black")) +
  labs(x = "Collection location", y = "Proportion of total variance of temperature", 
       fill = "Type of variation:") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))
              
sd <- ggplot(compare, aes(x = location, y = sd, fill = location)) +
  geom_col(width = .60) + 
  scale_fill_manual(values = c("red3", "black")) +
  labs(x = "Collection location", y = "Standard deviation in temperature") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none")

var <- ggplot(compare, aes(x = location, y = var, fill = location)) +
  geom_col(width = .60) + 
  scale_fill_manual(values = c("red3", "black")) +
  labs(x = "Collection location", y = "Total variance in temperature") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none") 
  
both <- grid.arrange(sd, var, ncol = 2)

ggsave(both, filename = "BerkEarth_600and3_sd-var.png", path = "./figures", dpi = 300, device = "png", width = 6.47*1.863265, height = 6.47)

ggsave(c, filename = "BerkEarth_600and3_variation-breakdown.png", path = "./figures", dpi = 300, device = "png", width = 6.47*1.625541, height = 6.47)

ggsave(c2, filename = "BerkEarth_600and3_variation-comparison.png", path = "./figures", dpi = 300, device = "png", width = 6.47*1.625541, height = 6.47)




## interesting figures to explore:
##  - latitude vs proportion of variation that is seasonal, annual, interannual 
##  - standard deviation/experienced variaation in temp values vs proportion of variation that is seasonal, annual, interannual
##  - ARR vs proportion of variation that is seasonal, annual, interannual 





## does proportion of seasonal variation predict ARR?
## calculate prop of seasonal variation:
total <- c()
seasonal <- c()
annual <- c()
halfyearly <- c()
interannual <- c()
  
i = 1
while (i < ncol(temps)) {
  loc <- temps[1:50770,i+1]
  
  if(length(which(is.na(loc))) > (51042/5) & length(which(!is.na(loc))) > 3) {
    loc <- na_seadec(loc, find_frequency = TRUE, algorithm = "interpolation")
    total <- append(total, get.total.SS(loc), after = length(total))
    seasonal <- append(seasonal, get.seasonal.SS(loc), after = length(seasonal))
    annual <- append(annual, get.annual.SS(loc), after = length(annual))
    halfyearly <- append(halfyearly, get.halfyearly.SS(loc), after = length(halfyearly))
    interannual <- append(interannual, get.interannual.SS(loc), after = length(interannual))
    
    i = i + 1
  }
  else if (length(which(is.na(loc))) > (51042/5) & length(which(!is.na(loc))) < 3) {
    total <- append(total, NA, after = length(total))
    seasonal <- append(seasonal, NA, after = length(seasonal))
    annual <- append(annual, NA, after = length(annual))
    halfyearly <- append(halfyearly, NA, after = length(halfyearly))
    interannual <- append(interannual, NA, after = length(interannual))
    
    i = i + 1
  }
  else {
    total <- append(total, get.total.SS(loc), after = length(total))
    seasonal <- append(seasonal, get.seasonal.SS(loc), after = length(seasonal))
    annual <- append(annual, get.annual.SS(loc), after = length(annual))
    halfyearly <- append(halfyearly, get.halfyearly.SS(loc), after = length(halfyearly))
    interannual <- append(interannual, get.interannual.SS(loc), after = length(interannual))
    
    i = i + 1
  }
}

output <- data.frame(population_id = colnames(temps)[2:ncol(temps)], seasonal, 
                     annual, halfyearly, interannual, total) %>%
  mutate(proportion_seasonal = seasonal/total) %>%
  mutate(proportion_annual = annual/total) %>%
  mutate(proportion_halfyearly = halfyearly/total) %>%
  mutate(proportion_interannual = interannual/total)
  
intratherm <- read.csv("~/Documents/arr/data-processed/arr_sliding-window-output-with-arr.csv")

arr <- intratherm %>%
  select(population_id, genus_species, latitude, class, realm_general2, ARR, season_when_away_100km, season_when_away_10km, season_inactive, lifespan_days)

merged <- left_join(output, arr) %>%
  filter(!duplicated(population_id))

merged %>%
  filter(realm_general2 == "Terrestrial") %>%
  ggplot(., aes(x = proportion_annual, y = ARR)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "none")


                         

### QUESTION: 
## are terrestrial species that experience more variation/seasonal variation in max temps more often behavioural thermoregulators?

## create logical column for behavioural thermoregulator (TRUE/FALSE) 
nku <- c("unk", "kunk", "none")
inactivity <- merged %>%
  filter(!duplicated(genus_species)) %>%
  mutate(is.nku1 = !as.character(season_when_away_100km) %in% nku & !is.na(season_when_away_100km)) %>%
  mutate(is.nku2 = !as.character(season_when_away_10km) %in% nku & !is.na(season_when_away_10km)) %>%
  mutate(is.nku3 = !as.character(season_inactive) %in% nku & !is.na(season_inactive)) %>%
  mutate(behaviourally_thermoregulates = ifelse(is.nku1 == TRUE, TRUE, 
                                                ifelse(is.nku2== TRUE, TRUE,
                                                       ifelse(is.nku3 == TRUE, TRUE, FALSE)))) %>%
  select(-is.nku1, -is.nku2, -is.nku3)

## variation in general:
behav1_1 <- inactivity %>%
  filter(realm_general2 == "Terrestrial") %>%
  ggplot(., aes(x = behaviourally_thermoregulates, y = total, colour = behaviourally_thermoregulates)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = "Total thermal variance at collection location (°C)", x = "Behaviourally thermoregulates") + 
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

t_totalvar <- t.test(x = inactivity$total, y = inactivity$behaviourally_thermoregulates)

## can we go further and ask whether terrestrial species that experience more annual variation (high summer highs, low winter lows) as opposed to more seasonal or interannual variation more often behavioural thermoregulators?
behav2_1 <- inactivity %>%
  filter(realm_general2 == "Terrestrial") %>%
  ggplot(., aes(x = behaviourally_thermoregulates, y = proportion_annual, colour = behaviourally_thermoregulates)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = "Proportion of total thermal variance accounted for by annual cycle", 
       x = "Behaviourally thermoregulates") + 
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

t_annualvar <- t.test(x = inactivity$proportion_annual, y = inactivity$behaviourally_thermoregulates)

## compare to seasonal and interannual
inactivity %>%
  filter(realm_general2 == "Terrestrial") %>%
  ggplot(., aes(x = behaviourally_thermoregulates, y = proportion_interannual, colour = behaviourally_thermoregulates)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = "Proportion of total thermal variance accounted for by interannual cycle", 
       x = "Behaviourally thermoregulates") + 
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

t_interannualvar <- t.test(x = inactivity$proportion_interannual, y = inactivity$behaviourally_thermoregulates)

inactivity %>%
  filter(realm_general2 == "Terrestrial") %>%
  ggplot(., aes(x = behaviourally_thermoregulates, y = proportion_seasonal, colour = behaviourally_thermoregulates)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = "Proportion of total thermal variance accounted for by seasonal cycle", 
       x = "Behaviourally thermoregulates") + 
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

t_seasonalvar <- t.test(x = inactivity$proportion_seasonal, y = inactivity$behaviourally_thermoregulates)


behav <- grid.arrange(behav1_1, behav2_1, ncol = 1)

ggsave(behav, width = 8.802532*0.7770992, height = 8.802532, filename = "BerkEarth_behaviouralthermo.png", path = "./figures", dpi = 300, device = "png")







## QUESTION:
## does changing the time period (ex. from 1880-2019 vs from 1990-2019) alter the proportion of variation that is annual vs seasonal?
## answer: yes a bit
total <- c()
seasonal <- c()
annual <- c()
halfyearly <- c()
interannual <- c()
index = which(temps$date > 1990.001)[1]

i = 1
while (i < ncol(temps)) {
  loc <- temps[index:50770,i+1]
  
  if(length(which(is.na(loc))) > (51042/5) & length(which(!is.na(loc))) > 3) {
    loc <- na_seadec(loc, find_frequency = TRUE, algorithm = "interpolation")
    total <- append(total, get.total.SS(loc), after = length(total))
    seasonal <- append(seasonal, get.seasonal.SS(loc), after = length(seasonal))
    annual <- append(annual, get.annual.SS(loc), after = length(annual))
    halfyearly <- append(halfyearly, get.halfyearly.SS(loc), after = length(halfyearly))
    interannual <- append(interannual, get.interannual.SS(loc), after = length(interannual))
    
    i = i + 1
  }
  else if (length(which(is.na(loc))) > (51042/5) & length(which(!is.na(loc))) < 3) {
    total <- append(total, NA, after = length(total))
    seasonal <- append(seasonal, NA, after = length(seasonal))
    annual <- append(annual, NA, after = length(annual))
    halfyearly <- append(halfyearly, NA, after = length(halfyearly))
    interannual <- append(interannual, NA, after = length(interannual))
    
    i = i + 1
  }
  else {
    total <- append(total, get.total.SS(loc), after = length(total))
    seasonal <- append(seasonal, get.seasonal.SS(loc), after = length(seasonal))
    annual <- append(annual, get.annual.SS(loc), after = length(annual))
    halfyearly <- append(halfyearly, get.halfyearly.SS(loc), after = length(halfyearly))
    interannual <- append(interannual, get.interannual.SS(loc), after = length(interannual))
    
    i = i + 1
  }
}

shortened <- data.frame(population_id = colnames(temps)[2:ncol(temps)], seasonal, 
                     annual, halfyearly, interannual, total) %>%
  mutate(proportion_seasonal = seasonal/total) %>%
  mutate(proportion_annual = annual/total) %>%
  mutate(proportion_halfyearly = halfyearly/total) %>%
  mutate(proportion_interannual = interannual/total)

short_merged <- left_join(shortened, arr) %>%
  filter(!duplicated(population_id))

nku <- c("unk", "kunk", "none")
inactivity2 <- short_merged %>%
  filter(!duplicated(genus_species)) %>%
  mutate(is.nku1 = !as.character(season_when_away_100km) %in% nku & !is.na(season_when_away_100km)) %>%
  mutate(is.nku2 = !as.character(season_when_away_10km) %in% nku & !is.na(season_when_away_10km)) %>%
  mutate(is.nku3 = !as.character(season_inactive) %in% nku & !is.na(season_inactive)) %>%
  mutate(behaviourally_thermoregulates = ifelse(is.nku1 == TRUE, TRUE, 
                                                ifelse(is.nku2== TRUE, TRUE,
                                                       ifelse(is.nku3 == TRUE, TRUE, FALSE)))) %>%
  select(-is.nku1, -is.nku2, -is.nku3)


behav1_2 <- inactivity2 %>%
  filter(realm_general2 == "Terrestrial") %>%
  ggplot(., aes(x = behaviourally_thermoregulates, y = total, colour = behaviourally_thermoregulates)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = "Total thermal variance over ~30 years \n at collection location (°C)", x = "Behaviourally thermoregulates") + 
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

behav2_2 <- inactivity2 %>%
  filter(realm_general2 == "Terrestrial") %>%
  ggplot(., aes(x = behaviourally_thermoregulates, y = proportion_annual, colour = behaviourally_thermoregulates)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = "Proportion of total thermal variance over ~30 years \n accounted for by annual cycle", 
       x = "Behaviourally thermoregulates") + 
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

behavplot2 <- grid.arrange(behav1_2, behav2_2, ncol = 1)

ggsave(behavplot2, width = 8.802532*0.7770992, height = 8.802532, filename = "BerkEarth_behaviouralthermo-30years.png", path = "./figures", dpi = 300, device = "png")


g <- grid.arrange(behav2_1, behav2_2, ncol = 1)
ggsave(g, width = 8.802532*0.7770992, height = 8.802532, filename = "BerkEarth_behaviouralthermo-comparison-annual.png", path = "./figures", dpi = 300, device = "png")
g <- grid.arrange(behav1_1, behav1_2, ncol = 1)
ggsave(g, width = 8.802532*0.7770992, height = 8.802532, filename = "BerkEarth_behaviouralthermo-comparison-total.png", path = "./figures", dpi = 300, device = "png")


L <- length(loc)

# Fourier transform the series: 
dft <- fft(loc)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2)	
amp <- amp[1:(L/2)]				

frequency <- 1:(L/2)/L ## sampling frequency of 1 day

plotdata_freq<- data.frame(frequency, amp)

g <- ggplot(plotdata_freq, aes(x=frequency, y=log(amp))) + 
  geom_line() +
  scale_colour_manual(values = c("red3", "black")) +
  labs(x = "Frequency (1/days)", y = "log amplitude", colour = "Collection location:") + 
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





## QUESTION: how much does the amplitude at 1 year change depending on the length of the time series used?
annual_all <- c()

loc <- temps[1:50770,600]

get.annual.power = function(series, i){
  if(length(which(is.na(series))) > 1) {
    series <- na_seadec(series, find_frequency = TRUE, algorithm = "interpolation")
  }
 
    L = length(series)
    dft<-  fft(series)/L
    amp<-sqrt(Im(dft[-1])^2+Re(dft[-1])^2)
    amp<-amp[1:(L/2)]	
    
    annual <- which(amp == max(amp))
    
    amp <- amp[annual] 
    return(amp)
  
}

i = 1
while (i < nrow(temps) + 1) {
  loc_short <- loc[i:length(loc)]
  annual_all <- append(annual_all, get.annual.power(loc_short, i), after = length(annual_all))
  
  i = i + 1
}


d <- data.frame(index = rep(1:length(annual_all)), amp = annual_all)
d <- readRDS("./data-processed/annual_all.rds")


g <- ggplot(d, aes(x = length(loc) - index, y = amp)) +
  geom_point() + 
  labs(x = "Number of days in time series", y = "Amplitude at 1 year") +
  scale_x_continuous(limits = c(2000,8000)) + 
  theme_bw() +
  theme(panel.grid.major = element_line(), panel.grid.minor = element_line())

ggsave(g, height = 8.802532/2, width = (8.802532*1.311111)/2, filename = "BerkEarth_ts-length-investigation.png", path = "./figures", dpi = 300, device = "png")






## trying with freshwater model data:

temps_fresh <- read_csv("~/Documents/intra-therm/data-processed/intratherm-freshwater-temp-data-daily.csv")
temps_fresh <- as.data.frame(temps_fresh)
total <- c()
seasonal <- c()
annual <- c()
halfyearly <- c()
interannual <- c()

i = 1
while (i < ncol(temps_fresh)) {
  loc <- temps_fresh[,i+1]
  
  if(length(which(is.na(loc))) > (16071/5) & length(which(!is.na(loc))) > 3) {
    loc <- na_seadec(loc, find_frequency = TRUE, algorithm = "interpolation")
    total <- append(total, get.total.SS(loc), after = length(total))
    seasonal <- append(seasonal, get.seasonal.SS(loc), after = length(seasonal))
    annual <- append(annual, get.annual.SS(loc), after = length(annual))
    halfyearly <- append(halfyearly, get.halfyearly.SS(loc), after = length(halfyearly))
    interannual <- append(interannual, get.interannual.SS(loc), after = length(interannual))
    
    i = i + 1
  }
  else if (length(which(is.na(loc))) > (16071/5) & length(which(!is.na(loc))) < 3) {
    total <- append(total, NA, after = length(total))
    seasonal <- append(seasonal, NA, after = length(seasonal))
    annual <- append(annual, NA, after = length(annual))
    halfyearly <- append(halfyearly, NA, after = length(halfyearly))
    interannual <- append(interannual, NA, after = length(interannual))
    
    i = i + 1
  }
  else {
    total <- append(total, get.total.SS(loc), after = length(total))
    seasonal <- append(seasonal, get.seasonal.SS(loc), after = length(seasonal))
    annual <- append(annual, get.annual.SS(loc), after = length(annual))
    halfyearly <- append(halfyearly, get.halfyearly.SS(loc), after = length(halfyearly))
    interannual <- append(interannual, get.interannual.SS(loc), after = length(interannual))
    
    i = i + 1
  }
}

output <- data.frame(population_id = colnames(temps_fresh)[2:ncol(temps_fresh)], seasonal, 
                     annual, halfyearly, interannual, total) %>%
  mutate(proportion_seasonal = seasonal/total) %>%
  mutate(proportion_annual = annual/total) %>%
  mutate(proportion_halfyearly = halfyearly/total) %>%
  mutate(proportion_interannual = interannual/total)

intratherm <- read.csv("~/Documents/arr/data-processed/arr_sliding-window-output-with-arr.csv")

intratherm$population_id <- paste(intratherm$genus_species, intratherm$latitude, intratherm$elevation, intratherm$longitude, sep = "_")

arr <- intratherm %>%
  select(population_id, genus_species, latitude, class, realm_general2, ARR, season_when_away_100km, season_when_away_10km, season_inactive, lifespan_days)

merged <- left_join(output, arr) %>%
  filter(!duplicated(population_id))

merged %>%
  ggplot(., aes(x = proportion_annual, y = ARR)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "none")


