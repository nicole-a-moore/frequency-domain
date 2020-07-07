# getting some species lifespans in an automated fashion for Loblolly Marsh
library(tidyverse)
library(taxize)

species <- read.csv("./data-raw/species-list_LoblollyMarsh.csv")

species <- species %>%
  mutate(genus_species = paste(genus, species, sep = " ")) 

## taxize names:
taxa <- data.frame(taxa = species$genus_species) ## create dataframe of names to check

tsn_search <- get_tsn(as.character(taxa$taxa), accepted = FALSE) ## find tsn for each unique taxa
tsn_search <- readRDS("./data-processed/tsn-search_loblolly.rds")
tsns <- data.frame(tsn_search)
tsns$taxa <- taxa$taxa
taxa <- tsns

found <- taxa %>%
  subset(match == "found") 

report <- lapply(found$ids, itis_acceptname)
report_df <- data.frame(matrix(unlist(report), nrow=78, byrow=T),stringsAsFactors=FALSE)
report_df <- readRDS("./data-processed/report-df_loblolly.rds")

found <- found %>%
  mutate(genus_species_corrected = report_df$X2)

## merge short unique list to long list of all taxa
merged_unique <- left_join(taxa, found)
merged <- left_join(taxa, merged_unique)
merged$taxa <- as.character(merged$taxa)


## if names found are not accepted names, then change to accepted name
i = 1
while (i < length(merged$taxa)) {
  if (!is.na(merged$genus_species_corrected[i])) {
    merged$taxa[i] <- merged$genus_species_corrected[i]
  }
  i = i+1
}

## create new genus and species columns, correct original dataset
split <- str_split_fixed(merged$taxa, pattern = " ", n = 2)
merged$genus <- split[,1]
merged$species <- split[,2]

## update the database :)
species <- species %>%
  ungroup()%>%
  mutate(genus = merged$genus) %>%
  mutate(species = merged$species) %>%
  mutate(genus_species = merged$taxa) 


amphibs <- species %>%
  filter(group == "Amphibia") %>%
  mutate(genus_species = paste(genus, species, sep = " ")) 

amphibio <- read.csv("~/Documents/intra-therm/AmphiBIO_v1/AmphiBIO_v1.csv")

found <- amphibs[which(amphibs$genus_species %in% amphibio$Species),] %>%
  droplevels()

lifespan <- amphibio %>%
  subset(Species %in% found$genus_species) %>%
  select(Species, Longevity_max_y) %>%
  rename(genus_species = Species) %>%
  mutate(source = ifelse(!is.na(Longevity_max_y), "Amphibio", ""))

species <- left_join(species, lifespan) %>%
  rename(lifespan = Longevity_max_y) %>%
  mutate(lifespan = lifespan*365)


# search intratherm traits for species
intratherm <- read.csv("~/Documents/intra-therm/data-processed/intratherm-may-2020-squeaky-clean.csv")

intratherm <- intratherm %>%
  select(genus_species, lifespan_days, lifespan_days_reference) %>%
  unique()

intratherm <- intratherm[which(intratherm$genus_species %in% species$genus_species),]

species <- species %>%
  left_join(., intratherm) %>%
  mutate(source = ifelse(is.na(lifespan), as.character(lifespan_days_reference), source)) %>%
  mutate(lifespan = ifelse(is.na(lifespan), as.character(lifespan_days), lifespan)) %>%
  select(-lifespan_days_reference, -lifespan_days)

## write out to file to fill in the rest  
write.csv(species, "./data-processed/species-list_LoblollyMarsh_after-script.csv", row.names = FALSE)
