library(countrycode)
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)
library(readxl)
library(writexl)
library(rotl)
library(openxlsx)
library(spatstat)
library(geosphere)
library(spThin)
library(sf)

########################################################################################
########################################################################################
#######R Script for Cleaning and Preparing Species Occurrence Records###################
########################################################################################
########################################################################################

#This script is designed to clean and prepare data on species occurrence. 
#The main output of the script is a list containing the coordinates, ott_id, and taxonomic order for each species occurrence record. 
#Subsequently, environmental variable values associated with each coordinate are extracted using ArcGIS.



########################################################################################
###########################     1) Coordinate cleaner     ##############################
########################################################################################


#In this part, the occurence coordinate are cleaned using Coordinate cleaner (https://cran.r-project.org/web/packages/CoordinateCleaner/index.html)
#Input: Original online database of species records (Here is an example for the SALVE online database)

dat<- read_xlsx("SALVE.xlsx") #Import

dat <- dat %>%filter(Ano >= 1990 & Ano <= 2022) #We only keep data between 1990 and 2022

#Series of coordinate check (e.g., coordinate system, is there a coordinate with the records...)
dat <- dat %>%filter(Datum=="WGS 84")
dat <- dat %>%filter(Precisao_da_coordenada=="Exata")
dat <- dat %>%
  dplyr::select(Nome_cientifico, Longitude, 
                Latitude, Country, Precisao_da_coordenada,
                Ano, Tipo_registro)
dat <- dat %>%
  filter(!is.na(Longitude)) %>%
  filter(!is.na(Latitude))

dat <- dat %>%
  filter(grepl("^-?\\d+(\\.\\d+)?$", Latitude) & grepl("^-?\\d+(\\.\\d+)?$", Longitude))

dat <- dat %>%
  mutate(
    Latitude = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  )

#Use of Coordinate cleaner to eliminate records with incorrect geo-referencing 
flags <- clean_coordinates(x = dat, 
                           lon = "Longitude", 
                           lat = "Latitude",
                           countries = "Country",
                           species = "Nome_cientifico",
                           tests = c("capitals", "centroids",
                                     "equal", "zeros", "countries")) # most test are on by default

plot(flags, lon = "Longitude", lat = "Latitude")
dat_cl <- dat[flags$.summary,]
dat_fl <- dat[!flags$.summary,]

#Export of the cleaned file 
write_xlsx(dat_cl,"SLAVE_CoordClean.xlsx")

########################################################################################
###########################     2) Taxonomic cleaning     ##############################
########################################################################################

#In this part, the occurence records are taxonimically cleaned using rotl (https://cran.r-project.org/web/packages/rotl/index.html)
#The inconsistency in standardization across online databases regarding taxonomic names could result in records of the same species being named differently. 
#To tackle this issue, the “rotl” R package was employed to align records with their taxonomic names as described in the Open Tree of Life (OTL), a digital phylogenetic tree including all organisms and their taxonomic details. 
#Records that could not be matched with an OTL taxonomic name were excluded. 
#Input:  Database of species records (Here is an example for the SALVE online database)



#Import of data
dat <- read_xlsx("SLAVE_CoordClean.xlsx")

#Extract Taxa name as provided by the online database
taxa <- as.list(dat[[2]])
taxa <- as.character(dat[[2]])
taxa <- unique(taxa)
taxa <- taxa[!is.na(taxa)]
taxa <- taxa[grepl("^[a-zA-Z ]+$", taxa)]
taxa <- trimws(taxa)

# Find the ott_id for each taxa
resolved_names <- tnrs_match_names(taxa, context_name = "Animals")
resolved_names <- resolved_names[!is.na(resolved_names$ott_id) & 
                                   resolved_names$number_matches <= 1, ]

#Function to get the taxonomic order
get_taxonomic_ranks <- function(ott_id) {
  # Attempt to retrieve taxonomic information
  tax_info <- tryCatch({
    taxonomy_taxon_info(ott_id, include_lineage = TRUE)
  }, error = function(e) {
    cat("Error retrieving taxonomic info for OTT ID:", ott_id, "\n")
    return(NULL)
  })
  
  if (is.null(tax_info) || is.null(tax_info[[1]]$lineage)) {
    cat("No lineage information for OTT ID:", ott_id, "\n")
    return(rep(NA, 9))  # Return NA if lineage is not available
  }
  
  lineage <- tax_info[[1]]$lineage
  
  ranks <- c("class" = NA, "subclass" = NA, "order" = NA, "suborder" = NA, 
             "superclass" = NA, "phylum" = NA, "subphylum" = NA, 
             "kingdom" = NA, "domain" = NA)
  
  for (taxon in lineage) {
    if (taxon$rank %in% names(ranks)) {
      ranks[taxon$rank] <- taxon$name
    }
  }
  
  return(ranks)
}

# Initialize columns for taxonomic ranks
resolved_names$class <- NA
resolved_names$subclass <- NA
resolved_names$order <- NA
resolved_names$suborder <- NA
resolved_names$superclass <- NA
resolved_names$phylum <- NA
resolved_names$subphylum <- NA
resolved_names$kingdom <- NA
resolved_names$domain <- NA

# Loop through each row and add taxonomic information
for (i in 1:nrow(resolved_names)) {
  ott_id <- resolved_names$ott_id[i]
  cat("Processing OTT ID:", ott_id, "\n")  # Debug print statement
  ranks <- get_taxonomic_ranks(ott_id)
  print(ranks)  # Debug print statement
  
  resolved_names$class[i] <- ranks["class"]
  resolved_names$subclass[i] <- ranks["subclass"]
  resolved_names$order[i] <- ranks["order"]
  resolved_names$suborder[i] <- ranks["suborder"]
  resolved_names$superclass[i] <- ranks["superclass"]
  resolved_names$phylum[i] <- ranks["phylum"]
  resolved_names$subphylum[i] <- ranks["subphylum"]
  resolved_names$kingdom[i] <- ranks["kingdom"]
  resolved_names$domain[i] <- ranks["domain"]
}

# Remove rows where 'class' is not 'Mammalia'
resolved_names <- resolved_names[resolved_names$class == "Mammalia", ]


# Merge resolved_names with dat based on resolved_names$unique_name and dat_species
merged_data <- merge( dat, resolved_names, by.x = "sciname", by.y = "unique_name", all.x = FALSE)

#Remove duplicate
merged_data <- merged_data[!duplicated(merged_data[, c("sciname", "Longitude", "Latitude","Ano")]), ]

#Export
write_xlsx(merged_data, "Salve_Final_rotl.xlsx")
write_xlsx(resolved_names,"PresenceSALVE_Species.xlsx")  


########################################################################################
###########################     3) Thinning     ########################################
########################################################################################


#To improve the model’s performance, spatial thinning was applied to the occurrence records using the “spThin” R package (https://cran.r-project.org/web/packages/spThin/index.html). 
#This tool randomly selects one record within a 10km radius (which corresponds to the resolution of the land use projections and environmental variables). 


#Import of data 
data <- read_xlsx("Salve_Final_rotl.xlsx")
data <- data[!duplicated(data[, c("sciname", "Longitude", "Latitude")]), ]
species_count <- data %>% group_by(Ott_id, sciname) %>% tally(name = "Occurrences")
filtered_data <- data %>% semi_join(species_count %>% filter(Occurrences > 10), by = c("Ott_id", "sciname"))

# Thinning process (Among all the random record selection, the subset maximizing the NNI is selected)
calculate_nni_spin <- function(data) {
  
  data <- data %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
  
  coords <- st_coordinates(data)
  bbox <- st_bbox(data)
  window <- owin(xrange = c(bbox["xmin"], bbox["xmax"]), yrange = c(bbox["ymin"], bbox["ymax"]))
  points_ppp <- ppp(coords[, 1], coords[, 2], window = window)
  
  nn_distances <- nndist(points_ppp)
  D_observed <- mean(nn_distances)
  lambda <- npoints(points_ppp) / area.owin(window)
  D_expected <- 1 / (2 * sqrt(lambda))
  NNI <- D_observed / D_expected
  
  return(NNI)
}  #Function to calculate Nearest Neighbor Index (NNI)

thinned_data <- filtered_data %>%
  group_by(Ott_id) %>%
  do({
    species_data <- .
    
    thinned_species <- thin(
      loc.data = species_data,
      lat.col = "Latitude",
      long.col = "Longitude",
      spec.col = "Ott_id",
      thin.par = 10,
      reps = 100,
      locs.thinned.list.return = TRUE,
      write.files = FALSE,
      write.log.file = FALSE
    )
    
    max_NNI_index <- which.max(lapply(thinned_species, calculate_nni_spin))
    Thin_data <- thinned_species[[max_NNI_index]]
    
    merged_df <- merge(species_data, Thin_data, by = c("Latitude", "Longitude"), all = FALSE)
    merged_df %>%
      group_by(Latitude, Longitude) %>%
      ungroup()
  }) %>%
  bind_rows()

spin_df <- thinned_data %>%
  group_by(Latitude, Longitude) %>%
  ungroup()

#ONLY >10

# Count occurrences per species↨
species_count_spin <- spin_df %>% group_by(Ott_id, sciname) %>% tally(name = "Occurrences")

# Filter species with more than 10 occurrences
filtered_spin <- spin_df %>% 
  semi_join(species_count_spin %>% filter(Occurrences > 10), by = c("Ott_id", "sciname"))


#Export as SHP

sf_object <- filtered_spin %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
st_write(sf_object, "SALVE_Final_thin.shp")



########################################################################################
##################     4) Extraction of environemental variables     ###################
########################################################################################


# Environmental Variable Extraction in ArcGIS: Using the generated SHP file in ArcGIS to extract environmental variable values for each species occurrence and absence record.
