## code to prepare `longitudinalData` dataset goes here

# We actually have to manually add the longitudinal data from the Williams paper
# based on the information contained in https://github.com/NickWilliamsSanger/mpn_phylogenies_and_evolution
# Specifically, we use the subdirectory found here https://github.com/NickWilliamsSanger/mpn_phylogenies_and_evolution/tree/master/supplementary_reports/supplementary_note_5_per_patient_reports
# Fortunately, only one of the clones has good longitudinal data, so we only have to manually add the data once.
# We include the .html file for PD9478 in the data-raw folder, just in case something changes.

# PD 9478
patient <-  "PD9478"
driver <- "JAK2_DNMT3A"
cellTypeVec <- "PB Gran"

ageVec <- c(54.8, 55.4, 61.7, 63.2, 64, 66.4, 67, 67.6, 68, 70, 70.5)
vafVec <- c(0, 0, 0.14, 0.14, 0.24, 0.33, 0.41, 0.41, 0.44, 0.45, 0.46)

pd9478 <- data.frame("Sample.ID" = patient, "Age" = ageVec, "VAF" = vafVec,
                     "Gene" = driver, "cellType" = cellTypeVec,
                     "Protein" = "Multiple",
                     "cloneName" = paste0(patient, "_clone1"))



# Load Fabre longitudinal data (downloaded from: )
# We also include the downloaded csv in data-raw
fabre_data <- read.csv("~/rotation_fall2021/fabre_longitudinal_vaf_all_01032023.csv", skip = 3)
# Remove empty col
fabre_data$X <- NULL
# Remove empty rows
fabre_data <- fabre_data[! is.na(fabre_data$Age),]

# Restrict cols of fabre data to the ones we need
all_long_tmp.df <- fabre_data[,c(1, 4, 10, 11, 12)]
all_long_tmp.df$cellType <- "PB (whole)"

subset_fabre <- all_long_tmp.df[all_long_tmp.df$Sample.ID %in% c("PD41276", "PD34493") & all_long_tmp.df$Protein == "p.K666N",]
subset_fabre$cloneName <- paste0(subset_fabre$Sample.ID, "_clone1")

# Combine the two dataframes, one df for both sources (fabre + williams)
longitudinalData <- rbind(subset_fabre, pd9478)

usethis::use_data(longitudinalData, overwrite = TRUE)
