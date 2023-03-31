## code to prepare `longitudinalData` dataset goes here

library(car)
library(ggplot2)

####################### Williams et al. ########################################
# We actually have to manually add the longitudinal data from the Williams paper
# based on the information contained in https://github.com/NickWilliamsSanger/mpn_phylogenies_and_evolution
# Specifically, we use the subdirectory found here https://github.com/NickWilliamsSanger/mpn_phylogenies_and_evolution/tree/master/supplementary_reports/supplementary_note_5_per_patient_reports
# Fortunately, only one of the clones has good longitudinal data, so we only have to manually add the data once.
# We include the .html file for PD9478 in the data-raw folder, just in case something changes.
# See Williams Ext. Data Fig. 7a for graphical view of longitudinal data.

# PD 9478
patient <- "PD9478"
driver <- "JAK2_DNMT3A"
cellTypeVec <- "PB Gran"

ageVec <- c(54.8, 55.4, 61.7, 63.2, 64, 66.4, 67, 67.6, 68, 70, 70.5)
vafVec <- c(0, 0, 0.14, 0.14, 0.24, 0.33, 0.41, 0.41, 0.44, 0.45, 0.46)

pd9478 <- data.frame(
  "Sample.ID" = patient, "Age" = ageVec, "VAF" = vafVec,
  "Gene" = driver, "cellType" = cellTypeVec,
  "Protein" = "Multiple",
  "cloneName" = paste0(patient, "_1_clone1")
)



######################### Fabre et al. #########################################
# Load Fabre longitudinal data (downloaded as Supp. Table 2 from Fabre et al.
# "The longitudinal dynamics and natural history of clonal haematopoiesis")
# We also include the downloaded csv in data-raw
fabre_data <- read.csv("~/rotation_fall2021/fabre_longitudinal_vaf_all_01032023.csv", skip = 3)
# Remove empty col
fabre_data$X <- NULL
# Remove empty rows
fabre_data <- fabre_data[!is.na(fabre_data$Age), ]

# Restrict cols of fabre data to the ones we need
all_long.df <- fabre_data[, c(1, 4, 10, 11, 12)]
all_long.df$cellType <- "PB (whole)"

# Split the csv by patient + gene + protein combo (unique ID for a clone)
all_long.df$uniqueClone <- paste0(all_long.df$Sample.ID, "+", all_long.df$Gene, "+", all_long.df$Protein)

# Make a list containing a df for each clone
cloneList <- split(all_long.df, all_long.df$uniqueClone)

# Add columns to track fit params
addCols.df <- data.frame(
  "t_m" = NA, "r" = NA, "K" = NA,
  "r_lb" = NA, "r_ub" = NA,
  "fitError" = NA, "stdError" = NA,
  "NumberSamples" = NA
)
results.vaf_data <- cbind(all_long.df, addCols.df)

# Define function to return fits and a string telling us if there's a fitting error
threeParamFn <- function(df) {
  startingParams <- coef(lm(logit(VAF / .5) ~ Age, data = df))
  tryCatch(
    {
      out <- nls(VAF ~ K / (1 + exp(-(phi2 + r * Age))),
        start = list(K = 0.2, phi2 = min(c(-0.000001, startingParams[1])), r = max(c(0.00001, startingParams[2]))), data = df, trace = F, algorithm = "port",
        lower = c(0, -500, 0.00001), upper = c(0.5, 0, 5)
      )
      list("No Error", out)
    },
    error = function(cond) {
      # Error in fitting occurred, we will fit anyway, but take note of it
      out <- nls(VAF ~ K / (1 + exp(-(phi2 + r * Age))),
        start = list(K = 0.2, phi2 = min(c(-0.000001, startingParams[1])), r = max(c(0.00001, startingParams[2]))), data = df, trace = F, algorithm = "port",
        lower = c(0, -500, 0.00001), upper = c(0.5, 0, 5), control = list(warnOnly = T)
      )
      list("Fit Error", out)
    }
  )
}



for (i in c(1:length(cloneList))) {
  df <- cloneList[[i]]

  # Remove data points that decrease (and anything after)
  df_check <- df
  for (j in c(2:nrow(df))) {
    maxPastVAF <- max(df$VAF[1:(j - 1)])
    presentVAF <- df$VAF[j]
    if (presentVAF < maxPastVAF * .8 & maxPastVAF > 0.05) {
      vafDrop <- T
      vafDropIndex <- j
      df_check <- df[c(1:(j - 1)), ]
      break
    }
  }

  # Skip if doesn't pass checks
  if (nrow(df_check) < 4 | max(df_check$VAF) - min(df_check$VAF) < .05 |
    nrow(df_check[df_check$VAF != 0 & df_check$VAF < 0.25, ]) < 2) {
    next
  } else {
    df <- df_check
  }

  # Run the logistic growth model fit
  fn_return_list <- threeParamFn(df)
  threeParamModel <- fn_return_list[[2]]
  fitError <- fn_return_list[[1]]

  # Set output equal to summary of model
  output <- summary(threeParamModel)

  # Assign fitted params
  r <- output$coefficients["r", 1]
  t_m <- -output$coefficients["phi2", 1] / r
  K <- output$coefficients["K", 1]

  # Get std. error from model
  stdError <- output$coefficients["r", 2]

  # Get 95% confidence intervals for model
  lb <- max(c(0, r - 1.96 * output$coefficients["r", 2]))
  ub <- min(c(4, r + 1.96 * output$coefficients["r", 2]))

  # Add to results df
  matching_index <- which(results.vaf_data$uniqueClone %in% df$uniqueClone[1])
  results.vaf_data$t_m[matching_index] <- t_m
  results.vaf_data$r[matching_index] <- r
  results.vaf_data$r_lb[matching_index] <- lb
  results.vaf_data$r_ub[matching_index] <- ub
  results.vaf_data$NumberSamples[matching_index] <- nrow(df)
  results.vaf_data$K[matching_index] <- K
  results.vaf_data$fitError[matching_index] <- fitError
  results.vaf_data$stdError[matching_index] <- stdError

  # Construct predict.df for plotting
  x <- c((min(df$Age) - 10):(max(df$Age) + 10)) # construct a range of x values bounded by the data
  y <- K / (1 + exp(-(x - t_m) * r)) # curve VAF (3 param model)
  predict.df <- data.frame("x" = x, "y" = y, "model" = rep(paste0("Three param. model, r: ", round(r, 3)), length(x)))

  # Set color and name for plot
  colorLong <- "#0072B2"
  plotName <- paste0(
    "~/rotation_fall2021/logistic_fits/individual_longitudinal/3paramOnly_fit", fitError, "_",
    df$uniqueClone[1], "_", Sys.Date(), ".pdf"
  )

  # Plot and manually review
  pdf(plotName, height = 5, width = 4)
  print(ggplot(df, aes(x = Age, y = VAF)) +
    theme_bw() +
    coord_cartesian(xlim = c(min(x), max(x)), ylim = c(-0.01, 0.52), expand = 0) +
    labs(x = "Person Age (yr)", y = "Variant Allele Frequency (VAF)") +
    ggtitle(gsub("_", " & ", gsub("\\+", " ", df$uniqueClone[1]))) +
    theme(
      axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.title = element_blank(),
      legend.text = element_text(size = 13), plot.title = element_text(size = 14)
    ) +
    geom_line(data = predict.df, aes(x = x, y = y), color = colorLong, size = 1, show.legend = T) +
    geom_point(color = "black", size = 3) +
    geom_vline(xintercept = t_m, linetype = "dotted", color = colorLong, linewidth = .6) +
    geom_hline(yintercept = K, linetype = "dotted", color = colorLong, linewidth = .6))
  dev.off()
}

# Get the results that passed checks
fitsOnly_results <- results.vaf_data[!is.na(results.vaf_data$r), ]

# Remove clones that have fit (convergence) errors or the one clone that is clearly being outcompeted by SRSF2 (manually reviewed)
fitsOnly_cut <- fitsOnly_results[fitsOnly_results$fitError == "No Error" & !fitsOnly_results$uniqueClone %in% c("PD41233+JAK2+p.V617F"), ]

all_fabre <- fitsOnly_cut[, c("Sample.ID", "Age", "VAF", "Gene", "Protein", "cellType")]

# If there's a matching single-cell clone, include the clone name as annotated in single cell data
all_fabre$cloneName <- "No matching single-cell clone"
all_fabre$cloneName[all_fabre$Sample.ID == "PD34493" & all_fabre$Gene == "SF3B1"] <- "PD34493_clone1"
all_fabre$cloneName[all_fabre$Sample.ID == "PD41276" & all_fabre$Gene == "SF3B1"] <- "PD41276_clone1"

# Only include the expansions that match a driver observed in any single-cell data.
# Only the SF3B1 clones with exact clone matches are included because these are
# the only single-cell clones with SF3B1 as a driver.
matching_phylo_driver_genes <- c("JAK2", "DNMT3A", "PPM1D", "CBL", "JAK2_DNMT3A")
matching_fabre_tmp <- all_fabre[all_fabre$Gene %in% matching_phylo_driver_genes, ]
matching_fabre <- rbind(matching_fabre_tmp, all_fabre[all_fabre$cloneName %in% c("PD34493_clone1", "PD41276_clone1"), ])

longitudinalData <- rbind(matching_fabre, pd9478)

usethis::use_data(longitudinalData, overwrite = TRUE)
