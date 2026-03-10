#####
#Species accumulation curve
#Last updated: 1/7/2026
#Author: Kelly Wu 
#Description: 
#Species accumulation curve for germination experiment. As the curve plateaus, 
#we are closer and closer to reaching total diversity; in this case, total 
# viable seed bank
###

#vegan package: ordination methods, diversity analysis and other functions for community and vegetation ecologists

# Load necessary packages
library(tidyverse) # organization; includes tidyr and dplyr
library(vegan) #community analysis

# Load your CSV: already in wide format with presence absence
##data was created using pivot table in excel 
##basic function built into R
dataframe <- read.csv("2_data_raw/speciessample_sacdec19.csv")

# Replace blanks/NA with 0 in presence absence data
dataframe[is.na(dataframe)] <- 0    # in case blanks are read as ""
nolastrowdf <- dataframe[, -c(1, ncol(dataframe))] # removes the first and last column

# Convert to matrix
##returns all values of a raster layer as a matrix
##mat <- as.matrix(dataframe) <-this makes everything strings

# Run species accumulation curve (by plots)
spec_accum <- specaccum(nolastrowdf, method = "random")

# Plot the curve
plot(spec_accum, ci.type = "bar", ci.col = "lightblue",
     xlab="Number of plots", ylab="Cumulative species")
grid(nx = NULL, ny = NULL, col = "lightgray", lwd = 1)



#descriptive statistics
str(spec_accum)

curve_stats <- data.frame(
  Plots = spec_accum$sites,
  Richness = spec_accum$richness,
  SD = spec_accum$sd,
  CI_lower = spec_accum$richness - 1.96*spec_accum$sd,
  CI_upper = spec_accum$richness + 1.96*spec_accum$sd,
  Gain = c(NA, diff(spec_accum$richness))
)
print(curve_stats)

mean_richness <- mean(spec_accum$richness)
max_richness  <- max(spec_accum$richness)
min_richness  <- min(spec_accum$richness)

