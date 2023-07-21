# Load the Rquake package

pacman::p_load(pacman, Rquake) 
pacman::p_load(pacman, rio) 
library(Rquake)
library(datasets)

# Load the earthquake dataset
rio_csv <- import("/Users/aatal/Downloads/Copy of maths project dataset (1).xlsx")
head(rio_csv)

mag1<-rio_csv$magnitude
# Set the minimum and maximum magnitudes to consider
min_mag <- min (mag, na.rm=T)
max_mag <- max(mag, na.rm = F)

mfd <- function(mag, min_mag = NULL, max_mag = NULL, bin_width = 0.1) {
  
  # Check input arguments
  stopifnot(is.numeric(mag))
  if (!is.null(min_mag)) stopifnot(is.numeric(min_mag))
  if (!is.null(max_mag)) stopifnot(is.numeric(max_mag))
  stopifnot(is.numeric(bin_width) && bin_width > 0)
  
  # Set default min_mag and max_mag values if not provided
  if (is.null(min_mag)) min_mag <- floor(min(mag))
  if (is.null(max_mag)) max_mag <- ceiling(max(mag))
  
  # Compute the number of bins
  n_bins <- ceiling((max_mag - min_mag) / bin_width)
  
  # Compute the bin edges
  bin_edges <- seq(min_mag, min_mag + n_bins * bin_width, by = bin_width)
  
  # Compute the bin centers
  bin_centers <- bin_edges[1:n_bins] + bin_width / 2
  
  # Bin the magnitudes
  bin_counts <- table(cut(mag, breaks = bin_edges))
  
  # Convert bin counts to a numeric vector
  bin_counts <- as.numeric(bin_counts)
  
  # Remove any bins with zero counts
  non_zero_bins <- bin_counts > 0
  bin_centers <- bin_centers[non_zero_bins]
  bin_counts <- bin_counts[non_zero_bins]
  
  # Return a list with the bin centers and counts
  return(list(x = bin_centers, y = bin_counts))
}
mbass <- function(mag, freq, plot = TRUE) {
  
  # Check input arguments
  stopifnot(is.numeric(mag) && length(mag) > 1)
  stopifnot(is.numeric(freq) && length(freq) == length(mag))
  
  # Compute the cumulative sum of frequencies
  cum_freq <- cumsum(freq)
  
  # Compute the slope of each segment
  slopes <- diff(log(cum_freq))
  
  # Compute the midpoints of each segment
  midpoints <- mag[-1] - diff(mag) / 2
  
  # Compute the median slope and its 95% confidence interval
  median_slope <- median(slopes, na.rm = TRUE)
  slope_std <- mad(slopes, na.rm = TRUE) * 1.4826
  slope_ci <- c(median_slope - 1.96 * slope_std, median_slope + 1.96 * slope_std)
  
  # Find the Mc estimate by finding the magnitude with the smallest difference between
  # the segment slope and the median slope
  diff_slope <- abs(slopes - median_slope)
  idx <- which.min(diff_slope)
  mc <- midpoints[idx]
  
  # Plot the segment slopes and the Mc estimate, if requested
  if (plot) {
    plot(midpoints, slopes, type = "l", xlab = "Magnitude", ylab = "Segment slope")
    abline(h = median_slope, col = "red")
    abline(v = mc, col = "blue")
    legend("topright", legend = c("Median slope", "Mc"), col = c("red", "blue"), lty = 1)
  }
  
  # Return the Mc estimate
  return(mc)
}


# Compute the magnitude-frequency distribution (MFD) of the earthquake catalog
mfd <- mfd(mag1, min_mag = min_mag, max_mag = max_mag, bin_width = 0.1)

# Compute the MBASS estimate of Mc
mbass_mc <- mbass(mfd$x, mfd$y, plot = T)

# Print the estimated Mc value
mbass_mc


# CLEAN UP #################################################

# Clear environment
rm(list = ls()) 

# Clear packages

p_unload(all)  # Remove all add-ons

# Clear console


cat("\014")  # ctrl+L

# Clear mind :)
