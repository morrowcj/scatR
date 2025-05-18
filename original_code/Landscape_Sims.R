# Load necessary libraries
library(terra)
library(NLMR)
library(tidyverse)

#------------------------------
# 1. Simulate Landscape
#------------------------------
#' Title
#'
#' @param ncol 
#' @param nrow 
#' @param fractal_dimension 
#' @param n_classes 
#'
#' @returns
#' @export
#'
#' @examples
simulate_landscape <- function(
    ncol = 1000, nrow = 1000, fractal_dimension = 1, n_classes = 10
  ) {
  raw_landscape <- nlm_mpd(ncol = ncol, nrow = nrow, roughness = 1)
  land_rast <- rast(raw_landscape)
  breaks <- quantile(
    values(land_rast)[,1], probs = seq(0, 1, length.out = n_classes + 1), 
    na.rm = TRUE
  )
  rcl <- matrix(ncol = 3, nrow = n_classes)
  for (i in 1:n_classes) {
    rcl[i, ] <- c(breaks[i], breaks[i + 1], i)
  }
  land_classified <- classify(
    land_rast, rcl = rcl, include.lowest = TRUE, 
    filename = tempfile(fileext = ".tif"), overwrite = TRUE
  )
  return(land_classified)
}

#------------------------------
# 2. Create Circular HR Mask
#------------------------------
create_hr_mask <- function(center, radius, template_rast) {
  grid <- as.data.frame(expand.grid(x = 1:ncol(template_rast), y = 1:nrow(template_rast)))
  distances <- sqrt((grid$x - center[1])^2 + (grid$y - center[2])^2)
  mask_vals <- ifelse(distances <= radius, 1, NA)
  mask_rast <- rast(ncol = ncol(template_rast), nrow = nrow(template_rast))
  ext(mask_rast) <- ext(template_rast)
  values(mask_rast) <- mask_vals
  crs(mask_rast) <- crs(template_rast)
  return(mask_rast)
}

#------------------------------
# 3. Check Overlap (Area-Based)
#------------------------------
check_overlap <- function(
    new_center, new_radius, existing_centers, existing_radii, max_overlap
  ) {
  if (nrow(existing_centers) == 0) return(TRUE)
  for (i in 1:nrow(existing_centers)) {
    dist <- sqrt(sum((new_center - existing_centers[i, ])^2))
    r1 <- new_radius
    r2 <- existing_radii[i]
    if (dist < (r1 + r2)) {
      if (dist <= abs(r1 - r2)) {
        overlap_area <- pi * min(r1, r2)^2
      } else {
        part1 <- r1^2 * acos((dist^2 + r1^2 - r2^2) / (2 * dist * r1))
        part2 <- r2^2 * acos((dist^2 + r2^2 - r1^2) / (2 * dist * r2))
        part3 <- 0.5 * sqrt(
          (-dist + r1 + r2)*(dist + r1 - r2)*(dist - r1 + r2)*(dist + r1 + r2)
        )
        overlap_area <- part1 + part2 - part3
      }
      new_area <- pi * r1^2
      if ((overlap_area / new_area) > max_overlap) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

#------------------------------
# 4. Generate Home Ranges with Consumer Strategies
#------------------------------
generate_home_ranges <- function(landscape, 
                                 num_hrs, 
                                 mean_radius, 
                                 sd_radius, 
                                 max_overlap, 
                                 strategy, 
                                 selection_strength_range
                                 ) {
  hr_centers <- matrix(ncol = 2, nrow = 0)
  hr_radii <- numeric(0)
  extracted_values <- list()
  available_values <- list()
  selection_strengths <- numeric(0)
  ncol_land <- ncol(landscape)
  nrow_land <- nrow(landscape)
  
  while (nrow(hr_centers) < num_hrs) {
    center <- runif(2, min = 1, max = c(ncol_land, nrow_land))
    radius <- abs(rnorm(1, mean = mean_radius, sd = sd_radius))
    
    if (all(center - radius >= 1) && all(center + radius <= c(ncol_land, nrow_land))) {
      if (check_overlap(center, radius, hr_centers, hr_radii, max_overlap)) {
        individual_strategy <- ifelse(strategy == "mixed", sample(c("uniform", "selective"), 1), strategy)
        individual_strength <- ifelse(individual_strategy == "selective", runif(1, selection_strength_range[1], selection_strength_range[2]), NA)
        hr_centers <- rbind(hr_centers, center)
        hr_radii <- c(hr_radii, radius)
        hr_mask <- create_hr_mask(center, radius, landscape)
        masked <- mask(landscape, hr_mask)
        values_in_hr <- values(masked)[!is.na(values(masked))]
        
        
        if (individual_strategy == "selective") {
          probs <- values_in_hr^individual_strength
          probs[is.na(probs) | probs <= 0] <- 1e-6
          sampled <- sample(values_in_hr, length(values_in_hr), replace = TRUE, prob = probs)
          extracted_values[[nrow(hr_centers)]] <- sampled
        } else {
          extracted_values[[nrow(hr_centers)]] <- values_in_hr
        }
        available_values[[nrow(hr_centers)]] <- values_in_hr
        selection_strengths[nrow(hr_centers)] <- ifelse(is.na(individual_strength), 0, individual_strength)
      }
    }
  }
  
  return(list(centers = hr_centers, radii = hr_radii, values = extracted_values, available = available_values, selection_strengths = selection_strengths))
}

#------------------------------
# 5. Calculate Overlap Matrix
#------------------------------
calculate_overlap_matrix <- function(centers, radii) {
  n <- nrow(centers)
  overlap_matrix <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:(n - 1)) {
    area_i <- pi * radii[i]^2
    for (j in (i + 1):n) {
      dist <- sqrt(sum((centers[i, ] - centers[j, ])^2))
      if (dist >= (radii[i] + radii[j])) {
        overlap_area <- 0
      } else if (dist <= abs(radii[i] - radii[j])) {
        overlap_area <- pi * min(radii[i], radii[j])^2
      } else {
        r1 <- radii[i]; r2 <- radii[j]
        part1 <- r1^2 * acos((dist^2 + r1^2 - r2^2) / (2 * dist * r1))
        part2 <- r2^2 * acos((dist^2 + r2^2 - r1^2) / (2 * dist * r2))
        part3 <- 0.5 * sqrt((-dist + r1 + r2)*(dist + r1 - r2)*(dist - r1 + r2)*(dist + r1 + r2))
        overlap_area <- part1 + part2 - part3
      }
      area_j <- pi * radii[j]^2
      overlap_matrix[i, j] <- (overlap_area / area_i) * 100
      overlap_matrix[j, i] <- (overlap_area / area_j) * 100
    }
  }
  return(overlap_matrix)
}

#------------------------------
# 6. Plot Landscape and HRs
#------------------------------
plot_landscape_and_hrs <- function(landscape, centers, radii) {
  plot(landscape, main = "Landscape with Home Ranges")
  for (i in 1:nrow(centers)) {
    theta <- seq(0, 2 * pi, length.out = 100)
    x <- centers[i, 1] + radii[i] * cos(theta)
    y <- centers[i, 2] + radii[i] * sin(theta)
    lines(x, y, col = "red", lwd = 2)
  }
}

#------------------------------
# Example Usage
#------------------------------
landscape <- simulate_landscape()

hr_data <- generate_home_ranges(
  landscape, num_hrs = 500, mean_radius = 5, sd_radius = 0.5, 
  max_overlap = 10, strategy = "mixed", 
  selection_strength_range = c(1, 5)
)

overlap_matrix <- calculate_overlap_matrix(hr_data$centers, hr_data$radii)

plot_landscape_and_hrs(landscape, hr_data$centers, hr_data$radii)

# Non-zero overlaps
overlap_tbl <- as.data.frame(as.table(overlap_matrix))

# Ensure column names are valid
names(overlap_tbl) <- c("hr1", "hr2", "OverlapPercentage")

# Filter non-zero and unique overlap pairs
overlap_df <- subset(overlap_tbl, OverlapPercentage > "0")

print(overlap_df)

summary(overlap_df$OverlapPercentage)

# Extracted values
hr_data$values

# Summary of use vs availability
summary_df <- tibble(
  hr_id = 1:length(hr_data$values),
  selection_strength = hr_data$selection_strengths,
  mean_used = sapply(hr_data$values, mean),
  mean_available = sapply(hr_data$available, mean),
  sd_used = sapply(hr_data$values, sd),
  sd_available = sapply(hr_data$available, sd)
)

print(summary_df)

summary(summary_df)
