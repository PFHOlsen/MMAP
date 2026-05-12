################################################################################
# Preprocessing of forest fire data ############################################
################################################################################

# Load tidyverse package
library(tidyverse)

# Forest fire downloaded from https://catalog.data.gov/dataset/odf-fire-occurrence-data-2000-2022
# on 25 May 2025.

# Load forest fire data - adjust path
fire <- read.csv2(file = '.../ODF_Fire_Occurrence_Data_2000-2022.csv', sep = ",")

# Set initial time point
time0 <- as.POSIXct("01/01/2018 12:00:00 AM", format = "%m/%d/%Y %I:%M:%S %p", tz = "UTC")

# Select 2018-2022 and handle multiple arrivals by assigning area code with 
# largest EstTotalAcres corresponding to the largest area affected
fire_clean <- fire %>%
  mutate(
    time = as.POSIXct(Ign_DateTime, format = "%m/%d/%Y %I:%M:%S %p", tz = "UTC"),
    year_var = year(time),
    month_var = month(time)
  ) %>%
  filter(year_var %in% 2018:2022) %>%
  mutate(
    arrival_time = as.numeric(difftime(time, time0, units = "days"))
  ) %>% 
  group_by(arrival_time) %>%
  slice_max(order_by = EstTotalAcres, n = 1, with_ties = FALSE) %>%
  ungroup()

# Helper function to construct area specific data sets
make_data <- function(Areas){
  fire_clean %>% dplyr::filter(Area %in% Areas) %>%
    dplyr::mutate(time = as.POSIXct(Ign_DateTime, format = "%m/%d/%Y %I:%M:%S %p", tz = "UTC")) %>% 
    dplyr::mutate(year_var = year(time), 
                  month_var = month(time)) %>%
    dplyr::filter(year_var %in% c(2018:2022)) %>%
    dplyr::arrange(time) %>%
    dplyr::mutate(arrival_time = as.numeric(difftime(time, time0, units = "days")),
                  interarrival_days = as.numeric(difftime(time, lag(time), units = "days"))) %>% 
    dplyr::mutate(interarrival_days = ifelse(is.na(interarrival_days), arrival_time,interarrival_days))
  }
  
# Construct arrival sequence data for all areas and EOA, NOA and SOA respectively.
fire_all <- make_data(Areas = c("EOA","NOA","SOA"))
fire_EOA <- make_data(Areas = c("EOA"))
fire_NOA <- make_data(Areas = c("NOA"))
fire_SOA <- make_data(Areas = c("SOA"))

# Select interarrival_days
x <- fire_all %>% pull(interarrival_days)
x <- x[!is.na(x)]
x_EOA <- fire_EOA %>% pull(interarrival_days)
x_EOA <- x_EOA[!is.na(x_EOA)]
x_NOA <- fire_NOA %>% pull(interarrival_days)
x_NOA <- x_NOA[!is.na(x_NOA)]
x_SOA <- fire_SOA %>% pull(interarrival_days)
x_SOA <- x_SOA[!is.na(x_SOA)]

# Length of each data set
N <- length(x)
N_EOA <- length(x_EOA)
N_NOA <- length(x_NOA)
N_SOA <- length(x_SOA)

# Uniform sequences for QQ-plots
U_sq_full <- (1:N)/N
U_sq_EOA <- (1:N_EOA)/N_EOA
U_sq_NOA <- (1:N_NOA)/N_NOA
U_sq_SOA <- (1:N_SOA)/N_SOA
