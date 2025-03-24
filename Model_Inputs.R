#-------------------------------------------------------------------#
# INITIALIZATION AND DATA INPUTS -----------------------------------#
#-------------------------------------------------------------------#


# This section defines the key demographic inputs required for the model, 
# including age groups, vaccination doses, and baseline risks for severe 
# health outcomes, intialisation for data, and the compartmental model setup.

#Packages######

pacman::p_load(dplyr, surveillance,  tidyverse,data.table)

#-------------------------------------------------------------------#
# DEMOGRAPHIC SETUP ---------------------------------------------
#-------------------------------------------------------------------#
# Define age groups and dose categories
ages <- paste0('a', 1:5)  # 5 broad age bands
ages_10bds <- c('0_5', '6_11', '12_17', '18_29', '30_39', 
                '40_49', '50_59', '60_69', '70_79', '80+')  # 10 detailed age bands
doses <- c(0, 1, 2, 3)  # Vaccination dose categories (0 = unvaccinated, 3 = booster)

#-------------------------------------------------------------------#
# POPULATION DATA LOADING  -----------------------------------------
#-------------------------------------------------------------------#

# Load census data for the population of Hong Kong (end of 2021)
census <- read.csv('temp_data/Clean_age_sex_10AgeBandsDemoHk_end21.csv')

# Process the census data
N <- census %>%
  # Add factor levels for age groups
  mutate(age_grp = factor(age_grp, 
                          levels = c('0-5', '6-11', '12-17', 
                                     '18-29', '30-39', '40-49', 
                                     '50-59', '60-69', '70-79', '80+'))) %>%
  arrange(age_grp) %>% 
  filter(Sex == 'Both sexes' & !is.na(age_grp)) %>%  # Keep only relevant rows
  # Convert population counts to thousands
  mutate(year.end.n21_1 = year.end.n21_1 * 1000) %>% 
  select(year.end.n21_1) %>%  # Extract the population column
  unlist()

# Assign names to the 10 detailed age bands
names(N) <- paste0('n.', ages_10bds)

# Aggregate population into 5 broad age groups
N.5bands <- c(sum(N[1:3]), sum(N[4:5]), sum(N[6:7]), sum(N[8:9]), N[10])
names(N.5bands) <- ages  # Assign names to the 5 broad age groups

# Check total population (sanity check, no values are dropped)
sum(N.5bands)

#-------------------------------------------------------------------#
# COMPARTMENTAL MODEL SETUP -------------------------------------
#-------------------------------------------------------------------#

# Define the number of age groups, doses, and compartments
n.ages <- length(ages)  # Number of broad age groups
n.ages10bds <- length(ages_10bds)  # Number of detailed age bands
n.doses <- length(doses)  # Number of dose categories
comparts <- n.ages * n.doses  # Total number of compartments

# Define IDs for compartments across broad and detailed age groups
id_v <- paste0(rep(ages, each = n.doses), paste0('d', rep(doses, times = n.ages)))
id_v10bds <- paste0(rep(ages_10bds, each = n.doses), paste0('d', rep(doses, times = n.ages10bds)))

#-------------------------------------------------------------------#
# FUNCTION FOR AGE GROUP CONVERSION -----------------------------
#-------------------------------------------------------------------#

# Function to map detailed age groups to broad age groups (10 bands to 5 bands)
muta_FTN_a10toa5 <- function(x) {
  x <- gsub("0_5d", "a1d", x)
  x <- gsub("6_11", "a1", x)
  x <- gsub("12_17", "a1", x)
  x <- gsub("18_29", "a2", x)
  x <- gsub("30_39", "a2", x)
  x <- gsub("40_49", "a3", x)
  x <- gsub("50_59", "a3", x)
  x <- gsub("60_69", "a4", x)
  x <- gsub("70_79", "a4", x)
  x <- gsub("80\\+", "a5", x)
  return(x)
}

#-------------------------------------------------------------------#
#  OBSERVED DATA AND STUDY DATES ---------------------------------
#-------------------------------------------------------------------#

studydates<-structure(c(19014, 19112), class = "Date")
# Define end time for the simulation
set.tau <- 8  # Set the time step (tau)
set.end <- 99 # Total number of days in the study period
end <- set.end  # End time for the simulation




#-------------------------------------------------------------------#
#  VACCINATION PROPORTIONS AND POPULATION SETUP ------------------
#-------------------------------------------------------------------#

N.vacp <- readRDS("temp_data/N_vacp.rds")

# Define indices for vaccination groups (unvaccinated, 1st dose, 2nd dose, 3rd dose)
uvid_v <- which(grepl('d0', id_v))
v1id_v <- which(grepl('d1', id_v))
v2id_v <- which(grepl('d2', id_v))
v3id_v <- which(grepl('d3', id_v))

uvid_v10bds <- which(grepl('d0', id_v10bds))
v1id_v10bds <- which(grepl('d1', id_v10bds))
v2id_v10bds <- which(grepl('d2', id_v10bds))
v3id_v10bds <- which(grepl('d3', id_v10bds))

# Map 10 age bands to 5 broader age groups
a5bds_10bds <- c(1, 1, 1, 2, 2, 3, 3, 4, 4, 5)

N.pop <- rep(N, each = n.doses) * N.vacp 
names(N.pop) <- id_v10bds

# Verify the total population remains consistent
sum(N.pop) == sum(N)

# Define severe health outcome indicators
Sev_indicators <- c('Hosp', 'SCC', 'Death')



#-------------------------------------------------------------------#
# 9. KEPPA VALUES --------------------------------------------------#
#-------------------------------------------------------------------#
# Define keppa values for different vaccination groups and age groups
# Keppa (asymptomatic proportion) values derived from systematic reviews and studies
IE<-1
keppa0 <- c(0.187, 0.1, 0.351)[IE]  # Unvaccinated
keppa1 <- c(0.187, 0.1, 0.351)[IE]  # 1st dose
keppa2 <- c(0.193, 0.123, 0.305)[IE]  # 2nd dose
keppa3 <- c(0.672, 0.573, 0.788)[IE]  # Booster dose

# Keppa by age group
keppa_a1 <- c(0.4375, 0.38, 0.4905)[IE]  # Age group 1 (0-5)
keppa_a2 <- c(0.3259, 0.2326, 0.4192)[IE]  # Age group 2 (6-11)
keppa_a3 <- c(0.30, 0.20, 0.35)[IE]  # Age group 3 (12-17)
keppa_a4 <- c(0.25, 0.18, 0.32)[IE]  # Age group 4 (18-29, etc.)
keppa_a5 <- c(0.20, 0.15, 0.25)[IE]  # Age group 5 (80+)

# Combine age group and dose-specific keppa values
keppa_a1a5 <- c(keppa_a1 * c(keppa0, keppa1, keppa2, keppa3), 
                keppa_a2 * c(keppa0, keppa1, keppa2, keppa3),
                keppa_a3 * c(keppa0, keppa1, keppa2, keppa3), 
                keppa_a4 * c(keppa0, keppa1, keppa2, keppa3),
                keppa_a5 * c(keppa0, keppa1, keppa2, keppa3))
names(keppa_a1a5) <- id_v

keppa <- c(keppa_a1 * c(keppa0, keppa1, keppa2, keppa3), 
           keppa_a1 * c(keppa0, keppa1, keppa2, keppa3),
           keppa_a1 * c(keppa0, keppa1, keppa2, keppa3), 
           keppa_a2 * c(keppa0, keppa1, keppa2, keppa3),
           keppa_a2 * c(keppa0, keppa1, keppa2, keppa3),
           keppa_a3 * c(keppa0, keppa1, keppa2, keppa3),
           keppa_a3 * c(keppa0, keppa1, keppa2, keppa3),
           keppa_a4 * c(keppa0, keppa1, keppa2, keppa3),
           keppa_a4 * c(keppa0, keppa1, keppa2, keppa3),
           keppa_a5 * c(keppa0, keppa1, keppa2, keppa3))
names(keppa) <- id_v10bds

# Extract specific keppa values for unvaccinated, 2nd dose, and booster groups
keppa_uv <- keppa[which(grepl('d0', names(keppa)))]
keppa_v2 <- keppa[which(grepl('d2', names(keppa)))]
keppa_v3 <- keppa[which(grepl('d3', names(keppa)))]

# Clean up intermediate variables
rm(keppa_a1, keppa_a2, keppa_a3, keppa_a4, keppa_a5, keppa0, keppa2, keppa3)

#-------------------------------------------------------------------#
# ALPHA VALUES -------------------------------------------------
#-------------------------------------------------------------------#
# Alpha: Age- and dose-specific susceptibility to infection
# Data from Tsang et al. and other sources
alphaC0 <- 1
alphaC2 <- 1
alphaC3 <- 1

alphaC2_12_17 <- c(1 - 0.12, 1 - 0.29, 1)[IE]  # Vaccine effectiveness for 12-17 years

alphaC3_12_17 <- c(1 - 0.25, 1 - 0.50, 1 - 0.10)[IE]  # Adjust for reduced effectiveness for Sinovac

alphaA0 <- 1
alphaA2 <- 1
alphaA3 <- c(1 - 0.35, 1 - 0.584, 1 - 0.10)[IE]

alphaS0 <- 1
alphaS2 <- 1
alphaS3 <- 1  # No protection for 60+ (Tsang et al.)

# Combine alpha values
alpha_uv <- c(alphaC0, alphaC0, alphaC0, alphaA0, alphaA0, alphaA0, alphaA0, alphaS0, alphaS0, alphaS0)
alpha_v2 <- c(alphaC2, alphaC2, alphaC2_12_17, alphaA2, alphaA2, alphaA2, alphaA2, alphaS2, alphaS2, alphaS2)
alpha_v3 <- c(alphaC3, alphaC3, alphaC3_12_17, alphaA3, alphaA3, alphaA3, alphaA3, alphaS3, alphaS3, alphaS3)

# Clean up intermediate variables
rm(alphaC0, alphaA0, alphaS0, alphaC2, alphaA2, alphaS2, alphaC3, alphaA3, alphaS3)

#-------------------------------------------------------------------#
#PHI VALUES ---------------------------------------------------
#-------------------------------------------------------------------#
# Define phi values for asymptomatic transmissibility
phi_asymp_v <- c(0.35, 0.35, 0.35, 0.3035, 0.2442, 0.2338, 0.2904, 0.5433, 0.60, 0.60)

# Add a 20% error margin for lower and upper bounds
phi_asymp_lower <- phi_asymp_v - (0.2 * phi_asymp_v)
phi_asymp_upper <- phi_asymp_v + (0.2 * phi_asymp_v)

# Use lower bound for modeling
phi_asymp_v <- phi_asymp_lower

# Define dose-specific reductions in transmissibility
phi_v1_reduc0 <- c(1, 0.9, 1)[IE]
phi_v2_reduc0 <- c(0.8, 0.7, 0.9)[IE]
phi_v3_reduc0 <- c(0.7, 0.6, 0.9)[IE]





#-------------------------------------------------------------------#
# INITIAL INFECTIONS AND PARAMETERS ----------------------------
#-------------------------------------------------------------------#

# Latent period (time from exposure to infectiousness)
latent <- c(3.1, 2.8, 3.5)[IE]  # Bounds for latent period (source: studies)

#-------------------------------------------------------------------#
# DISTRIBUTIONS AND DELAY TIMES --------------------------------
#-------------------------------------------------------------------#
# PMFs for exposure to severe event 
# For peak periods
Peak_ProwtNDays_mat_list<-readRDS("temp_data/Peak_pmf_mat.rds")
# For non-peak periods
NonPeak_ProwtNDays_mat_list<-readRDS("temp_data/NonPeak_pmf_mat.rds")
# aggregated hospitalisation rate for SEIR model specifically
med_hosp_days_peak<-readRDS("temp_data/med_hosp_days_peak.rds")





#-------------------------------------------------------------------#
# VACCINATION RATE CHANGE DATES --------------------------------
#-------------------------------------------------------------------#

#splits between before and after change in vac rate, age specific, vac per day. 
pi_v1_0<-readRDS('temp_data/pi_v1_0.rds')
pi_v2_0<-readRDS('temp_data/pi_v2_0.rds')
pi_v3_0<-readRDS('temp_data/pi_v3_0.rds')

#vac rates
# 
#d1 change dates
incrVR_d1_u5<-as.Date("2022-02-28")
incrVR_d1_6_11<-as.Date("2022-02-28")

incrVR_d1_12_17<-as.Date("2022-03-28")
incrVR_d1_18_59<-as.Date("2022-03-15")
incrVR_d1_60pls<-as.Date("2022-03-15")

#d2 change dates
incrVR_d2_u5<-as.Date("2022-03-28")
incrVR_d2_6_11<-as.Date("2022-03-28")

incrVR_d2_12_17<-as.Date("2022-03-28")
incrVR_d2_18_59<-as.Date("2022-02-18")
incrVR_d2_60pls<-as.Date("2022-02-18")

#d3 change dates
incrVR_d3_u17<-as.Date("2022-03-25")
incrVR_d3_18_59<-as.Date("2022-02-17")
incrVR_d3_60_69<-as.Date("2022-02-17")
incrVR_d3_70pls<-as.Date("2022-03-14")


#calculate Day for each age group
Day_incrVR_d1_u5<-as.numeric(incrVR_d1_u5 - studydates[1])
Day_incrVR_d1_6_11<-as.numeric(incrVR_d1_6_11 - studydates[1])
Day_incrVR_d1_12_17<-as.numeric(incrVR_d1_12_17 - studydates[1])
Day_incrVR_d1_18_59<-as.numeric(incrVR_d1_18_59 - studydates[1])
Day_incrVR_d1_60pls<-as.numeric(incrVR_d1_60pls - studydates[1])

Day_incrVR_d2_u5<-as.numeric(incrVR_d2_u5 - studydates[1])
Day_incrVR_d2_6_11<-as.numeric(incrVR_d2_6_11 - studydates[1])
Day_incrVR_d2_12_17<-as.numeric(incrVR_d2_12_17 - studydates[1])
Day_incrVR_d2_18_59<-as.numeric(incrVR_d2_18_59 - studydates[1])
Day_incrVR_d2_60pls<-as.numeric(incrVR_d2_60pls - studydates[1])

Day_incrVR_d3_u17<-as.numeric(incrVR_d3_u17 - studydates[1])
Day_incrVR_d3_18_59<-as.numeric(incrVR_d3_18_59 - studydates[1])
Day_incrVR_d3_60_69<-as.numeric(incrVR_d3_60_69 - studydates[1])
Day_incrVR_d3_70pls<-as.numeric(incrVR_d3_70pls - studydates[1])
# 

#-------------------------------------------------------------------#
# INITIALISING PARAMS-------------------------------
#-------------------------------------------------------------------#

# Load sample data for severe health probabilities
# samplevers_27 <- read_csv("temp_data/samplevers.27.2.csv")
# # Median of calibrated parameters. 
# sample_med_init <- lapply(samplevers_27, median) %>% unlist()
# saveRDS(sample_med_init, "temp_data/sample_med_init.rds")


# INITIAL IOTA VALUES ----------------------------------------------#
iota_values0 <- c(0.5, 0.5, 0.5)

# lrr place holders
v1lrr<-c(1L, 1L, 4L, 7L, 10L)
v2lrr<-c(2L, 2L, 5L, 8L, 11L)
v3lrr<-c(3L, 3L, 6L, 9L, 12L)




# SOCIAL MIXING MATRICES -------------------------------------------
# Load the social mixing survey data from the HK 2020 survey
contact_school<-readRDS("temp_data/contact_school.rds")
contact_work<-readRDS("temp_data/contact_work.rds")
contact_home<-readRDS("temp_data/contact_home.rds")
contact_transport<-readRDS("temp_data/contact_transport.rds")
contact_other<-readRDS("temp_data/contact_other.rds")
prim_school <- 1:2
sec_school  <- 3
teacher     <- 4:10

# Load the time-varying iota data
time_varying_iota <- readRDS("temp_data/time_varying_iota.rds")
# Convert to a matrix and remove the date column
time_varying_iota_0 <- time_varying_iota %>%
  select(-date) %>%
  as.matrix()

# Remove the original time-varying iota data to free memory
rm(time_varying_iota)







#-------------------------------------------------------------------#
# PEAK SEVERE HEALTH INCREASE --------------------------------------
#-------------------------------------------------------------------#
# Define peak period start and end dates
peak_tstart_dt <- as.Date("2022-02-16")  # Peak start date
peak_tend_dt <- as.Date("2022-03-30")    # Peak end date

# Convert peak dates into numeric days relative to study start date
peak_tstart <- (peak_tstart_dt - studydates[1]) %>% as.numeric()  # Peak start day
peak_all <- (as.Date("2022-03-10") - studydates[1]) %>% as.numeric()  # Peak midpoint day
peak_tend <- (peak_tend_dt - studydates[1]) %>% as.numeric()  # Peak end day

# Calculate the duration of the peak rise and fall
lngthup <- (peak_all - peak_tstart) %>% round(0)  # Length of peak rise
lngthdn <- (peak_tend - peak_all) %>% floor()  # Length of peak decline

