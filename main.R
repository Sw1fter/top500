
###
#   Preamble
###

library("dplyr")
library("rjags")
library("coda")
library("reshape2")
library("ggplot2")
library("ggfan")
library("lubridate")



###
#   Section 1 - modelling the development of new species per year for Sweden.
###

# Get the data on when people twitched what
top100 <- read.csv("top100.csv", header=TRUE, sep=",")

# Get all the data over first records for Sweden
swe <- read.csv("first records.csv", header=TRUE)

# Format date column
top100$date <- as.Date(top100$date)

# Add a column with the first year a birder has a twitch with date known.
top100 <- top100 %>% 
  group_by(member_id) %>%
  mutate(start_year = year(min(date)))

# Defining constants for which time window we care about
max_yr <- max(swe$year, na.rm = TRUE)
min_yr <- 1950

# Subset first records so we only have those from time window
swe_subset <- subset(swe, year >= min_yr)

# Add one observation for each year so that all years are in the data.
# Remove one obs later. This is a trick to not miss the years with zero new species.
row <- swe_subset[1, ]
row[, -which(names(swe_subset) == "year")] <- NA
year_column_index <- which(names(swe_subset) == "year")
swe_subset_add_one <- lapply(min_yr:max_yr, function(y) {
  row <- swe[1, ]
  row[, year_column_index] <- NA
  row$year <- y
  return(row)
})
swe_subset_add_one <- do.call("rbind", swe_subset_add_one)
swe_subset <- rbind(swe_subset, swe_subset_add_one)

# Turn new species per year into a time series.
swe_timeseries <- swe_subset %>% 
  group_by(year) %>% 
  summarise(n = n() - 1) %>%
  arrange(year)

# Visualize the time series.
ggplot(data = swe_timeseries, aes(x=year, y=n)) +
  geom_bar(stat="identity", color="black", fill="dark orange") + ylab("Antal nya arter") + xlab("År")


# The model of development of future new species per year for JAGS
model_string <- "model {
  for (year in 1:n_years) {
    N_birds[year] ~ dpois(lambda[year])
    lambda[year] <- exp(alpha + beta*year)
  }
  
  alpha ~ dnorm(0, 0.0625)
  beta ~ dnorm(0, 0.0625)
  
  }
"

# Create a list with the data for JAGS
data_list <- list(N_birds = swe_timeseries$n,
                  n_years = nrow(swe_timeseries))

# Specifying MCMC iterations.
n_adapt <- 5000
n_burnin <- 5000
n_mcmc_iterations <- 10000

# Compiling model
m <- jags.model(textConnection(model_string), 
                data = data_list,
                n.chains = 3, 
                n.adapt = n_adapt)

# Running burnin
update(m, n_burnin)

# Generating MCMC samples
s <- coda.samples(m, 
                  variable.names = c("alpha", "beta"), 
                  n.iter = n_mcmc_iterations, 
                  thin = 1)

# Plot the mcmc samples to assess mixing and posterior distribution estimate
plot(s)

# Merging the three MCMC chains into one matrix
ms <- as.matrix(s)

# Simulate trajectories of lambda based on samples of alpha and beta
n_years_pred <- 50
time <- matrix(1:n_years_pred, ncol = n_years_pred)
alpha <- ms[, 1]
beta <- ms[, 2]
lambda_trajs <- exp(sweep(beta %*% time, 1, alpha, "+"))

# Make a fanplot of the lambda trajectories
data.frame(time = t(time), t(lambda_trajs)) %>% 
  reshape2::melt(id = c("time")) %>%
  subset(variable == "X1") %>%
  ggplot(aes(x = time, y = value)) +
  geom_line(alpha=0.4) +
  # geom_interval() 
  scale_linetype_manual(values=c("solid", "dashed", "dotted"))

# Simulate trajectories of future birds
N_birds_trajs <- rpois(length(lambda_trajs), lambda_trajs)
dim(N_birds_trajs) <- dim(lambda_trajs)

# Making path example
path_ex <- temp[51:100,]
path_ex$year <- 2018:2067
combo <- rbind.data.frame(swe_timeseries, path_ex[,c(1,4)])
combo$Typ <- c(rep("Historiskt",nrow(swe_timeseries)), rep("Predikterat", n_years_pred))

# Plotting path example
ggplot(data = combo, aes(x=year, y=n, fill=Typ)) +
  geom_bar(stat="identity", color="black") + 
  ylab("Antal nya arter") + 
  xlab("År") + 
  geom_segment(aes(x=2018, y=0), yend=c(0-10), xend=2018, size=2 )



###
#   Section 2 - modelling the development of peoples' national lists.
###

# Read data on survival probabilities and clean/format it
surv_prob <- read.table(file = "Mx_1x1.txt", skip = 2, header = TRUE, stringsAsFactors = FALSE) %>%
  filter(Year == 2016)
tempmat <- as.data.frame(matrix(data=as.numeric(surv_prob[105, ]), ncol=5, nrow=26, byrow=T))
colnames(tempmat) <- colnames(surv_prob)
surv_prob <- rbind.data.frame(surv_prob[1:105, ], tempmat)
surv_prob$Age <- 0:130

# Read data on twitchers ages
ages <- read.csv(file = "ages.csv", stringsAsFactors = FALSE)

# Add sex information (females have higher survival probability)
ages$sex <- "Male"
ages$sex[c(51, 65) ] <- "Female"

# Removing those without public ages
ages <- na.omit(ages)

# Generate ID-name key
member_ids <- unique(ages$member_id)
ids <- 1:length(member_ids)
names(member_ids) <- ids
names(ids) <- member_ids

# Subset twitch data on those with age data
top100_with_age <- subset(top100, member_id %in% unique(ages$member_id))

# Subset the new species based on twitcher data above
swe_subset_top100 <- subset(swe, year >= year(min(top100_with_age$date)))

# Extract everyones twitches. First create placeholder
new_birds_seen_df <- data.frame(member_id = c(), seen = c())

# Loop over all new species since 2001 and check who has seen what.
for (i in 1:nrow(swe_subset_top100)) {
  # Take out year and species of row
  year <- swe_subset_top100$year[i]
  species <- swe_subset_top100$species[i]
  
  # Check who are the active bird watchers
  active_member_ids <- unique(top100_with_age$member_id[top100_with_age$start_year <= year])
  
  # Check who has seen the bird and who has not during the year
  member_ids_seen <- unique(top100_with_age$member_id[which(as.character(species) == as.character(top100_with_age$bird))])
  member_ids_not_seen <- active_member_ids[!active_member_ids %in% member_ids_seen]
  
  # Make a data.frame of the above
  one_bird_seen_df <- data.frame(member_id = c(member_ids_seen, member_ids_not_seen),
                                 seen = c(rep(1, length(member_ids_seen)), 
                                          rep(0, length(member_ids_not_seen))),
                                 year = year,
                                 species = species)
  
  # Add to the whole
  new_birds_seen_df <- rbind(new_birds_seen_df, one_bird_seen_df)
}

# Add a column of ids starting at 1 to data.frame
new_birds_seen_df$ids <- ids[as.character(new_birds_seen_df$member_id)]

# Make the (very simple) model of twitching skill for JAGS
model_string2 <- "model {
  for (i in 1:n) {
    bird_seen[i] ~ dbern(p[ids[i]])
  }
  for (i in 1:n_ids) {
    p[i] ~ dbeta(1, 1)
  }
}
"

# Make a list of data for JAGS
data_list2 <- list(n = nrow(new_birds_seen_df),
                   n_ids = length(unique(new_birds_seen_df$ids)),
                   bird_seen = new_birds_seen_df$seen,
                   ids = new_birds_seen_df$ids)

# Compile model
m2 <- jags.model(textConnection(model_string2), 
                 data = data_list2,
                 n.chains = 3, 
                 n.adapt = n_adapt)
# Burnin
update(m2, n_burnin)

# Generating MCMC samples
s2 <- coda.samples(m2, 
                   variable.names = c("p"), 
                   n.iter = n_mcmc_iterations, 
                   thin = 1)

# Plot the mcmc samples - or don't, it is massive.
# plot(s2)

# Merging the three MCMC chains into one matrix and naming it
p <- as.matrix(s2)
colnames(p) <- names(ids)
colnames(p) <- ages$member_name[ match(colnames(p), ages$member_id) ]

# Simulate all the twitchers future list developments.
N_member_trajs_incr <- vector(mode = "list", length = length(ids))
for (id in ids) {
  
  # Select twitcher
  i <- ids[id]
  
  # Get data for that twitcher
  age_row <- which(ages$member_id == as.numeric(names(i)))
  
  # Simulate list based on development of new birds
  p_traj <- sweep(1 + 0*N_birds_trajs, 1, p[, i], "*")
  
  # Catch if there for some reason isn't data
  if (length(age_row) == 0) {
    p_traj[] <- 0
  } else {
    
    # Compute survival probability trajectory
    age <- ages$age[age_row]
    surv_prob_traj <- surv_prob %>% 
      filter(Age %in% age:(age + n_years_pred - 1)) %>%
      arrange(Age)
    
    # Get appropriate survival probabilities based on sex.
    if (ages$sex[age_row] == "Male") {
      surv_prob_traj <- as.numeric(surv_prob_traj$Male)
    } else {
      surv_prob_traj <- as.numeric(surv_prob_traj$Female)
    }
    
    # Compute actual survival of twitcher
    for (j in 1:nrow(p_traj)) {
      p_traj[j, ] <- p_traj[j, ] * cumprod(rbinom(n_years_pred, 1, 1 - surv_prob_traj))
    }
  }
  
  # Compoute trajectory of new twitches per year
  traj <- rbinom(length(N_birds_trajs), N_birds_trajs, p_traj)
  dim(traj) <- dim(N_birds_trajs)
  
  # Compute accumulation of new species.
  N_member_trajs_incr[[i]] <- t(apply(traj, 1, cumsum))
}

# Reformat data to compute individual trajectories since 2001
df <- top100_with_age %>% 
  mutate(year = year(date)) %>% 
  arrange(member_id, date) %>%
  group_by(member_id) %>%
  mutate(n = n_birds_seen - (n() - 1):0) %>%
  group_by(member_id, year) %>%
  summarise(n = max(n))

# Extract current number of species in Sweden
number_birds_seen_df <- top100_with_age %>% 
  group_by(member_id) %>% 
  summarise(n = unique(n_birds_seen))
n_birds_seen <- number_birds_seen_df$n
names(n_birds_seen) <- number_birds_seen_df$member_id

# Adding the simulated trajectories
N_member_trajs <- vector(mode = "list", length = length(N_member_trajs_incr))
for (i in ids) {
  starting_value <- n_birds_seen[as.character(member_ids[i])]
  N_member_trajs[[i]] <- N_member_trajs_incr[[i]] + starting_value
}

# Computing how many years into the future a person gets 500 species
first_500 <- lapply(ids, function(i) {
  apply(N_member_trajs[[i]], 1, function(vec) {
    indexes <- which(vec >= 500)
    if (length(indexes) == 0) {
      return(Inf)
    } else {
      return(min(indexes))
    }
  })
})
first_500 <- do.call("cbind", first_500)

# Computing who wins which year
winner <-  apply(first_500, 1, which.min)
winner <- member_ids[as.character(winner)]
first_500_year <- first_500 + 2017

# Computing probability of ever reaching 500
prob_500 <- apply(first_500_year, 2, function(x) {
  temp <- table(x)/sum(table(x))
  return(1 - temp[length(temp)])
})
names(prob_500) <- ages$member_name[ match(names(prob_500), ages$member_id) ]

# Computing when a person reached 500 with various probabilities
when_500 <- apply(first_500_year, 2, function(x) {
  # keep <- x[which(x<Inf)]
  return(quantile(x=x, probs=c(1/10, 1/2, 0.75)))
})

# 10% 50% and 90% probabilities of someone having reached 500 species.
quantile(apply(first_500_year, 1, min), c(1/10, 1/2, 9/10))

# The probability of nobody of the analysed people ever reaching 500
nobody <- apply(first_500_year, 1, function(r) {
  return(all(r==Inf))
})
table(nobody)

# Winning percentages
winners <- table(winner)/sum(table(winner))*100
names(winners) <- ages$member_name[ match(names(winners), ages$member_id) ]
winners

