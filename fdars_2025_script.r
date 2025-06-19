###################################################################################################################################################

### FDAR Models and Related Content ###

## Loading R packages
library("devtools")
#devtools::install_github("txm676/sars")
library("sars")
library("ggplot2")
library("dplyr")
library("UpSetR")
library("corrplot")

#################################################################

# Initial procedures

## Directory.
setwd("C:\\Users\\Rafael\\Desktop\\fdars_2025\\fdars_2025")

set.seed(999)

#################################################################

## SARs Modeling - Bulrush ##

# Import the dataset
bulrush_SES_RESULTS <- read.table(
  file = "bulrush_SES_RESULTS.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

##########

# Build dataframes for SARs
bulrush_fric_ab <- bulrush_SES_RESULTS[, c("area", "Fric_obs")]
colnames(bulrush_fric_ab) <- c("a", "s")

bulrush_fric_pres_abs <- bulrush_SES_RESULTS[, c("area", "PA_Fric_obs")]
colnames(bulrush_fric_pres_abs) <- c("a", "s")

bulrush_fric_ses <- bulrush_SES_RESULTS[, c("area", "Fric_SES")]
colnames(bulrush_fric_ses) <- c("a", "s")

bulrush_fdiv <- bulrush_SES_RESULTS[, c("area", "Disp_obs")]
colnames(bulrush_fdiv) <- c("a", "s")

bulrush_fdiv_ses <- bulrush_SES_RESULTS[, c("area", "Disp_SES")]
colnames(bulrush_fdiv_ses) <- c("a", "s")

bulrush_freg <- bulrush_SES_RESULTS[, c("area", "Eve_obs")]
colnames(bulrush_freg) <- c("a", "s")

bulrush_freg_ses <- bulrush_SES_RESULTS[, c("area", "Eve_SES")]
colnames(bulrush_freg_ses) <- c("a", "s")

# Preview dataframes
head(bulrush_fric_ab)
head(bulrush_fric_ses)
head(bulrush_fric_pres_abs)
head(bulrush_fdiv)
head(bulrush_fdiv_ses)
head(bulrush_freg)
head(bulrush_freg_ses)

##########

# Functions for SARs models
# Filter best SAR models (ΔAICc ≤ 2 & weight ≥ 0.1)
filter_sar_models <- function(fit_obj, delta = 2, weight_min = 0.1) {
  rank_tbl <- summary(fit_obj)$Model_table
  best_aic <- min(rank_tbl$AICc, na.rm = TRUE)
  
  keep <- (rank_tbl$AICc - best_aic <= delta) & (rank_tbl$Weight >= weight_min)
  keep <- keep | rank_tbl$AICc == best_aic
  
  rank_tbl[keep, , drop = FALSE]
}

# Refit selected SAR models with their specific sar_<model> functions
refit_top_sar_models <- function(fit_obj, data, delta = 2, weight_min = 0.1) {
  top_models <- filter_sar_models(fit_obj, delta, weight_min)$Model
  message("Refitting models: ", paste(top_models, collapse = ", "))
  
  refits <- lapply(top_models, function(model_name) {
    fn <- get(paste0("sar_", model_name), mode = "function")
    
    if (model_name == "linear") {
      # Linear model requires fewer arguments
      fn(
        data      = data,
        normaTest = "shapiro",
        homoTest  = "cor.fitted"
      )
    } else {
      # All other models
      fn(
        data       = data,
        normaTest  = "shapiro",
        homoTest   = "cor.fitted",
        grid_start = "exhaustive",
        grid_n     = 1000
      )
    }
  })
  
  names(refits) <- top_models
  refits
}

####################

# FRic observed - abundance

# Run the sar_average function using all models
fit_bulrush_fric_ab <- sar_average(data = bulrush_fric_ab,
                                     #normaTest ="shapiro",
                                     #homoTest = "cor.fitted",
                                     neg_check = FALSE,
                                     confInt = FALSE,
                                     ciN = 100,
                                     grid_start = "exhaustive",
                                     grid_n = 1000)
summary(fit_bulrush_fric_ab)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_bulrush_fric_ab,
  data    = bulrush_fric_ab
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_bulrush_frics_ab)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "fric_ab_bul.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Bulrush -", best_model_name),
     ylab     = "FRic (abundance)",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

####################

# FRic observed - presence absence

# Run the sar_average function using all models
fit_bulrush_fric_pres_abs <- sar_average(data = bulrush_fric_pres_abs,
                                              #normaTest ="shapiro",
                                              #homoTest = "cor.fitted",
                                              neg_check = FALSE,
                                              confInt = FALSE,
                                              ciN = 100,
                                              grid_start = "exhaustive",
                                              grid_n = 1000)
summary(fit_bulrush_fric_pres_abs)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_bulrush_fric_pres_abs,
  data    = bulrush_fric_pres_abs
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_bulrush_fric_pres_abs)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "fric_pa_bul.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Bulrush -", best_model_name),
     ylab     = "FRic (presence-absence)",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

####################

# FDiv observed

# Run the sar_average function using all models
fit_bulrush_fdiv <- sar_average(data = bulrush_fdiv,
                                             #normaTest ="shapiro",
                                             #homoTest = "cor.fitted",
                                             neg_check = FALSE,
                                             confInt = FALSE,
                                             ciN = 100,
                                             grid_start = "exhaustive",
                                             grid_n = 1000)
summary(fit_bulrush_fdiv)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_bulrush_fdiv,
  data    = bulrush_fdiv
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_bulrush_fdiv)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "fdiv_bul.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Bulrush -", best_model_name),
     ylab     = "FDiv",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

####################

# FReg observed

# Run the sar_average function using all models
fit_bulrush_freg <- sar_average(data = bulrush_freg,
                                #normaTest ="shapiro",
                                #homoTest = "cor.fitted",
                                neg_check = FALSE,
                                confInt = FALSE,
                                ciN = 100,
                                grid_start = "exhaustive",
                                grid_n = 1000)
summary(fit_bulrush_freg)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_bulrush_freg,
  data    = bulrush_freg
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_bulrush_freg)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "freg_bul.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Bulrush -", best_model_name),
     ylab     = "FReg",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

########################################

# Log-log power models

# List of data frames and y-axis labels
data_list <- list(
  bulrush_fric_ab       = list(df = bulrush_fric_ab,
                               ylab = "Log(FRic abundance)"),
  bulrush_fric_pres_abs = list(df = bulrush_fric_pres_abs,
                               ylab = "Log(FRic presence–absence)"),
  bulrush_fdiv          = list(df = bulrush_fdiv,
                               ylab = "Log(FDiv abundance)"),
  bulrush_freg          = list(df = bulrush_freg,
                               ylab = "Log(FReg abundance)")
)

# Loop over datasets, fit model, save summary and plot
for (name in names(data_list)) {
  
  dat   <- data_list[[name]]$df
  y_lab <- data_list[[name]]$ylab
  
  # Fit log-log power model
  fit_fdar_log <- lin_pow(
    data       = dat,
    normaTest  = "shapiro",
    homoTest   = "cor.fitted",
    compare    = TRUE,
    con        = 1
  )
  
  # Print summary
  message("\n----- ", name, " -----")
  print(summary(fit_fdar_log))
  
  # Save PNG plot
  png(filename = paste0(name, "_loglog_power.png"),
      width    = 944,    
      height   = 566,    
      units    = "px",
      res      = 300)
  
  plot(
    fit_fdar_log,
    ylab     = y_lab,
    xlab     = "Log(Area)",
    ModTitle = paste("Bulrush –", name, "log-log power"),
    cex.lab  = 1.2,
    cex.axis = 1
  )
  
  dev.off()
}

########################################

# Threshold SARs for SES metrics

# Assemble datasets and labels
ses_datasets <- list(
  Fric_SES = list(df = bulrush_fric_ses, ylab = "FRic SES"),
  Disp_SES = list(df = bulrush_fdiv_ses, ylab = "FDiv SES"),
  Eve_SES  = list(df = bulrush_freg_ses, ylab = "FReg SES")
)

# Global parameters
intz_set  <- 0.001   # breakpoint search interval
n_starts  <- 5       # number of initial start lines
use_cores <- 4       # parallel cores

# Helper: fit breakpoint models (ContOne + ZslopeOne)
fit_thr_models <- function(df_A_S) {
  sar_threshold(
    data      = df_A_S,
    mod       = c("ContOne", "ZslopeOne"),
    interval  = intz_set,
    logAxes   = "area",   # semi-log (log10 Area, unlogged S)
    con       = 0.1,
    nisl      = n_starts,
    logT      = log10,
    parallel  = TRUE,
    cores     = use_cores
  )
}

# Main loop over SES metrics
bulrush_thr_results <- lapply(names(ses_datasets), function(metric) {
  
  # Prepare data: rename to A (area) and S (response)
  df_raw <- ses_datasets[[metric]]$df
  ylab   <- ses_datasets[[metric]]$ylab
  df_AS  <- df_raw %>% rename(A = a, S = s)
  
  # Fit breakpoint models
  fit_obj <- fit_thr_models(df_AS)
  
  # Extract and rename ranking table
  rank_tbl <- summary(fit_obj)[[2]]
  bp_rank <- rank_tbl %>%
    mutate(Model = rownames(rank_tbl)) %>%
    select(Model, AICc, R2, Th1, Th2, seg1, seg2, seg3) %>%
    arrange(AICc)
  
  list(
    ranking   = bp_rank,
    fit_bp    = fit_obj,
    ylab      = ylab
  )
})

# Set names of results list
names(bulrush_thr_results) <- names(ses_datasets)

# Results
bulrush_thr_results$Fric_SES$ranking
bulrush_thr_results$Disp_SES$ranking
bulrush_thr_results$Eve_SES$ranking

########################################

# SES patterns for FRic, FDiv, FReg

# Create a named list of SES dataframes
ses_datasets <- list(
  Fric_SES = list(
    df = bulrush_SES_RESULTS %>%
      transmute(a = area, s = Fric_SES, p = Fric_SES_P),
    ylab = "FRic SES"
  ),
  Disp_SES = list(
    df = bulrush_SES_RESULTS %>%
      transmute(a = area, s = Disp_SES, p = Disp_SES_P),
    ylab = "FDiv SES"
  ),
  Eve_SES = list(
    df = bulrush_SES_RESULTS %>%
      transmute(a = area, s = Eve_SES, p = Eve_SES_P),
    ylab = "FReg SES"
  )
)

# Define colors for each pattern
colors <- c(random = "dodgerblue",
            clustering = "green",
            overdispersion = "red")

# Loop over the three metrics
for (met in names(ses_datasets)) {
  
  df <- ses_datasets[[met]]$df
  ylab_text <- ses_datasets[[met]]$ylab
  
  # Add pattern classification
  use_p <- "p" %in% names(df)
  
  df <- df %>%
    mutate(pattern = case_when(
      s > 0  & (!use_p | p < 0.05) ~ "overdispersion",
      s < 0  & (!use_p | p < 0.05) ~ "clustering",
      TRUE                         ~ "random"
    ))
  
  # Print frequency table
  cat("\n=== ", met, " ===\n")
  print(table(df$pattern))
  
  # Null (intercept-only) model if its line is needed later
  null_mod <- lm(s ~ 1, data = df)
  
  # Make and save plot
  png(filename = paste0(tolower(met), "_bulrush.png"),
      width    = 1600, height = 1200, units = "px", res = 300)
  
  plot(s ~ log10(a),
       data  = df,
       col   = colors[df$pattern],
       pch   = 19,
       xlab  = "log10(Area)",
       ylab  = ylab_text,
       main  = paste("Bulrush –", ylab_text))
  
  # Uncomment next line to add the null (horizontal) line
  # abline(h = coef(null_mod), col = "black", lwd = 2)
  
  dev.off()
}

#################################################################

## SARs Modeling - Grassland ##

# Import the dataset
grass_SES_RESULTS <- read.table(
  file = "grass_SES_RESULTS.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

##########

# Build dataframes for SARs
grass_fric_ab <- grass_SES_RESULTS[, c("area", "Fric_obs")]
colnames(grass_fric_ab) <- c("a", "s")

grass_fric_pres_abs <- grass_SES_RESULTS[, c("area", "PA_Fric_obs")]
colnames(grass_fric_pres_abs) <- c("a", "s")

grass_fric_ses <- grass_SES_RESULTS[, c("area", "Fric_SES")]
colnames(grass_fric_ses) <- c("a", "s")

grass_fdiv <- grass_SES_RESULTS[, c("area", "Disp_obs")]
colnames(grass_fdiv) <- c("a", "s")

grass_fdiv_ses <- grass_SES_RESULTS[, c("area", "Disp_SES")]
colnames(grass_fdiv_ses) <- c("a", "s")

grass_freg <- grass_SES_RESULTS[, c("area", "Eve_obs")]
colnames(grass_freg) <- c("a", "s")

grass_freg_ses <- grass_SES_RESULTS[, c("area", "Eve_SES")]
colnames(grass_freg_ses) <- c("a", "s")

# Preview dataframes
head(grass_fric_ab)
head(grass_fric_ses)
head(grass_fric_pres_abs)
head(grass_fdiv)
head(grass_fdiv_ses)
head(grass_freg)
head(grass_freg_ses)

##########

# Functions for SARs models
# Filter best SAR models (ΔAICc ≤ 2 & weight ≥ 0.1)
filter_sar_models <- function(fit_obj, delta = 2, weight_min = 0.1) {
  rank_tbl <- summary(fit_obj)$Model_table
  best_aic <- min(rank_tbl$AICc, na.rm = TRUE)
  
  keep <- (rank_tbl$AICc - best_aic <= delta) & (rank_tbl$Weight >= weight_min)
  keep <- keep | rank_tbl$AICc == best_aic
  
  rank_tbl[keep, , drop = FALSE]
}

# Refit selected SAR models with their specific sar_<model> functions
refit_top_sar_models <- function(fit_obj, data, delta = 2, weight_min = 0.1) {
  top_models <- filter_sar_models(fit_obj, delta, weight_min)$Model
  message("Refitting models: ", paste(top_models, collapse = ", "))
  
  refits <- lapply(top_models, function(model_name) {
    fn <- get(paste0("sar_", model_name), mode = "function")
    
    if (model_name == "linear") {
      # Linear model requires fewer arguments
      fn(
        data      = data,
        normaTest = "shapiro",
        homoTest  = "cor.fitted"
      )
    } else {
      # All other models
      fn(
        data       = data,
        normaTest  = "shapiro",
        homoTest   = "cor.fitted",
        grid_start = "exhaustive",
        grid_n     = 1000
      )
    }
  })
  
  names(refits) <- top_models
  refits
}

####################

# FRic observed - abundance

# Run the sar_average function using all models
fit_grass_fric_ab <- sar_average(data = grass_fric_ab,
                                   #normaTest ="shapiro",
                                   #homoTest = "cor.fitted",
                                   neg_check = FALSE,
                                   confInt = FALSE,
                                   ciN = 100,
                                   grid_start = "exhaustive",
                                   grid_n = 1000)
summary(fit_grass_fric_ab)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_grass_fric_ab,
  data    = grass_fric_ab
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_grass_fric_ab)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "fric_ab_grass.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Grassland -", best_model_name),
     ylab     = "FRic (abundance)",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

####################

# FRic observed - presence absence

# Run the sar_average function using all models
fit_grass_fric_pres_abs <- sar_average(data = grass_fric_pres_abs,
                                         #normaTest ="shapiro",
                                         #homoTest = "cor.fitted",
                                         neg_check = FALSE,
                                         confInt = FALSE,
                                         ciN = 100,
                                         grid_start = "exhaustive",
                                         grid_n = 1000)
summary(fit_grass_fric_pres_abs)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_grass_fric_pres_abs,
  data    = grass_fric_pres_abs
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_grass_fric_pres_abs)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "fric_pa_grass.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Grassland -", best_model_name),
     ylab     = "FRic (presence-absence)",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

####################

# FDiv observed

# Run the sar_average function using all models
fit_grass_fdiv <- sar_average(data = grass_fdiv,
                                #normaTest ="shapiro",
                                #homoTest = "cor.fitted",
                                neg_check = FALSE,
                                confInt = FALSE,
                                ciN = 100,
                                grid_start = "exhaustive",
                                grid_n = 1000)
summary(fit_grass_fdiv)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_grass_fdiv,
  data    = grass_fdiv
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_grass_fdiv)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "fdiv_grass.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Grassland -", best_model_name),
     ylab     = "FDiv",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

####################

# FReg observed

# Run the sar_average function using all models
fit_grass_freg <- sar_average(data = grass_freg,
                                #normaTest ="shapiro",
                                #homoTest = "cor.fitted",
                                neg_check = FALSE,
                                confInt = FALSE,
                                ciN = 100,
                                grid_start = "exhaustive",
                                grid_n = 1000)
summary(fit_grass_freg)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_grass_freg,
  data    = grass_freg
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_grass_freg)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "freg_grass.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Grassland -", best_model_name),
     ylab     = "FReg",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

########################################

# Log-log power models

# List of data frames and y-axis labels
data_list <- list(
  grass_fric_ab       = list(df = grass_fric_ab,
                               ylab = "Log(FRic abundance)"),
  grass_fric_pres_abs = list(df = grass_fric_pres_abs,
                               ylab = "Log(FRic presence–absence)"),
  grass_fdiv          = list(df = grass_fdiv,
                               ylab = "Log(FDiv abundance)"),
  grass_freg          = list(df = grass_freg,
                               ylab = "Log(FReg abundance)")
)

# Loop over datasets, fit model, save summary and plot
for (name in names(data_list)) {
  
  dat   <- data_list[[name]]$df
  y_lab <- data_list[[name]]$ylab
  
  # Fit log-log power model
  fit_fdar_log <- lin_pow(
    data       = dat,
    normaTest  = "shapiro",
    homoTest   = "cor.fitted",
    compare    = TRUE,
    con        = 1
  )
  
  # Print summary
  message("\n----- ", name, " -----")
  print(summary(fit_fdar_log))
  
  # Save PNG plot
  png(filename = paste0(name, "_loglog_power.png"),
      width    = 944,    
      height   = 566,    
      units    = "px",
      res      = 300)
  
  plot(
    fit_fdar_log,
    ylab     = y_lab,
    xlab     = "Log(Area)",
    ModTitle = paste("Grassland –", name, "log-log power"),
    cex.lab  = 1.2,
    cex.axis = 1
  )
  
  dev.off()
}

########################################

# Threshold SARs for SES metrics

# Assemble datasets and labels
ses_datasets <- list(
  Fric_SES = list(df = grass_fric_ses, ylab = "FRic SES"),
  Disp_SES = list(df = grass_fdiv_ses, ylab = "FDiv SES"),
  Eve_SES  = list(df = grass_freg_ses, ylab = "FReg SES")
)

# Global parameters
intz_set  <- 0.001   # breakpoint search interval
n_starts  <- 5       # number of initial start lines
use_cores <- 4       # parallel cores

# Helper: fit breakpoint models (ContOne + ZslopeOne)
fit_thr_models <- function(df_A_S) {
  sar_threshold(
    data      = df_A_S,
    mod       = c("ContOne", "ZslopeOne"),
    interval  = intz_set,
    logAxes   = "area",   # semi-log (log10 Area, unlogged S)
    con       = 0.1,
    nisl      = n_starts,
    logT      = log10,
    parallel  = TRUE,
    cores     = use_cores
  )
}

# Main loop over SES metrics
grass_thr_results <- lapply(names(ses_datasets), function(metric) {
  
  # Prepare data: rename to A (area) and S (response)
  df_raw <- ses_datasets[[metric]]$df
  ylab   <- ses_datasets[[metric]]$ylab
  df_AS  <- df_raw %>% rename(A = a, S = s)
  
  # Fit breakpoint models
  fit_obj <- fit_thr_models(df_AS)
  
  # Extract and rename ranking table
  rank_tbl <- summary(fit_obj)[[2]]
  bp_rank <- rank_tbl %>%
    mutate(Model = rownames(rank_tbl)) %>%
    select(Model, AICc, R2, Th1, Th2, seg1, seg2, seg3) %>%
    arrange(AICc)
  
  list(
    ranking   = bp_rank,
    fit_bp    = fit_obj,
    ylab      = ylab
  )
})

# Set names of results list
names(grass_thr_results) <- names(ses_datasets)

# Results
grass_thr_results$Fric_SES$ranking
grass_thr_results$Disp_SES$ranking
grass_thr_results$Eve_SES$ranking

########################################

# SES patterns for FRic, FDiv, FReg

# Create a named list of SES dataframes
ses_datasets <- list(
  Fric_SES = list(
    df = grass_SES_RESULTS %>%
      transmute(a = area, s = Fric_SES, p = Fric_SES_P),
    ylab = "FRic SES"
  ),
  Disp_SES = list(
    df = grass_SES_RESULTS %>%
      transmute(a = area, s = Disp_SES, p = Disp_SES_P),
    ylab = "FDiv SES"
  ),
  Eve_SES = list(
    df = grass_SES_RESULTS %>%
      transmute(a = area, s = Eve_SES, p = Eve_SES_P),
    ylab = "FReg SES"
  )
)

# Define colors for each pattern
colors <- c(random = "dodgerblue",
            clustering = "green",
            overdispersion = "red")

# Loop over the three metrics
for (met in names(ses_datasets)) {
  
  df <- ses_datasets[[met]]$df
  ylab_text <- ses_datasets[[met]]$ylab
  
  # Add pattern classification
  use_p <- "p" %in% names(df)
  
  df <- df %>%
    mutate(pattern = case_when(
      s > 0  & (!use_p | p < 0.05) ~ "overdispersion",
      s < 0  & (!use_p | p < 0.05) ~ "clustering",
      TRUE                         ~ "random"
    ))
  
  # Print frequency table
  cat("\n=== ", met, " ===\n")
  print(table(df$pattern))
  
  # Null (intercept-only) model if its line is needed later
  null_mod <- lm(s ~ 1, data = df)
  
  # Make and save plot
  png(filename = paste0(tolower(met), "_grass.png"),
      width    = 1600, height = 1200, units = "px", res = 300)
  
  plot(s ~ log10(a),
       data  = df,
       col   = colors[df$pattern],
       pch   = 19,
       xlab  = "log10(Area)",
       ylab  = ylab_text,
       main  = paste("Grassland –", ylab_text))
  
  # Uncomment next line to add the null (horizontal) line
  # abline(h = coef(null_mod), col = "black", lwd = 2)
  
  dev.off()
}

#################################################################
#################################################################

## SARs Modeling - Washout ##

# Import the dataset
wash_SES_RESULTS <- read.table(
  file = "wash_SES_RESULTS.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

##########

# Build dataframes for SARs
wash_fric_ab <- wash_SES_RESULTS[, c("area", "Fric_obs")]
colnames(wash_fric_ab) <- c("a", "s")

wash_fric_pres_abs <- wash_SES_RESULTS[, c("area", "PA_Fric_obs")]
colnames(wash_fric_pres_abs) <- c("a", "s")

wash_fric_ses <- wash_SES_RESULTS[, c("area", "Fric_SES")]
colnames(wash_fric_ses) <- c("a", "s")

wash_fdiv <- wash_SES_RESULTS[, c("area", "Disp_obs")]
colnames(grass_fdiv) <- c("a", "s")

wash_fdiv_ses <- wash_SES_RESULTS[, c("area", "Disp_SES")]
colnames(wash_fdiv_ses) <- c("a", "s")

wash_freg <- wash_SES_RESULTS[, c("area", "Eve_obs")]
colnames(wash_freg) <- c("a", "s")

wash_freg_ses <- wash_SES_RESULTS[, c("area", "Eve_SES")]
colnames(wash_freg_ses) <- c("a", "s")

# Preview dataframes
head(wash_fric_ab)
head(wash_fric_ses)
head(wash_fric_pres_abs)
head(wash_fdiv)
head(wash_fdiv_ses)
head(wash_freg)
head(wash_freg_ses)

##########

# Functions for SARs models
# Filter best SAR models (ΔAICc ≤ 2 & weight ≥ 0.1)
filter_sar_models <- function(fit_obj, delta = 2, weight_min = 0.1) {
  rank_tbl <- summary(fit_obj)$Model_table
  best_aic <- min(rank_tbl$AICc, na.rm = TRUE)
  
  keep <- (rank_tbl$AICc - best_aic <= delta) & (rank_tbl$Weight >= weight_min)
  keep <- keep | rank_tbl$AICc == best_aic
  
  rank_tbl[keep, , drop = FALSE]
}

# Refit selected SAR models with their specific sar_<model> functions
refit_top_sar_models <- function(fit_obj, data, delta = 2, weight_min = 0.1) {
  top_models <- filter_sar_models(fit_obj, delta, weight_min)$Model
  message("Refitting models: ", paste(top_models, collapse = ", "))
  
  refits <- lapply(top_models, function(model_name) {
    fn <- get(paste0("sar_", model_name), mode = "function")
    
    if (model_name == "linear") {
      # Linear model requires fewer arguments
      fn(
        data      = data,
        normaTest = "shapiro",
        homoTest  = "cor.fitted"
      )
    } else {
      # All other models
      fn(
        data       = data,
        normaTest  = "shapiro",
        homoTest   = "cor.fitted",
        grid_start = "exhaustive",
        grid_n     = 1000
      )
    }
  })
  
  names(refits) <- top_models
  refits
}

####################

# FRic observed - abundance

# Run the sar_average function using all models
fit_wash_fric_ab <- sar_average(data = wash_fric_ab,
                                 #normaTest ="shapiro",
                                 #homoTest = "cor.fitted",
                                 neg_check = FALSE,
                                 confInt = FALSE,
                                 ciN = 100,
                                 grid_start = "exhaustive",
                                 grid_n = 1000)
summary(fit_wash_fric_ab)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_wash_fric_ab,
  data    = wash_fric_ab
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_wash_fric_ab)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "fric_ab_wash.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Washout -", best_model_name),
     ylab     = "FRic (abundance)",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

####################

# FRic observed - presence absence

# Run the sar_average function using all models
fit_wash_fric_pres_abs <- sar_average(data = wash_fric_pres_abs,
                                       #normaTest ="shapiro",
                                       #homoTest = "cor.fitted",
                                       neg_check = FALSE,
                                       confInt = FALSE,
                                       ciN = 100,
                                       grid_start = "exhaustive",
                                       grid_n = 1000)
summary(fit_wash_fric_pres_abs)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_wash_fric_pres_abs,
  data    = wash_fric_pres_abs
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_wash_fric_pres_abs)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "fric_pa_wash.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Washout -", best_model_name),
     ylab     = "FRic (presence-absence)",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

####################

# FDiv observed

# Run the sar_average function using all models
fit_wash_fdiv <- sar_average(data = wash_fdiv,
                              #normaTest ="shapiro",
                              #homoTest = "cor.fitted",
                              neg_check = FALSE,
                              confInt = FALSE,
                              ciN = 100,
                              grid_start = "exhaustive",
                              grid_n = 1000)
summary(fit_wash_fdiv)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_wash_fdiv,
  data    = wash_fdiv
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_wash_fdiv)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "fdiv_wash.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Washout -", best_model_name),
     ylab     = "FDiv",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

####################

# FReg observed

# Run the sar_average function using all models
fit_wash_freg <- sar_average(data = wash_freg,
                              #normaTest ="shapiro",
                              #homoTest = "cor.fitted",
                              neg_check = FALSE,
                              confInt = FALSE,
                              ciN = 100,
                              grid_start = "exhaustive",
                              grid_n = 1000)
summary(fit_wash_freg)

####

# Refit best models
refitted_models <- refit_top_sar_models(
  fit_obj = fit_wash_freg,
  data    = wash_freg
)

# Summarize all refitted models
lapply(refitted_models, summary)

# Identify the best model (lowest AICc)
best_model_name <- filter_sar_models(fit_wash_freg)[1, "Model"]
best_model_obj  <- refitted_models[[best_model_name]]

# Save plot as a PNG file
png(filename = "freg_wash.png",
    width    = 1600,
    height   = 1200,
    units    = "px",
    res      = 300)

# Plot the best model
plot(best_model_obj,
     ModTitle = paste("Washout -", best_model_name),
     ylab     = "FReg",
     xlab     = "Area",
     cex.lab  = 1.2,
     cex.axis = 1)

# Close graphics device
dev.off()

########################################

# Log-log power models

# List of data frames and y-axis labels
data_list <- list(
  wash_fric_ab       = list(df = wash_fric_ab,
                             ylab = "Log(FRic abundance)"),
  wash_fric_pres_abs = list(df = wash_fric_pres_abs,
                             ylab = "Log(FRic presence–absence)"),
  wash_fdiv          = list(df = wash_fdiv,
                             ylab = "Log(FDiv abundance)"),
  wash_freg          = list(df = wash_freg,
                             ylab = "Log(FReg abundance)")
)

# Loop over datasets, fit model, save summary and plot
for (name in names(data_list)) {
  
  dat   <- data_list[[name]]$df
  y_lab <- data_list[[name]]$ylab
  
  # Fit log-log power model
  fit_fdar_log <- lin_pow(
    data       = dat,
    normaTest  = "shapiro",
    homoTest   = "cor.fitted",
    compare    = TRUE,
    con        = 1
  )
  
  # Print summary
  message("\n----- ", name, " -----")
  print(summary(fit_fdar_log))
  
  # Save PNG plot
  png(filename = paste0(name, "_loglog_power.png"),
      width    = 944,    
      height   = 566,    
      units    = "px",
      res      = 300)
  
  plot(
    fit_fdar_log,
    ylab     = y_lab,
    xlab     = "Log(Area)",
    ModTitle = paste("Washout –", name, "log-log power"),
    cex.lab  = 1.2,
    cex.axis = 1
  )
  
  dev.off()
}

########################################

# Threshold SARs for SES metrics

# Assemble datasets and labels
ses_datasets <- list(
  Fric_SES = list(df = wash_fric_ses, ylab = "FRic SES"),
  Disp_SES = list(df = wash_fdiv_ses, ylab = "FDiv SES"),
  Eve_SES  = list(df = wash_freg_ses, ylab = "FReg SES")
)

# Global parameters
intz_set  <- 0.001   # breakpoint search interval
n_starts  <- 5       # number of initial start lines
use_cores <- 4       # parallel cores

# Helper: fit breakpoint models (ContOne + ZslopeOne)
fit_thr_models <- function(df_A_S) {
  sar_threshold(
    data      = df_A_S,
    mod       = c("ContOne", "ZslopeOne"),
    interval  = intz_set,
    logAxes   = "area",   # semi-log (log10 Area, unlogged S)
    con       = 0.1,
    nisl      = n_starts,
    logT      = log10,
    parallel  = TRUE,
    cores     = use_cores
  )
}

# Main loop over SES metrics
wash_thr_results <- lapply(names(ses_datasets), function(metric) {
  
  # Prepare data: rename to A (area) and S (response)
  df_raw <- ses_datasets[[metric]]$df
  ylab   <- ses_datasets[[metric]]$ylab
  df_AS  <- df_raw %>% rename(A = a, S = s)
  
  # Fit breakpoint models
  fit_obj <- fit_thr_models(df_AS)
  
  # Extract and rename ranking table
  rank_tbl <- summary(fit_obj)[[2]]
  bp_rank <- rank_tbl %>%
    mutate(Model = rownames(rank_tbl)) %>%
    select(Model, AICc, R2, Th1, Th2, seg1, seg2, seg3) %>%
    arrange(AICc)
  
  list(
    ranking   = bp_rank,
    fit_bp    = fit_obj,
    ylab      = ylab
  )
})

# Set names of results list
names(wash_thr_results) <- names(ses_datasets)

# Results
wash_thr_results$Fric_SES$ranking
wash_thr_results$Disp_SES$ranking
wash_thr_results$Eve_SES$ranking

########################################

# SES patterns for FRic, FDiv, FReg

# Create a named list of SES dataframes
ses_datasets <- list(
  Fric_SES = list(
    df = wash_SES_RESULTS %>%
      transmute(a = area, s = Fric_SES, p = Fric_SES_P),
    ylab = "FRic SES"
  ),
  Disp_SES = list(
    df = wash_SES_RESULTS %>%
      transmute(a = area, s = Disp_SES, p = Disp_SES_P),
    ylab = "FDiv SES"
  ),
  Eve_SES = list(
    df = wash_SES_RESULTS %>%
      transmute(a = area, s = Eve_SES, p = Eve_SES_P),
    ylab = "FReg SES"
  )
)

# Define colors for each pattern
colors <- c(random = "dodgerblue",
            clustering = "green",
            overdispersion = "red")

# Loop over the three metrics
for (met in names(ses_datasets)) {
  
  df <- ses_datasets[[met]]$df
  ylab_text <- ses_datasets[[met]]$ylab
  
  # Add pattern classification
  use_p <- "p" %in% names(df)
  
  df <- df %>%
    mutate(pattern = case_when(
      s > 0  & (!use_p | p < 0.05) ~ "overdispersion",
      s < 0  & (!use_p | p < 0.05) ~ "clustering",
      TRUE                         ~ "random"
    ))
  
  # Print frequency table
  cat("\n=== ", met, " ===\n")
  print(table(df$pattern))
  
  # Null (intercept-only) model if its line is needed later
  null_mod <- lm(s ~ 1, data = df)
  
  # Make and save plot
  png(filename = paste0(tolower(met), "_wash.png"),
      width    = 1600, height = 1200, units = "px", res = 300)
  
  plot(s ~ log10(a),
       data  = df,
       col   = colors[df$pattern],
       pch   = 19,
       xlab  = "log10(Area)",
       ylab  = ylab_text,
       main  = paste("Washout –", ylab_text))
  
  # Uncomment next line to add the null (horizontal) line
  # abline(h = coef(null_mod), col = "black", lwd = 2)
  
  dev.off()
}

#################################################################
#################################################################

## Upset plots - SES patterns for FRic, FDiv, FReg ##

# Grassland

# Create the dataset
island_data <- data.frame(
  Island = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t"),
  FRic = c("Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", 
           "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", 
           "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "-"),
  FDiv = c("Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", 
           "Clustered", "Clustered", "Clustered", "Clustered", "-", "Clustered", "-", "Clustered", 
           "Clustered", "-", "Clustered", "Clustered", "Clustered"),
  FReg = c("Overdispersed", "Overdispersed", "-", "-", "Overdispersed", "-", "-", "Overdispersed", 
           "Overdispersed", "-", "-", "Overdispersed", "Overdispersed", "-", "Overdispersed", "-", "-", "-", "-","Overdispersed")
)

# Convert to a binary matrix where 1 represents the presence of a pattern
island_data_binary <- data.frame(
  FRic = ifelse(island_data$FRic == "Clustered", 1, 0),
  FDiv = ifelse(island_data$FDiv == "Clustered", 1, 0),
  FReg = ifelse(island_data$FReg == "Overdispersed", 1, 0)
)

# Add row names (Islands) to the binary data
rownames(island_data_binary) <- island_data$Island

# Create the upset plot
upset_plot_grass <- upset(island_data_binary, 
                          sets = c("FRic", "FDiv", "FReg"), 
                          order.by = "freq", 
                          main.bar.color = "skyblue", 
                          sets.bar.color = "darkgreen", 
                          set_size.show = TRUE)

upset_plot_grass

####

# Bulrush

# Create the dataset
island_data <- data.frame(
  Island = c("J4", "J2", "J154", "J91", "J25", "J19", "J1", "J31", "J169", "J70", "J69"),
  FRic = c("Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "-", "-", "-"),
  FDiv = c("-", "-", "-", "-", "Clustered", "Clustered", "Clustered", "-", "Clustered", "-", "-"),
  FReg = c("-", "Overdispersed", "-", "Overdispersed", "-", "Overdispersed", "Overdispersed", "Overdispersed", "-", "Overdispersed", "Overdispersed")
)

# Convert to a binary matrix where 1 represents the presence of a pattern
island_data_binary <- data.frame(
  FRic = ifelse(island_data$FRic == "Clustered", 1, 0),
  FDiv = ifelse(island_data$FDiv == "Clustered", 1, 0),
  FReg = ifelse(island_data$FReg == "Overdispersed", 1, 0)
)

# Add row names (Islands) to the binary data
rownames(island_data_binary) <- island_data$Island

# Create the upset plot
upset_plot_bul <- upset(island_data_binary, 
                        sets = c("FRic", "FDiv", "FReg"), 
                        order.by = "freq", 
                        main.bar.color = "skyblue", 
                        sets.bar.color = "darkgreen", 
                        set_size.show = TRUE)

upset_plot_bul

####

# Washout

# Create the dataset
island_data <- data.frame(
  Island = c("A4", "A8", "A7", "A5", "A15", "A9", "A6"),
  FRic = c("Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "Clustered", "-"),
  FDiv = c("Clustered", "Clustered", "-", "Clustered", "Clustered", "-", "-"),
  FReg = c("-", "-", "Overdispersed", "-", "-", "-", "Overdispersed")
)

# Convert the data to a binary matrix
island_data_binary <- data.frame(
  FRic = ifelse(island_data$FRic == "Clustered", 1, 0),
  FDiv = ifelse(island_data$FDiv == "Clustered", 1, 0),
  FReg = ifelse(island_data$FReg == "Overdispersed", 1, 0)
)

# Add row names (Islands) to the binary data
rownames(island_data_binary) <- island_data$Island

# Create the upset plot
upset_plot_wash <- upset(island_data_binary, 
                         sets = c("FRic", "FDiv", "FReg"), 
                         order.by = "freq", 
                         main.bar.color = "skyblue", 
                         sets.bar.color = "darkgreen", 
                         set_size.show = TRUE)
upset_plot_wash

#################################################################
#################################################################

## Correlations ##

# Create a list of datasets
ses_data_list <- list(
  bulrush = bul_SES_RESULTS,
  grass = grass_SES_RESULTS,
  wash = wash_SES_RESULTS
)

# Define a general function to plot correlations
create_corr_plots <- function(data_list, prefix) {
  
  # Define variable sets
  basic_vars <- c("species_richness", "abundance", "Fric_obs", "PA_Fric_obs", "Disp_obs", "Eve_obs")
  ses_vs_es_vars <- c("Fric_SES", "Disp_SES", "Eve_SES", "Fric_ES", "Disp_ES", "Eve_ES")
  fd_3ax_vars <- c("Fric_obs", "Disp_obs", "Eve_obs", "Rich_3ax", "Disp_3ax", "Even_3ax")
  
  rename_labels <- list(
    basic = c("Species richness", "Abundance", "FRic (presence-absence)", "FRic (abundance)", "FDiv", "FReg"),
    ses_vs_es = c("FRic (SES)", "FDiv (SES)", "FReg (SES)", "FRic (ES)", "FDiv (ES)", "FReg (ES)"),
    fd_3ax = c("FRic (Multiple Axes)", "FDiv (Multiple Axes)", "FReg (Multiple Axes)", "FRic (3 Axes)", "FDiv (3 Axes)", "FReg (3 Axes)")
  )
  
  plot_titles <- c(wash = "Washout", bulrush = "Bulrush", grass = "Grassland")
  
  for (type in c("basic", "ses_vs_es", "fd_3ax")) {
    vars <- get(paste0(type, "_vars"))
    labels <- rename_labels[[type]]
    
    for (name in names(data_list)) {
      df <- data_list[[name]]
      correlation_data <- df %>%
        select(all_of(vars)) %>%
        setNames(labels)
      
      correlation_matrix <- cor(correlation_data, use = "pairwise.complete.obs")
      
      # Save table
      table_file <- paste0("correlation_table_", name, "_", type, ".csv")
      write.csv(correlation_matrix, file = table_file, row.names = TRUE)
      
      # Save plot
      png(filename = paste0("correlation_plot_", name, "_", type, ".png"),
          width = 1600, height = 1200, units = "px", res = 300)
      corrplot(correlation_matrix,
               method = "color",
               type = "upper",
               tl.col = "black",
               tl.srt = 45,
               addCoef.col = "black",
               number.cex = 0.7,
               col = colorRampPalette(c("red", "white", "blue"))(200),
               title = plot_titles[[name]],
               mar = c(0, 0, 2, 0),
               cl.pos = "n")
      dev.off()
    }
  }
}

# Run plotting
create_corr_plots(ses_data_list, "SES_results")

#################################################################
#################################################################

## Summary Statistics ##

# Put the three data frames into a named list
ses_tables <- list(
  bulrush = bul_SES_RESULTS,
  grass   = grass_SES_RESULTS,
  wash    = wash_SES_RESULTS
)

# Columns to summarize
columns_to_summarize <- c(
  "abundance", "species_richness",
  "Fric_obs", "PA_Fric_obs",
  "Disp_obs", "Eve_obs"
)

# Function to compute Total, Min, Max, Mean, SD for each column
calculate_summary_stats <- function(df) {
  stats <- sapply(df[columns_to_summarize], function(x) {
    c(
      Total = sum(x, na.rm = TRUE),
      Min   = round(min(x, na.rm = TRUE), 2),
      Max   = round(max(x, na.rm = TRUE), 2),
      Mean  = round(mean(x, na.rm = TRUE), 2),
      SD    = round(sd(x, na.rm = TRUE), 2)
    )
  })
  as.data.frame(t(stats))
}

# Apply to each dataset
summary_stats_list <- lapply(ses_tables, calculate_summary_stats)

# Name the results and print
names(summary_stats_list) <- names(ses_tables)
summary_stats_list

#################################################################