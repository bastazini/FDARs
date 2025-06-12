##############################################################
###Main Script for Calculating FRic, FDisp and Feve using ####
###Hypervolumes, as well as running the null models###########
##############################################################

library("FD") #version 1.1-12
library("hypervolume") #version 2.0.11
library("BAT")
library("dplyr")
library(foreach)
library(doParallel)
library(cluster)

#Set up parallel processing
cores = 12
cl = makeCluster(cores)
registerDoParallel(cl)
i = 1 #Dummy line for RStudio warnings

tr1 <- read.delim("traits_bulrush.txt")
tr2 <- read.delim("traits_grass.txt")
tr3 <- read.delim("traits_wash.txt")
#merge and remove duplicate rows
trAll <- rbind(tr1, tr2, tr3) %>%
  distinct()
if (length(unique(trAll$X)) != nrow(trAll)){
  stop("Pulse")
}
trAll$mass <- log(trAll$mass)
# Standardize:
trAll <- trAll %>% 
  mutate_at(c("mass", "bill_lenght", "bill_angle", 
              "claw_angle", "foot", "wing"),
            ~(scale(.) %>% as.vector))
trAll2 <- trAll
trAll <- select(trAll, -X)

##Use same bandwidth for all hypervolumes (have done this by
#pooling traits across all 3 datasets).
BAND <- hypervolume::estimate_bandwidth(as.matrix(trAll)) 

#Dataset names
ds <- c("bulrush", "grass", "wash")

##LOOP across each dataset, and calculate the observed
#and null model values for each one.
#NOTE: if you want to run the analyses using presence-absence data,
#you will need to manually change the 'abund argument in kernel.build
for (j in 1:length(ds)){

## Loading the databases and preparing the data
DATASET <- ds[j]
traits_dataset <- paste0("traits_", DATASET, ".txt")
comm_dataset <- paste0("commun_", DATASET, ".txt")

trait <- read.delim(traits_dataset, row.names=1)

if (!all(rownames(trait) %in% trAll2$X)){
  stop("Boc")
}
trAll3 <- filter(trAll2, X %in% rownames(trait))
rownames(trAll3) <- trAll3$X
trAll3 <- select(trAll3, !X)
trAll3 <- trAll3[order(rownames(trAll3)),]
if(!identical(rownames(trait), rownames(trAll3))){
  stop("Boc2")
}
trait <- trAll3

# Load data:
comm <- read.delim(comm_dataset, row.names=1)
comm <- as.matrix(comm)

if (!identical(colnames(comm), rownames(trait))){
  stop("sp names do not match")
}

## Build the hypervolume
hv <- kernel.build(comm = as.matrix(comm), 
                   trait = as.matrix(trait), 
                   distance = "euclidean",
                  abund = TRUE, #samples are abundances
     #              abund = FALSE, #samples are presence-absence
                   method = "gaussian", 
                   axes = 0, #If 0, no transformation of data is done.
                   kde.bandwidth = BAND) #bandwidths

##Calculate Fric, Disp and Feve. Note sometimes kernel.evenness errors,
#so tryCatch has been used
richness <- kernel.alpha(comm = hv) 
results_rich <- data.frame(Richness = richness)
dispersion_divergence <- kernel.dispersion(comm = hv,frac=0.1, 
                                           func= 'divergence')
results_disp <- data.frame(Dispersion=dispersion_divergence,
                           Dispersion_div=dispersion_divergence)
evenness <- tryCatch(kernel.evenness(comm = hv),
                     error = function(e) NA)
results_even <- data.frame(Evenness = evenness)

observed <- list(results_rich, results_disp,
                 results_even)

#################################################################
# Null models
# Calculated by randomizing the lines, so that combinations of
# traits that do not make biological sense are not allowed. This
# null model tests specifically if the observed value differs
# from random functional compositions for the communities.
#########################################################

nreps <- 999 # number of replicates of the null model
reps <- replicate(nreps, sample(1:nrow(trait), 
                                size = nrow(trait), 
                                replace = FALSE))
list_null_trait <- 
  apply(reps, 2, function(x){
    trait_null <- trait[x,]
    rownames(trait_null) <- rownames(trait)
    return(trait_null)
  }) # list with matrices of null traits

##Uses a parallel for loop, as the cores argument inside
#kernel.build does not seem to work.
list_df_null = foreach(i=seq(from=1, 
                             to=length(list_null_trait), 
                             by=1), 
                       .inorder = TRUE)  %dopar% { 
    library(BAT)
  
  null_kernel <- 
    kernel.build(comm = as.matrix(comm), 
                 trait = list_null_trait[[i]], 
                 distance = "euclidean",
                 abund = TRUE, #samples are presences/absences
                 method = "gaussian", 
                 axes = 0, #If 0, no transformation of data is done. If 0 < axes <= 1 a PCoA is done with Gower/euclidean distances and as many axes as needed to achieve this proportion of variance explained are selected; if axes > 1 these many axes are selected
                 cores = 1,
                 kde.bandwidth = BAND)
  
  richness <- kernel.alpha(comm = null_kernel) # vector with richness value for each subplot
  dispersion_divergence <- kernel.dispersion(comm = null_kernel,
                                             frac=0.1, 
                                             func= 'divergence')
  evenness <- tryCatch(kernel.evenness(comm = null_kernel),
                       error = function(e) NA)
  df_null <- data.frame(Fric = richness, 
                        Fdis = dispersion_divergence, 
                        Feve = evenness)
  df_null
}#eo foreach                         

res <- list("observed" = observed, "null" = list_df_null)
fn <- paste0(DATASET,"_results.rds")
saveRDS(res, file = fn)

}#eo j

#############################################
#######CREATE THE SES VALUE FILES#################
#######################################################

##Set working directory to location of _results.rds files

for (j in 1:length(ds)){

DATASET <- ds[j]

res <- readRDS(paste0(DATASET, "_results.rds"))

OBS <- res[[1]]
#convert res output into observed and null
list_df_null <- res[[2]]

#check if any null runs have NAs for each variable (only affected Feve)
sunshine <- sapply(list_df_null, function(x){
apply(x, 2, anyNA)
})
if(any(apply(sunshine, 1, any))){
  stop("Koffing")
  }

# Grouping null results for each dimension (Fric, Fdis e Feve) in matrices
null_mat <- do.call(cbind, list_df_null)
null_mat_rich <- null_mat[, 
                          grep("Fric", colnames(null_mat))] # values from the null model for Fric
null_mat_dis <- null_mat[, 
                         grep("Fdis", colnames(null_mat))] # values from the null model for Fdis
null_mat_eve <- null_mat[, 
                         grep("Feve", colnames(null_mat))] # values from the null model for Feve

#################################################################

# Calculating SES for all matrices

comm <- OBS[[1]]

ses_res_rich <- matrix(NA, nrow = nrow(comm), 
                       ncol = 5, 
                       dimnames = list(rownames(comm), 
                                       c("Obs","ses", "ses_p-value",
                                         "es", "es_P-value")))
ses_res_disp <- matrix(NA, nrow = nrow(comm), 
                       ncol = 5, 
                       dimnames = list(rownames(comm), 
                                       c("Obs","ses", "ses_p-value",
                                         "es", "es_P-value")))
ses_res_eve <- matrix(NA, nrow = nrow(comm),
                      ncol = 5,
                      dimnames = list(rownames(comm), 
                                      c("Obs","ses", "ses_p-value",
                                        "es", "es_P-value")))

if (!identical(rownames(OBS[[3]]), rownames(null_mat_rich))){
  stop("41")
}

for(i in 1:nrow(comm)){
  ses_res_rich[i, 1] <- OBS[[1]]$Richness[i]
  ses_res_rich[i, 2:3] <- BAT::ses(obs = OBS[[1]]$Richness[i], 
                                est = unlist(null_mat_rich[i, ]), 
                                param = TRUE, p = TRUE)
  ses_res_rich[i, 4:5] <- BAT::ses(obs = OBS[[1]]$Richness[i], 
                                   est = unlist(null_mat_rich[i, ]), 
                                   param = FALSE, p = TRUE)
  
  ses_res_disp[i, 1] <- OBS[[2]]$Dispersion[i]
  ses_res_disp[i, 2:3] <- BAT::ses(obs = OBS[[2]]$Dispersion[i], 
                                est = unlist(null_mat_dis[i, ]), 
                                param = TRUE, p = TRUE)
  ses_res_disp[i, 4:5] <- BAT::ses(obs = OBS[[2]]$Dispersion[i], 
                                   est = unlist(null_mat_dis[i, ]), 
                                   param = FALSE, p = TRUE)
  
  ses_res_eve[i, 1] <- OBS[[3]]$Evenness[i]
  ses_res_eve[i, 2:3] <- BAT::ses(obs = OBS[[3]]$Evenness[i],
                               est = unlist(null_mat_eve[i, ]), 
                               param = TRUE, p = TRUE)
  ses_res_eve[i, 4:5] <- BAT::ses(obs = OBS[[3]]$Evenness[i],
                                        est = unlist(null_mat_eve[i, ]), 
                                        param = FALSE, p = TRUE)
}

if(!identical(rownames(PA2[[j]]), rownames(OBS[[2]]))){
  stop("your_song")
}

SES_RES <- list("ses_res_rich" = ses_res_rich, 
                "ses_res_disp" = ses_res_disp, 
                "ses_res_eve" = ses_res_eve,
                "presence_absence_obs" = PA2[[j]])

SES_RES <- lapply(SES_RES, function(x){round(x,3)})
saveRDS(SES_RES, file = paste0(DATASET,"_SES_RESULTS.rds"))

}#eo for j

#####################################################
#########Convert SES Res Files into Tables############
####################################################

#Set working directory to location of SES_RESULTS files
lf2 <- c("bulrush_SES_RESULTS.rds", "grass_SES_RESULTS.rds",
         "wash_SES_RESULTS.rds")
ldf2 <- lapply(lf2, readRDS)

CN <- c("Fric_obs", "Fric_SES", "Fric_SES_P",
        "Fric_ES", "Fric_ES_P",
        "Disp_obs", "Disp_SES", "Disp_SES_P",
        "Disp_ES", "Disp_ES_P",
        "Eve_obs", "Eve_SES", "Eve_SES_P",
        "Eve_ES", "Eve_ES_P",
        "PA_Fric_obs", "PA_Disp_obs", "PA_Eve_obs")

FN <- gsub("_.*$", "", lf2)

for (j in 1:length(ldf2)){
  ldf3 <- do.call(cbind, ldf2[[j]])
  colnames(ldf3) <- CN
  write.csv(ldf3, 
            file = paste0(FN[j], "_MainTable.csv"))
}

###################################################################################################################################################
###################################################################################################################################################