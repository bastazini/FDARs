#######################################################
########Run the Jackknife-type cross validation#########
########sensitivity test###############################
####################################################

library(sars)
library(dplyr)
library(foreach)
library(doParallel)
library(cluster)

#Set up parallel processing
cores = 4
cl = makeCluster(cores)
registerDoParallel(cl)
i = 1 #Dummy line for RStudio warnings

#load in datasets (in this order)
dat1 <- read.csv("grass_MainTable.csv")
dat2 <- read.csv("bulrush_MainTable.csv")
dat3 <- read.csv("wash_MainTable.csv")

datL <- list(dat1, dat2, dat3)

##function to remove each of the n largest islands
#in turn, and re-run the sar model comparison
#note - first column in dat must be area
jackknife <- function(dat, n = 3, best){
  
  dat2 <- dat[order(dat[,1], decreasing = TRUE),]
  
  if (n > nrow(dat2)) stop("Idle")
  
  BM <- matrix(nrow = n, ncol = 4)
  
  for (i in 1:n){
    dum <- dat2[-i,]
    sfr2 <- sar_average(data = dum, 
                        grid_start = "exhaustive",
                        grid_n = 1000,
                        verb = FALSE, 
                        display = FALSE)
    sfr3 <- summary(sfr2)
    BM[i,1] <- sfr3$Model_table$Model[1]
    if (!best %in% sfr3$Model_table$Model) stop("Debra")
    BM[i,2] <- which(sfr3$Model_table$Model == best)
    BM[i,3] <- sfr3$Model_table$Shape[1]
    BM[i,4] <- sfr3$Model_table$Asymptote[1]
  }
  colnames(BM) <- c("Jacknife_BM", "Orig_BM_rank",
                    "Jacknife_BM_shape",
                    "Jacknife_BM_asymp")
  return(BM)
}

##Function to iterate the main jackknife function across
#each dataset and FD metric, and format the output into
#a table for the SI.
#fig1bm = a vector of best model fit names, matching the 
#order of the fitting process (dataset followed by ys inside
#the function).
sar_jackknife <- function(datL,
                     dataset = c("Grassland", "Bulrush",
                        "Washout"),
                     n = 3,
                     fig1bm = fig1bm){
  
  ys <- c("PA_Fric_obs", "Fric_obs",
          "Disp_obs", "Eve_obs")
  
  la = foreach(i=seq(from=1, to=length(1:3), 
                          by=1), 
                         .inorder = TRUE,
               .export = c("jackknife"))  %dopar% { 
  library(sars)
  library(dplyr)
    x <- i
    k <- fig1bm[[i]]
    names(k) <- ys
    
    lapply(ys, function(y){
      df <- data.frame(datL[[x]]$area, datL[[x]][,y])
      jn <- jackknife(df, n = 3, best = as.vector(k[y]))
      jn <- as.data.frame(jn)
      jn$Resp <- y
      jn <- relocate(jn, Resp, .before = Jacknife_BM)
     # write.csv(1, file = paste0(i,"_",y,".csv"))
      jn
    })
  }#eo for each
  
  ##Convert to results table
  names(la) <- dataset
  
  la2 <- lapply(la, function(x){
    do.call(rbind.data.frame,x)
  })
  la3 <- do.call(rbind.data.frame,la2)
  la3$System <- c(rep("Grassland", 12),
                  rep("Bulrush", 12),
                  rep("Washout", 12))
  la3$Point_rem <- rep(1:3, 12)
  la3 <- relocate(la3, c(System, Point_rem), 
                  .before = Resp)
  return(la3)
}

##The best model fits, as presented in Table 2 in the
#main paper
fig1bm <- list("Grassland" = c("monod", "logistic", 
                               "monod", "power"),
               "Bulrush" = c("power", "power", 
                             "epm2", "loga"),
               "Washout" = c("logistic", "logistic", 
                             "epm1", "p1"))

#Note this takes a few hours to run
tt <- sar_jackknife(datL, n = 3, fig1bm = fig1bm)
