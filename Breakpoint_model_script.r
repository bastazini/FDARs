
library(sars)
library(dplyr)

####################################################
#########Fit Breakpoint / piecewise models to############
#########the SES-area relationship data################
#####################################################

#Parameters to use within sar_threshold
intz_set = 0.001 #interval arg
Ncore <- 4 #number of cores in parallel processes
LOG <- log10 #log function to use

##function to run sar_threshold and return the model fits
#Only fits the  "ContOne" and "ZslopeOne" models
fit_thr_int2 <- function(datz, logAxesz = "none",
                        logBase,
                        conz = 0.1,  nizl = NULL,
                        parallelz = FALSE,
                        coresz = NULL,
                        intz = NULL){
  
  if (ncol(datz) != 2) stop("incorrect col N")
  
  modz <- c("ContOne", "ZslopeOne")
  
  s <- sar_threshold(data = datz,
                     mod = modz,
                     interval = intz,
                     logAxes = logAxesz,
                     con = conz,
                     nisl = nizl,
                     logT = logBase,
                     parallel = parallelz,
                     cores = coresz)
  return(s)
}

##Main function, runs the above internal function and then 
#formats the output to match the results table in the main ms.
#Orders model fits by AICc
fit_thr2 <- function(datz,
                    logBase,
                    conz = 0.1,
                    nizl = NULL,
                    parallelz = FALSE,
                    coresz = NULL,
                    intz = NULL){

  m1 <- fit_thr_int2(datz, logAxesz = "area",
                    logBase = logBase,
                    conz = conz,
                    nizl = nizl,
                    parallelz = parallelz,
                    coresz = coresz,
                    intz = intz)
  
  sm1 <- summary(m1)
  
  rsm1 <- sm1[[2]] %>%
    select(AICc, R2, Th1, Th2, seg1, seg2, seg3)
  
  rsm1 <- rsm1[order(rsm1$AICc),]
  
  snow <- list(list(rsm1), m1)
   return(snow)

}

##Read in the three results tables with SES values (see Main_
##hypervolume_script) as a list called ldf
lf <- c("bulrush_MainTable.csv", "grass_MainTable.csv",
        "wash_MainTable.csv")
ldf <- lapply(lf, read.csv)

metr <- c("Fric_SES",
          "Disp_SES",
          "Eve_SES")

##Iterate across each of the results tables, 
##running the fit_thr2 function. This fits 
##the models for each of Fric, Disp and eve
##separately
th_res_all <- lapply(ldf, function(y){
  
da <- y

th_res_int <- lapply(metr, function(x){

datAll <- select(da, area, 
                 "native_count" = all_of(x))

ResAll <- fit_thr2(datz = datAll,
                  logBase = LOG,
                  parallelz = TRUE,
                  coresz = Ncore,
                  nizl = 5,
                  intz = intz_set)

ResAll

})#eo lapply int

names(th_res_int) <- metr
th_res_int

})#eo main lapply


names(th_res_all) <- sub("\\_.*", "", lf)

All_res <- lapply(th_res_all, 
       function(x) lapply(x,function(y) y[[1]][[1]]))

All_res
