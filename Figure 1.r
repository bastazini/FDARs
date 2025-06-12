
library(sars)

dat1 <- read.csv("grass_MainTable.csv")
dat2 <- read.csv("bulrush_MainTable.csv")
dat3 <- read.csv("wash_MainTable.csv")

datL <- list(dat1, dat2, dat3)

Nams <- sars::display_sars_models()

##function to fit sar_average and build a plot
#of the best model
sar_fit <- function(dat, dataset, xaxis = NULL,
                    yaxis = NULL, Nams = Nams){
  sfr <- sar_average(data = dat, 
                     grid_start = "none",
                     verb = FALSE, 
                     display = FALSE)
  s2 <- summary(sfr)
  wb <- s2$Model_table$Model[1]
  bm <- sfr$details$fits[wb][[1]]
  MF <-  paste0("sar_", wb, "()")
  MT <- Nams$Model[which(Nams$Function.name == MF)]
  if (MT == "Logistic (standard)"){
    MT <- "Logistic"
  }
  MT2 <- paste0(dataset, " - ", MT)
  if (is.null(xaxis)) xaxis <- ""
  if (is.null(yaxis)) yaxis <- ""
  par(mar=c(5,5,4,1)+.1)
  plot(bm, xlab = xaxis, ylab = yaxis,
       ModTitle = MT2,
       cex.main = 2.05,
       cex.lab = 1.9,
       cex.axis = 1.5)
}

##function to iterate across all datasets and
#metrics in correct order and run sar_fit, filling
#the plotting window

sar_plot <- function(datL, 
                     dataset = c("Grassland", "Bulrush",
                                 "Washout"),
                     xaxis = NULL,
                     yaxis = NULL,
                     Nams = Nams){
  
  ys <- c("PA_Fric_obs", "Fric_obs",
          "Disp_obs", "Eve_obs")
  
  #set plotting window
  par(mfcol = c(4,3))
  
  lapply(1:3, function(x){
    ar <- datL[[x]]$area
    k <- 1
    if (x == 1){
      yns <- c("FRic (presence-absence)",
               "FRic abundance",
               "FDiv",
               "FReg")
    } else {
      yns <- NULL
    }
    lapply(ys, function(y){
      if (k == 4){
        xns <- "Area"
      } else {
        xns <- NULL
      }
      df <- data.frame(ar, datL[[x]][,y])
      sar_fit(dat = df, dataset = dataset[x],
              xaxis = xns,
              yaxis = yns[k],
              Nams = Nams)
      k <<- k + 1
    })
  })

}

jpeg("Fig1_test.jpeg", width = 34, height = 36,
     units = "cm", res = 300)
sar_plot(datL,    dataset = c("Grassland", "Bulrush",
                              "Washout"),
         Nams = Nams)
dev.off()
