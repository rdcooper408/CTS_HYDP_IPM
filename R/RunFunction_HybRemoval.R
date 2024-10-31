# Use this script to run the demographic model
# Note: there is an R Shiny app that is more user-friendly here: <https://rdcooper408.shinyapps.io/cts_ipm_shiny1/>

setwd("/Users/rdcooper/Google Drive/Manuscripts/HydroperiodIPM/GitHub_CTS_HYDP_IPM/") # Set to main github directory

# Load packages
library(foreach)
library(doParallel)
source("R/demographicFunctions.R") # Load demographic functions
load("data/data.Rdata") # Load Data


# This function runs any of three models
    # model.sel = "densDepEnvtStoch" # Density-Dependent Model with Environmental Stochasticity (PVA)
    # model.sel = "densDep" # Density-Dependent Model (K)
    # model.sel = "densInd" # Density-Independent Model (Lambda)


# Build input file. The model will loop over each row of this file. 
# You can run the model with a large input file explores the range of all params of interest (e.g., hydp = 80:120, p.hyb = 0:1, iter = 0:500)
input.df <- data.frame(
    years = 100, # edit to change n.years (affects compute time)
    nat.his = 0.05, # set native HIS
    hyb.his = 0.75, # set hybrid HIS
    n.meta = 1000, # set starting number of metamorphs
    n.adult = 1000, # set starting number of adults
    pond.volume = 1000, # set pond volume: Cubic meters of pond volume (for density-dependent model)
    pop.limit = 1e6, # set the maximum population size allowed in the model. To limit runaway simulations
    model.sel = "densDep", # Select model to run, options: "densInd", "densDep", "densDepEnvtStoch"
    iter = 1, # default 1 iteration, can select additional iterations below
    hydp = 100, # Hydroperiod: Manually set pond hydroperiod
    p.hyb = 0.5, # Proportion of Hybrids: The starting proportion of hybrids
    pop.limit = 1e5 # max population size, above which random animals are removed
)




run.IPM <- function(input.df=input.df, 
                    param.file = cts,
                    n.threads=4, # Set maximum threads for parallel computing
                    adult.age=3, # Age of starting adults (for initialization of population)
                    file.path="modelOutput/", # path to output folder (file names are auto generated)
                    ){
    
        dir.create(file.path("modelOutput/"), showWarnings = FALSE)
  
        # # Load previous runs to determine K
        #       # Load RData files from model.sel = "densDep" to determine starting values for model.sel = "densDepEnvtStoch"
        #     if(pond.volume == 1000){load(file = paste(path_to,"ind.DenDep_Y200_PV1000M_1kA.M.Rdata", sep =""))}
        #     if(pond.volume == 10000){load(file = paste(path_to,"ind.DenDep_Y100_PV10000M_1kA.M.Rdata", sep =""))}
        #     
        #     d.dep.nt <- nt.par.out
        #     d.dep.sim <- input.df
        #     burn.in=round(nrow(d.dep.sim)*0.75, 0)
        #     d.dep.sim$k_adult <- apply(X = d.dep.nt[burn.in:nrow(d.dep.sim), "n.adult", ], MARGIN = 2, FUN = median)
        #     d.dep.sim$k_meta <- apply(X = d.dep.nt[burn.in:nrow(d.dep.sim), "n.meta", ], MARGIN = 2, FUN = median)
        #     
     
        # Generate file name
        file.name = paste(input.df$model.sel[1],"Iter_y", input.df$years[1], "_PondV",pond.volume, "_Iter",min(input.df$iter),".",max(input.df$iter), 
                          "_Hydp",min(input.df$hydp),".",max(input.df$hydp), 
                          "_Phyb", min(input.df$p.hyb), ".", max(input.df$p.hyb)  ,"_nA",
                          input.df$n.adult[1],"_",Sys.Date(),"_",format(Sys.time(),"%H.%M"), sep = "") # create file name based on input params
        log.file=paste(file.path, file.name, ".log",sep="") # output logfile
        out.file=paste(file.path, file.name, ".Rdata",sep="")# output Rdata file 
        out.file
        
        #setup parallel cluster to use multiple threads
        cl <- makeCluster(n.threads) # set max threads
        registerDoParallel(cl) # register cluster
        start.t.all <- Sys.time() # Set starting system time
        
        # Begin parallel loop
        nt.par.out <- foreach(ii=1:nrow(input.df), .combine='cbind', .multicombine=TRUE) %dopar% {
        # Use these commands if you do not want to run in parallel
            #rm(nt.par.out) # For debugging not in parallel
            #for(ii in 1:nrow(input.df)){ # For debugging not in parallel
            
            # Set initial parameters
            start.t <- Sys.time() # set iteration start time
            n.meta <- round(input.df$n.meta[ii], digits = 0) # number of metamorphs
            n.hyb.m <- round(n.meta*input.df$p.hyb[ii], digits = 0) # number of hybrid metamorphs
            n.nat.m <- round(n.meta-n.hyb.m, digits = 0) # number of native metamorphs
            n.adult <- round(input.df$n.adult[ii], digits = 0) # number of adults
            dd <- input.df$iter[ii] # select parameter
            
            #Setup Starting population
            
            # Metamorphs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            hyb.m <- recruits.ind.m(n = n.hyb.m,
                                    mean = log(logmass.bayes(hydp = input.df$hydp[ii], his = input.df$hyb.his[ii], params = cts[,dd])),
                                    params = cts[,dd]) # find hybrid metamorph mass using hydp and HIS
            
            nat.m <- recruits.ind.m(n = n.nat.m,
                                    mean = log(logmass.bayes(hydp = input.df$hydp[ii], his = input.df$nat.his[ii], params = cts[,dd])),
                                    params = cts[,dd]) # find native metamorph mass using hydp and HIS
            
            his <- c(rep(input.df$hyb.his[ii],n.hyb.m), 
                     rep(input.df$nat.his[ii], n.nat.m)) # create vector of HIS
            
            meta.df <- data.frame(ID=paste(0, 1:n.meta, sep ="."), 
                                  HYDP=rep(input.df$hydp[ii],n.meta),
                                  HIS=his, MASS=exp(c(hyb.m, nat.m))) # Create starting metamorph data.frame
            
            # Adults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            n.hyb.a <- round(n.adult*input.df$p.hyb[ii], digits = 0) # number of hybrid adults
            n.nat.a <- round(n.adult-n.hyb.a, digits = 0) # number of native adults
            
            hyb.a <- recruits.ind.m(n = n.hyb.a,
                                    mean = log(logmass.bayes(hydp = input.df$hydp[ii], his = input.df$hyb.his[ii], params = cts[,dd])), 
                                    params = cts[,dd]) # find hybrid adult mass using hydp and HIS
            
            nat.a <- recruits.ind.m(n = n.nat.a,
                                    mean = log(logmass.bayes(hydp = input.df$hydp[ii], his = input.df$nat.his[ii], params = cts[,dd])), 
                                    params = cts[,dd]) # find native adult mass using hydp and HIS
            
            his <- c(rep(input.df$hyb.his[ii],n.hyb.a), 
                     rep(input.df$nat.his[ii], n.nat.a)) # for some reason i had to round here to get the correct number of repitions
            
            adult.df <- data.frame(ID=paste(0, 1:n.adult, sep ="."), 
                                   HYDP=rep(input.df$hydp[ii],n.adult),
                                   HIS=his,
                                   MASS=exp(c(hyb.a, nat.a))) # Create starting adult data.frame
            
            # Loop to calculate adult mass starting from metam. using growth equations
            adult.df$MASS <- exp(vapply(X = log(adult.df$MASS), FUN = gxy.m.ind, FUN.VALUE = 1.0, y=y, params = cts[,1])) # find 1st year adult mass
            for(gg in 1:adult.age){
                adult.df$MASS <- exp(vapply(X = log(adult.df$MASS), FUN = gxy.a.ind, FUN.VALUE = 1.0, y=y, params = cts[,1])) # find growth each year 
            }
            
            # Run demographic model
            nt.df <- ind.fullBayesIPM(model.sel = input.df$model.sel[ii], hydp = input.df$hydp[ii], 
                                      years = input.df$years[ii], meta.df=meta.df, adult.df=adult.df, pond.volume = pond.volume, 
                                      params = cts, param.sel = input.df$iter[ii], clim = clim, pop.limit = input.df$pop.limit[ii])
            
            # Write update to log file
            cat(paste("Finished Iteration =", input.df$iter[ii], "| ii = ", ii,"| HYDP =",input.df$hydp[ii],
                      " | p.hyb =",input.df$p.hyb[ii]," | N. Adult Clipped = ", sum(nt.df$n.a.clip)," | N. Metam Clipped = ", sum(nt.df$n.m.clip), " | \n"), file = log.file, append = T)
            cat(paste("N = ", nt.df$n.adult[nrow(nt.df)], "\n"), file = log.file, append = T)
            cat(paste("Iteration Time",round(difftime(Sys.time(), start.t, units = "min"), digits = 2), "min \n"), 
                file = log.file, sep = ".", append = T)
            cat(paste("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", "\n"), 
                file = log.file, append = T)
            
            # return the model output
            return(unlist(nt.df))
            # nt.par.out <- ifelse(exists("nt.par.out"), nt.par.out <- rbind(nt.par.out, unlist(nt.df)), nt.par.out <- unlist(nt.df)) # for debugging not in parallel
        }
        
        
        # Write final update to log file
        cat(paste("Program completed in",round(difftime(Sys.time(), start.t.all, units = "min"), digits = 2), "min \n"), file = log.file, append = T)
        
        # Stop cluster
        stopCluster(cl) 
        
        # Cleanup output file
        dim(nt.par.out) <- c(input.df$years[1], 17, nrow(input.df)) # set array structure for output file
        dimnames(nt.par.out) <- list(1:input.df$years[1], c("year","n.meta","his.meta", "n.adult", "his.adult", 
                                                                "n.hyb.m","n.hyb.a","n.nat.m","n.nat.a","climate.year", "clim.OJ", "clim.DJ",
                                                                "clim.prop.fem.breed", "clim.meta.coeff", "n.a.clipped","n.m.clipped" ,"n.all"), 
                                     paste(input.df$hydp, input.df$p.hyb, input.df$iter, sep="_")) # set dimnames for output array
        
        # Save output and input files
        save(nt.par.out, input.df, file = out.file)
    }


# # Visualize output
    # matplot(nt.par.out[,"n.all",], type = "l", main = "Population Size", ylab = "Population Size", xlab = "Years")
    # nt <- nt.par.out[,"n.all",]
    # nt.m <- reshape2::melt(nt)
    # colnames(nt.m) <- c("Year", "ModelRun", "Nt")
    # library(ggplot2)
    # t_col <- function(color, percent = 50, name = NULL) {
    #     rgb.val <- col2rgb(color)
    #     t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
    #                  max = 255,
    #                  alpha = (100 - percent) * 255 / 100,
    #                  names = name)
    #     invisible(t.col)
    # }
    # table(nt.m$ModelRun)
    # nt.m$p.hyb <- sapply(strsplit(as.character(nt.m$ModelRun), split = "_"), "[", 2 )
    # table(nt.m$p.hyb)
    # ggplot(data = nt.m, aes(x = Year, y = Nt, col = p.hyb, group = p.hyb))+geom_point()+scale_color_manual(values = c(t_col("blue", 99),t_col("red", 99),t_col("green", 99) ))+geom_smooth(size = 5)+ylim(0,2000)



input.df <- data.frame(
    years = 100, # edit to change n.years (affects compute time)
    nat.his = 0.05, # set native HIS
    hyb.his = 0.75, # set hybrid HIS
    n.meta = 1000, # set starting number of metamophs
    n.adult = 1000, # set starting number of adults
    pond.volume = 1000, # set pond volume
    pop.limit = 1e6, # set the maximum population size allowed in the model. To limit runaway simulations
    model.sel = "hybridRemove", # Select model to run, options: "densInd", "densDep", "densDepEnvtStoch", "hybridRemove"
    iter = 1, # default 1 iteration, can select additional iterations below
    hydp = 100, # Hydroperiod: Manually set pond hydroperiod
    p.hyb = 0.5, # Proportion of Hybrids: The starting proportion of hybrids
    # Parameters for hybridRemove only
    hybrid.detect.his = 0.2, # Minimum HIS Detection (Precision): The lowest % BTS an individual can be while still being identified as hybrid (depends on # SNPs in assay)
    prob.hybrid.capture = 0.95, # Probability of Adult Capture: what is the probability of capturing a breeding adult in a pond (based on survey effort)
    max.sal.screen = 200, # Maximum Number of Salamanders Screened: how many salamanders you can capture and assay each year
    n.years.hybrid.removal = 20, # Number of Years Hybrids are Screened and Removed: starting year 1 extending to selected year
    removeHybridOverHIS = 0.5, # HIS Threshold for Removal (Cutoff): remove all assayed hybrids with HIS greater than this value
    dup = 1)


run.hrIPM <- function(input.df=input.df, param.file = cts, nr.iter = 5){
    # Increment the progress bar, and update the detail text.
    # incProgress(0, detail = paste("Initializing Simulation"))
    # prog = 1/((nr.iter*2)+2)
    adult.age = 3
    ii <- 1
    pond.volume=input.df$pond.volume[ii]
    cts.raw <- cts
    if(input.df$model.sel[ii] == "hybridRemove"){cts=as.matrix(rowMeans(cts.raw))} # for hybrid removal, use average param values
    
    
    # Set initial parameters
    #start.t <- Sys.time() # set iteration start time
    n.meta <- max(round(input.df$n.meta[ii], digits = 0), 1) # number of metamorphs; If no meta at K, then add 1 so model doesnt break
    n.hyb.m <- round(n.meta*input.df$p.hyb[ii], digits = 0) # number of hybrid metamorphs
    n.nat.m <- round(n.meta-n.hyb.m, digits = 0) # number of native metamorphs 
    n.adult <- max(round(input.df$n.adult[ii], digits = 0), 1) # number of adults; If no adults at K, then add 1 so model doesnt break
    dd <- input.df$iter[ii] # select parameter
    
    #Setup Starting population
    
    
    # Metamorphs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bb =  b.fun(a = 10, mu = input.df$hyb.his[ii]) # calculate beta.binom shape parameter Beta from mean hybrid HIS
    hyb.m.rbeta.his <- rbeta(n = n.hyb.m, shape1=10, shape2 = bb) # Draw a distribtuion of hybrid HIS values
    # hist(hyb.m.rbeta.his)
    meta.df <- data.frame(ID=paste(0, 1:n.meta, sep ="."), 
                          HYDP=rep(input.df$hydp[ii],n.meta),
                          HIS=c(hyb.m.rbeta.his, rep(input.df$nat.his[ii], n.nat.m)), 
                          MASS=NA)
    
    meta.df$MASS <- exp(log(logmass.bayes(hydp = meta.df$HYDP, his = meta.df$HIS, params = cts[,dd])) + rnorm(n = nrow(meta.df),mean=0,sd=cts[19,dd])) # determine meta. sizes based on hydp equation and adding variance from Searcy et al
    #plot(MASS ~ HIS, data = meta.df)
    
    
    # Adults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    n.hyb.a <- round(n.adult*input.df$p.hyb[ii], digits = 0) # number of hybrid adults
    n.nat.a <- round(n.adult-n.hyb.a, digits = 0) # number of native adults
    
    hyb.a.rbeta.his <- rbeta(n = n.hyb.a, shape1=10, shape2 = bb) # Draw a distribtuion of hybrid HIS values
    
    adult.df <- data.frame(ID=paste(0, 1:n.adult, sep ="."), 
                           HYDP=rep(input.df$hydp[ii],n.adult),
                           HIS=c(hyb.a.rbeta.his, rep(input.df$nat.his[ii], n.nat.a)), 
                           MASS=NA)
    
    adult.df$MASS <- exp(log(logmass.bayes(hydp = adult.df$HYDP, his = adult.df$HIS, params = cts[,dd])) + rnorm(n = nrow(adult.df),mean=0,sd=cts[19,dd])) # determine meta. sizes based on hydp equation and adding variance from Searcy et al
    # plot(MASS ~ HIS, data = adult.df)
    
    
    # Loop to calculate adult mass starting from metam. using growth equations
    adult.df$MASS <- exp(vapply(X = log(adult.df$MASS), FUN = gxy.m.ind, FUN.VALUE = 1.0, y=y, params = cts[,1])) # find 1st year adult mass
    for(gg in 1:adult.age){
        adult.df$MASS <- exp(vapply(X = log(adult.df$MASS), FUN = gxy.a.ind, FUN.VALUE = 1.0, y=y, params = cts[,1])) # find growth each year 
    }
    
    # Increment the progress bar, and update the detail text.
    # incProgress(prog, detail = paste("Running non-removal baseline demography"))
    
    for(kk in 1:nr.iter){
        # incProgress(prog, detail = paste("Running non-removal baseline demography | iter =", kk))
    nt.df.tmp <- ind.fullBayesIPM(model.sel = input.df$model.sel[ii], hydp = input.df$hydp[ii], 
                                 years = input.df$years[ii], meta.df=meta.df, adult.df=adult.df, pond.volume = input.df$pond.volume[ii], 
                                 params = cts, param.sel = input.df$iter[ii], clim = clim, pop.limit = input.df$pop.limit[ii], 
                                 hybrid.detect.his = input.df$hybrid.detect.his[ii], prob.hybrid.capture = input.df$prob.hybrid.capture[ii], 
                                 max.sal.screen = input.df$max.sal.screen[ii], removeHybridOverHIS = input.df$removeHybridOverHIS[ii], 
                                 n.years.hybrid.removal = 0)
    nt.df.tmp <- data.frame(nt.df.tmp, hybRemove = "No Removal", iter = kk)
    if(kk<2){nt.df <- nt.df.tmp} else nt.df <- rbind(nt.df, nt.df.tmp)
    }
    
    # Increment the progress bar, and update the detail text.
    # incProgress(prog, detail = paste("Running hybrid removal"))
    
    # Run demographic model
    for(kk in 1:nr.iter){
        # incProgress(prog, detail = paste("Running hybrid removal | iter =", kk))
    nt.df.tmp <- ind.fullBayesIPM(model.sel = input.df$model.sel[ii], hydp = input.df$hydp[ii], 
                              years = input.df$years[ii], meta.df=meta.df, adult.df=adult.df, pond.volume = input.df$pond.volume[ii], 
                              params = cts, param.sel = input.df$iter[ii], clim = clim, pop.limit = input.df$pop.limit[ii], 
                              hybrid.detect.his = input.df$hybrid.detect.his[ii], prob.hybrid.capture = input.df$prob.hybrid.capture[ii], 
                              max.sal.screen = input.df$max.sal.screen[ii], removeHybridOverHIS = input.df$removeHybridOverHIS[ii], 
                              n.years.hybrid.removal = input.df$n.years.hybrid.removal[ii])
    nt.df.tmp <- data.frame(nt.df.tmp, hybRemove = "Yes Removal", iter = kk)
    nt.df <- rbind(nt.df, nt.df.tmp)
    }
    # Increment the progress bar, and update the detail text.
    # incProgress(prog, detail = paste("Building output data"))

    return(nt.df)

}


# test.out <- run.hrIPM(input.df = input.df, param.file = cts)



