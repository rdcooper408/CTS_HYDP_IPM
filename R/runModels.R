# Use this script to run the demographic model

    # Note: there is an R Shiny app that for running the hybrid removal simulations  
        # URL: <https://rdcooper408.shinyapps.io/cts_ipm_shiny1/>

# setwd("CTS_HYDP_IPM/") # Set to main github directory

# Load packages
library(foreach)
library(doParallel)
source("R/demographicFunctions.R") # Load demographic functions
load("data/data.Rdata") # Load Data


# This function runs any of three models
    # model.sel = "densDepEnvtStoch" # Density-Dependent Model with Environmental Stochasticity (PVA)
    # model.sel = "densDep" # Density-Dependent Model (K)
    # model.sel = "densInd" # Density-Independent Model (Lambda)
    # model.sel = "hybridRemove" # PVA with hybrid removal strategies, see below for specific input.df

# The model will output two files: 
    # 1) a log file that is continually updated as the function runs, you can see the progress by running "grep -c 'N = ' PATH/TO/OUTPUT.log" in the terminal and compare to nrow(input.df). 
    # 2) once the full run has completed it will save the output in an RData file

# Build input file. The model will loop over each row of this file. 
# You can run the model with a large input file explores the range of all params of interest (e.g., hydp = 80:120, p.hyb = 0:1, iter = 0:500)
    # To generate a large input file that runs every combination of input params
        # Note: these models can be computationally intensive and may take considerable time to complete. 
        input.df <- expand.grid(
            years = 100, # edit to change n.years (affects compute time)
            nat.his = 0.05, # set native HIS
            hyb.his = 0.75, # set hybrid HIS
            n.meta = 1000, # set starting number of metamorphs
            n.adult = 1000, # set starting number of adults
            pond.volume = 1000, # set pond volume: Cubic meters of pond volume (for density-dependent model)
            pop.limit = 1e6, # set the maximum population size allowed in the model. To limit runaway simulations
            model.sel = "densDep", # Select model to run, options: "densInd", "densDep", "densDepEnvtStoch"
            iter = seq(1,10,1), # How many iterations to run, each iteration uses a different draw from bayesian models used to build IPM (1 to 500)
            hydp = seq(90,110, 10), # Hydroperiod: Manually set pond hydroperiod
            p.hyb = seq(0.25, 0.75, 0.25), # Proportion of Hybrids: The starting proportion of hybrids
            pop.limit = 1e5 # max population size, above which random animals are removed
        )

# Run the function
test.out <- run.IPM(input.df = input.df, param.file = cts, n.threads = 8, file.path = "modelOutput/")


# Visualize output
    nt.par.out <- test.out
    # load("modelOutput/ExampleData_densDepIter_y100_PondV1000_Iter1.5_Hydp90.90_Phyb0.1_nA500_2024-10-31_10.45.Rdata") # alternatively load "nt.par.out" from the saved output file
    
    # Plot all individual runs
    matplot(nt.par.out[,"n.all",], type = "l", main = "Population Size", ylab = "Population Size", xlab = "Years") 
    
    # Summarize the runs by group
    nt <- nt.par.out[,"n.all",]
    nt.m <- reshape2::melt(nt)
    colnames(nt.m) <- c("Year", "ModelRun", "Nt")
    nt.m$p.hyb <- sapply(strsplit(as.character(nt.m$ModelRun), split = "_"), "[", 2 ) # extract proportion hybrids from model run info (also available in input.df file)
    nt.m$hydp <- sapply(strsplit(as.character(nt.m$ModelRun), split = "_"), "[", 1 ) # extract hydroperiod from model run info (also available in input.df file)
    
    library(ggplot2)
    ggplot2::ggplot(data = nt.m, aes(x = Year, y = Nt, col = p.hyb, group = p.hyb))+geom_point(alpha = 0.1)+scale_color_brewer()+geom_smooth(linewidth = 2)+ylim(0,2000)+theme_bw()
    ggplot2::ggplot(data = nt.m, aes(x = Year, y = Nt, col = hydp, group = hydp))+geom_point(alpha = 0.1)+scale_color_brewer()+geom_smooth(linewidth = 2)+ylim(0,2000)+theme_bw()


#################################################################################
                # Hybrid Removal Simulations
##################################################################################
    
    # Build hybrid removal simulation input dataframe. Note that this input.df can also be used in the parallelized "runIPM" function if raw data is required. 
input.df <- data.frame(
    years = 100, # edit to change n.years (affects compute time)
    nat.his = 0.05, # set native HIS
    hyb.his = 0.75, # set hybrid HIS
    n.meta = 1000, # set starting number of metamorphs
    n.adult = 1000, # set starting number of adults
    pond.volume = 1000, # set pond volume
    pop.limit = 1e6, # set the maximum population size allowed in the model. To limit runaway simulations
    model.sel = "hybridRemove", # Select model to run, options: "densInd", "densDep", "densDepEnvtStoch", "hybridRemove"
    iter = 1, # default 1 iteration, can select additional iterations using the nr.iter in function
    hydp = 100, # Hydroperiod: Manually set pond hydroperiod
    p.hyb = 0.50, # Proportion of Hybrids: The starting proportion of hybrids
    # Parameters for hybridRemove only
    hybrid.detect.his = 0.2, # Minimum HIS Detection (Precision): The lowest % BTS an individual can be while still being identified as hybrid (depends on # SNPs in assay)
    prob.hybrid.capture = 0.95, # Probability of Adult Capture: what is the probability of capturing a breeding adult in a pond (based on survey effort)
    max.sal.screen = 200, # Maximum Number of Salamanders Screened: how many salamanders you can capture and assay each year
    n.years.hybrid.removal = 20, # Number of Years Hybrids are Screened and Removed: starting year 1 extending to selected year
    removeHybridOverHIS = 0.5, # HIS Threshold for Removal (Cutoff): remove all assayed hybrids with HIS greater than this value
    dup = 1)


# Run the function
    nt.df  <- run.hrIPM(input.df = input.df, param.file = cts, nr.iter = 2)



# Summary Table
mean.25y <- function(x){mean(x[75:100])}
n.adult.df <- reshape2::dcast(data = nt.df, iter ~ hybRemove, value.var = "n.adult", fun.aggregate = mean.25y)
n.adult.df$iter <- as.character(n.adult.df$iter) 
n.adult.df

# Total reduction in Hybrid Index Score (HIS)
avg.df <- reshape2::dcast(data = nt.df, year ~ hybRemove, fun.aggregate = mean, value.var = "his.adult", na.rm = T)
avg.df$diff <-avg.df$`No Removal` -  avg.df$`Yes Removal`
ggplot2::ggplot(data = avg.df, aes(x = year, y = diff)) + geom_point(col = "darkgreen") + geom_line(size = 1.2, col = "darkgreen") + theme_bw(base_size = 14) + 
    labs(title = "Reduction in Hybrid Index Score (HIS)", subtitle = "No Removal - Yes Removal", x = "Year", y = "Reduction in HIS") + 
    geom_text(mapping = aes(x = 80, y = max(avg.df$diff)*0.9), label = paste("Max Reduction HIS:", round(max(avg.df$diff), digits = 2)), size = 5)

# Percent reduction in HIS
avg.df <- reshape2::dcast(data = nt.df, year ~ hybRemove, fun.aggregate = mean, value.var = "his.adult", na.rm = T)
avg.df$perc <-((avg.df$`No Removal` -  avg.df$`Yes Removal`)/avg.df$`No Removal`)*100
ggplot2::ggplot(data = avg.df, aes(x = year, y = perc)) + geom_point(col = "darkblue") + geom_line(size = 1.2, col = "darkblue") + theme_bw(base_size = 14) + 
    labs(title = "Percent Reduction in Hybrid Index Score (HIS)", subtitle = "(No Removal - Yes Removal) / (No Removal)", x = "Year", y = "Percent Reduction HIS") + 
    geom_text(mapping = aes(x = 80, y = max(avg.df$perc)*0.9), label = paste("Max Reduction HIS: ", round(max(avg.df$perc), digits = 1), "%", sep =""), size = 5)

