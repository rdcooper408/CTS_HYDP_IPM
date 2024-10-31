


# Functions taken directly from Messerman et al 2022, 2023
    # Messerman A, Clause A, Gray L, Krkošek M, Rollins H, Trenham P, Shaffer B, Searcy C. 2022, October 21. Data: Applying stochastic and Bayesian integral projection modeling to amphibian population viability analysis. Dryad. Available from https://datadryad.org/stash/dataset/doi:10.5061/dryad.59zw3r291 (accessed February 6, 2024).
    # Messerman AF, Clause AG, Gray LN, Krkošek M, Rollins HB, Trenham PC, Shaffer HB, Searcy CA. 2023. Applying stochastic and Bayesian integral projection modeling to amphibian population viability analysis. Ecological Applications 33:e2783.



##To calculate inverse logit functions
invlogit<-function(x){exp(x)/(1+exp(x))}

###########################################################################
##Vital rate functions. Parameter indices are hard-coded infunctions 
## and must correspond to rows of MCMC matrix as defined in IPM script
###########################################################################
##Metamorph growth from size x to size y
gxy.m<-function(x,y,params){
    xb=pmin(pmax(x,params[16]),params[17]) #Transforms all values below/above limits in min/max size
    xg=params[1] + params[2]*(xb-2.39) #Creates one mean growth value per body size bin per iteration of posterior draws, centered
    return(dnorm(y,mean=xg, sd=params[3]))
}

##Adult growth from size x to size y
gxy.a<-function(x,y,params){
    xb=pmin(pmax(x,params[16]),params[17]) 
    xg=params[4] + params[5]*(xb-2.95) 
    return(dnorm(y,mean=xg,sd=params[6]))
}

##Metamorph survival at size x
sx.m<-function(x,params){
    xb=pmin(pmax(x,params[16]),params[17])
    smeta<-(invlogit(params[7] + params[8]*(xb-2.28))/params[10])
    smeta<-ifelse(smeta>1, 1, smeta)
    return(smeta)
}

##Adult/juvenile survival at size x
sx.a<-function(x,params){
    xb=pmin(pmax(x,params[16]),params[17])
    sadj<-(invlogit(params[9] + params[8]*(xb-2.49))/params[11])
    sadj<-ifelse(sadj>1, 1, sadj)
    return(sadj)
}

#Metamorph survival*growth
pxy.m<-function(x,y,params){
    xb=pmin(pmax(x,params[16]),params[17])
    sx.m(xb,params)*gxy.m(xb,y,params)
}

#Adult survival*growth
pxy.a<-function(x,y,params){
    xb=pmin(pmax(x,params[16]),params[17])
    sx.a(xb,params)*gxy.a(xb,y,params)
}

#Fecundity (fx) is the product of maturity (mx), fertility (fer), percent females breeding (34%), 
# percent population that is female (50%), and maximum possible embryonic/larval survival
# representing the number of juveniles produced by a female of a size class (y)
mx<-function(y,params){
    yb=pmin(pmax(y,params[16]),params[17])
    invlogit(params[12] + params[13]*(yb-2.45))
}

fer<-function(y,params){
    yb=pmin(pmax(y,params[16]),params[17])
    return(params[14] + params[15]*(yb-3.48)) 
}

fx<-function(y,params){
    yb=pmin(pmax(y,params[16]),params[17])
    nfem<-params[20]
    nfem.b<-params[21]
    sx.l<-params[22]
    return(mx(yb, params)*fer(yb, params)*nfem*nfem.b*sx.l)
}

#Size distribution of new metamorphs under low densities
recruits<-function(y,params){
    yb=pmin(pmax(y,params[16]),params[17])
    dnorm(x=yb,mean=params[18],sd=params[19])
}


###Fecundity and body size at metamorphosis under density dependence
#Fecundity without embryonic/larval survival accounted for
fx1<-function(y,params){
    yb=pmin(pmax(y,params[16]),params[17])
    nfem<-params[20]
    nfem.b<-params[21]
    return(mx(yb, params)*fer(yb, params)*nfem*nfem.b)
}

#Function for metamorph body size given egg density
met<-function(d,params){
    return(params[24] + params[23]*(d-1.23)) #Mean-centered on Searcy et al 2015 data
}

#Function for larval survival to metamorphosis given egg density
lar<-function(d,params){
    return(params[25] + params[26]*(d-3.38)) #Mean-centered on observed data
}

#Metamorph survival*growth applied to current population size distribution (m)
pxy.m1<-function(x,y,m,params){
    xb=pmin(pmax(x,params[16]),params[17])
    (sx.m(xb,params)*m)*gxy.m(xb,y,params)
}

#Adult survival*growth applied to current population size distribution (a)
pxy.a1<-function(x,y,a,params){
    xb=pmin(pmax(x,params[16]),params[17])
    (sx.a(xb,params)*a)*gxy.a(xb,y,params)
}
#Size distribution of new metamorphs given egg densities (d)
recruits2<-function(y,d,params){ 
    yb=pmin(pmax(y,params[16]),params[17])
    return(dnorm(x=yb, mean = (params[24] + params[23]*(d-1.23)), sd=0.227))#0.227 set variance, mean-centered from Searcy et al 2015 data
}

###Fecundity under environmental and density dependence
#Fecundity without embryonic/larval survival accounted for
#Where percent females breeding is now predicted by Dec-Jan rainfall
nfem.b1<-function(e,params){
    return(min(exp(params[27] + params[28]*(e-221.09)), 1)) 
    #Where 221.09 mm is mean Dec-Jan precip from the study, and output is back-transformed from log and constrained to a maximum probability of 1
}

fx2<-function(y,e,params){
    yb=pmin(pmax(y,params[16]),params[17])
    nfem<-params[20]
    return(mx(yb, params)*fer(yb, params)*nfem*nfem.b1(e, params))
}




# Functions modified by Cooper et al. (submitted)
# "Building genotype-specific demographic models to guide management of invasive hybrids"

# Color transparency function for plotting
t_col <- function(color, percent = 50, name = NULL) {
    rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    invisible(t.col)
}

# Larval Survival
# Density-Independent larval survival
            hydp.surv.bayes <- function(hydp, his, params,
                                        m.HYDP = 98.60286, # mean HYDP for centering (dont change)
                                        m.HIS = 0.6956714){ # mean HIS for centering (dont change)

                # Params
                bHYDP = params[33] # mean = 0.07877391, # slope: effect of HYDP on survival
                bHIS = params[34] #  mean = 0.67888270,  # slope: effect of HIS on survival
                b0 = params[35] #  mean = -2.68017337,  # Intercept: survival at mean HYDP and mean HIS

                HYDP_C <- hydp - m.HYDP # center HYDP
                HIS_C <- his - m.HIS # Center HIS
                log.s <- b0 + (bHYDP * HYDP_C) + (bHIS * HIS_C) # calculate the logit of survival
                prob.s <- plogis(log.s) # convery logit to probability
                return(prob.s)
            }


            hydp.surv.eggdens.bayes <- function(hydp, his, params, egg.d,
                                                      mean.egg.d = 0.8967111){ # estimated from Messerman et al. file# mean HIS for centering (dont change))

                      aa = log(hydp.surv.bayes(hydp=hydp, his=his, params = params)) - l.egg(mean.egg.d) # This takes the differnce my calculated survival
                      surv = (params[25] + params[26]*(egg.d-3.38))+aa
                      return(surv)
                  }


 # Metamorph Mass
            
    # Density-Independent metamorph mass function
        # bayesian function which uses the log(mass) ~ hydp^2 + hydp + HIS function
        ## This returns the exp(log(mass)) = mass
        logmass.bayes <- function(hydp, his, params,
                                  m.HYDP2 = 11223.33, # mean HYDP^2 for centering (dont change)
                                  m.HYDP = 105.5263, # mean HYDP for centering (dont change)
                                  m.HIS =0.7506176){ # mean HIS for centering (dont change))
            # params
            bHYDP2 =  params[29] # mean = -0.001114589, # slope: effect of HYDP^2 on survival
            bHYDP = params[30] # mean = 0.220949884, # slope: effect of HYDP on survival
            bHIS = params[31] # mean =  0.651188858,  # slope: effect of HIS on survival
            b0 = params[32] # mean =  2.262043081,  # Intercept: survival at mean HYDP and mean HIS


            # Center the predictors
            hydp2_C = (hydp^2) - m.HYDP2
            hydp_C = hydp - m.HYDP
            his_C = his - m.HIS

            # calculate mass in grams
            mass = exp((b0) + (bHYDP2*hydp2_C) + (bHYDP*hydp_C) + (bHIS*his_C))
            return(mass)
        }

            


    # Density-Dependent Metamorph mass function
        logmass.eggdens.bayes <- function(hydp, his, params, egg.d, # Egg density (calculated as: log((sum(fem.df$fer)/pond.volume)+1e-20))
                                          m.HYDP2 = 11223.33, # mean HYDP^2 for centering (dont change)
                                          m.HYDP = 105.5263, # mean HYDP for centering (dont change)
                                          m.HIS =0.7506176){ # mean HIS for centering (dont change))

            # calculate mass in grams
            mass = logmass.bayes(hydp = hydp, his=his, params = params) # from dens-ind function
            d.mass = exp(log(mass) + (params[23]*(egg.d-1.23)))
            return(d.mass)
        }

            
                    
        # Determine beta binomial shape parameter Beta, from Mu and alpha; Used in generating distribution of HIS values to start sim.
        b.fun <- function(a, mu){(a*(1-mu))/mu}

# Other modified functions from Messerman / Serarcy IPM
        
        # l.egg function
        m.egg <- function(x){2.28 + (-0.272*(x-1.23))}
        l.egg <- function(d){-4.578732 + -0.6335772*(d-3.38)}

        
        #Function for metamorph body size given egg density
        met<-function(d,params){
            return(params[24] + params[23]*(d-1.23)) #Mean-centered on Searcy et al 2015 data
        }
        
        #Function for larval survival to metamorphosis given egg density
        lar<-function(d,params){
            return(params[25] + params[26]*(d-3.38)) #Mean-centered on observed data
        }
        
        
        # rep HIS (hybrid index score) for mapply in model
        rep.his <- function(his, n.rep){rep(his, n.rep)}
        
        # Modified IPM functions: draw single value from density dist.
        ## Metamorph Growth
        gxy.m.ind<-function(x,y,params){
            xb=pmin(pmax(x,params[16]),params[17]) #Transforms all values below/above limits in min/max size
            xg=params[1] + params[2]*(xb-2.39) #Creates one mean growth value per body size bin per iteration of posterior draws, centered
            return(rnorm(1,mean=xg, sd=params[3]))
        }
        
        ## Adult Growth from size x to size y
        gxy.a.ind<-function(x,y,params){
            xb=pmin(pmax(x,params[16]),params[17]) 
            xg=params[4] + params[5]*(xb-2.95) 
            return(rnorm(1,mean=xg,sd=params[6]))
        }
        
        #Size distribution of new metamorphs under low densities
        recruits.ind<-function(n, params){
            #yb=pmin(pmax(y,params[16]),params[17])
            rnorm(n = n,mean=params[18],sd=params[19])
        }
        
        recruits.ind.m<-function(n,mean, params){
            #yb=pmin(pmax(y,params[16]),params[17])
            rnorm(n = n,mean=mean,sd=params[19])
        }
        
        # metamorph survival coefficient (what fraction of metam. survive based on Oct-Jun rainfall)
        clim.meta.coeff <- function(clim.OJ, params){
            rep.infl.none = params[36] # mean = 404.5 # changed to match new param file after HydpGen pub; Changed from 29 -> 36
            rep.infl.full = params[37] # mean = 674.5 # changed to match new param file after HydpGen pub; Changed from 30 -> 37
            if(clim.OJ < rep.infl.none){meta.coeff <- 0} # if rainfall less than lower inflection point, 0 metam survive
            if(clim.OJ > rep.infl.full){meta.coeff <- 1} # if rainfall greater than upper inflection point, all metam survive
            if(clim.OJ >= rep.infl.none & clim.OJ <= rep.infl.full){meta.coeff <- ((clim.OJ-rep.infl.none)/(rep.infl.full-rep.infl.none))} # between infl points, follows linear model
            return(meta.coeff) # return coefficient
        }
        
        
        
        ind.fullBayesIPM <- function(hydp, years, meta.df, adult.df, model.sel = "densDep", 
                                     pond.volume = 1000, params = cts, param.sel = 1 , clim = clim, pop.limit = 1e5, hybrid.detect.his = 0.10,
                                     prob.hybrid.capture = 0.90, max.sal.screen = 100,  n.years.hybrid.removal = 10, removeHybridOverHIS = 0.3){
            
            # Pond volume is used to determine egg density
            # Olcott pond volume = 101000 m^3
            # This has a huge K, so lets use volume = 1000m^3  which is ~ 100x smaller
            # Model Select
            # This determines which model to use. Options = c("densInd", "densDep", "densDepEnvtStoch")
            
            dd <- param.sel # This selects which column in the parameter file to use
            
            # Build input metamorph df
            ind.m <- meta.df
            ind.m$ln.mass <- log(ind.m$MASS)
            ind.m$SEX <- sample(c("M", "F"), size = nrow(ind.m), replace = T)
            
            # Build input adult df
            ind.a <- adult.df
            ind.a$ln.mass <- log(ind.a$MASS)
            ind.a$SEX <- sample(c("M", "F"), size = nrow(ind.a), replace = T)
            
            # Set up df variables
            ind.a$HYDP <- ind.m$HYDP <- hydp # set hydroperiod
            ind.a$p.surv <- ind.m$p.surv <- NA # set prob. of survival
            ind.a$surv <- ind.m$surv <- NA # set survival (binary: 0 = dead, 1 = survive)
            ind.a$grow <- ind.m$grow <- NA # set new mass (after growth)
            ind.a$fec <- ind.m$fec <- NA # set fecundity
            ind.a$mat <- ind.m$mat <- NA # set maturity
            ind.a$fer <- ind.m$fer <- NA # set fertility
            ind.a$breed <- ind.m$breed <- NA # set breeding (binary: 0 = not breeding, 1 = breeding this year)
            ind.a$age <- 3 # set starting adult age
            ind.m$age <- 1 # set starting metam age
            
            # keep track of number of times population is "clipped" (limited)
            n.a.clips <- n.m.clips <- 0
            
            # Select climate-year
            clim.year <- sample(x = 1:nrow(clim), size = years, replace = T) # randomly select climate data for each sim year
            
            # Build output df (will be updated each model year)
            nt.df <- data.frame(year=1:years, n.meta=NA, his.meta=NA, n.adult=NA, his.adult=NA, n.hyb.m=NA, n.hyb.a=NA, n.nat.m=NA, n.nat.a=NA,
                                climate.year= clim$Year[clim.year], clim.OJ=clim$Oct.Jun[clim.year], clim.DJ=clim$Dec.Jan[clim.year], clim.prop.fem.breed=NA, clim.meta.coeff=NA, n.a.clipped=0, n.m.clipped=0)
            
            # Estimate influence of climate on demog.
            nt.df$clim.prop.fem.breed=vapply(X = nt.df$clim.DJ, FUN = nfem.b1, FUN.VALUE = 1.0, params=params[,dd]) # estimate prop females breeding based on climate
            nt.df$clim.meta.coeff= vapply(X = nt.df$clim.OJ, FUN = clim.meta.coeff, FUN.VALUE = 1.0, params=params[,dd] ) # estimate prop of metamorphs surviving based on climate
            
            # Begin model iteration
            
            
            for (yy in 1:years){
                
                # check to make sure there are individuals left
                if(is.null(ind.m)&nrow(ind.a)<1){break} # Check if pop is extinct (if metamorph dataframe is null and if there are no more adults)
                
                # Update output df (calculate values if there are metamorphs / adults)
                nt.df$n.meta[yy]<-ifelse(!is.null(ind.m),nrow(ind.m),0) 
                nt.df$his.meta[yy] <- ifelse(!is.null(ind.m),mean(ind.m$HIS),NA)
                nt.df$n.adult[yy]<-ifelse(!is.null(ind.a),nrow(ind.a),0)
                nt.df$his.adult[yy] <- ifelse(!is.null(ind.a),mean(ind.a$HIS),NA)
                nt.df$n.hyb.m[yy] <- ifelse(!is.null(ind.m),sum(ind.m$HIS>0.05),NA)
                nt.df$n.hyb.a[yy] <- ifelse(!is.null(ind.a),sum(ind.a$HIS>0.05),NA)
                nt.df$n.nat.m[yy] <- ifelse(!is.null(ind.m),sum(ind.m$HIS<=0.05),NA)
                nt.df$n.nat.a[yy] <- ifelse(!is.null(ind.a),sum(ind.a$HIS<=0.05),NA)
                
                
                # Begin individual demography
                
                # METAMORPHS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if(!is.null(ind.m)){ # check if there are metamorphs at start of year
                    for(ii in 1:nrow(ind.m)){ # Loop through each metamorph
                        
                        # metamorph survival
                        ind.m$p.surv[ii] <- sx.m(x = ind.m$ln.mass[ii], params = params[,dd]) # calc prob metam survived
                        ind.m$surv[ii] <- rbinom(n = 1, size = 1, prob = ind.m$p.surv[ii]) # draw binary survival (0=died, 1=lived)
                        
                        # metamorph growth
                        ind.m$grow[ii] <- gxy.m.ind(x = ind.m$ln.mass[ii], y=y, params = params[,dd]) # calc new metam ln(mass)
                        
                        # metamorph maturity
                        if(ind.m$mat[ii]==0 | is.na(ind.m$mat[ii])){ # Check if metam is immature
                            ind.m$mat[ii] <- rbinom(n = 1, size = 1, prob = mx(ind.m$grow[ii], params = params[,dd])) # Use the metam new mass to determine if it matures
                        } 
                        
                        # metamorph breeding
                        ind.m$fer[ii] <- fer(ind.m$grow[ii], params = params[,dd]) # what is metam fertility
                        if(ind.m$fer[ii] < 0){ind.m$fer[ii] <- 0} # set fertility to 0 if < 0
                        ind.m$breed[ii] <- rbinom(n = 1, size = 1, prob = params[21,dd]) # did metam breed this year
                    }#ii
                    
                    # Create next years adults (metam -> adults)
                    if(sum(ind.m$surv)>0){ # if metam survive
                        ind.m2a <- ind.m[which(ind.m$surv == 1), ] # create metam to adult df
                        ind.m2a$MASS <- exp(ind.m2a$grow) # transform to real mass
                        ind.m2a$ln.mass <- ind.m2a$grow # record ln.mass
                        ind.m2a$p.surv <- ind.m2a$surv <- ind.m2a$grow <- NA # reset other variables to NA
                    }
                    if(sum(ind.m$surv)<1){ind.m2a <- ind.m[0,]} # if no metam survive, make ind.m2a null
                }# end metam loop
                
                
                # Adults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if(nrow(ind.a)>0){ # check if there are adults
                    for(aa in 1:nrow(ind.a)){ # loop over adults
                        
                        # adult survival
                        ind.a$p.surv[aa] <- sx.a(x = ind.a$ln.mass[aa], params = params[,dd]) # calc prob. of adult survival
                        ind.a$surv[aa] <- rbinom(n = 1, size = 1, prob = ind.a$p.surv[aa]) # draw binary survival (0=died, 1=lived)
                        if(ind.a$age[aa] > 15){ind.a$surv[aa] <- 0} # add a maximum age of 15 years, kill when age > 15 years
                        
                        # adult growth
                        ind.a$grow[aa] <- gxy.a.ind(x = ind.a$ln.mass[aa], y=y, params = params[,dd]) # ln(mass) next year
                        
                        # adult maturity
                        if(ind.a$mat[aa]==0 | is.na(ind.a$mat[aa])){ # check if adult is immature
                            ind.a$mat[aa] <- rbinom(n = 1, size = 1, prob = mx(ind.a$grow[aa], params = params[,dd])) # does adult mature?
                        } 
                        
                        # adult fertility
                        ind.a$fer[aa] <- fer(ind.a$grow[aa], params = params[,dd]) # calc fertility based on NEW size
                        if(ind.a$fer[aa] < 0){ind.a$fer[aa] <- 0} # set fertility to 0 if < 0
                        if(model.sel != "densDepEnvtStoch"){ind.a$breed[aa] <- rbinom(n = 1, size = 1, prob = params[21,dd])} # did adult breed this year?
                        if(model.sel == "densDepEnvtStoch"){ind.a$breed[aa] <- rbinom(n = 1, size = 1, prob = nt.df$clim.prop.fem.breed[yy])} # did adult breed this year based on climate data?
                        
                    }#aa
                    
                    # Create next years adults (adults -> adults)
                    if(sum(ind.a$surv)>0){ # if  adults survive
                        ind.a2a <- ind.a[which(ind.a$surv == 1), ] # create adult to adult df
                        ind.a2a$MASS <- exp(ind.a2a$grow) # transform to real mass
                        ind.a2a$ln.mass <- ind.a2a$grow # record ln.mass
                        ind.a2a$p.surv <- ind.a2a$surv <- ind.a2a$grow <- NA # reset other variables to NA
                    }
                    if(sum(ind.a$surv)<1){ind.a2a <- ind.a[0,]} # if no adults survive, make ind.a2a null
                    
                } # end adult loop
                
                
                # Set up matrices for next year
                
                # Set up breeding crosses
                # Select breeders
                # fem.df <- rbind(ind.m2a[ind.m2a$mat>0 & ind.m2a$breed>0 & ind.m2a$SEX == "F", ], 
                #                 ind.a2a[ind.a2a$mat>0 & ind.a2a$breed>0 & ind.a2a$SEX == "F", ]) # combine breeding females from metam and adults
                # male.df <- rbind(ind.m2a[ind.m2a$mat>0 & ind.m2a$breed>0 & ind.m2a$SEX == "M", ], 
                #                  ind.a2a[ind.a2a$mat>0 & ind.a2a$breed>0 & ind.a2a$SEX == "M", ]) # combine breeding males from metam and adults
                
                breed.df <- rbind(ind.m2a[ind.m2a$mat>0 & ind.m2a$breed>0, ], 
                                  ind.a2a[ind.a2a$mat>0 & ind.a2a$breed>0, ]) # combine breeding metam and adults
                
                #### Remove Hybrids ######
                
                # print(paste("HIS before: ", round(mean(breed.df$HIS), 3)))
                if(nrow(breed.df)>0){breed.df$remove <- 0}
                if(model.sel == "hybridRemove" & yy < n.years.hybrid.removal){
                    if(nrow(breed.df)>0){
                        breed.df$capture <- NA
                        breed.df$capture[sample(nrow(breed.df), min(nrow(breed.df)*prob.hybrid.capture, max.sal.screen), replace = F)] <- 1 # which individuals are captured to test HIS, only screen the max num. allowed 
                        # index <- which(breed.df$HIS > hybrid.detect.his & breed.df$capture == 1) # of the individuals captured, which pass HIS detection threshold
                        prop.g <- breed.df$HIS/runif(n = nrow(breed.df), min = 1, max = 2) # what is the average zygosity of loci (1= all heterozygous, 2 = all homozygous)
                        #index <- which(prop.g > hybrid.detect.his & breed.df$capture == 1) # of the individuals captured, which pass HIS proportion of genome detection threshold
                        index <- which(prop.g > hybrid.detect.his & breed.df$capture == 1 & breed.df$HIS > removeHybridOverHIS) # of the individuals captured, which pass HIS proportion of genome detection threshold and have HIS above removal threshold
                        breed.df$remove[sample(index, size = min(max.sal.screen, length(index)), replace = F)] <- 1 # remove random eligible hybrids
                    }
                }
                #print(paste("HIS after: ", round(mean(breed.df$HIS[breed.df$remove == 0]), 3)))
                
                fem.df <- breed.df[which(breed.df$SEX == "F" & breed.df$remove == 0), ]
                male.df <- breed.df[which(breed.df$SEX == "M" & breed.df$remove == 0), ]

                # make offspring
                if(nrow(fem.df)>0 & nrow(male.df)>0){ # if there are males and females
                    fem.df$mate.his <- male.df$HIS[sample(1:nrow(male.df), size = nrow(fem.df), replace = T)] # Randomly pick male mates
                    fem.df$off.his <- apply(fem.df[,c("HIS", "mate.his")], MARGIN = 1, FUN = mean) # Take the average of the HIS (could make more complex)
                    
                    # Larval Survival: density-independent
                    if(model.sel == "densInd"){
                        fem.df$larv.surv <- hydp.surv.bayes(hydp = fem.df$HYDP, his = fem.df$off.his, params = params[,dd] ) # calc number of larv that survive based on hydp exp
                    }
                    
                    # Larval Survival: density-dependent
                    if(model.sel != "densInd"){
                        egg.d <- log((sum(fem.df$fer)/pond.volume)+1e-20) # calculate egg density from female fertility
                        fem.df$larv.surv <- pmin(exp(hydp.surv.eggdens.bayes(hydp = fem.df$HYDP, 
                                                                             his = fem.df$off.his, 
                                                                             params = params[,dd], 
                                                                             egg.d = egg.d)), 1) # Density Dependent larval survival function using hydp and HIS
                    }
                    
                    if(model.sel != "densDepEnvtStoch"){fem.df$fec <- fem.df$fer*fem.df$larv.surv} # multiply survival prob. and fertility to get num of larvae that survive
                    if(model.sel == "densDepEnvtStoch"){fem.df$fec <- fem.df$fer*fem.df$larv.surv*nt.df$clim.meta.coeff[yy]} # multiply survival prob., fertility, and climate mortality to get num of larvae that survive
                }
                
                # create new metam
                n.off <- sum(round(fem.df$fec, digits = 0)) # round the number of new larvae
                if(is.null(n.off) | is.na(n.off)){n.off <- 0} # set to zero if no offspring
                if(n.off > 0){ # if there are offspring
                    ind.new.m <- data.frame(ID=rep(NA, n.off),HYDP=NA,HIS=NA,MASS=NA,ln.mass=NA,SEX=NA, 
                                            p.surv=NA, surv=NA, grow=NA, fec=NA, mat=NA, fer=NA, breed=NA) # create new blank metam dataframe
                    ind.new.m$ID <- paste(yy, rownames(ind.new.m), sep = ".") # give each a unique ID "year.ID"
                    ind.new.m$SEX <- sample(c("M", "F"), size = n.off, replace = T) # determine offspring sex
                    ind.new.m$HYDP <- hydp # set pond hydroperiod (can make this more complex)
                    
                    # Assign HIS (hybrid index score)
                    ind.new.m$HIS <- as.vector(unlist(mapply(rep.his, fem.df$off.his, round(fem.df$fec, digits = 0)))) # assign to all metam. using rep HIS function
                    
                    # Calculate Mass: density-independent
                    if(model.sel == "densInd"){for(ii in 1:nrow(ind.new.m)){
                        ind.new.m$MASS[ii] <- logmass.bayes(ind.new.m$HYDP[ii], ind.new.m$HIS[ii], params = params[,dd])}
                    } # calculate mass using hydp equation
                    
                    # Calculate Mass: density-dependent
                    if(model.sel != "densInd"){for(ii in 1:nrow(ind.new.m)){
                        ind.new.m$MASS[ii] <- logmass.eggdens.bayes(ind.new.m$HYDP[ii], ind.new.m$HIS[ii], params = params[,dd], egg.d)}
                    } # calculate mass using density-dep hydp equation
                    
                    ind.new.m$ln.mass <- log(ind.new.m$MASS) # log-transform mass
                } # end creating new metamorphs
                
                
                if(n.off == 0){ind.new.m <- NULL} # make null if no new offspring
                
                ind.m <- ind.new.m # create next years metamorphs
                if(!is.null(ind.m)){ind.m$age <- 1} # set metam age
                

                ind.a <- rbind(ind.m2a[ind.m2a$mat == 0 | ind.m2a$breed == 0, ], # add metamorphs that didnt breed
                               ind.a2a[ind.a2a$mat == 0 | ind.a2a$breed == 0, ], # add adults that didnt breed
                               breed.df[breed.df$remove == 0, !(colnames(breed.df) %in% c("remove", "capture"))]) # add breeders that were not removed as hybrids
                ind.a$age <- ind.a$age + 1 # add 1 year to adult age
                
               
                # Extinction check
                if(is.null(ind.m) & is.null(ind.m2a) & nrow(ind.a2a)<1){ # if no metamorphs, or adults
                    nt.df$n.meta[yy:nrow(nt.df)] <- 0 # set metam number to 0
                    nt.df$n.adult[yy:nrow(nt.df)] <- 0 # set metam number to 0
                    print("XXXXXXXXXX Extinction! XXXXXXXXXX")
                    message("XXXXXXXXXX Population Extinction! XXXXXXXXXX")}
                if(is.null(ind.m) & is.null(ind.m2a) & nrow(ind.a2a)<1){break}
                ind.m2a <- ind.a2a <- NULL # clear the transition df
                print(paste("Finished Year = ", yy, " | Pop Size = ", (max(nrow(ind.m),0) + max(nrow(ind.a), 0))))
                
                # Make a limit to pop size, nescessary for runaway computation
                if(!is.null(ind.a)){
                    if(nrow(ind.a)>=pop.limit){
                    n.a.clips <- n.a.clips + 1
                    nt.df$n.a.clipped[yy] <- nrow(ind.a) - pop.limit
                    ind.a <- ind.a[sample(1:nrow(ind.a), size = pop.limit), ]
                    #print(paste("Population Clipped to", pop.limit))
                    }
                }

                if(!is.null(ind.m)){
                    if(nrow(ind.m)>=pop.limit){
                    n.m.clips <- n.m.clips + 1
                    nt.df$n.m.clipped[yy] <- nrow(ind.m) - pop.limit
                    ind.m <- ind.m[sample(1:nrow(ind.m), size = pop.limit), ]
                    #print(paste("Population Clipped to", pop.limit))
                    }
                }

                
            }
            #if(n.clips > 0){print(paste("Pop was clipped", n.clips, "times this iter:", iter))}
            nt.df$n.all <- nt.df$n.adult+nt.df$n.meta # find total population size
            
            
            return(nt.df) # return total population size
            
        }
        
    
        
        
    run.IPM <- function(input.df=input.df, 
                        param.file = cts,
                        n.threads=4, # Set maximum threads for parallel computing
                        adult.age=3, # Age of starting adults (for initialization of population)
                        file.path="modelOutput/"){ # path to output folder (file names are auto generated)
        
        
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
            
            # Source functions for each iteration
            source("R/demographicFunctions.R")
            load("data/data.Rdata")
            
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
                                      years = input.df$years[ii], meta.df=meta.df, adult.df=adult.df, pond.volume =input.df$pond.volume[ii], 
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
        return(nt.par.out)
    }
    
        
    
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
    
    
        
        
        
        
        
        
        