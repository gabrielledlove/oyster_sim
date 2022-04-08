### Cursory examination of pilot field experiment and power analysis on oyster mortality at different
### controlled densities (high vs low) and structural orientations (structured vs unstructured)
### for determining final necessary experimental sample size

### See Read_me file


####### Data Preparation ####### 

rm(list = ls())

library(MASS)
library(dplyr)

data2 <- read.csv('PilotDataFull.csv')            # data from pilot experiment with additional simulated data appended

data2 <- filter(data2, density != "M") %>%        # remove medium density reps; were useful for equipment testing but irrelevant to my analysis
  select(-c("notes", "ret_per")) %>%              # remove notes and tray retention percentage (not needed for this analysis)
  filter(total > 5)                               # total = number of combined dead and alive oysters remaining in place upon tray retrieval; remove trays with <5 oysters
colnames(data2)[1] <- "tray_id"                   # fix name issue for tray id column

data2$mortality <- data2$dead/data2$total         # add mortality rate 
data2$survival <- data2$live/data2$total          # add survival rate
#


####### Sample Experimental Design for Data Simulation ####### 

hd_oysters <- 90          # n oysters per entire high density tray
ld_oysters <- 16          # low density tray n
hd_half <- hd_oysters/2   # oysters per with-tray treatment at high density
ld_half <- ld_oysters/2   # oysters per with-tray treatment at low density

tot_trays <- 38           # tot_trays here = treatments 
trays_half <- tot_trays/2 # 2 treatments per tray, so trays_half = number of physical trays needed 

base_mort <- 0.15         # baseline natural mortality (between 0.1-0.2 (Eggleston 1990, La Peyre et al 2013))
hd_effect <- 1.5          # effect of high density on mortality (numbers from Eggleston 1990)
ld_effect <- 0.9          # conservative, less certain value, but low density is expected to decrease mortality
str_effect <- 1           # structure is like a "natural default" so it won't necessarily reduce or add any mortality effect
unstr_effect <- 1.4       # based off Grabowski 2004 (rough estimate)

HDU_effect <- 2           # simulate interaction effects
LDS_effect <- 0.9         # based on effect*effect values

den_fac <- as.factor(c(rep("HD", times = trays_half), rep("HD", times = trays_half),   # treatment levels: density
                       rep("LD", times = trays_half), rep("LD", times = trays_half)))        
str_fac <- as.factor(c(rep("S", times = trays_half), rep("U", times = trays_half),     # treatment levels: structure
                       rep("S", times = trays_half), rep("U", times = trays_half)))  
tray_oysters <- c(rep(hd_half, times = trays_half), rep(hd_half, times = trays_half),  # oysters counts by tray half to match treatments
                  rep(ld_half, times = trays_half), rep(ld_half, times = trays_half))


## Interaction effects ##

mat <- as.matrix(model.matrix(~den_fac*str_fac)) ; mat ; colnames(mat)    # Intercept = HDS; LDS; HDU; LDU
        

                           
####### Simple Models of Effects on Mortality ####### 

hist(data2$dead)                  # check distribution: negative binomial
table(data2$dead)                 # confirm that it isn't a zero-inflated Poisson
hist(data2$mortality)             # alternatively, could model mortality rate as a beta-distributed model

mean(data2$dead); sd(data2$dead)
mean(data2$loose[-(3:5)]/data2$oysters_start[-(3:5)])   # average loss rate (minus outliers) = 0.14


mod.nb <- glm.nb(formula = dead ~ density*structure + offset(log(total)), data = data2) 
summary(mod.nb)                   # negative binomial mod with offset to account for oysters lost throughout sampling
anova(mod.nb, test = "Chisq")     # but depending on error structure of real data, Poisson model may also be appropriate
#



####### Power Analysis ####### 

ns <- seq(from = 4, to = 100, by = 2)                 # value to adjust to determine necessary sample size; by 2 because each tray equals 2 reps
n_sims <- 200                                         # number of simulations
alpha <- 0.05                                         # significance check
pval <- matrix(NA, nrow = n_sims, ncol = 3)           # matrix to store pvalues from model
powers <- matrix(NA, nrow = length(ns), ncol = 3)     # matrix to store power values
 

for(j in 1:length(ns)) { 
  n <- ns[j]
  ######################################
  tot_trays <- n                                       # will change as ns changes
  trays_half <- tot_trays/2       
  hd_oysters <- 90                
  ld_oysters <- 16              
  hd_half <- hd_oysters/2         
  ld_half <- ld_oysters/2         
  
  base_mort <- 0.15         
  hd_effect <- 1.5          
  ld_effect <- 0.9          
  str_effect <- 1           
  unstr_effect <- 1.4
  
  HDU_effect <- 2       # from hd_effect*unstr_effect          
  LDS_effect <- 0.9     # from ld_effect*str_effect
  
  den_fac <- as.factor(rep(c("HD", "HD", "LD", "LD"), each = trays_half)) # create treatment factor levels
  str_fac <- as.factor(rep(c("S", "U", "S", "U"), each = trays_half))
  tray_oysters <- rep(c(hd_half, hd_half, ld_half, ld_half), each = trays_half)
  
  oysters_lost_HD <- rbinom(tot_trays, rep(c(hd_half, hd_half)), 0.14)    # simulate losing oysters during sampling based on avg in pilot data
  oysters_lost_LD <- rbinom(tot_trays, rep(c(ld_half, ld_half)), 0.14)
  lost <- c(oysters_lost_HD, oysters_lost_LD)
  
  
  for(i in 1:n_sims){               # generating mortality for simulated trays/treatments; as number dead per tray treatment
    fakeHDS <- rbinom(trays_half, hd_half, base_mort*hd_effect*str_effect)
    fakeHDU <- rbinom(trays_half, hd_half, base_mort*hd_effect*unstr_effect*HDU_effect)
    fakeLDS <- rbinom(trays_half, ld_half, base_mort*ld_effect*str_effect*LDS_effect)   
    fakeLDU <- rbinom(trays_half, ld_half, base_mort*ld_effect*unstr_effect)
    mortality <- c(fakeHDS, fakeHDU, fakeLDS, fakeLDU)


    dat_tmp <- data.frame(          # create data frame for analysis with generated tray data
      trays = tot_trays, 
      den_treatment = den_fac, 
      str_treatment = str_fac,
      oysters <-  tray_oysters,     # oysters per treatment (half tray)
      mortality = mortality,
      lost = lost,
      end_oysters = oysters - lost  # create total number of oysters remaining for model offset
    )
    
    
    mod.po <- glm(mortality ~ den_treatment*str_treatment + offset(log(end_oysters)), data = dat_tmp, family = poisson)        
    summary(mod.po)                 # using tray treatments as defined above 

    #####################################
    
    pval[i, 1] <- summary(mod.po)$coefficients[14] # LD; [13] = HDS = intercept
    pval[i, 2] <- summary(mod.po)$coefficients[15] # U
    pval[i, 3] <- summary(mod.po)$coefficients[16] # LD*U
    
  }

  powers[j, 1] <- sum(pval[, 1] < alpha)/n_sims    # effect detection power as num trays increases
  powers[j, 2] <- sum(pval[, 2] < alpha)/n_sims
  powers[j, 3] <- sum(pval[, 3] < alpha)/n_sims
  
  #####################################
}


par(mfrow = c(1,3))

plot(ns, powers[, 1], ylim = c(0, 1), main = "Density Effect"); abline(h = 0.8)  #graphically show how many trays needed to reach
plot(ns, powers[, 2], ylim = c(0, 1), main = "S/U Effect"); abline(h = 0.8)      # 80% probability of detecting significant effects
plot(ns, powers[, 3], ylim = c(0, 1), main = "D*S/U Effect"); abline(h = 0.8)

# since each tray is divided in half (2 treatments), actual trays needed = ns/2




