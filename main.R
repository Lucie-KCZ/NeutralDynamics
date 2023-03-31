############################################
### Biodiversity time series are biased  ###
### towards increasing species richness  ###
###       in changing environments       ###
############################################

# Authors list
# Lucie Kuczynski, Vicente J. Ontiveros, Helmut Hillebrand

# Script infos
# L. Kuczynski
# lucie.kuczynski@hotmail.com
# Mar, 2022
# last edit: Dec, 2022

# dev.off() 
for(i in 1:10) gc(reset = T) ; rm(list = ls())

dataset.name <- 'RivFish'
# dataset.name <- 'bbs'
# dataset.name <- 'mzb'

dir.create(paste0('./out/', dataset.name))

# misc.  internal functions -----------------------------------------------
# # An mc-version of the sapply function.
# mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
#   FUN <- match.fun(FUN)
#   answer <- parallel::mclapply(X = X, FUN = FUN, ...)
#   if (USE.NAMES && is.character(X) && is.null(names(answer))) 
#     names(answer) <- X
#   if (!isFALSE(simplify) && length(answer)) 
#     simplify2array(answer, higher = (simplify == "array"))
#   else answer
# }

# import the data ---------------------------------------------------------
### a. OBSERVED DATA
if(dataset.name == 'RivFish') load('~/Dropbox/Neutral trends/data/RivFishTIME/RivFish.RData', v = T)
if(dataset.name == 'bbs') load('~/Documents/Work/Data/BBS/outputs/BBS.RData', v = T)
if(dataset.name == 'mzb') load('~/Dropbox/Neutral trends/data/mzb.RData', v = T)

if(dataset.name == 'RivFish') Obs.list <- RivFish.list ; rm(RivFish.list)
if(dataset.name == 'bbs') Obs.list <- BBS.list ; rm(BBS.list)

Obs.infos <- do.call(rbind, lapply(Obs.list, function(x) x$infos))
SR.distribution <- unlist(lapply(Obs.list, function(x) rowSums(ifelse(x$com.matrix > 0, 1, 0))))

if(dataset.name == 'RivFish') {
  # only European sites
  Obs.Europe <- Obs.infos[Obs.infos$BIOME_MAP == 'Palearctic', ]
  
  # only 10 sampled years
  Obs.Europe10 <- Obs.Europe[Obs.Europe$NB_SAMP > 9, ]
  
  # final data
  Obs.list <- Obs.list[Obs.Europe10$STUDY_ID]  
}

### b. RIVFISHTIME - RANDOMISATION
# 9 randomisations per observed TS
randomisation <- function(samp, k){
  
  rownames(samp$com.matrix) <- sample(rownames(samp$com.matrix))
  samp$com.matrix <- samp$com.matrix[order(rownames(samp$com.matrix)), ]
  
  samp$infos$STUDY_ID <- paste(samp$infos$STUDY_ID, k, sep = '_')
  
  random.samp <- list(infos = samp$infos, com.matrix = samp$com.matrix)
  return(random.samp)
  
}

library(parallel)
cl <- makeCluster(detectCores() - 3)
clusterExport(cl, "randomisation")

system.time(null.list <-
              parLapply(cl, Obs.list, function(y) sapply(1:9, function(x) randomisation(samp = y, k = x), simplify = F)))
# about 10 sec. / 15 sec

stopCluster(cl) ; rm(randomisation)

null.list <- unlist(null.list, recursive = F)

### c. SIMULATED DATA
library(tidyverse)

# SR and TS length obs distributions
# getting the Obs data formatted
# source(textConnection(readLines('main.R')[18:29]))

# getting the distribution of TS length and NB Samples
INFOS <- Obs.infos[, 1:4]

# simulations
files <- list.files('~/Dropbox/Neutral trends/analysis/functions island package/R/', full.names = T)
for(i in 1:17) source(files[i]) ; rm(files)

simulating.ts <- function(i, constraints = INFOS, autocor.factor = 1, rate.c, rate.e, SR.dist = SR.distribution) {
  
  print(i)
  
  length.sim <- sample(constraints$TS_LENGTH, size = 1)
  nb.samp.sim <- sample(constraints$NB_SAMP[constraints$NB_SAMP <= length.sim], size = 1)
  nb.sp.sim <- sample(SR.dist, size = 1)
  
  sp.pool <- 0
  while(sp.pool < nb.sp.sim){
    sp.pool <- sample(x = constraints$NB.TOT.SP[constraints$NB.TOT.SP >= nb.sp.sim], 1)
    # print(paste('settings:', nb.sp.sim, sp.pool))
  }
  
  
  initial <- matrix(c(rep(0, sp.pool - nb.sp.sim), rep(1, nb.sp.sim)), ncol = 1)
  
  # Set colonization and extinction rates
  col <- rate.c / autocor.factor
  ext <- rate.e / autocor.factor
  1 - sum(cetotrans(col, ext)) # This is the temporal autocorrelation. See Ontiveros et al 2021 appendix.
  
  # running the simulation to get the com dynamics
  sim.ts <- t(PA_simulation(initial, 1, cetotrans(col, ext), length.sim))
  
  # having the 'years'
  rownames(sim.ts) <- 1:nrow(sim.ts)
  
  # selecting the sampled years based on the nb of samples
  sim.ts <- sim.ts[sort(sample(1:length.sim, nb.samp.sim)), ]
  
  # format
  list.sim <- list(infos = data.frame(STUDY_ID = paste0('sim_', i),
                                      NB_SAMP = nb.samp.sim,
                                      TS_LENGTH = length.sim,
                                      NB.TOT.SP = ncol(sim.ts),
                                      LATITUDE = NA, LONGITUDE = NA, ID_LUCIE = NA,
                                      PROTECTED_AREA = NA, BIOME.MAP = NA, GRAIN_SQ_KM = NA,
                                      AREA_SQ_KM = NA, SUMMARY_METHODS = 0,
                                      ABUNDANCE_TYPE = 0, BIOMASS = NA),
                   com.matrix = sim.ts)
  
  return(list.sim)
  
}

# getting the rate
# For a first time, you have to run the NICE function first on the obs data
# as the rates are from the estimated COL and EXT events (i.e. by OLEs)
load(paste0('./out/', dataset.name, '/NICE.RData'))

rates <- do.call(rbind, lapply(NICE$observed, function(x) x[ncol(x), 3:4] / ncol(x)))
rates <- rates[-which(apply(rates, 1, sum) == 0), ]
rates <- apply(rates, 2, mean)
rm(NICE)

system.time(sim.list <- mapply(FUN = simulating.ts, i = 1:9999, 
                               autocor.factor = 1, rate.c = rates[1], rate.e = rates[1], 
                               SIMPLIFY = F)) 

names(sim.list) <- paste('autocor.', 1)
rm(rates)

### MERGE DATA
out <- list(observed = Obs.list, 
            simulated = sim.list, 
            null = null.list)
save(out, file = paste0('./out/', dataset.name, '/complete.datasets.RData'))

# SR trends ---------------------------------------------------------------
load(paste0('./out/', dataset.name, '/complete.datasets.RData'))
data.list <- out ; rm(list = ls()[-which(ls() %in% c('out', 'data.list', 'dataset.name'))])

sr.over.time <- function(samp, out = 'trend'){
  # out can be 'trend' or 'raw'
  #   'trend': the sample trend in SR
  #   'raw': the df w SR for each year and the STUDY_ID
  years <- as.numeric(as.character(rownames(samp$com.matrix)))
  sr <- log(apply(ifelse(samp$com.matrix > 0, 1, 0), 1, sum, na.rm = T) + 1)
  
  if(out == 'trend') {
    output <- data.frame(t(summary(lm(log(sr+1) ~ years))$coefficients[2, c(1, 4)]), TS_L = samp$infos$TS_LENGTH)
  } else {
    output <- data.frame(STUDY_ID = samp$infos$STUDY_ID,
                         YEARS = years,
                         SR = exp(sr))
  }
}

system.time(sr.trends <- lapply(data.list, function(x) do.call(rbind, lapply(x, sr.over.time)))) # ab 40 sec / 100 for bbs
system.time(sr.df <- lapply(data.list, function(x) do.call(rbind, lapply(x, sr.over.time, 'raw')))) # ab 20 sec / 70 for birds
rm(sr.over.time)
save(sr.trends, sr.df, file = paste0('./out/', dataset.name, '/log.sr.RData'))

# 1. Descriptive stats
# min, mean, sd and max
lapply(sr.trends, function(x) round(c(min(x[, 1]), mean(x[, 1]), sd(x[, 1], na.rm = T), max(x[, 1], na.rm = T)), 5))

# hw many sign. 
length(sr.trends$observed[sr.trends$observed$Pr...t.. < .05, 1])
length(sr.trends$simulated[sr.trends$simulated$Pr...t.. < .05, 1])
length(sr.trends$null[sr.trends$null$Pr...t.. < .05, 1])

# about the sign. trends
table(sign(sr.trends$observed[sr.trends$observed$Pr...t.. < .05, 1]))
table(sign(sr.trends$simulated[sr.trends$simulated$Pr...t.. < .05, 1]))
table(sign(sr.trends$null[sr.trends$null$Pr...t.. < .05, 1]))

# 2. General trends over time
load(paste0('./out/', dataset.name, '/log.sr.RData'))

# model settings
library(nlme)
ctrl <- lmeControl(optimizer = 'optim', optCtrl = list(maxfun = 200000), msMaxIter = 250, msMaxEval = 250, singular.ok = T, returnObject = T)
system.time(global.mod <-
              lapply(sr.df, function(x) lme(fixed = log(SR+1) ~ scale(YEARS), random = ~ 1 | STUDY_ID, method = 'ML', control = ctrl, data = x)))
# ca 6 sec
rm(ctrl)

# getting results
lapply(global.mod, summary)

# getting the r squared
library(MuMIn)
lapply(global.mod, r.squaredGLMM)

# 3. Effect of TS_length
library(gamlss)

for(i in 1:length(sr.trends)){
  # no na in the data
  trends.na <- na.omit(sr.trends[[i]])
  
  # gamglss
  Model <- gamlss(trends.na[, 1] ~ trends.na[, 3],
                  sigma.formula = ~ pb(as.numeric(trends.na[, 3])),
                  data = trends.na,
                  family = NO(),
                  method = RS())
  
  # printing the output
  print(paste(names(sr.trends)[i], paste0(rep('_', 80), collapse = '')))
  print(summary(Model))
  
  print('--- time to get a zero slope:')
  print(round(- Model$mu.coefficients[1] / Model$mu.coefficients[2]))
  
  print('--- R2:')
  print(Rsq(Model))
  
  # cleaning up
  rm(trends.na, Model)
  
} ; rm(i)


# 4. Comparison btw datasets
t.test(sr.trends$observed$Estimate, sr.trends$null$Estimate)
t.test(sr.trends$observed$Estimate, sr.trends$simulated$Estimate)

# NICE ---------------------------------------------------------------------
source('~/Dropbox/Neutral trends/analysis/sExtinct/R/OLE.R')
source('~/Dropbox/Neutral trends//analysis/sExtinct/R/OLE.fun.R')

NICE_function <- function(samp) {
  
  # print(samp$infos$STUDY_ID)
  
  samp.com <- samp$com.matrix
  
  # needs at least three occurrences
  if(any(which(apply(ifelse(samp.com > 0, 1, 0), 2, sum) > 2))) {
    samp.com <- as.matrix(samp.com[, which(apply(ifelse(samp.com > 0, 1, 0), 2, sum) > 2)])
  } else { return(NULL) }
  
  if(ncol(samp.com) < 2) return(NULL)
  
  if(is.null(colnames(samp.com))) colnames(samp.com) <- paste('sp', 1:ncol(samp.com), sep = '.')
  
  # str of the ts
  years <- as.numeric(rownames(samp.com))
  
  # extinction timing estimation
  est.EXT <- apply(samp.com, 2, function(x) OLE(sightingdata = cbind(years, x), alpha = .05))
  est.EXT <- data.frame(Sp = colnames(samp.com), do.call(rbind, lapply(est.EXT, function(x) cbind(x[1], x[2], x[3]))))
  colnames(est.EXT) <- c('Sp', 'Est.E', 'low.E', 'up.E')
  
  # colonisation timing estimation
  samp.COL <- cbind(abs(years - max(years)), samp.com)
  samp.COL <- samp.COL[order(samp.COL[, 1]), ]
  
  est.COL <- apply(samp.COL[, -1], 2, function(x) lapply(OLE(
    sightingdata = cbind(years, x), alpha = .05),
    function(z) abs(z - max(samp.com[, 1]))))
  est.COL <- data.frame(Sp = colnames(samp.COL[, -1]), do.call(rbind, lapply(est.COL, function(x) cbind(x[1], x[2], x[3]))))
  colnames(est.COL) <- c('Sp', 'Est.C', 'low.C', 'up.C')
  
  # combining EXT and COL
  sp.infos <- merge(est.COL, est.EXT, by = 'Sp', all = T) ; rm(est.COL, est.EXT)
  sp.infos[, -1] <- apply(sp.infos[, -1], 2, function(x) as.numeric(as.character(substr(x, 1, 4))))
  
  # counting events for each year
  COL <- apply(table(sp.infos$Sp, sp.infos$Est.C), 2, sum)
  COL <- data.frame(YEAR = as.numeric(names(COL)), COL)
  EXT <- apply(table(sp.infos$Sp, sp.infos$Est.E), 2, sum)
  EXT <- data.frame(YEAR = as.numeric(names(EXT)), EXT)
  events <- merge(COL, EXT, all = T) ; rm(COL, EXT)
  
  # removing years from outside the sampling window
  events <- merge(matrix(years, dimnames = list(NULL, 'YEAR')), events, all.x = T, all.y = F)
  
  # adding the zero when no event
  events[is.na(events)] <- 0
  
  # cumulative events over time
  events.cum <- cbind(events$YEAR, as.matrix(apply(events[, -1], 2, cumsum))) ; rm(events)
  colnames(events.cum)[1] <- 'YEAR'
  
  # if(nrow(events.cum) < 5) {
  #   out <- NULL
  # } else {
  NICE <- (events.cum[, 'COL'] - events.cum[, 'EXT']) / (events.cum[, 'COL'] + events.cum[, 'EXT'])

  NICE <- data.frame(YEAR = events.cum[, 1], NICE = NICE)
  NICE[NICE == Inf] <- NA
  NICE[NICE == -Inf] <- NA
  NICE[is.nan(NICE[, 2]), 2] <- NA
  
  out <- merge(events.cum, NICE) ; rm(events.cum, NICE)
  SR <- data.frame(YEAR = years, SR = apply(ifelse(samp.com > 0, 1, 0), 1, sum))
  out <- merge(out, SR) ; rm(SR)
  
  out <- data.frame(STUDY_ID = samp$infos$STUDY_ID,
                    out)
  
  return(out)
}

library(parallel)
cl <- makeCluster(detectCores() - 3)
clusterExport(cl, list("NICE_function", 'OLE', 'OLE.fun'))
system.time(NICE <- lapply(data.list, function(x) parLapply(cl, x, NICE_function))) # ab 1300 sec
stopCluster(cl)

save(NICE, file = paste0('./out/', dataset.name, '/NICE.RData'))

for(i in 1:10) gc(reset = T) ;  rm(list = ls()[-which(ls() %in% c('out', 'data.list', 'dataset.name'))])
load(paste0('./out/', dataset.name, '/NICE.RData'), v = T)


# 1. Descriptive metrics
# making a list to make it easier
NICE <- list()
for(i in 1:length(NICE)) NICE[[i]] <- do.call(rbind, NICE[[i]])
rm(i)
names(NICE) <- names(NICE)

lapply(NICE, function(x) round(c(min(x[, 5]), mean(x[, 5]), sd(x[, 5], na.rm = T), max(x[, 5], na.rm = T)), 5))
lapply(NICE, function(x) round(c(mean(x[, 5], na.rm = T), sd(x[, 5], na.rm = T)), 2))

lapply(NICE, summary)

NICE$observed <- na.omit(NICE$observed)
NICE$simulated <- na.omit(NICE$simulated)
NICE$null <- na.omit(NICE$null)

# 2. General trends over time
# model settings
library(nlme)
ctrl <- lmeControl(optimizer = 'optim', optCtrl = list(maxfun = 200000), msMaxIter = 250, msMaxEval = 250, singular.ok = T, returnObject = T)
system.time(global.mod <-
              lapply(NICE, function(x) lme(fixed = NICE ~ YEAR, random = ~ 1 | STUDY_ID, method = 'ML', control = ctrl, data = x)))
rm(ctrl)
# a little bit more than a min

# getting the results
lapply(global.mod, summary)

# getting the time to get a balanced value
print(round(- global.mod$observed$coefficients$fixed[1] / global.mod$observed$coefficients$fixed[2]) - 2022)
print(round(- global.mod$simulated$coefficients$fixed[1] / global.mod$simulated$coefficients$fixed[2]))

# getting the r squared
library(sjstats)
lapply(global.mod, r2)

# 3. Comparison btw datasets
t.test(NICE$observed$NICE, mu = 0)
t.test(NICE$observed$NICE, NICE$null$NICE)
t.test(NICE$observed$NICE, NICE$simulated$NICE)

 