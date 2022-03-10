#########################
### LRR and SR trends ###
#########################

# L. Kuczynski
# lucie.kuczynski@hotmail.com
# Mar 8, 2022
# last edit: Mar 9, 2022

dev.off() 
for(i in 1:10) gc(reset = T) ; rm(list = ls())

load('sr.RData')
load('NICE.RData')

# import the data ---------------------------------------------------------
### a. RIVFISHTIME - OBSERVED DATA
load('~/Dropbox/Neutral trends/data/RivFishTIME/RivFish.RData', v = T)

# only European sites
RivFish.infos <- do.call(rbind, lapply(RivFish.list, function(x) x$infos))
RivFish.Europe <- RivFish.infos[RivFish.infos$BIOME_MAP == 'Palearctic', ]

# only 10 sampled years
RivFish.Europe10 <- RivFish.Europe[RivFish.Europe$NB_SAMP > 9, ]

# final data
RivFish.list <- RivFish.list[RivFish.Europe10$STUDY_ID]

### b. RIVFISHTIME - RANDOMISATION
# 99 randomisations per observed TS
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
              parLapply(cl, RivFish.list, function(y) sapply(1:99, function(x) randomisation(samp = y, k = x), simplify = F)))
# about 15 sec.

stopCluster(cl) ; rm(randomisation)

null.list <- unlist(null.list, recursive = F)

### c. SIMULATED DATA
# Getting 50,000 simulations
# Subset 10,000 for which we can compute LRR
data <- read.csv('~/Dropbox/Neutral trends/data/simulated data/simcom_2.0.csv')

# getting the years
data$year <- as.numeric(substr(data$year, 1, 2))
data.y <- as.data.frame(do.call(rbind, as.list(by(data = data,
                                                  INDICES = data$year,
                                                  FUN = function(x) apply(x, 2, mean)))))
data.y$day <- NULL

# getting the distribution of TS length and NB Samples
INFOS <- do.call(rbind, lapply(RivFish.list, function(x) x$infos))

subset.sim <- function(i, simulations = data.y, constraints = INFOS) {
  length.sim <- sample(constraints$TS_LENGTH, size = 1)
  nb.samp.sim <- sample(constraints$NB_SAMP[constraints$NB_SAMP <= length.sim], size = 1)
  last.potential.starting.y <- max(unique(simulations$year) - length.sim)
  start.sim <- sample(x = 0:last.potential.starting.y, size = 1)
  years.sim <- sort(sample(x = start.sim:(start.sim+length.sim), size = nb.samp.sim, replace = F))

  data.sim <- simulations[simulations$year %in% years.sim, ]
  rownames(data.sim) <- data.sim$year ; data.sim$year <- NULL
  # data.sim <- data.sim[, -which(apply(data.sim, 2, sum) == 0)]


  nb.sp.sim <- sample(constraints[constraints$NB_SAMP <= nb.samp.sim & constraints$TS_LENGTH <= length.sim, 'NB.TOT.SP'], 1)

  data.sim <- data.sim[, sample(1:ncol(data.sim), size = nb.sp.sim, replace = F)]


  list.sim <- list(infos = data.frame(STUDY_ID = paste0('sim_', i),
                                      NB_SAMP = nb.samp.sim,
                                      TS_LENGTH = length.sim,
                                      NB.TOT.SP = ncol(data.sim),
                                      LATITUDE = NA, LONGITUDE = NA, ID_LUCIE = NA,
                                      PROTECTED_AREA = NA, BIOME.MAP = NA, GRAIN_SQ_KM = NA,
                                      AREA_SQ_KM = NA, SUMMARY_METHODS = 0,
                                      ABUNDANCE_TYPE = 0, BIOMASS = NA),
                   com.matrix = data.sim)


  return(list.sim)
}

system.time(sim.list <- mapply(FUN = subset.sim, i = 1:100000, SIMPLIFY = F)) # about 205 sec

data <- list(observed = RivFish.list,
             null = null.list,
             simulated = sim.list)

### saving the three final datasets
# save(data, file = '3datasets.RData')

# SR trends ---------------------------------------------------------------

sr.over.time <- function(samp, out = 'trend'){
  # out can be 'trend' or 'raw'
  #   'trend': the sample trend in SR
  #   'raw': the df w SR for each year and the STUDY_ID
  years <- as.numeric(as.character(rownames(samp$com.matrix)))
  sr <- log(apply(ifelse(samp$com.matrix > 0, 1, 0), 1, sum, na.rm = T) + 1)

  if(out == 'trend') {
    output <- data.frame(t(summary(lm(sr ~ years))$coefficients[2, c(1, 4)]), TS_L = nrow(samp$com.matrix))
  } else {
    output <- data.frame(STUDY_ID = samp$infos$STUDY_ID,
                         YEARS = years,
                         SR = exp(sr))
  }
}

system.time(sr.trends <- lapply(data, function(x) do.call(rbind, lapply(x, sr.over.time)))) # ab 500 sec
system.time(sr.df <- lapply(data, function(x) do.call(rbind, lapply(x, sr.over.time, 'raw')))) # ab 200 sec
rm(sr.over.time)

lapply(sr.trends, function(x) round(c(min(x[, 1]), mean(x[, 1]), sd(x[, 1], na.rm = T), max(x[, 1], na.rm = T)), 2))

library(nlme)
ctrl <- lmeControl(optimizer = 'optim', optCtrl = list(maxfun = 200000), msMaxIter = 250, msMaxEval = 250, singular.ok = T, returnObject = T)
system.time(global.mod <-
              lapply(sr.df, function(x) lme(fixed = log(SR+1) ~ YEARS, random = ~ 1 | STUDY_ID, method = 'ML', control = ctrl, data = x)))
# a little bit more than a min
rm(ctrl)

lapply(global.mod, summary)

# save(sr.trends, sr.df, file = 'sr.RData')

# Figure 1 (ie SR trends & SR trends x TS_L) ------------------------------
par(mfrow = c(3, 1), cex.axis = 2.5, mar = c(3, 8, 2, 2))
hist(sr.trends$observed[, 1], border = 'black', col = 'gray90', las = 1, xlab = '', ylab = '', main = '', breaks = 15)
hist(sr.trends$null[, 1], border = 'black', col = 'gray50', las = 1, xlab = '', ylab = '', main = '', breaks = 10)
hist(sr.trends$simulated[, 1], border = 'white', col = 'black', las = 1, xlab = '', ylab = '', main = '', breaks = 10, xlim = c(-.03, .03))

library(gamlss)
trends.na <- na.omit(sr.trends$observed)
trends.na <- na.omit(sr.trends$null)
trends.na <- na.omit(sr.trends$simulated)
Model <- gamlss(trends.na[, 1] ~ trends.na[, 3],
                sigma.formula = ~ pb(as.numeric(trends.na[, 3])),
                data = trends.na,
                family = NO(),
                method = RS())

par(mfrow = c(1, 1))
centiles(Model, xvar = trends.na[, 3], cent = c(5, 10, 25, 50, 75, 90, 95, 99),
         xlab = '', ylab = '', main = '', pch = 20, cex = 1.5, xlim = c(6, 55),
         points = T, col.centiles = 'black', lwd.centiles = 3, col = 'grey', axes = F)
         # points = T, col.centiles = 'grey', lwd.centiles = 3, col = 'black', axes = F)
axis(side = 1, tcl = -.5, lwd = 2, las = 1, at = seq(10, 50, by = 10), cex.axis = 2)
axis(side = 2, tcl = -.5, lwd = 2, las = 1, at = seq(-.15, .15, by = .1), cex.axis = 2)
abline(h = 0, lwd = 3)

summary(Model)

# LRR ---------------------------------------------------------------------
source('~/Dropbox/debts_local/analysis/sExtinct/R/OLE.R')
source('~/Dropbox/debts_local/analysis/sExtinct/R/OLE.fun.R')

LRR.CE <- function(samp) {

  print(samp$infos$STUDY_ID)

  samp.com <- samp$com.matrix

  # needs at least three occurrences
  if(any(which(apply(ifelse(samp.com > 0, 1, 0), 2, sum) > 2))) {
    samp.com <- as.matrix(samp.com[, which(apply(ifelse(samp.com > 0, 1, 0), 2, sum) > 2)])
  } else { return(NULL) }

  if(ncol(samp.com) < 2) return(NULL)

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
  LRR <- log(events.cum[, 'COL'] / events.cum[, 'EXT'])

  LRRs <- data.frame(YEAR = events.cum[, 1], LRR = LRR)
  LRRs[LRRs == Inf] <- NA
  LRRs[LRRs == -Inf] <- NA
  LRRs[is.nan(LRRs[, 2]), 2] <- NA

  out <- merge(events.cum, LRRs) ; rm(events.cum, LRR, LRRs)
  SR <- data.frame(YEAR = years, SR = apply(ifelse(samp.com > 0, 1, 0), 1, sum))
  out <- merge(out, SR) ; rm(SR)

  out <- data.frame(STUDY_ID = samp$infos$STUDY_ID,
                    out)

  return(out)
}

# system.time(NICE <- lapply(data, function(x) lapply(x, LRR.CE))) # ab 8hrs
# system.time(NICE$simulated <- lapply(data$simulated, LRR.CE))
# save(NICE, file = 'NICE.RData')

# histograms
LRR.o <- do.call(rbind, NICE$observed)
LRR.n <- do.call(rbind, NICE$null)
LRR.s <- do.call(rbind, NICE$simulated)

round(c(min(LRR.o$LRR, na.rm = T), mean(LRR.o$LRR, na.rm = T), sd(LRR.o$LRR, na.rm = T), max(LRR.o$LRR, na.rm = T)), 3)
round(c(min(LRR.n$LRR, na.rm = T), mean(LRR.n$LRR, na.rm = T), sd(LRR.n$LRR, na.rm = T), max(LRR.n$LRR, na.rm = T)), 3)
round(c(min(LRR.s$LRR, na.rm = T), mean(LRR.s$LRR, na.rm = T), sd(LRR.s$LRR, na.rm = T), max(LRR.s$LRR, na.rm = T)), 3)

t.test(LRR.o$LRR, LRR.n$LRR)
t.test(LRR.o$LRR, LRR.s$LRR)

hist(LRR.o$LRR)
hist(LRR.n$LRR)
hist(LRR.s$LRR)


par(mfrow = c(3, 1), mar = c(5, 7, 2, 2), oma = c(3, 3, 3, 3), cex.axis = 2.5)
hist(LRR.o$LRR, col = 'gray75', border = 'black',
     las = 1, main = '', ylab = '', xlab = '',
     xlim = c(-3.5, 4), ylim = c(0, 1300))
lines(y = c(1200, 1200), x = c(mean(LRR.o$LRR, na.rm = T) + sd(LRR.o$LRR, na.rm = T),
                               mean(LRR.o$LRR, na.rm = T) - sd(LRR.o$LRR, na.rm = T)),
      col = 'gray75', lwd = 4)
points(x = mean(LRR.o$LRR, na.rm = T), y = 1200, pch = 20, col = 'gray75', cex = 4)

hist(LRR.n$LRR, col = 'gray40', border = 'black',
     las = 1, main = '', ylab = '', xlab = '',
     xlim = c(-3.5, 4), ylim = c(0, 60000))
lines(y = c(59000, 59000), x = c(mean(LRR.n$LRR, na.rm = T) + sd(LRR.n$LRR, na.rm = T),
                                 mean(LRR.n$LRR, na.rm = T) - sd(LRR.n$LRR, na.rm = T)),
      col = 'gray40', lwd = 4)
points(x = mean(LRR.n$LRR, na.rm = T), y = 59000, pch = 20, col = 'gray40', cex = 4)

hist(LRR.s$LRR, col = 'black', border = 'white',
     las = 1, main = '', ylab = '', xlab = '',
     xlim = c(-3.5, 4), ylim = c(0, 160))
lines(y = c(150, 150), x = c(mean(LRR.s$LRR, na.rm = T) + sd(LRR.s$LRR, na.rm = T),
                             mean(LRR.s$LRR, na.rm = T) - sd(LRR.s$LRR, na.rm = T)),
      col = 'black', lwd = 4)
points(x = mean(LRR.s$LRR, na.rm = T), y = 150, pch = 20, col = 'black', cex = 4)

# # Figure 2 (ie COL, EXT & SR over time) -----------------------------------
# 
# LRR.o <- do.call(rbind, NICE$observed)
# LRR.n <- do.call(rbind, NICE$null)
# LRR.s <- do.call(rbind, NICE$simulated)
# 
# # # getting all TS to finish the same year (ie 2019)
# # LRR.o <- do.call(rbind, as.list(by(LRR.o, LRR.o$STUDY_ID, function(x)
# #   data.frame(x$STUDY_ID, x$YEAR + (2019 - max(x$YEAR)), x$COL, x$EXT, x$LRR, x$SR))))
# # LRR.n <- do.call(rbind, as.list(by(LRR.n, LRR.n$STUDY_ID, function(x)
# #   data.frame(x$STUDY_ID, x$YEAR + (2019 - max(x$YEAR)), x$COL, x$EXT, x$LRR, x$SR))))
# # LRR.s <- do.call(rbind, as.list(by(LRR.s, LRR.s$STUDY_ID, function(x)
# #   data.frame(x$STUDY_ID, x$YEAR + (2019 - max(x$YEAR)), x$COL, x$EXT, x$LRR, x$SR))))
# # 
# # colnames(LRR.o) <- colnames(LRR.n) <- colnames(LRR.s) <-
# #   c('STUDY_ID', 'YEAR', 'COL', 'EXT', 'LRR', 'SR')
# 
# par(mfrow = c(2, 1), mar = c(5, 5, 2, 2), oma = c(3, 3, 3, 3), cex.axis = 2.5)
# 
# ## COLONISATIONS
# plot(x = 1, y = 1, xlim = c(1950, 2020), ylim = c(0, 40), type = 'n', axes = F, xlab = '', ylab = '')
# axis(side = 1, tcl = -.5, lwd = 2, las = 1, at = seq(1950, 2020, by = 10))
# axis(side = 2, tcl = -.5, lwd = 2, las = 1, at = seq(0, 40, by = 5))
# 
# by(LRR.o, LRR.o$STUDY_ID, function(x) lines(x$COL ~ x$YEAR, col = 'gray75', lwd = 1, lty = 1))
# 
# # null means
# MEANS <- do.call(rbind, as.list(by(data = LRR.n, INDICES = LRR.n$YEAR, FUN = function(x)
#   data.frame(MEAN = mean(x$COL), t(t.test(x$COL)$conf.int)))))
# MEANS <- data.frame(YEAR = as.numeric(as.character(rownames(MEANS))), MEANS)
# 
# for(i in 1:nrow(MEANS)){
#   lines(y = c(MEANS[i, 3], MEANS[i, 4]), x = c(MEANS[i, 1], MEANS[i, 1]), col = 'chocolate' ,lwd = 3)
#   points(x = MEANS[i, 1], y = MEANS[i, 2], pch = 21, bg = 'chocolate', cex = 2)
# } ; rm(i, MEANS)
# 
# # simulated means
# MEANS <- do.call(rbind, as.list(by(data = LRR.s, INDICES = LRR.s$YEAR, FUN = function(x)
#   data.frame(MEAN = mean(x$COL), t(t.test(x$COL)$conf.int)))))
# MEANS <- data.frame(YEAR = as.numeric(as.character(rownames(MEANS))), MEANS)
# 
# for(i in 1:nrow(MEANS)){
#   lines(y = c(MEANS[i, 3], MEANS[i, 4]), x = c(MEANS[i, 1], MEANS[i, 1]), col = 'chocolate4' ,lwd = 3)
#   points(x = MEANS[i, 1], y = MEANS[i, 2], pch = 21, bg = 'chocolate4', cex = 2)
# } ; rm(i, MEANS)
# 
# # observed means
# MEANbefre74 <- do.call(rbind, as.list(by(data = LRR.o, INDICES = LRR.o$YEAR, FUN =
#                                            function(x) data.frame(MEAN = mean(x$COL)))))
# 
# LRR.o.recent <- LRR.o[LRR.o$YEAR > 1970, ]
# MEANS <- do.call(rbind, as.list(by(data = LRR.o.recent, INDICES = LRR.o.recent$YEAR, FUN = function(x)
#   data.frame(MEAN = mean(x$COL), t(t.test(x$COL)$conf.int)))))
# 
# MEANbefre74 <- data.frame(YEAR = as.numeric(as.character(rownames(MEANbefre74))), MEANbefre74)
# MEANS <- data.frame(YEAR = as.numeric(as.character(rownames(MEANS))), MEANS)
# 
# MEANS <- merge(MEANS, MEANbefre74, by = 'YEAR', all = T)[, -2]
# rm(MEANbefre74, LRR.o.recent)
# 
# colnames(MEANS) <- c('YEAR', 'LOW', 'UP', 'MEAN')
# 
# for(i in 1:nrow(MEANS)){
#   lines(y = c(MEANS[i, 2], MEANS[i, 3]), x = c(MEANS[i, 1], MEANS[i, 1]), col = 'chocolate1' ,lwd = 3)
#   points(x = MEANS[i, 1], y = MEANS[i, 4], pch = 21, bg = 'chocolate1', cex = 2)
# } ; rm(i, MEANS)
# 
# 
# ## EXTINCTIONS
# plot(x = 1, y = 1, xlim = c(1950, 2020), ylim = c(0, 40), type = 'n', axes = F, xlab = '', ylab = '')
# axis(side = 1, tcl = -.5, lwd = 2, las = 1, at = seq(1950, 2020, by = 10))
# axis(side = 2, tcl = -.5, lwd = 2, las = 1, at = seq(0, 40, by = 5))
# 
# by(LRR.o, LRR.o$STUDY_ID, function(x) lines(x$EXT ~ x$YEAR, col = 'gray75', lwd = 1, lty = 1))
# 
# # null means
# MEANS <- do.call(rbind, as.list(by(data = LRR.n, INDICES = LRR.n$YEAR, FUN = function(x)
#   data.frame(MEAN = mean(x$EXT), t(t.test(x$EXT)$conf.int)))))
# MEANS <- data.frame(YEAR = as.numeric(as.character(rownames(MEANS))), MEANS)
# 
# for(i in 1:nrow(MEANS)){
#   lines(y = c(MEANS[i, 3], MEANS[i, 4]), x = c(MEANS[i, 1], MEANS[i, 1]), col = 'darkslategray3' ,lwd = 3)
#   points(x = MEANS[i, 1], y = MEANS[i, 2], pch = 21, bg = 'darkslategray3', cex = 2)
# } ; rm(i, MEANS)
# 
# # simulated means
# MEANS <- do.call(rbind, as.list(by(data = LRR.s, INDICES = LRR.s$YEAR, FUN = function(x)
#   data.frame(MEAN = mean(x$EXT), t(t.test(x$EXT)$conf.int)))))
# MEANS <- data.frame(YEAR = as.numeric(as.character(rownames(MEANS))), MEANS)
# 
# for(i in 1:nrow(MEANS)){
#   lines(y = c(MEANS[i, 3], MEANS[i, 4]), x = c(MEANS[i, 1], MEANS[i, 1]), col = 'darkslategray4' ,lwd = 3)
#   points(x = MEANS[i, 1], y = MEANS[i, 2], pch = 21, bg = 'darkslategray4', cex = 2)
# } ; rm(i, MEANS)
# 
# # observed means
# MEANbefre74 <- do.call(rbind, as.list(by(data = LRR.o, INDICES = LRR.o$YEAR, FUN =
#                                            function(x) data.frame(MEAN = mean(x$EXT)))))
# 
# LRR.o.recent <- LRR.o[LRR.o$YEAR > 1970, ]
# MEANS <- do.call(rbind, as.list(by(data = LRR.o.recent, INDICES = LRR.o.recent$YEAR, FUN = function(x)
#   data.frame(MEAN = mean(x$EXT), t(t.test(x$EXT)$conf.int)))))
# 
# MEANbefre74 <- data.frame(YEAR = as.numeric(as.character(rownames(MEANbefre74))), MEANbefre74)
# MEANS <- data.frame(YEAR = as.numeric(as.character(rownames(MEANS))), MEANS)
# 
# MEANS <- merge(MEANS, MEANbefre74, by = 'YEAR', all = T)[, -2]
# rm(MEANbefre74, LRR.o.recent)
# 
# colnames(MEANS) <- c('YEAR', 'LOW', 'UP', 'MEAN')
# 
# for(i in 1:nrow(MEANS)){
#   lines(y = c(MEANS[i, 2], MEANS[i, 3]), x = c(MEANS[i, 1], MEANS[i, 1]), col = 'darkslategray1' ,lwd = 3)
#   points(x = MEANS[i, 1], y = MEANS[i, 4], pch = 21, bg = 'darkslategray1', cex = 2)
# } ; rm(i, MEANS)
# 
# 
# 
# # SR
# par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), oma = c(3, 3, 3, 3), cex.axis = 2.5)
# 
# sr.df <- lapply(sr.df, function(x)
#   do.call(rbind, as.list(by(x, x$STUDY_ID, function(y)
#     data.frame(y$STUDY_ID, y$YEAR + (2019 - max(y$YEAR)), y$SR)))))
# colnames(sr.df$observed) <- colnames(sr.df$null) <- colnames(sr.df$simulated) <-
#   c('STUDY_ID', 'YEAR', 'SR')
# 
# # Linear models for obs data
# library(nlme)
# model <- lme(fixed = SR ~ YEAR,
#              random = ~ 1 | STUDY_ID, method = "ML",
#              data = sr.df$observed)
# 
# library(ggeffects)
# pre <- as.data.frame(ggpredict(model, "YEAR", nsim = 999, type = 'sim'))
# 
# plot(x = 1, y = 1, xlim = c(1950, 2020), ylim = c(0, 30),
#      type = 'n', axes = F, xlab = '', ylab = '')
# axis(side = 1, tcl = -.5, lwd = 2, las = 1, at = seq(1950, 2020, by = 10))
# axis(side = 2, tcl = -.5, lwd = 2, las = 1, at = seq(0, 30, by = 5))
# 
# by(sr.df$observed, sr.df$observed$STUDY_ID, function(x) lines(x$SR ~ x$YEAR, col = 'gray75', lwd = 1, lty = 1))
# 
# polygon(y = c(pre$conf.low, rev(pre$conf.high), pre$conf.low[1]),
#         x = c(pre$x, rev(pre$x), pre$x[1]), col = "gray60", border = 'gray60')
# lines(pre$x, pre$predicted, col = 'black', lwd = 2)
# 
# # for null model
# library(nlme)
# ctrl <- lmeControl(optimizer = 'optim', optCtrl = list(maxfun = 200000), msMaxIter = 250, msMaxEval = 250, singular.ok = T, returnObject = T)
# model <- lme(fixed = SR ~ YEAR,
#              random = ~ 1 | STUDY_ID, method = "ML",
#              control = ctrl, data = sr.df$null)
# 
# library(ggeffects)
# pre <- as.data.frame(ggpredict(model, "YEAR", nsim = 99, type = 'sim'))
# 
# polygon(y = c(pre$conf.low, rev(pre$conf.high), pre$conf.low[1]),
#         x = c(pre$x, rev(pre$x), pre$x[1]), col = "gray45", border = 'gray45')
# lines(pre$x, pre$predicted, col = 'black', lwd = 1)
# 
# # simulated means
# MEANS <- do.call(rbind, as.list(by(data = sr.df$simulated, INDICES = sr.df$simulated$YEAR, FUN = function(x)
#   data.frame(MEAN = mean(x$SR), t(t.test(x$SR)$conf.int)))))
# MEANS <- data.frame(YEAR = as.numeric(as.character(rownames(MEANS))), MEANS)
# 
# for(i in 1:nrow(MEANS)){
#   lines(y = c(MEANS[i, 3], MEANS[i, 4]), x = c(MEANS[i, 1], MEANS[i, 1]), col = 'black' ,lwd = 3)
#   points(x = MEANS[i, 1], y = MEANS[i, 2], pch = 20, bg = 'black', cex = 2)
# } ; rm(i, MEANS)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
