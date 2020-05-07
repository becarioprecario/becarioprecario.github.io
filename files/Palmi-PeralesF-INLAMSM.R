### R code from vignette source 'Palmi-PeralesF-INLAMSM.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Palmi-PeralesF-INLAMSM.Rnw:3-12
###################################################
# JSS style
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)

# This is set to reduce the length of the output
options(digits = 3)

# Color for border of areas in maps
#bordercolor <- NA # No border
bordercolor <- "gray"


###################################################
### code chunk number 2: Palmi-PeralesF-INLAMSM.Rnw:782-788
###################################################
library("rgdal")

#Load SIDS data
nc.sids <- readOGR(system.file("shapes/sids.shp", package = "spData")[1],
  verbose = FALSE)
proj4string(nc.sids) <- CRS("+proj=longlat +ellps=clrk66")


###################################################
### code chunk number 3: Palmi-PeralesF-INLAMSM.Rnw:797-801
###################################################
library("spdep")
#Compute adjacency matrix, as nb object 'adj' and sparse matrix 'W'
adj <- poly2nb(nc.sids)
W <- as(nb2mat(adj, style = "B"), "Matrix")


###################################################
### code chunk number 4: Palmi-PeralesF-INLAMSM.Rnw:831-848
###################################################
# First time period
# Compute expected cases
r74 <- sum(nc.sids$SID74) / sum(nc.sids$BIR74)
nc.sids$EXP74 <- r74 * nc.sids$BIR74
# SMR
nc.sids$SMR74 <- nc.sids$SID74 / nc.sids$EXP74
# Proportion of non-white births
nc.sids$NWPROP74 <- nc.sids$NWBIR74 / nc.sids$BIR74

# Second time period
# Compute expected cases
r79 <- sum(nc.sids$SID79) / sum(nc.sids$BIR79)
nc.sids$EXP79 <- r79 * nc.sids$BIR79
# SMR
nc.sids$SMR79 <- nc.sids$SID79 / nc.sids$EXP79
# Proportion of non-white births
nc.sids$NWPROP79 <- nc.sids$NWBIR79 / nc.sids$BIR79


###################################################
### code chunk number 5: Palmi-PeralesF-INLAMSM.Rnw:856-862
###################################################
d <- data.frame(OBS = c(nc.sids$SID74, nc.sids$SID79),
  PERIOD = c(rep("74", 100), rep("79", 100)), 
  NWPROP = c(nc.sids$NWPROP74, nc.sids$NWPROP79),
  EXP = c(nc.sids$EXP74, nc.sids$EXP79))
# County-period index
d$idx <- 1:length(d$OBS)


###################################################
### code chunk number 6: Palmi-PeralesF-INLAMSM.Rnw:885-893
###################################################
library("INLAMSM")
library("INLA")

# Number of variables (i.e., periods)
k <- 2
# Define bivariate latent effect
model.indimcar <- inla.rgeneric.define(inla.rgeneric.indep.IMCAR.model,
  list(k = k, W = W))


###################################################
### code chunk number 7: Palmi-PeralesF-INLAMSM.Rnw:905-906
###################################################
model.indimcar <- inla.INDIMCAR.model(k = k, W = W)


###################################################
### code chunk number 8: Palmi-PeralesF-INLAMSM.Rnw:917-919
###################################################
A <- kronecker(Diagonal(k, 1), Matrix(1, ncol = nrow(W), nrow = 1))
e  = rep(0, k)


###################################################
### code chunk number 9: Palmi-PeralesF-INLAMSM.Rnw:928-934
###################################################
#Fit model
IIMCAR <- inla(OBS ~ 0 + PERIOD + f(idx, model = model.indimcar,
    extraconstr = list(A = as.matrix(A), e = e)),
  data = d,
  E = EXP, family = "poisson", control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE))


###################################################
### code chunk number 10: Palmi-PeralesF-INLAMSM.Rnw:949-950
###################################################
summary(IIMCAR)


###################################################
### code chunk number 11: Palmi-PeralesF-INLAMSM.Rnw:967-982
###################################################
#Get fitted data, i.e., relative risk
nc.sids$FITTED74 <- IIMCAR$summary.fitted.values[1:100, "mean"]
nc.sids$FITTED79 <- IIMCAR$summary.fitted.values[100 + 1:100, "mean"]

# Set legend
library("RColorBrewer")
ubrks <- c(0, 0.5, 0.9, 1, 1.1, 1.5, 5.1)
colorkey <- list()
colorkey$labels <- as.character(trunc(1000 * ubrks) / 1000)
colorkey$at <- ubrks[1] + ubrks[7] * 0:6 / 6

#Display fitted relative risks
spplot(nc.sids, c("SMR74", "FITTED74", "SMR79", "FITTED79"),
  col = bordercolor, lwd = 0.1, at = ubrks, colorkey = colorkey,
  col.regions = brewer.pal(6, "Oranges"))


###################################################
### code chunk number 12: Palmi-PeralesF-INLAMSM.Rnw:999-1010
###################################################
# Independent Proper MCAR model
# Define range for the autocorrelation parameter
alpha.min <- 0
alpha.max <- 1
model.indmcar <- inla.INDMCAR.model(k = k, W = W, alpha.min = alpha.min,
  alpha.max = alpha.max)

#Fit model
IPMCAR <- inla(OBS ~ 0 + PERIOD + f(idx, model = model.indmcar), data = d,
  E = EXP, family = "poisson", control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE))


###################################################
### code chunk number 13: Palmi-PeralesF-INLAMSM.Rnw:1017-1026
###################################################
# IMCAR model
model.imcar <- inla.IMCAR.model(k = k, W = W)

#Fit model
IMCAR <- inla(OBS ~ 0 + PERIOD + f(idx, model = model.imcar,
    extraconstr = list(A = as.matrix(A), e = e)),
  data = d,
  E = EXP, family = "poisson", control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE))


###################################################
### code chunk number 14: Palmi-PeralesF-INLAMSM.Rnw:1033-1041
###################################################
# Proper MCAR model
model.mcar <- inla.MCAR.model(k = k, W = W, alpha.min = alpha.min,
  alpha.max = alpha.max)

#Fit model
PMCAR <- inla(OBS ~ 0 + PERIOD + f(idx, model = model.mcar), data = d,
  E = EXP, family = "poisson", control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE))


###################################################
### code chunk number 15: Palmi-PeralesF-INLAMSM.Rnw:1044-1052
###################################################
# M-model
model.mmodel <- inla.Mmodel.model(k = k, W = W, alpha.min = alpha.min,
  alpha.max = alpha.max)

# Fit model
Mmodel <- inla(OBS ~ 0 + PERIOD + f(idx, model = model.mmodel), data = d,
  E = EXP, family = "poisson", control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE))


###################################################
### code chunk number 16: Palmi-PeralesF-INLAMSM.Rnw:1060-1075
###################################################
# Add results to nc.sids
nc.sids$IIMCAR74 <- IIMCAR$summary.fitted[1:100, "mean"]
nc.sids$IIMCAR79 <- IIMCAR$summary.fitted[100 + 1:100, "mean"]

nc.sids$IPMCAR74 <- IPMCAR$summary.fitted[1:100, "mean"]
nc.sids$IPMCAR79 <- IPMCAR$summary.fitted[100 + 1:100, "mean"]

nc.sids$IMCAR74 <- IMCAR$summary.fitted[1:100, "mean"]
nc.sids$IMCAR79 <- IMCAR$summary.fitted[100 + 1:100, "mean"]

nc.sids$PMCAR74 <- PMCAR$summary.fitted[1:100, "mean"]
nc.sids$PMCAR79 <- PMCAR$summary.fitted[100 + 1:100, "mean"]

nc.sids$Mmodel74 <- Mmodel$summary.fitted[1:100, "mean"]
nc.sids$Mmodel79 <- Mmodel$summary.fitted[100 + 1:100, "mean"]


###################################################
### code chunk number 17: Palmi-PeralesF-INLAMSM.Rnw:1081-1086
###################################################
spplot(nc.sids, c("IIMCAR74", "IIMCAR79", "IPMCAR74", "IPMCAR79",
  "IMCAR74", "IMCAR79", "PMCAR74", "PMCAR79", "Mmodel74", "Mmodel79"),
  col = bordercolor, lwd = 0.1, at = ubrks, colorkey = colorkey,
  col.regions = brewer.pal(6, "Oranges")
)


###################################################
### code chunk number 18: Palmi-PeralesF-INLAMSM.Rnw:1109-1116
###################################################
# New need a palette with 10 colors below
nw.pal <- colorRampPalette(brewer.pal(6, "Blues")[-1])
# PLot maps
spplot(nc.sids, c("NWPROP74", "NWPROP79"),
  col = bordercolor, lwd = 0.1, at = 0:10 / 10,
  col.regions = nw.pal(10)
)


###################################################
### code chunk number 19: Palmi-PeralesF-INLAMSM.Rnw:1127-1132
###################################################
# Number of areas
n <- nrow(W)
NWPROP <- matrix(NA, ncol = 2, nrow = 2 * n)
NWPROP[1:n, 1] <- nc.sids$NWPROP74
NWPROP[n + 1:n, 2] <- nc.sids$NWPROP79


###################################################
### code chunk number 20: Palmi-PeralesF-INLAMSM.Rnw:1138-1140
###################################################
d <- as.list(d)
d$NWPROP <- NWPROP


###################################################
### code chunk number 21: Palmi-PeralesF-INLAMSM.Rnw:1146-1152
###################################################
# Independent ICAR model
IIMCAR2 <- inla(OBS ~ 0 + PERIOD + NWPROP + f(idx, model = model.indimcar,
    extraconstr = list(A = as.matrix(A), e = e)),
  data = d,
  E = EXP, family = "poisson", control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE))


###################################################
### code chunk number 22: Palmi-PeralesF-INLAMSM.Rnw:1154-1159
###################################################
# Independent proper CAR model
IPMCAR2 <- inla(OBS ~ 0 + PERIOD + NWPROP + f(idx, model = model.indmcar),
  data = d,
  E = EXP, family = "poisson", control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE))


###################################################
### code chunk number 23: Palmi-PeralesF-INLAMSM.Rnw:1161-1167
###################################################
# Instrinsic MCAR model
IMCAR2 <- inla(OBS ~ 0 + PERIOD + NWPROP + f(idx, model = model.imcar,
    extraconstr = list(A = as.matrix(A), e = e)),
  data = d,
  E = EXP, family = "poisson", control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE))


###################################################
### code chunk number 24: Palmi-PeralesF-INLAMSM.Rnw:1169-1174
###################################################
# Proper MCAR model
PMCAR2 <- inla(OBS ~ 0 + PERIOD + NWPROP + f(idx, model = model.mcar),
  data = d, E = EXP, family = "poisson", 
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE))


###################################################
### code chunk number 25: Palmi-PeralesF-INLAMSM.Rnw:1176-1181
###################################################
# M-model
Mmodel2 <- inla(OBS ~ 0 + PERIOD + NWPROP + f(idx, model = model.mmodel),
  data = d, E = EXP, family = "poisson",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE))


###################################################
### code chunk number 26: Palmi-PeralesF-INLAMSM.Rnw:1193-1194
###################################################
summary(IIMCAR2)


###################################################
### code chunk number 27: Palmi-PeralesF-INLAMSM.Rnw:1206-1229
###################################################
# Add results to nc.sids
nc.sids$IIMCAR274 <- IIMCAR2$summary.fitted[1:100, "mean"]
nc.sids$IIMCAR279 <- IIMCAR2$summary.fitted[100 + 1:100, "mean"]

nc.sids$IPMCAR274 <- IPMCAR2$summary.fitted[1:100, "mean"]
nc.sids$IPMCAR279 <- IPMCAR2$summary.fitted[100 + 1:100, "mean"]

nc.sids$IMCAR274 <- IMCAR2$summary.fitted[1:100, "mean"]
nc.sids$IMCAR279 <- IMCAR2$summary.fitted[100 + 1:100, "mean"]

nc.sids$PMCAR274 <- PMCAR2$summary.fitted[1:100, "mean"]
nc.sids$PMCAR279 <- PMCAR2$summary.fitted[100 + 1:100, "mean"]

nc.sids$Mmodel274 <- Mmodel2$summary.fitted[1:100, "mean"]
nc.sids$Mmodel279 <- Mmodel2$summary.fitted[100 + 1:100, "mean"]

# PLot maps
spplot(nc.sids, c("IIMCAR274", "IIMCAR279", "IPMCAR274", "IPMCAR279",
  "IMCAR274", "IMCAR279", "PMCAR274", "PMCAR279", "Mmodel274", "Mmodel279"),
  col = bordercolor, lwd = 0.1, at = ubrks, colorkey = colorkey,
  col.regions = brewer.pal(6, "Oranges")
)



###################################################
### code chunk number 28: Palmi-PeralesF-INLAMSM.Rnw:1310-1312
###################################################
#Load data
data(CV)


###################################################
### code chunk number 29: Palmi-PeralesF-INLAMSM.Rnw:1318-1320
###################################################
#Compute sparse adjacency matrix W
W <- as(nb2mat(CV.nb, style = "B"), "Matrix")


###################################################
### code chunk number 30: Palmi-PeralesF-INLAMSM.Rnw:1327-1335
###################################################
#Data
d <- data.frame(OBS = c(CV$Obs.Cirrhosis, CV$Obs.Lung, CV$Obs.Oral),
  EXP = c(CV$Exp.Cirrhosis, CV$Exp.Lung, CV$Exp.Oral)
)
# Add disease specific intercept
d$Intercept <- rep(c("Cirrhosis", "Lung", "Oral"), each = nrow(W))
# Index for latent effect
d$idx <- 1:length(d$OBS)


###################################################
### code chunk number 31: Palmi-PeralesF-INLAMSM.Rnw:1342-1347
###################################################
#Number of diseases
k <- 3 
# Range of autocorrelation parameter
alpha.min <- 0
alpha.max <- 1


###################################################
### code chunk number 32: Palmi-PeralesF-INLAMSM.Rnw:1354-1356
###################################################
A <- kronecker(Diagonal(k, 1), Matrix(1, ncol = nrow(W), nrow = 1))
e = rep(0, k)


###################################################
### code chunk number 33: Palmi-PeralesF-INLAMSM.Rnw:1361-1393
###################################################
# Define latent IMCAR model
model <- inla.IMCAR.model(k = k, W = W)
# FIT IMCAR model
IMCAR.cval <- inla(OBS ~ 0 + Intercept + f(idx, model = model,
    extraconstr = list(A = as.matrix(A), e = e)),
  data = d, E = EXP, family = "poisson",
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  # Increase starting values of marginal precisions for model fitting
  control.mode = list(theta = c(1, 1, 1, 0, 0, 0), restart = TRUE),
  # Required to obtain more accurate gradients for model fitting
  control.inla = list(h = 0.001),
  control.predictor = list(compute = TRUE))

# Define latent PMCAR model
model <- inla.MCAR.model(k = k, W = W, alpha.min = alpha.min,
  alpha.max = alpha.max)
# Fit PMCAR model
PMCAR.cval <-  inla(OBS ~ 0 + Intercept + f(idx, model = model), data = d,
  E = EXP, family = "poisson",
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE))

# Define latent M-model
model <- inla.Mmodel.model(k = k, W = W, alpha.min = alpha.min,
  alpha.max = alpha.max)

# Fit M-model
Mmodel.cval <- inla(OBS ~ 0 + Intercept + f(idx, model = model), data = d,
  E = EXP, family = "poisson",
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE)
)


###################################################
### code chunk number 34: Palmi-PeralesF-INLAMSM.Rnw:1411-1447
###################################################
# New dataset
d2 <- subset(d,  Intercept != "Oral")
d2$idx <- 1:nrow(d2)
k <- 2

# Linear constraint for the new models
A <- kronecker(Diagonal(k, 1), Matrix(1, ncol = nrow(W), nrow = 1))
e = rep(0, k)

model <- inla.IMCAR.model(k = k, W = W)
# FIT IMCAR model
IMCAR.cval2 <- inla(OBS ~ 0 + Intercept + f(idx, model = model,
    extraconstr = list(A = as.matrix(A), e = e)),
  data = d2, E = EXP, family = "poisson",
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE))

# Define latent PMCAR model
model <- inla.MCAR.model(k = k, W = W, alpha.min = alpha.min,
  alpha.max = alpha.max)
# Fit PMCAR model
PMCAR.cval2 <-  inla(OBS ~ 0 + Intercept + f(idx, model = model),
  data = d2, E = EXP, family = "poisson",
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE))

# Define latent M-model
model <- inla.Mmodel.model(k = k, W = W, alpha.min = alpha.min,
  alpha.max = alpha.max)

# Fit M-model
Mmodel.cval2 <- inla(OBS ~ 0 + Intercept + f(idx, model = model),
  data = d2, E = EXP, family = "poisson",
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE)
)


###################################################
### code chunk number 35: Palmi-PeralesF-INLAMSM.Rnw:1485-1498
###################################################
n <- nrow(W)

CV$IMCAR.CIR <- IMCAR.cval$summary.fitted[1:n, "mean"] 
CV$IMCAR.LUN <- IMCAR.cval$summary.fitted[n + 1:n, "mean"] 
CV$IMCAR.ORA <- IMCAR.cval$summary.fitted[2 * n + 1:n, "mean"] 

CV$PMCAR.CIR <- PMCAR.cval$summary.fitted[1:n, "mean"] 
CV$PMCAR.LUN <- PMCAR.cval$summary.fitted[n + 1:n, "mean"] 
CV$PMCAR.ORA <- PMCAR.cval$summary.fitted[2 * n + 1:n, "mean"] 

CV$Mmodel.CIR <- Mmodel.cval$summary.fitted[1:n, "mean"] 
CV$Mmodel.LUN <- Mmodel.cval$summary.fitted[n + 1:n, "mean"] 
CV$Mmodel.ORA <- Mmodel.cval$summary.fitted[2 * n + 1:n, "mean"] 


###################################################
### code chunk number 36: Palmi-PeralesF-INLAMSM.Rnw:1503-1515
###################################################
# Set legend
library("RColorBrewer")
ubrks <- c(0, 0.5, 0.9, 1, 1.1, 1.5, 1.8)
colorkey <- list()
colorkey$labels <- as.character(trunc(1000 * ubrks) / 1000)
colorkey$at <- ubrks[1] + ubrks[7] * 0:6 / 6

spplot(CV, c("IMCAR.CIR", "PMCAR.CIR", "Mmodel.CIR", "IMCAR.LUN",
  "PMCAR.LUN", "Mmodel.LUN", "IMCAR.ORA", "PMCAR.ORA", "Mmodel.ORA"),
  col = bordercolor, lwd = 0.005, at = ubrks, colorkey = colorkey,
  col.regions = brewer.pal(6, "Oranges")
)


###################################################
### code chunk number 37: Palmi-PeralesF-INLAMSM.Rnw:1531-1536
###################################################
hyper.imcar <- inla.MCAR.transform(IMCAR.cval, 3)
hyper.pmcar <- inla.MCAR.transform(PMCAR.cval, 3, model = "PMCAR",
  alpha.min = alpha.min, alpha.max = alpha.max)
hyper.mmodel <- inla.Mmodel.transform(Mmodel.cval, 3,
  alpha.min = alpha.min, alpha.max = alpha.max)


###################################################
### code chunk number 38: Palmi-PeralesF-INLAMSM.Rnw:1549-1557
###################################################
# IMCAR
hyper.imcar$summary.hyperpar

# PMCAR
hyper.pmcar$summary.hyperpar

#M-model
hyper.mmodel$summary.hyperpar


###################################################
### code chunk number 39: Palmi-PeralesF-INLAMSM.Rnw:1595-1600
###################################################
# IMCAR
hyper.imcar$VAR.m

# PMCAR 
hyper.pmcar$VAR.m


###################################################
### code chunk number 40: Palmi-PeralesF-INLAMSM.Rnw:1607-1609
###################################################
# M-model
hyper.mmodel$VAR.m


