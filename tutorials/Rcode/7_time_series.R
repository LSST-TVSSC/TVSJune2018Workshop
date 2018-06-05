# Time series analysis  
# based on Chapter 11, Modern Statistical Methods for Astronomy (Feigelson & Babu 2012 & 2nd ed)
# Eric Feigelson
# NAOJ tutorials, February 2017

setwd('/Users/ericfeigelson/Desktop/NAOJ_2017')

###
### HD 3651 radial velocities: observed and folded time series 
###

rv <- read.table('HD3651_rv.dat')[1:3]
names(rv) <- c('Day', 'Vel', 'sigVel')

install.packages('Hmisc')  ;  library(Hmisc)
plot(rv[,1], rv[,2], pch=20)
errbar(rv[,1], rv[,2], (rv[,2]+rv[,3]), (rv[,2]-rv[,3]), xlab='JD - 244000', 
	ylab='Velocity (m/s)', cap=0.01)
text(10800, 30, 'HD 3651    P=62.257d   e=0.63', pos=4)

install.packages('locfit')  ;  library(locfit)
locfit_phase <- locfit(rv[,2] ~ lp(rv[,1] %% 62.257, nn=0.5), weights=1/rv[,3]^2)
plot(locfit_phase, ylim=c(-30,20), band='local', 
	xlab='JD - 2440000 mod P=62.257d', ylab='Velocity (m/s)')
points(rv[,1] %% 62.257, rv[,2], pch=20)

# Lomb-Scargle periodogram for HD 3651

install.packages('lomb')  ;  library(lomb)
lsp(rv[,1:2], from=0.00, to=0.05, ylab='Lomb-Scargle', alpha=NULL, main='')
arrows(1/64.257, 35, 1/64.257, 30, lwd=2, length=0.1)

# Various astronomical periodograms for HD 3651

install.packages('RobPer')  ; library(RobPer)
RobPer_per <- seq(from=1, to=1000, length.out=10000) 
RobPer_freq <- 1/RobPer_per

par(mfrow=c(2,2))  ;  par(mar=c(5,4,1,1)) ; par(mgp=c(2,1,0))
HD3651_AOV <- RobPer(rv[,1:3], weighting=TRUE,  regression='L2', 
	model='step', steps=10, periods=RobPer_per)
plot(RobPer_freq, HD3651_AOV, type='l', xlab='frequency', ylab='AOV', 
	xlim=c(0,0.03), ylim=c(0,0.85))
arrows(1/64.257, 0.85, 1/64.257, 0.75, length=0.05)

HD3651_PDM <- RobPer(rv[,1:3], weighting=TRUE,  regression='L2', 
	model='2step', steps=10, periods=RobPer_per)
plot(RobPer_freq, HD3651_PDM, type='l', xlab='frequency', ylab='PDM', 
	xlim=c(0,0.03), ylim=c(0,0.80))
arrows(1/64.257, 0.80, 1/64.257, 0.70, length=0.05)

HD3651_L1 <- RobPer(rv[,1:3], weighting=TRUE,  regression='L1', 
	model='sine', steps=10, periods=RobPer_per)
plot(RobPer_freq, HD3651_L1, type='l', xlab='frequency', ylab='L1', 
	xlim=c(0,0.03), ylim=c(0,0.4))
arrows(1/64.257, 0.40, 1/64.257, 0.35, length=0.05)

HD3651_GCV <- RobPer(rv[,1:3], weighting=TRUE,  regression='L2', 
	model='splines', steps=10, periods=RobPer_per)
plot(RobPer_freq, HD3651_GCV, type='l', xlab='frequency', ylab='GCV', 
	xlim=c(0,0.03), ylim=c(0,0.70))
arrows(1/64.257, 0.70, 1/64.257, 0.63, length=0.05)

###
### Ingest and plot light curves for two Kepler stars
###

Kepler1 <- read.table('Kepler1.dat')[[1]]  # KIC 007596240
Kepler2 <- read.table('Kepler2.dat')[[1]]  # KIC 007609553

length(which(is.na(Kepler1))) / length(Kepler1)  # 16% NAs
length(which(is.na(Kepler2))) / length(Kepler2)  # 29% NAs

par(mfrow=c(2,1))  ;  par(mar=c(5,4,1,2))
plot(Kepler1, type='l', xlab='Time')
plot(Kepler2, type='l', xlab='Time')


###
### Properties of the Kepler 1 lightcurve (KIC 007596240)
###

median(Kepler1, na.rm=TRUE)  ;  IQR(Kepler1, na.rm=TRUE)

par(mfrow=c(3,1))  ;  par(mar=c(5,4,1,2))
acf(Kepler1, na.action=na.pass, ylim=c(-0.05, 0.2),
	 xlab='Kepler 1 lag', main='',  ci.col='black')
hist(Kepler1, freq=FALSE, breaks=200, main='', xlim=c(-20,20), 
	xlab='Kepler 1 values')

library(MASS)  # Comparison with normal distribution
Kepler1_noNA <- na.omit(Kepler1)
Kep1_mn <- fitdistr(Kepler1_noNA, 'normal')[[1]][1]  
Kep1_sd <- fitdistr(Kepler1_noNA, 'normal')[[1]][2]
curve(dnorm(x, Kep1_mn, Kep1_sd), -20, 20, add=TRUE)
library(nortest)  ;  ad.test(Kepler1)

cor.test(1:length(Kepler1), Kepler1, method='kendall')  # Tests for trend

Kep1_arima <- arima(Kepler1, order=c(2,1,2))  #  Autoregressive model
str(Kep1_arima)
acf(Kep1_arima$residuals, na.action=na.pass, ylim=c(-0.05, 0.2), 
	xlab='Kepler 1 lag', ci.col='black')

# Diagnostics tests for Kepler 1 ARIMA residuals

install.packages('nortest')  ;  library(nortest) 
ad.test(Kep1_arima$residuals)	# test for normality

install.packages('imputeTS')  ;  library(imputeTS)
arima_resids <- na.kalman(Kep1_arima$residuals)

library(car)
x <- 1:length(arima_resids) ; y <- rnorm(length(arima_resids)) 
lmobject <- lm(y ~ x)
lmobject$residuals <- arima_resids
durbinWatsonTest(lmobject)	# test for serial correlation

Box.test(arima_resids)  		#  test for independence

install.packages('tseries') ; library(tseries)
adf.test(arima_resids)			# tests for stationarity (CPU intensive)
bds.test(arima_resids)			# test for i.i.d.

###
### Properties  of Kepler 2 lightcurve (KIC 007609553)
###

# Imputation of missing values

install.packages("imputeTS") ; library(imputeTS)
Kepler2_impute <- na.kalman(Kepler2)   # impute NAs with Kalman smoother
par(mfrow=c(2,1))  ;  par(mar=c(5,4,1,2))
plot(Kepler2, type='l', xlim=c(60000,70000), xlab='Time', ylab='Kepler 2')
plot(Kepler2_impute, type='l', xlim=c(60000,70000), xlab='Time', 
	ylab='Kepler 2 imputed')

# Three periodograms: Fourier, Lomb-Scargle, epoch folding

par(mfrow=c(3,1)) ;  par(mar=c(5,4,1,2))
spec.pgram(Kepler2_impute, xlim=c(0,0.005), spans=5, taper=0.0, 
	main='', ylab='Fourier', sub='')

install.packages('lomb')  ;  library(lomb)
lsp(Kepler2, from=0.00, to=0.005, ylab='Lomb-Scargle', main='')

install.packages('RobPer')  ; library(RobPer)
Kepler2_temp <- cbind(1:length(Kepler2), Kepler2)
Kepler2_irreg <- Kepler2_temp[!is.na(Kepler2_temp[,2]),]
PDM_per <- seq(from=1, to=10001, length.out=1000)
Kepler_PDM <- RobPer(Kepler2_irreg, weighting=FALSE,  regression='L2', 
	model='step', steps=10, periods=PDM_per)	# (CPU intensive)
PDM_freq <- 1/PDM_per
plot(PDM_freq, Kepler_PDM, type='l', xlab='frequency', ylab='PDM', 
	xlim=c(0,0.005))

# Discrete wavelet transform for Kepler 3 (KIC 010191257)

install.packages('waveslim') ; library(waveslim) 
Kepler2_impute <- na.kalman(Kepler2)   # impute NAs with Kalman smoother
Kepler2_wavdat <- Kepler2_impute[1:2^16]
Kepler2_dwt <-  dwt(Kepler2_wavdat,n.levels=10) 
par(mfrow=c(3,1))
plot.ts(up.sample(Kepler2_dwt[[5]],2^5), type='h', ylab='') ; abline(h=0)
plot.ts(up.sample(Kepler2_dwt[[7]],2^7), type='h', ylab='', lwd=2) ; abline(h=0)
plot.ts(up.sample(Kepler2_dwt[[9]],2^{9}), type='h', ylab='', lwd=2) ; abline(h=0)

# Wavelet denoising

install.packages('wavethresh')  ;  library(wavethresh)
Kepler2_wd <- wd(Kepler2_wavdat)
Kepler2_wdth <- threshold(Kepler2_wd, policy='universal')		
Kepler2_thresh <- wr(Kepler2_wdth)		
par(mfrow=c(2,1))  ;  par(mar=c(4,4,1,1))
plot(Kepler2_wavdat, type='l', xlim=c(30000,40000), ylim=c(-200,200), ylab='Kepler 2')
plot(Kepler2_thresh, type='l', xlim=c(30000,40000), ylim=c(-200,200), ylab='Kepler 2')
