# FINANCIAL TIME SERIES ANALYSIS: GAME STOP 
rm(list=ls())

library(tseries)  
library(sandwich)
library(lmtest)
library(urca)     ## For unit root
library(rugarch)  ## For GARCH models
library(FinTS)    ## For ArchTest (download from RForge)
library(car)
library(forecast) 
library(xts)      ## For time stamps
library(quantmod) 
source("TSA-Predict-Student-Functions.R")
source("TSA-Finance-Functions.R")

########################################
# Symbol legend
# spec1.1 s_GARCH on reduced scale
# fit1.1 s_GARCH fitted on reduced scale

#########################################
symbol<-"GME"
data <- getSymbols(Symbols = symbol, src = "yahoo", auto.assign = FALSE,
                   from = "2002-01-03")#2002-01-03
#View(data)
#head(data)
#class(data)
colnames(data) <- gsub(x = colnames(data), pattern = paste0(symbol, "."),
                       replacement = "")
data <- data.frame(Date = index(data), data, check.names = FALSE)
#head(data)
#data[1,]
#data[NROW(data),]
head(data)
data <- data.frame(data, 
                   cc.ret = c(NA, diff(log(data$Adjusted))), 
                   gkVol = .garmanklass(data = data, sd = TRUE),
                   check.names = TRUE)# Add log returns and gkVol

data <- data[-1, , drop = FALSE]# Remove 1 row of NA

data<-data[1:5482,]# Truncated from 2023-11-22

time  <- as.Date(x = data$Date)
yc    <- data$Close # Close prices
yclog <- log(yc) # Log of close prices
y     <- data$Adjusted # Adjusted for dividends
ylog  <- log(y) # Log of adjusted prices

# Number of observations
nobs <- NROW(y)
par(mfrow=c(1,1))
# PRELIMINARY ANALYSIS
# 0. Plots of yc, yclog, y, ylog
par(mfrow = c(2,2))
plot(x = time, y = yc,    main = "Close",        xlab = "", ylab = "Dollars", type = "l")
plot(x = time, y = yclog, main = "Ln(close)",    xlab = "", ylab = "Dollars", type = "l")
plot(x = time, y = y,     main = "AdjClose",     xlab = "", ylab = "Dollars", type = "l")
plot(x = time, y = ylog,  main = "Ln(AdjClose)", xlab = "", ylab = "Dollars", type = "l")

# 1., 2. Seasonality and Stationarity 
# ACF on ylog
par(mfrow = c(2,1))
Acf(x = ylog, lag.max = 100, type = "correlation", main = "Price")
Acf(x = ylog, lag.max = 100, type = "partial", main = "")
# Apparent non-stationarity, absence of seasonal component, as expected

# Time series log returns in percentage
yret <- xts(x = 100 * data$cc.ret, order.by = time)
#yret
par(mfrow = c(1,1))
plot(x = time, y = yret, main = "Returns", 
     xlab = "", ylab = "", type = "l", lwd=1)
# OBS: volatility clustering in 2021

# 3. Serial correlation of returns
par(mfrow = c(2,1))
Acf(x = yret, lag.max = 100, type = "correlation", main = "Returns")
Acf(x = yret, lag.max = 100, type = "partial", main = "")

# 4. Heteroskedasticity
# 4.1 ACF log returns
par(mfrow = c(3,1))
Acf(x = yret, lag.max = 100, type = "correlation", main = "Returns")
Acf(x = abs(yret), lag.max = 100, type = "correlation", main = "|Returns|")
Acf(x = yret^2, lag.max = 100, type = "correlation", main = expression(Returns^2))

# 4.2 Ljung Box Test
npar <- 0
lag <- c(2, 5, 10, 15, 20, 30, 50) + npar
lb <- mapply(FUN = Box.test, lag = lag, 
             MoreArgs = list(x = yret, type = "Ljung-Box", fitdf = npar))[1:3,]
print(rbind(lag = lag, lb))
# Reject for all lags

# 4.3 ARCH test for heteroskedasticity
lag <- c(4, 8, 12, 16)
at <- mapply(FUN = ArchTest, lags = lag, 
             MoreArgs = list(x = yret, demean = TRUE))
print(rbind(lag = lag, at[1:3,]))
# Reject for every lag

# Unconditional distribution
par(mfrow = c(1,2))
.hist(x = yret, xlim = c(-10, 10), n = 200, breaks = 200, main = "Returns")
#xx<-seq(-10, 10, by=.01)
#lines(xx, dt(xx, df=1))
qqnorm(y = scale(yret))
abline(a = 0, b = 1, col = "red", lwd=2) # Leptokurtic distribution

##########################################################################
# MODELING USING ARMA
# ARMA(0,0)-std Akaike 5.1403 Bayes 5.1439
# ARMA(1,0)-std Akaike  5.1415 Bayes 5.1463
# ARMA(0,1)-std Akaike  5.1415 Bayes 5.1463
# ARMA(1,1)-std Akaike  5.1398 Bayes 5.1458 
# sGARCH(1,1)-std Akaike 4.988567 Bayes 4.994596

spec0<-arfimaspec(mean.model = list(armaOrder = c(1, 1), include.mean = TRUE, 
                                    external.regressors = NULL), distribution.model = "std") 
fit0 <- arfimafit(spec = spec0, data = yret, 
                  solver = "solnp")
#class(fit0)
np0 <- NROW(fit0@fit$coef)
#np0
print(infocriteria(fit0))
print(fit0@fit$matcoef)
print(fit0@fit$robust.matcoef)

res <- as.numeric(residuals(fit0))
par(mfrow = c(1,1))
Acf(x = res, lag.max = 100, type = "correlation", main = "Returns")
Acf(x = abs(res), lag.max = 100, type = "correlation", main = "|res|")
Acf(x = res^2, lag.max = 100, type = "correlation", main = expression(res^2))
# OBSERVATIONS:
# ACF of residuals: correlation present but not particularly strong, although 
# some values exit the bands, they don't exceed |.10|
# Strong correlation for absolute value of residuals, values reaching 0.4
# Poor correlation for squared residuals

par(mfrow = c(1,2))
xlim <- c(-5, 5)
.hist.fit(fit=fit0, xlim = xlim, ylim = c(0,0.9), n = 200, breaks = 100, 
          plot.norm = TRUE, main = "")
legend(x="topleft", legend=c("norm","std"), col=c("blue", "red"), lty=c(1,1), bty="n")
.qqplot.fit(fit = fit0)
# The std distribution seems suitable for these residuals, more than norm

# Modeling with GARCH
# Simple GARCH + constant
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
                    submodel = NULL, external.regressors = NULL, variance.targeting = FALSE), 
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,  
                    external.regressors = NULL), distribution.model = "std")
fit1 <- ugarchfit(spec = spec1, data = yret, solver = "solnp")
np1 <- NROW(fit1@fit$coef)
print(infocriteria(fit1));print(fit1@fit$matcoef);print(fit1@fit$robust.matcoef)

# Diagnostics: using (rt-mu)/sigma_t
fit <- fit1
# 1. ACF on residuals
par(mfrow = c(3,1))
Acf(x = fit@fit$z,      lag.max = 100, type = "correlation", main = "z")
Acf(x = abs(fit@fit$z), lag.max = 100, type = "correlation", main = "|z|")
Acf(x = fit@fit$z^2,    lag.max = 100, type = "correlation", main = expression(z^2))

# 2. Ljung Box test for residuals, abs(residuals) and residuals^2
lag1 <- np1 + c(1, 2, 5, 10, 15, 20)
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = fit@fit$z, type = "Ljung-Box", fitdf = np1) )
print(rbind(lag = lag1, lb1[1:3,]))
# No significance from lag 10 onwards
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = abs(fit@fit$z), type = "Ljung-Box", fitdf = np1) )
print(rbind(lag = lag1, lb1[1:3,]))
# Highly significant at all lags
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = fit@fit$z^2, type = "Ljung-Box", fitdf = np1) )
print(rbind(lag = lag1, lb1[1:3,]))
# Not significant for any lag

# 3. ARCH test
lag <- c(4, 8, 12, 16)
at <- mapply(FUN = ArchTest, lags = lag, 
             MoreArgs = list(x = fit1@fit$z, demean = TRUE))
print(at[1:3,])
# Not significant for any lag

par(mfrow = c(1,2))
xlim <- c(-4, 4)
.hist.fit(fit = fit1, xlim = xlim, ylim = c(0,0.55), n = 200, breaks = 100, 
          plot.norm = TRUE, main = "")
legend(x="topleft", legend=c("norm","std"), col=c("blue", "red"), lty=c(1,1), bty="n")
.qqplot.fit(fit = fit1)
#########################################
# Leverage effect
print(signbias(fit1))
# Positive Sign Bias 2.6481340 0.00811695 ***
# Joint Effect       8.9555465 0.02988780  **

# Moving to GJR-GARCH
# GJR-GARCH model
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), 
                    submodel = NULL, external.regressors = NULL, variance.targeting = FALSE), 
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE, 
                    external.regressors = NULL), distribution.model = "std")
fit2 <- ugarchfit(spec = spec2, data = yret, solver = "solnp")
print(infocriteria(fit2));print( fit2@fit$matcoef ); print(fit2@fit$robust.matcoef)

# T-GARCH model
spec3<-ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(1, 1), 
                                         submodel ="TGARCH", external.regressors = NULL, variance.targeting = FALSE), 
                   mean.model = list(armaOrder = c(0, 0), include.mean = TRUE, 
                                     external.regressors = NULL), distribution.model = "std")
fit3 <- ugarchfit(spec = spec3, data = yret, solver = "solnp")
print(infocriteria(fit3));print( fit3@fit$matcoef )

##############################################

# Parameter stability
print(nyblom(fit1))# Stability is rejected, shortening the series
names(data) 
data2<-data[data$Date>"2008-01-01",]# New shortened series
time2<-data2$Date
yret2<-xts(x = 100 * data2$cc.ret, order.by = time2)
# Simple-GARCH on reduced series
spec1.1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
                                          submodel = NULL, external.regressors = NULL, variance.targeting = TRUE), 
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,  
                                      external.regressors = NULL), distribution.model = "std")
fit1.1 <- ugarchfit(spec = spec1.1, data = yret2, solver = "solnp")
print(infocriteria(fit1.1));print(fit1.1@fit$matcoef);print(fit1.1@fit$robust.matcoef)
print(nyblom(fit1.1))
np1.1<-NROW(fit1.1@fit$coef)
# Joint Nyblom test; Year <2023-23-11 
# 2.76                   >2004
# 2.95                   >2006
# 2.69                   >2008 --> alpha1 (most significant) very unstable 
# 3.49                   >2010
# 4.02                   >2012
# 4.31                   >2014
# 5.40                   >2016
# 3.73                   >2018 --> individual parameter stability rejected
signbias(fit1.1)
spec1.2<-ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1,1), 
                                          submodel = NULL, external.regressors = NULL, variance.targeting =TRUE), 
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,  
                                      external.regressors = NULL), distribution.model = "std")
fit1.2 <- ugarchfit(spec = spec1.2, data = yret2, solver = "solnp")
print(infocriteria(fit1.2)); print(fit1.2@fit$matcoef); print(fit1.2@fit$robust.matcoef)
np1.2<-NROW(fit1.2@fit$coef)
# Diagnostics shortened sGARCH (fit1.1) and GJR (fit1.2)
fit <- fit1.1 
fit<-fit1.2
np<-np1.1  
np<-np1.2
# 1. ACF on residuals
par(mfrow = c(3,1))
Acf(x = fit@fit$z,      lag.max = 100, type = "correlation", main = "z")
Acf(x = abs(fit@fit$z), lag.max = 100, type = "correlation", main = "|z|")
Acf(x = fit@fit$z^2,    lag.max = 100, type = "correlation", main = expression(z^2))

# 2. Ljung Box test for residuals, abs(residuals) and residuals^2
lag1 <- np + c(1, 2, 5, 10, 15, 20)
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = fit@fit$z, type = "Ljung-Box", fitdf = np) )
print(rbind(lag = lag1, lb1[1:3,]))
#
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = abs(fit@fit$z), type = "Ljung-Box", fitdf = np) )
print(rbind(lag = lag1, lb1[1:3,]))
#
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = fit@fit$z^2, type = "Ljung-Box", fitdf = np) )
print(rbind(lag = lag1, lb1[1:3,]))

par(mfrow = c(1,2))
xlim <- c(-4, 4)
.hist.fit(fit = fit, xlim = xlim, ylim = c(0,0.55), n = 200, breaks = 100, 
          plot.norm = TRUE, main = "")
legend(x="topleft", legend=c("norm","std"), col=c("blue", "red"), lty=c(1,1), bty="n")
.qqplot.fit(fit = fit)

###########################################################################
# TGARCH
spec1.3<-ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
                                         submodel = "TGARCH", external.regressors = NULL, variance.targeting = TRUE), 
                   mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,  
                                     external.regressors = NULL), distribution.model = "std")
fit1.3<-ugarchfit(spec = spec1.3, data = yret2, solver = "solnp")
print(infocriteria(fit1.3)); print(fit1.3@fit$matcoef); print(fit1.3@fit$robust.matcoef)

###############################################################################

# News Impact Curve (NIC)
ni1 <- newsimpact(z = NULL, fit1.1)
ni2 <- newsimpact(z = NULL, fit1.2)
ni3 <- newsimpact(z = NULL, fit1.3)
legend <- c("Simple-GARCH", "GJR-GARCH", "T-GARCH")
col  <- c("black", "red", "blue")
ylim <- range( ni1$zy, ni2$zy, ni3$zy )
par(mfrow = c(1,1), mar = c(4, 4.5, 3, 1) + 0.1, lwd = 2)
plot(x = ni1$zx, y = ni1$zy, ylab = ni1$yexpr, xlab = ni1$xexpr, type = "l", 
     ylim = ylim, main = "News Impact Curve", col = col[1])
lines(x = ni2$zx, y = ni2$zy, col = col[2], lwd = 2)
lines(x = ni3$zx, y = ni3$zy, col = col[3], lwd=2)
legend(x = "topright", y = NULL, legend = legend, border = FALSE, col = col, 
       lty = 1, text.col = col, bty="n")
#abline(h=0)
#abline(v=0)
par(mfrow=c(1,1), lwd=1)

# Outlier insertion
# 1 outlier search through yret2 quantiles
#quantile(yret2, probs=seq(0,1, by=.025))
#q1<-quantile(yret2, prob=.01)
#q2<-quantile(yret2, prob=.99)
#hist(yret2, breaks=200)
#dummy0<-ifelse(yret2 < q1, 1, 0)
#dummy1<-ifelse(yret2 > q2, 1, 0)
#dummy<-as.matrix(cbind(dummy0, dummy1))
#dummy<-ifelse(data2$Date>"2021-09-01",1,0)

#class(dummy)
#spec1.2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
#                                           submodel = NULL, external.regressors=dummy,
#                                           variance.targeting = FALSE),
#                                           mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,  
#                                           external.regressors = dummy), distribution.model = "std")
#fit1.2 <- ugarchfit(spec = spec1.2, data = yret2, solver = "solnp")
#print(infocriteria(fit1.2));print(fit1.2@fit$matcoef);print(fit1.2@fit$robust.matcoef)
#print(nyblom(fit1.2))
#summary(yret2[data2$Date>="2021-01-01"])
#summary(yret2[data2$Date<"2021-01-01"])

# Residuals
fit <- fit1.1
par(mfrow = c(3,1), lty=1)
Acf(x = fit@fit$z, lag.max = 100, type = "correlation", main = "z")
Acf(x = abs(fit@fit$z), lag.max = 100, type = "correlation", main = "|z|")
Acf(x = fit@fit$z^2, lag.max = 100, type = "correlation", main = expression(z^2))
# Ljung-Box statistics on z residuals
lag1 <- np1.1 + c(1, 2, 5, 10, 15, 20)
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = fit@fit$z, type = "Ljung-Box", fitdf = np1.1) )
print(rbind(lag = lag1, lb1[1:3,]))
# Ljung-Box statistics on |z residuals|
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = abs(fit@fit$z), type = "Ljung-Box", fitdf = np1.1) )
print(rbind(lag = lag1, lb1[1:3,]))
# Ljung-Box statistics on (z residuals)^2
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = fit@fit$z^2, type = "Ljung-Box", fitdf = np1.1) )
print(rbind(lag = lag1, lb1[1:3,]))

###################################################################
# FORECASTING
# Benchmark: Garman Klass
gk<-data2$gkVol*100
par(mfrow=c(2,1),lwd=1)
plot(x = time2, y = gk, type = "l", main = "Garman-Klass", ylab = "")
plot(x = time2, y = sqrt(2/pi)*abs(yret2), type = "l", main = "Absolute returns", ylab = "")

par(mfrow = c(1,1), lwd = 1)
plot(x = time2, y = gk, type = "l", xlab="",ylab = "", ylim=c(0,40))
lines(x = time2, y = fit1.1@fit$sigma, col = "red", lwd=1.5)   ## sGARCH
lines(x=time2, y=fit1.2@fit$sigma, col="blue", lwd=1.5)
legend(x="topleft", legend=c("Simple-GARCH","GJR-GARCH"), col=c("red", "blue") ,lwd=1, bty="n")
# gjrGARCH
#lines(x=time2, y=fit1.3@fit$sigma, col="Orange", lwd=1.5) #tGARCH

# Naive benchmark
naive.vol<-sd(yret2)
naive.var <- naive.vol^2

# Error measures
ErrorMeas <- data.frame(
  measure = c("Volatility", "Volatility", "Variance", "Variance"), 
  model = c("GARCH", "gjrGARCH", "GARCH", "gjrGARCH"), 
  rbind( 
    .ErrorMeasures(y = gk,   fit = fit1.1@fit$sigma,   naive = naive.vol),
    .ErrorMeasures(y = gk,   fit = fit1.2@fit$sigma,   naive = naive.vol),
    .ErrorMeasures(y = gk^2, fit = fit1.1@fit$sigma^2, naive = naive.var),
    .ErrorMeasures(y = gk^2, fit = fit1.2@fit$sigma^2, naive = naive.var)))
print( ErrorMeas )
#require(plm)
#require(sandwich)
x1 <- .MincerZarnowitz(y = gk, fit = fit1.1@fit$sigma, msg = "GARCH\n")
x1 <- .MincerZarnowitz(y = gk, fit = fit1.2@fit$sigma, msg = "GJR-GARCH\n")

h<-5; msg <- "GARCH vs GJR-GARCH ->"
x1 <- .DieboldMariano(y = gk^2, f1 = fit1.1@fit$sigma^2, f2 = fit1.2@fit$sigma^2, h = h, loss = "SE", msg = msg)
x1 <- .DieboldMariano(y = gk^2, f1 = fit1.1@fit$sigma^2, f2 = fit1.2@fit$sigma^2, h = h, loss = "AE", msg = msg)
x1 <- .DieboldMariano(y = gk^2, f1 = fit1.1@fit$sigma^2, f2 = fit1.2@fit$sigma^2, h = h, loss = "LLE", msg = msg)

H<-10
forc1 <- ugarchforecast(fitORspec = fit1.2, n.ahead = H, 
                        data = NULL, out.sample =0, n.roll =0)
plot(1:10, forc1@forecast$sigmaFor, type="l")
#names(forc1)
#?ugarchforecast
plot(forc1)

spec1x <- getspec(fit1.2)
setfixed(spec1x) <- as.list(coef(fit1.2))
forc2a <- ugarchforecast(fitORspec = spec1x, n.ahead = H, 
                         data = yret2[1:(NROW(yret2)-H)], n.roll = 0)
forc2a
forc2 <- ugarchforecast(fitORspec = spec1x, n.ahead = H, 
                        data = yret2, out.sample = 30, n.roll = 0)
forc2
?ugarchforecast

J<-20
nobs2<-NROW(yret2)
t1 <- nobs2 - J + 1
 
pred1 <- .GARCH.predict(object = fit1.1, n.ahead = 1, t = t1, fixed.n.ahead = TRUE)
pred1
par(mfrow=c(1,1))
plot(1:20, gk[t1:NROW(yret2)], type="l", xlab="",ylab="volatility")
lines(1:20, (pred1$pred[,3]), col="red")
legend(x="topleft", legend=c("Garman-Klass","Ex-post"), col=c("black", "red"),
       lwd=1, bty="n")