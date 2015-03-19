pkgname <- "MASS"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MASS')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Insurance")
### * Insurance

flush(stderr()); flush(stdout())

### Name: Insurance
### Title: Numbers of Car Insurance claims
### Aliases: Insurance
### Keywords: datasets

### ** Examples

## main-effects fit as Poisson GLM with offset
glm(Claims ~ District + Group + Age + offset(log(Holders)),
    data = Insurance, family = poisson)

# same via loglm
loglm(Claims ~ District + Group + Age + offset(log(Holders)),
      data = Insurance)



cleanEx()
nameEx("Null")
### * Null

flush(stderr()); flush(stdout())

### Name: Null
### Title: Null Spaces of Matrices
### Aliases: Null
### Keywords: algebra

### ** Examples

# The function is currently defined as
function(M)
{
    tmp <- qr(M)
    set <- if(tmp$rank == 0L) seq_len(ncol(M)) else  -seq_len(tmp$rank)
    qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}



cleanEx()
nameEx("OME")
### * OME

flush(stderr()); flush(stdout())

### Name: OME
### Title: Tests of Auditory Perception in Children with OME
### Aliases: OME
### Keywords: datasets

### ** Examples

# Fit logistic curve from p = 0.5 to p = 1.0
fp1 <- deriv(~ 0.5 + 0.5/(1 + exp(-(x-L75)/scal)),
             c("L75", "scal"),
             function(x,L75,scal)NULL)
nls(Correct/Trials ~ fp1(Loud, L75, scal), data = OME,
    start = c(L75=45, scal=3))
nls(Correct/Trials ~ fp1(Loud, L75, scal),
    data = OME[OME$Noise == "coherent",],
    start=c(L75=45, scal=3))
nls(Correct/Trials ~ fp1(Loud, L75, scal),
    data = OME[OME$Noise == "incoherent",],
    start = c(L75=45, scal=3))

# individual fits for each experiment

aa <- factor(OME$Age)
ab <- 10*OME$ID + unclass(aa)
ac <- unclass(factor(ab))
OME$UID <- as.vector(ac)
OME$UIDn <- OME$UID + 0.1*(OME$Noise == "incoherent")
rm(aa, ab, ac)
OMEi <- OME

library(nlme)
fp2 <- deriv(~ 0.5 + 0.5/(1 + exp(-(x-L75)/2)),
            "L75", function(x,L75) NULL)
dec <- getOption("OutDec")
options(show.error.messages = FALSE, OutDec=".")
OMEi.nls <- nlsList(Correct/Trials ~ fp2(Loud, L75) | UIDn,
   data = OMEi, start = list(L75=45), control = list(maxiter=100))
options(show.error.messages = TRUE, OutDec=dec)
tmp <- sapply(OMEi.nls, function(X)
              {if(is.null(X)) NA else as.vector(coef(X))})
OMEif <- data.frame(UID = round(as.numeric((names(tmp)))),
         Noise = rep(c("coherent", "incoherent"), 110),
         L75 = as.vector(tmp), stringsAsFactors = TRUE)
OMEif$Age <- OME$Age[match(OMEif$UID, OME$UID)]
OMEif$OME <- OME$OME[match(OMEif$UID, OME$UID)]
OMEif <- OMEif[OMEif$L75 > 30,]
summary(lm(L75 ~ Noise/Age, data = OMEif, na.action = na.omit))
summary(lm(L75 ~ Noise/(Age + OME), data = OMEif,
           subset = (Age >= 30 & Age <= 60),
           na.action = na.omit), cor = FALSE)

# Or fit by weighted least squares
fpl75 <- deriv(~ sqrt(n)*(r/n - 0.5 - 0.5/(1 + exp(-(x-L75)/scal))),
               c("L75", "scal"),
               function(r,n,x,L75,scal) NULL)
nls(0 ~ fpl75(Correct, Trials, Loud, L75, scal),
    data = OME[OME$Noise == "coherent",],
    start = c(L75=45, scal=3))
nls(0 ~ fpl75(Correct, Trials, Loud, L75, scal),
    data = OME[OME$Noise == "incoherent",],
    start = c(L75=45, scal=3))

# Test to see if the curves shift with age
fpl75age <- deriv(~sqrt(n)*(r/n -  0.5 - 0.5/(1 +
                  exp(-(x-L75-slope*age)/scal))),
                  c("L75", "slope", "scal"),
                  function(r,n,x,age,L75,slope,scal) NULL)
OME.nls1 <-
nls(0 ~ fpl75age(Correct, Trials, Loud, Age, L75, slope, scal),
    data = OME[OME$Noise == "coherent",],
    start = c(L75=45, slope=0, scal=2))
sqrt(diag(vcov(OME.nls1)))

OME.nls2 <-
nls(0 ~ fpl75age(Correct, Trials, Loud, Age, L75, slope, scal),
    data = OME[OME$Noise == "incoherent",],
    start = c(L75=45, slope=0, scal=2))
sqrt(diag(vcov(OME.nls2)))

# Now allow random effects by using NLME
OMEf <- OME[rep(1:nrow(OME), OME$Trials),]
OMEf$Resp <- with(OME, rep(rep(c(1,0), length(Trials)),
                          t(cbind(Correct, Trials-Correct))))
OMEf <- OMEf[, -match(c("Correct", "Trials"), names(OMEf))]

## Not run: 
##D ## these fail in R on most platforms
##D fp2 <- deriv(~ 0.5 + 0.5/(1 + exp(-(x-L75)/exp(lsc))),
##D              c("L75", "lsc"),
##D              function(x, L75, lsc) NULL)
##D try(summary(nlme(Resp ~ fp2(Loud, L75, lsc),
##D      fixed = list(L75 ~ Age, lsc ~ 1),
##D      random = L75 + lsc ~ 1 | UID,
##D      data = OMEf[OMEf$Noise == "coherent",], method = "ML",
##D      start = list(fixed=c(L75=c(48.7, -0.03), lsc=0.24)), verbose = TRUE)))
##D 
##D try(summary(nlme(Resp ~ fp2(Loud, L75, lsc),
##D      fixed = list(L75 ~ Age, lsc ~ 1),
##D      random = L75 + lsc ~ 1 | UID,
##D      data = OMEf[OMEf$Noise == "incoherent",], method = "ML",
##D      start = list(fixed=c(L75=c(41.5, -0.1), lsc=0)), verbose = TRUE)))
## End(Not run)


cleanEx()
nameEx("Skye")
### * Skye

flush(stderr()); flush(stdout())

### Name: Skye
### Title: AFM Compositions of Aphyric Skye Lavas
### Aliases: Skye
### Keywords: datasets

### ** Examples

# ternary() is from the on-line answers.
ternary <- function(X, pch = par("pch"), lcex = 1,
                    add = FALSE, ord = 1:3, ...)
{
  X <- as.matrix(X)
  if(any(X < 0)) stop("X must be non-negative")
  s <- drop(X %*% rep(1, ncol(X)))
  if(any(s<=0)) stop("each row of X must have a positive sum")
  if(max(abs(s-1)) > 1e-6) {
    warning("row(s) of X will be rescaled")
    X <- X / s
  }
  X <- X[, ord]
  s3 <- sqrt(1/3)
  if(!add)
  {
    oldpty <- par("pty")
    on.exit(par(pty=oldpty))
    par(pty="s")
    plot(c(-s3, s3), c(0.5-s3, 0.5+s3), type="n", axes=FALSE,
         xlab="", ylab="")
    polygon(c(0, -s3, s3), c(1, 0, 0), density=0)
    lab <- NULL
    if(!is.null(dn <- dimnames(X))) lab <- dn[[2]]
    if(length(lab) < 3) lab <- as.character(1:3)
    eps <- 0.05 * lcex
    text(c(0, s3+eps*0.7, -s3-eps*0.7),
         c(1+eps, -0.1*eps, -0.1*eps), lab, cex=lcex)
  }
  points((X[,2] - X[,3])*s3, X[,1], ...)
}

ternary(Skye/100, ord=c(1,3,2))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("addterm")
### * addterm

flush(stderr()); flush(stdout())

### Name: addterm
### Title: Try All One-Term Additions to a Model
### Aliases: addterm addterm.default addterm.glm addterm.lm
### Keywords: models

### ** Examples

quine.hi <- aov(log(Days + 2.5) ~ .^4, quine)
quine.lo <- aov(log(Days+2.5) ~ 1, quine)
addterm(quine.lo, quine.hi, test="F")

house.glm0 <- glm(Freq ~ Infl*Type*Cont + Sat, family=poisson,
                   data=housing)
addterm(house.glm0, ~. + Sat:(Infl+Type+Cont), test="Chisq")
house.glm1 <- update(house.glm0, . ~ . + Sat*(Infl+Type+Cont))
addterm(house.glm1, ~. + Sat:(Infl+Type+Cont)^2, test = "Chisq")



cleanEx()
nameEx("anova.negbin")
### * anova.negbin

flush(stderr()); flush(stdout())

### Name: anova.negbin
### Title: Likelihood Ratio Tests for Negative Binomial GLMs
### Aliases: anova.negbin
### Keywords: regression

### ** Examples

m1 <- glm.nb(Days ~ Eth*Age*Lrn*Sex, quine, link = log)
m2 <- update(m1, . ~ . - Eth:Age:Lrn:Sex)
anova(m2, m1)
anova(m2)



cleanEx()
nameEx("area")
### * area

flush(stderr()); flush(stdout())

### Name: area
### Title: Adaptive Numerical Integration
### Aliases: area
### Keywords: nonlinear

### ** Examples

area(sin, 0, pi)  # integrate the sin function from 0 to pi.



cleanEx()
nameEx("bacteria")
### * bacteria

flush(stderr()); flush(stdout())

### Name: bacteria
### Title: Presence of Bacteria after Drug Treatments
### Aliases: bacteria
### Keywords: datasets

### ** Examples

contrasts(bacteria$trt) <- structure(contr.sdif(3),
     dimnames = list(NULL, c("drug", "encourage")))
## fixed effects analyses
summary(glm(y ~ trt * week, binomial, data = bacteria))
summary(glm(y ~ trt + week, binomial, data = bacteria))
summary(glm(y ~ trt + I(week > 2), binomial, data = bacteria))

# conditional random-effects analysis
library(survival)
bacteria$Time <- rep(1, nrow(bacteria))
coxph(Surv(Time, unclass(y)) ~ week + strata(ID),
      data = bacteria, method = "exact")
coxph(Surv(Time, unclass(y)) ~ factor(week) + strata(ID),
      data = bacteria, method = "exact")
coxph(Surv(Time, unclass(y)) ~ I(week > 2) + strata(ID),
      data = bacteria, method = "exact")

# PQL glmm analysis
library(nlme)
summary(glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
                family = binomial, data = bacteria))



cleanEx()
nameEx("bandwidth.nrd")
### * bandwidth.nrd

flush(stderr()); flush(stdout())

### Name: bandwidth.nrd
### Title: Bandwidth for density() via Normal Reference Distribution
### Aliases: bandwidth.nrd
### Keywords: dplot

### ** Examples

# The function is currently defined as
function(x)
{
    r <- quantile(x, c(0.25, 0.75))
    h <- (r[2] - r[1])/1.34
    4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
}



cleanEx()
nameEx("bcv")
### * bcv

flush(stderr()); flush(stdout())

### Name: bcv
### Title: Biased Cross-Validation for Bandwidth Selection
### Aliases: bcv
### Keywords: dplot

### ** Examples

bcv(geyser$duration)



cleanEx()
nameEx("beav1")
### * beav1

flush(stderr()); flush(stdout())

### Name: beav1
### Title: Body Temperature Series of Beaver 1
### Aliases: beav1
### Keywords: datasets

### ** Examples

beav1 <- within(beav1,
               hours <- 24*(day-346) + trunc(time/100) + (time%%100)/60)
plot(beav1$hours, beav1$temp, type="l", xlab="time",
   ylab="temperature", main="Beaver 1")
usr <- par("usr"); usr[3:4] <- c(-0.2, 8); par(usr=usr)
lines(beav1$hours, beav1$activ, type="s", lty=2)
temp <- ts(c(beav1$temp[1:82], NA, beav1$temp[83:114]),
           start = 9.5, frequency = 6)
activ <- ts(c(beav1$activ[1:82], NA, beav1$activ[83:114]),
            start = 9.5, frequency = 6)

acf(temp[1:53])
acf(temp[1:53], type = "partial")
ar(temp[1:53])
act <- c(rep(0, 10), activ)
X <- cbind(1, act = act[11:125], act1 = act[10:124],
          act2 = act[9:123], act3 = act[8:122])
alpha <- 0.80
stemp <- as.vector(temp - alpha*lag(temp, -1))
sX <- X[-1, ] - alpha * X[-115,]
beav1.ls <- lm(stemp ~ -1 + sX, na.action = na.omit)
summary(beav1.ls, cor = FALSE)
rm(temp, activ)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("beav2")
### * beav2

flush(stderr()); flush(stdout())

### Name: beav2
### Title: Body Temperature Series of Beaver 2
### Aliases: beav2
### Keywords: datasets

### ** Examples

attach(beav2)
beav2$hours <- 24*(day-307) + trunc(time/100) + (time%%100)/60
plot(beav2$hours, beav2$temp, type = "l", xlab = "time",
   ylab = "temperature", main = "Beaver 2")
usr <- par("usr"); usr[3:4] <- c(-0.2, 8); par(usr = usr)
lines(beav2$hours, beav2$activ, type = "s", lty = 2)

temp <- ts(temp, start = 8+2/3, frequency = 6)
activ <- ts(activ, start = 8+2/3, frequency = 6)
acf(temp[activ == 0]); acf(temp[activ == 1]) # also look at PACFs
ar(temp[activ == 0]); ar(temp[activ == 1])

arima(temp, order = c(1,0,0), xreg = activ)
dreg <- cbind(sin = sin(2*pi*beav2$hours/24), cos = cos(2*pi*beav2$hours/24))
arima(temp, order = c(1,0,0), xreg = cbind(active=activ, dreg))

library(nlme) # for gls and corAR1
beav2.gls <- gls(temp ~ activ, data = beav2, corr = corAR1(0.8),
                 method = "ML")
summary(beav2.gls)
summary(update(beav2.gls, subset = 6:100))
detach("beav2"); rm(temp, activ)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("birthwt")
### * birthwt

flush(stderr()); flush(stdout())

### Name: birthwt
### Title: Risk Factors Associated with Low Infant Birth Weight
### Aliases: birthwt
### Keywords: datasets

### ** Examples

bwt <- with(birthwt, {
race <- factor(race, labels = c("white", "black", "other"))
ptd <- factor(ptl > 0)
ftv <- factor(ftv)
levels(ftv)[-(1:2)] <- "2+"
data.frame(low = factor(low), age, lwt, race, smoke = (smoke > 0),
           ptd, ht = (ht > 0), ui = (ui > 0), ftv)
})
options(contrasts = c("contr.treatment", "contr.poly"))
glm(low ~ ., binomial, bwt)



base::options(contrasts = c(unordered = "contr.treatment",ordered = "contr.poly"))
cleanEx()
nameEx("boxcox")
### * boxcox

flush(stderr()); flush(stdout())

### Name: boxcox
### Title: Box-Cox Transformations for Linear Models
### Aliases: boxcox boxcox.default boxcox.formula boxcox.lm
### Keywords: regression models hplot

### ** Examples

boxcox(Volume ~ log(Height) + log(Girth), data = trees,
       lambda = seq(-0.25, 0.25, length = 10))

boxcox(Days+1 ~ Eth*Sex*Age*Lrn, data = quine,
       lambda = seq(-0.05, 0.45, len = 20))



cleanEx()
nameEx("caith")
### * caith

flush(stderr()); flush(stdout())

### Name: caith
### Title: Colours of Eyes and Hair of People in Caithness
### Aliases: caith
### Keywords: datasets

### ** Examples

corresp(caith)
dimnames(caith)[[2]] <- c("F", "R", "M", "D", "B")
par(mfcol=c(1,3))
plot(corresp(caith, nf=2)); title("symmetric")
plot(corresp(caith, nf=2), type="rows"); title("rows")
plot(corresp(caith, nf=2), type="col"); title("columns")
par(mfrow=c(1,1))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("cement")
### * cement

flush(stderr()); flush(stdout())

### Name: cement
### Title: Heat Evolved by Setting Cements
### Aliases: cement
### Keywords: datasets

### ** Examples

lm(y ~ x1 + x2 + x3 + x4, cement)



cleanEx()
nameEx("confint")
### * confint

flush(stderr()); flush(stdout())

### Name: confint-MASS
### Title: Confidence Intervals for Model Parameters
### Aliases: confint.glm confint.nls confint.profile.glm
###   confint.profile.nls
### Keywords: models

### ** Examples

expn1 <- deriv(y ~ b0 + b1 * 2^(-x/th), c("b0", "b1", "th"),
               function(b0, b1, th, x) {})

wtloss.gr <- nls(Weight ~ expn1(b0, b1, th, Days),
   data = wtloss, start = c(b0=90, b1=95, th=120))

expn2 <- deriv(~b0 + b1*((w0 - b0)/b1)^(x/d0),
         c("b0","b1","d0"), function(b0, b1, d0, x, w0) {})

wtloss.init <- function(obj, w0) {
  p <- coef(obj)
  d0 <-  - log((w0 - p["b0"])/p["b1"])/log(2) * p["th"]
  c(p[c("b0", "b1")], d0 = as.vector(d0))
}

out <- NULL
w0s <- c(110, 100, 90)
for(w0 in w0s) {
    fm <- nls(Weight ~ expn2(b0, b1, d0, Days, w0),
              wtloss, start = wtloss.init(wtloss.gr, w0))
    out <- rbind(out, c(coef(fm)["d0"], confint(fm, "d0")))
  }
dimnames(out) <- list(paste(w0s, "kg:"),  c("d0", "low", "high"))
out

ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive = 20 - numdead)
budworm.lg0 <- glm(SF ~ sex + ldose - 1, family = binomial)
confint(budworm.lg0)
confint(budworm.lg0, "ldose")



cleanEx()
nameEx("contr.sdif")
### * contr.sdif

flush(stderr()); flush(stdout())

### Name: contr.sdif
### Title: Successive Differences Contrast Coding
### Aliases: contr.sdif
### Keywords: models

### ** Examples

(A <- contr.sdif(6))
zapsmall(ginv(A))



cleanEx()
nameEx("corresp")
### * corresp

flush(stderr()); flush(stdout())

### Name: corresp
### Title: Simple Correspondence Analysis
### Aliases: corresp corresp.xtabs corresp.data.frame corresp.default
###   corresp.factor corresp.formula corresp.matrix
### Keywords: category multivariate

### ** Examples

(ct <- corresp(~ Age + Eth, data = quine))
plot(ct)

corresp(caith)
biplot(corresp(caith, nf = 2))



cleanEx()
nameEx("cov.rob")
### * cov.rob

flush(stderr()); flush(stdout())

### Name: cov.rob
### Title: Resistant Estimation of Multivariate Location and Scatter
### Aliases: cov.rob cov.mve cov.mcd
### Keywords: robust multivariate

### ** Examples

set.seed(123)
cov.rob(stackloss)
cov.rob(stack.x, method = "mcd", nsamp = "exact")



cleanEx()
nameEx("cov.trob")
### * cov.trob

flush(stderr()); flush(stdout())

### Name: cov.trob
### Title: Covariance Estimation for Multivariate t Distribution
### Aliases: cov.trob
### Keywords: multivariate

### ** Examples

cov.trob(stackloss)



cleanEx()
nameEx("denumerate")
### * denumerate

flush(stderr()); flush(stdout())

### Name: denumerate
### Title: Transform an Allowable Formula for 'loglm' into one for 'terms'
### Aliases: denumerate denumerate.formula
### Keywords: models

### ** Examples

denumerate(~(1+2+3)^3 + a/b)
## which gives ~ (.v1 + .v2 + .v3)^3 + a/b



cleanEx()
nameEx("dose.p")
### * dose.p

flush(stderr()); flush(stdout())

### Name: dose.p
### Title: Predict Doses for Binomial Assay model
### Aliases: dose.p print.glm.dose
### Keywords: regression models

### ** Examples

ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive = 20 - numdead)
budworm.lg0 <- glm(SF ~ sex + ldose - 1, family = binomial)

dose.p(budworm.lg0, cf = c(1,3), p = 1:3/4)
dose.p(update(budworm.lg0, family = binomial(link=probit)),
       cf = c(1,3), p = 1:3/4)



cleanEx()
nameEx("dropterm")
### * dropterm

flush(stderr()); flush(stdout())

### Name: dropterm
### Title: Try All One-Term Deletions from a Model
### Aliases: dropterm dropterm.default dropterm.glm dropterm.lm
### Keywords: models

### ** Examples

quine.hi <- aov(log(Days + 2.5) ~ .^4, quine)
quine.nxt <- update(quine.hi, . ~ . - Eth:Sex:Age:Lrn)
dropterm(quine.nxt, test=  "F")
quine.stp <- stepAIC(quine.nxt,
    scope = list(upper = ~Eth*Sex*Age*Lrn, lower = ~1),
    trace = FALSE)
dropterm(quine.stp, test = "F")
quine.3 <- update(quine.stp, . ~ . - Eth:Age:Lrn)
dropterm(quine.3, test = "F")
quine.4 <- update(quine.3, . ~ . - Eth:Age)
dropterm(quine.4, test = "F")
quine.5 <- update(quine.4, . ~ . - Age:Lrn)
dropterm(quine.5, test = "F")

house.glm0 <- glm(Freq ~ Infl*Type*Cont + Sat, family=poisson,
                   data = housing)
house.glm1 <- update(house.glm0, . ~ . + Sat*(Infl+Type+Cont))
dropterm(house.glm1, test = "Chisq")



cleanEx()
nameEx("eagles")
### * eagles

flush(stderr()); flush(stdout())

### Name: eagles
### Title: Foraging Ecology of Bald Eagles
### Aliases: eagles
### Keywords: datasets

### ** Examples

eagles.glm <- glm(cbind(y, n - y) ~ P*A + V, data = eagles,
                  family = binomial)
dropterm(eagles.glm)
prof <- profile(eagles.glm)
plot(prof)
pairs(prof)



cleanEx()
nameEx("epil")
### * epil

flush(stderr()); flush(stdout())

### Name: epil
### Title: Seizure Counts for Epileptics
### Aliases: epil
### Keywords: datasets

### ** Examples

summary(glm(y ~ lbase*trt + lage + V4, family = poisson,
            data = epil), cor = FALSE)
epil2 <- epil[epil$period == 1, ]
epil2["period"] <- rep(0, 59); epil2["y"] <- epil2["base"]
epil["time"] <- 1; epil2["time"] <- 4
epil2 <- rbind(epil, epil2)
epil2$pred <- unclass(epil2$trt) * (epil2$period > 0)
epil2$subject <- factor(epil2$subject)
epil3 <- aggregate(epil2, list(epil2$subject, epil2$period > 0),
   function(x) if(is.numeric(x)) sum(x) else x[1])
epil3$pred <- factor(epil3$pred,
   labels = c("base", "placebo", "drug"))

contrasts(epil3$pred) <- structure(contr.sdif(3),
    dimnames = list(NULL, c("placebo-base", "drug-placebo")))
summary(glm(y ~ pred + factor(subject) + offset(log(time)),
            family = poisson, data = epil3), cor = FALSE)

summary(glmmPQL(y ~ lbase*trt + lage + V4,
                random = ~ 1 | subject,
                family = poisson, data = epil))
summary(glmmPQL(y ~ pred, random = ~1 | subject,
                family = poisson, data = epil3))



cleanEx()
nameEx("farms")
### * farms

flush(stderr()); flush(stdout())

### Name: farms
### Title: Ecological Factors in Farm Management
### Aliases: farms
### Keywords: datasets

### ** Examples

farms.mca <- mca(farms, abbrev = TRUE)  # Use levels as names
eqscplot(farms.mca$cs, type = "n")
text(farms.mca$rs, cex = 0.7)
text(farms.mca$cs, labels = dimnames(farms.mca$cs)[[1]], cex = 0.7)



cleanEx()
nameEx("fitdistr")
### * fitdistr

flush(stderr()); flush(stdout())

### Name: fitdistr
### Title: Maximum-likelihood Fitting of Univariate Distributions
### Aliases: fitdistr
### Keywords: distribution htest

### ** Examples

## avoid spurious accuracy
op <- options(digits = 3)
set.seed(123)
x <- rgamma(100, shape = 5, rate = 0.1)
fitdistr(x, "gamma")
## now do this directly with more control.
fitdistr(x, dgamma, list(shape = 1, rate = 0.1), lower = 0.001)

set.seed(123)
x2 <- rt(250, df = 9)
fitdistr(x2, "t", df = 9)
## allow df to vary: not a very good idea!
fitdistr(x2, "t")
## now do fixed-df fit directly with more control.
mydt <- function(x, m, s, df) dt((x-m)/s, df)/s
fitdistr(x2, mydt, list(m = 0, s = 1), df = 9, lower = c(-Inf, 0))

set.seed(123)
x3 <- rweibull(100, shape = 4, scale = 100)
fitdistr(x3, "weibull")

set.seed(123)
x4 <- rnegbin(500, mu = 5, theta = 4)
fitdistr(x4, "Negative Binomial")
options(op)



cleanEx()
nameEx("fractions")
### * fractions

flush(stderr()); flush(stdout())

### Name: fractions
### Title: Rational Approximation
### Aliases: fractions Math.fractions Ops.fractions Summary.fractions
###   [.fractions [<-.fractions as.character.fractions as.fractions
###   is.fractions print.fractions t.fractions
### Keywords: math

### ** Examples

X <- matrix(runif(25), 5, 5)
zapsmall(solve(X, X/5)) # print near-zeroes as zero
fractions(solve(X, X/5))
fractions(solve(X, X/5)) + 1



cleanEx()
nameEx("galaxies")
### * galaxies

flush(stderr()); flush(stdout())

### Name: galaxies
### Title: Velocities for 82 Galaxies
### Aliases: galaxies
### Keywords: datasets

### ** Examples

gal <- galaxies/1000
c(width.SJ(gal, method = "dpi"), width.SJ(gal))
plot(x = c(0, 40), y = c(0, 0.3), type = "n", bty = "l",
     xlab = "velocity of galaxy (1000km/s)", ylab = "density")
rug(gal)
lines(density(gal, width = 3.25, n = 200), lty = 1)
lines(density(gal, width = 2.56, n = 200), lty = 3)



cleanEx()
nameEx("gamma.shape.glm")
### * gamma.shape.glm

flush(stderr()); flush(stdout())

### Name: gamma.shape
### Title: Estimate the Shape Parameter of the Gamma Distribution in a GLM
###   Fit
### Aliases: gamma.shape gamma.shape.glm print.gamma.shape
### Keywords: models

### ** Examples

clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100),
    lot1 = c(118,58,42,35,27,25,21,19,18),
    lot2 = c(69,35,26,21,18,16,13,12,12))
clot1 <- glm(lot1 ~ log(u), data = clotting, family = Gamma)
gamma.shape(clot1)

gm <- glm(Days + 0.1 ~ Age*Eth*Sex*Lrn,
          quasi(link=log, variance="mu^2"), quine,
          start = c(3, rep(0,31)))
gamma.shape(gm, verbose = TRUE)
summary(gm, dispersion = gamma.dispersion(gm))  # better summary



cleanEx()
nameEx("gehan")
### * gehan

flush(stderr()); flush(stdout())

### Name: gehan
### Title: Remission Times of Leukaemia Patients
### Aliases: gehan
### Keywords: datasets

### ** Examples

library(survival)
gehan.surv <- survfit(Surv(time, cens) ~ treat, data = gehan,
     conf.type = "log-log")
summary(gehan.surv)
survreg(Surv(time, cens) ~ factor(pair) + treat, gehan, dist = "exponential")
summary(survreg(Surv(time, cens) ~ treat, gehan, dist = "exponential"))
summary(survreg(Surv(time, cens) ~ treat, gehan))
gehan.cox <- coxph(Surv(time, cens) ~ treat, gehan)
summary(gehan.cox)



cleanEx()
nameEx("ginv")
### * ginv

flush(stderr()); flush(stdout())

### Name: ginv
### Title: Generalized Inverse of a Matrix
### Aliases: ginv
### Keywords: algebra

### ** Examples

## Not run: 
##D # The function is currently defined as
##D function(X, tol = sqrt(.Machine$double.eps))
##D {
##D ## Generalized Inverse of a Matrix
##D   dnx <- dimnames(X)
##D   if(is.null(dnx)) dnx <- vector("list", 2)
##D   s <- svd(X)
##D   nz <- s$d > tol * s$d[1]
##D   structure(
##D     if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
##D     dimnames = dnx[2:1])
##D }
## End(Not run)


cleanEx()
nameEx("glm.convert")
### * glm.convert

flush(stderr()); flush(stdout())

### Name: glm.convert
### Title: Change a Negative Binomial fit to a GLM fit
### Aliases: glm.convert
### Keywords: regression models

### ** Examples

quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
quine.nbA <- glm.convert(quine.nb1)
quine.nbB <- update(quine.nb1, . ~ . + Sex:Age:Lrn)
anova(quine.nbA, quine.nbB)



cleanEx()
nameEx("glm.nb")
### * glm.nb

flush(stderr()); flush(stdout())

### Name: glm.nb
### Title: Fit a Negative Binomial Generalized Linear Model
### Aliases: glm.nb family.negbin logLik.negbin
### Keywords: regression models

### ** Examples

quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
quine.nb2 <- update(quine.nb1, . ~ . + Sex:Age:Lrn)
quine.nb3 <- update(quine.nb2, Days ~ .^4)
anova(quine.nb1, quine.nb2, quine.nb3)
## Don't show: 
## PR#1695
y <- c(7, 5, 4, 7, 5, 2, 11, 5, 5, 4, 2, 3, 4, 3, 5, 9, 6, 7, 10, 6, 12,
6, 3, 5, 3, 9, 13, 0, 6, 1, 2, 0, 1, 0, 0, 4, 5, 1, 5, 3, 3, 4)

lag1 <- c(0, 7, 5, 4, 7, 5, 2, 11, 5, 5, 4, 2, 3, 4, 3, 5, 9, 6, 7, 10,
6, 12, 6, 3, 5, 3, 9, 13, 0, 6, 1, 2, 0, 1, 0, 0, 4, 5, 1, 5, 3, 3)

lag2 <- c(0, 0, 7, 5, 4, 7, 5, 2, 11, 5, 5, 4, 2, 3, 4, 3, 5, 9, 6, 7,
10, 6, 12, 6, 3, 5, 3, 9, 13, 0, 6, 1, 2, 0, 1, 0, 0, 4, 5, 1, 5, 3)

lag3 <- c(0, 0, 0, 7, 5, 4, 7, 5, 2, 11, 5, 5, 4, 2, 3, 4, 3, 5, 9, 6,
7, 10, 6, 12, 6, 3, 5, 3, 9, 13, 0, 6, 1, 2, 0, 1, 0, 0, 4, 5, 1, 5)

(fit <- glm(y ~ lag1+lag2+lag3, family=poisson(link=identity),
            start=c(2, 0.1, 0.1, 0.1)))
try(glm.nb(y ~ lag1+lag2+lag3, link=identity))
glm.nb(y ~ lag1+lag2+lag3, link=identity,  start=c(2, 0.1, 0.1, 0.1))
glm.nb(y ~ lag1+lag2+lag3, link=identity,  start=coef(fit))
glm.nb(y ~ lag1+lag2+lag3, link=identity, etastart=rep(5, 42))
## End(Don't show)


cleanEx()
nameEx("glmmPQL")
### * glmmPQL

flush(stderr()); flush(stdout())

### Name: glmmPQL
### Title: Fit Generalized Linear Mixed Models via PQL
### Aliases: glmmPQL
### Keywords: models

### ** Examples

library(nlme) # will be loaded automatically if omitted
summary(glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
                family = binomial, data = bacteria))
## Don't show: 
 # an example of offset
summary(glmmPQL(y ~ trt + week, random = ~ 1 | ID,
                family = binomial, data = bacteria))
summary(glmmPQL(y ~ trt + week + offset(week), random = ~ 1 | ID,
                family = binomial, data = bacteria))
## End(Don't show)


cleanEx()
nameEx("housing")
### * housing

flush(stderr()); flush(stdout())

### Name: housing
### Title: Frequency Table from a Copenhagen Housing Conditions Survey
### Aliases: housing
### Keywords: datasets

### ** Examples

options(contrasts = c("contr.treatment", "contr.poly"))

# Surrogate Poisson models
house.glm0 <- glm(Freq ~ Infl*Type*Cont + Sat, family = poisson,
                  data = housing)
summary(house.glm0, cor = FALSE)

addterm(house.glm0, ~. + Sat:(Infl+Type+Cont), test = "Chisq")

house.glm1 <- update(house.glm0, . ~ . + Sat*(Infl+Type+Cont))
summary(house.glm1, cor = FALSE)

1 - pchisq(deviance(house.glm1), house.glm1$df.residual)

dropterm(house.glm1, test = "Chisq")

addterm(house.glm1, ~. + Sat:(Infl+Type+Cont)^2, test  =  "Chisq")

hnames <- lapply(housing[, -5], levels) # omit Freq
newData <- expand.grid(hnames)
newData$Sat <- ordered(newData$Sat)
house.pm <- predict(house.glm1, newData,
                    type = "response")  # poisson means
house.pm <- matrix(house.pm, ncol = 3, byrow = TRUE,
                   dimnames = list(NULL, hnames[[1]]))
house.pr <- house.pm/drop(house.pm %*% rep(1, 3))
cbind(expand.grid(hnames[-1]), round(house.pr, 2))

# Iterative proportional scaling
loglm(Freq ~ Infl*Type*Cont + Sat*(Infl+Type+Cont), data = housing)


# multinomial model
library(nnet)
(house.mult<- multinom(Sat ~ Infl + Type + Cont, weights = Freq,
                       data = housing))
house.mult2 <- multinom(Sat ~ Infl*Type*Cont, weights = Freq,
                        data = housing)
anova(house.mult, house.mult2)

house.pm <- predict(house.mult, expand.grid(hnames[-1]), type = "probs")
cbind(expand.grid(hnames[-1]), round(house.pm, 2))

# proportional odds model
house.cpr <- apply(house.pr, 1, cumsum)
logit <- function(x) log(x/(1-x))
house.ld <- logit(house.cpr[2, ]) - logit(house.cpr[1, ])
(ratio <- sort(drop(house.ld)))
mean(ratio)

(house.plr <- polr(Sat ~ Infl + Type + Cont,
                   data = housing, weights = Freq))

house.pr1 <- predict(house.plr, expand.grid(hnames[-1]), type = "probs")
cbind(expand.grid(hnames[-1]), round(house.pr1, 2))

Fr <- matrix(housing$Freq, ncol  =  3, byrow = TRUE)
2*sum(Fr*log(house.pr/house.pr1))

house.plr2 <- stepAIC(house.plr, ~.^2)
house.plr2$anova



base::options(contrasts = c(unordered = "contr.treatment",ordered = "contr.poly"))
cleanEx()
nameEx("huber")
### * huber

flush(stderr()); flush(stdout())

### Name: huber
### Title: Huber M-estimator of Location with MAD Scale
### Aliases: huber
### Keywords: robust

### ** Examples

huber(chem)



cleanEx()
nameEx("hubers")
### * hubers

flush(stderr()); flush(stdout())

### Name: hubers
### Title: Huber Proposal 2 Robust Estimator of Location and/or Scale
### Aliases: hubers
### Keywords: robust

### ** Examples

hubers(chem)
hubers(chem, mu=3.68)



cleanEx()
nameEx("immer")
### * immer

flush(stderr()); flush(stdout())

### Name: immer
### Title: Yields from a Barley Field Trial
### Aliases: immer
### Keywords: datasets

### ** Examples

immer.aov <- aov(cbind(Y1,Y2) ~ Loc + Var, data = immer)
summary(immer.aov)

immer.aov <- aov((Y1+Y2)/2 ~ Var + Loc, data = immer)
summary(immer.aov)
model.tables(immer.aov, type = "means", se = TRUE, cterms = "Var")



cleanEx()
nameEx("isoMDS")
### * isoMDS

flush(stderr()); flush(stdout())

### Name: isoMDS
### Title: Kruskal's Non-metric Multidimensional Scaling
### Aliases: isoMDS Shepard
### Keywords: multivariate

### ** Examples

swiss.x <- as.matrix(swiss[, -1])
swiss.dist <- dist(swiss.x)
swiss.mds <- isoMDS(swiss.dist)
plot(swiss.mds$points, type = "n")
text(swiss.mds$points, labels = as.character(1:nrow(swiss.x)))
swiss.sh <- Shepard(swiss.dist, swiss.mds$points)
plot(swiss.sh, pch = ".")
lines(swiss.sh$x, swiss.sh$yf, type = "S")



cleanEx()
nameEx("kde2d")
### * kde2d

flush(stderr()); flush(stdout())

### Name: kde2d
### Title: Two-Dimensional Kernel Density Estimation
### Aliases: kde2d
### Keywords: dplot

### ** Examples

attach(geyser)
plot(duration, waiting, xlim = c(0.5,6), ylim = c(40,100))
f1 <- kde2d(duration, waiting, n = 50, lims = c(0.5, 6, 40, 100))
image(f1, zlim = c(0, 0.05))
f2 <- kde2d(duration, waiting, n = 50, lims = c(0.5, 6, 40, 100),
            h = c(width.SJ(duration), width.SJ(waiting)) )
image(f2, zlim = c(0, 0.05))
persp(f2, phi = 30, theta = 20, d = 5)

plot(duration[-272], duration[-1], xlim = c(0.5, 6),
     ylim = c(1, 6),xlab = "previous duration", ylab = "duration")
f1 <- kde2d(duration[-272], duration[-1],
            h = rep(1.5, 2), n = 50, lims = c(0.5, 6, 0.5, 6))
contour(f1, xlab = "previous duration",
        ylab = "duration", levels  =  c(0.05, 0.1, 0.2, 0.4) )
f1 <- kde2d(duration[-272], duration[-1],
            h = rep(0.6, 2), n = 50, lims = c(0.5, 6, 0.5, 6))
contour(f1, xlab = "previous duration",
        ylab = "duration", levels  =  c(0.05, 0.1, 0.2, 0.4) )
f1 <- kde2d(duration[-272], duration[-1],
            h = rep(0.4, 2), n = 50, lims = c(0.5, 6, 0.5, 6))
contour(f1, xlab = "previous duration",
        ylab = "duration", levels  =  c(0.05, 0.1, 0.2, 0.4) )
detach("geyser")



cleanEx()
nameEx("lda")
### * lda

flush(stderr()); flush(stdout())

### Name: lda
### Title: Linear Discriminant Analysis
### Aliases: lda lda.default lda.data.frame lda.formula lda.matrix
###   model.frame.lda print.lda coef.lda
### Keywords: multivariate

### ** Examples

Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
                   Sp = rep(c("s","c","v"), rep(50,3)))
train <- sample(1:150, 75)
table(Iris$Sp[train])
## your answer may differ
##  c  s  v
## 22 23 30
z <- lda(Sp ~ ., Iris, prior = c(1,1,1)/3, subset = train)
predict(z, Iris[-train, ])$class
##  [1] s s s s s s s s s s s s s s s s s s s s s s s s s s s c c c
## [31] c c c c c c c v c c c c v c c c c c c c c c c c c v v v v v
## [61] v v v v v v v v v v v v v v v
(z1 <- update(z, . ~ . - Petal.W.))



cleanEx()
nameEx("leuk")
### * leuk

flush(stderr()); flush(stdout())

### Name: leuk
### Title: Survival Times and White Blood Counts for Leukaemia Patients
### Aliases: leuk
### Keywords: datasets

### ** Examples

library(survival)
plot(survfit(Surv(time) ~ ag, data = leuk), lty = 2:3, col = 2:3)

# now Cox models
leuk.cox <- coxph(Surv(time) ~ ag + log(wbc), leuk)
summary(leuk.cox)



cleanEx()
nameEx("lm.ridge")
### * lm.ridge

flush(stderr()); flush(stdout())

### Name: lm.ridge
### Title: Ridge Regression
### Aliases: lm.ridge plot.ridgelm print.ridgelm select select.ridgelm
### Keywords: models

### ** Examples

longley # not the same as the S-PLUS dataset
names(longley)[1] <- "y"
lm.ridge(y ~ ., longley)
plot(lm.ridge(y ~ ., longley,
              lambda = seq(0,0.1,0.001)))
select(lm.ridge(y ~ ., longley,
               lambda = seq(0,0.1,0.0001)))



cleanEx()
nameEx("loglm")
### * loglm

flush(stderr()); flush(stdout())

### Name: loglm
### Title: Fit Log-Linear Models by Iterative Proportional Scaling
### Aliases: loglm
### Keywords: category models

### ** Examples

# The data frames  Cars93, minn38 and quine are available
# in the MASS package.

# Case 1: frequencies specified as an array.
sapply(minn38, function(x) length(levels(x)))
## hs phs fol sex f
##  3   4   7   2 0
##minn38a <- array(0, c(3,4,7,2), lapply(minn38[, -5], levels))
##minn38a[data.matrix(minn38[,-5])] <- minn38$f

## or more simply
minn38a <- xtabs(f ~ ., minn38)

fm <- loglm(~ 1 + 2 + 3 + 4, minn38a)  # numerals as names.
deviance(fm)
## [1] 3711.9
fm1 <- update(fm, .~.^2)
fm2 <- update(fm, .~.^3, print = TRUE)
## 5 iterations: deviation 0.075
anova(fm, fm1, fm2)

# Case 1. An array generated with xtabs.

loglm(~ Type + Origin, xtabs(~ Type + Origin, Cars93))

# Case 2.  Frequencies given as a vector in a data frame
names(quine)
## [1] "Eth"  "Sex"  "Age"  "Lrn"  "Days"
fm <- loglm(Days ~ .^2, quine)
gm <- glm(Days ~ .^2, poisson, quine)  # check glm.
c(deviance(fm), deviance(gm))          # deviances agree
## [1] 1368.7 1368.7
c(fm$df, gm$df)                        # resid df do not!
c(fm$df, gm$df.residual)               # resid df do not!
## [1] 127 128
# The loglm residual degrees of freedom is wrong because of
# a non-detectable redundancy in the model matrix.



cleanEx()
nameEx("logtrans")
### * logtrans

flush(stderr()); flush(stdout())

### Name: logtrans
### Title: Estimate log Transformation Parameter
### Aliases: logtrans logtrans.formula logtrans.lm logtrans.default
### Keywords: regression models hplot

### ** Examples

logtrans(Days ~ Age*Sex*Eth*Lrn, data = quine,
         alpha = seq(0.75, 6.5, len=20))



cleanEx()
nameEx("lqs")
### * lqs

flush(stderr()); flush(stdout())

### Name: lqs
### Title: Resistant Regression
### Aliases: lqs lqs.formula lqs.default lmsreg ltsreg
### Keywords: models robust

### ** Examples

set.seed(123) # make reproducible
lqs(stack.loss ~ ., data = stackloss)
lqs(stack.loss ~ ., data = stackloss, method = "S", nsamp = "exact")



cleanEx()
nameEx("mca")
### * mca

flush(stderr()); flush(stdout())

### Name: mca
### Title: Multiple Correspondence Analysis
### Aliases: mca print.mca
### Keywords: category multivariate

### ** Examples

farms.mca <- mca(farms, abbrev=TRUE)
farms.mca
plot(farms.mca)



cleanEx()
nameEx("menarche")
### * menarche

flush(stderr()); flush(stdout())

### Name: menarche
### Title: Age of Menarche in Warsaw
### Aliases: menarche
### Keywords: datasets

### ** Examples

mprob <- glm(cbind(Menarche, Total - Menarche) ~ Age,
             binomial(link = probit), data = menarche)



cleanEx()
nameEx("motors")
### * motors

flush(stderr()); flush(stdout())

### Name: motors
### Title: Accelerated Life Testing of Motorettes
### Aliases: motors
### Keywords: datasets

### ** Examples

library(survival)
plot(survfit(Surv(time, cens) ~ factor(temp), motors), conf.int = FALSE)
# fit Weibull model
motor.wei <- survreg(Surv(time, cens) ~ temp, motors)
summary(motor.wei)
# and predict at 130C
unlist(predict(motor.wei, data.frame(temp=130), se.fit = TRUE))

motor.cox <- coxph(Surv(time, cens) ~ temp, motors)
summary(motor.cox)
# predict at temperature 200
plot(survfit(motor.cox, newdata = data.frame(temp=200),
     conf.type = "log-log"))
summary( survfit(motor.cox, newdata = data.frame(temp=130)) )



cleanEx()
nameEx("muscle")
### * muscle

flush(stderr()); flush(stdout())

### Name: muscle
### Title: Effect of Calcium Chloride on Muscle Contraction in Rat Hearts
### Aliases: muscle
### Keywords: datasets

### ** Examples

A <- model.matrix(~ Strip - 1, data=muscle)
rats.nls1 <- nls(log(Length) ~ cbind(A, rho^Conc),
   data = muscle, start = c(rho=0.1), algorithm="plinear")
(B <- coef(rats.nls1))

st <- list(alpha = B[2:22], beta = B[23], rho = B[1])
(rats.nls2 <- nls(log(Length) ~ alpha[Strip] + beta*rho^Conc,
                  data = muscle, start = st))

Muscle <- with(muscle, {
Muscle <- expand.grid(Conc = sort(unique(Conc)), Strip = levels(Strip))
Muscle$Yhat <- predict(rats.nls2, Muscle)
Muscle <- cbind(Muscle, logLength = rep(as.numeric(NA), 126))
ind <- match(paste(Strip, Conc),
            paste(Muscle$Strip, Muscle$Conc))
Muscle$logLength[ind] <- log(Length)
Muscle})

lattice::xyplot(Yhat ~ Conc | Strip, Muscle, as.table = TRUE,
   ylim = range(c(Muscle$Yhat, Muscle$logLength), na.rm = TRUE),
   subscripts = TRUE, xlab = "Calcium Chloride concentration (mM)",
   ylab = "log(Length in mm)", panel =
   function(x, y, subscripts, ...) {
      panel.xyplot(x, Muscle$logLength[subscripts], ...)
      llines(spline(x, y))
   })



cleanEx()
nameEx("mvrnorm")
### * mvrnorm

flush(stderr()); flush(stdout())

### Name: mvrnorm
### Title: Simulate from a Multivariate Normal Distribution
### Aliases: mvrnorm
### Keywords: distribution multivariate

### ** Examples

Sigma <- matrix(c(10,3,3,2),2,2)
Sigma
var(mvrnorm(n = 1000, rep(0, 2), Sigma))
var(mvrnorm(n = 1000, rep(0, 2), Sigma, empirical = TRUE))



cleanEx()
nameEx("negative.binomial")
### * negative.binomial

flush(stderr()); flush(stdout())

### Name: negative.binomial
### Title: Family function for Negative Binomial GLMs
### Aliases: negative.binomial
### Keywords: regression models

### ** Examples

# Fitting a Negative Binomial model to the quine data
#   with theta = 2 assumed known.
#
glm(Days ~ .^4, family = negative.binomial(2), data = quine)



cleanEx()
nameEx("nlschools")
### * nlschools

flush(stderr()); flush(stdout())

### Name: nlschools
### Title: Eighth-Grade Pupils in the Netherlands
### Aliases: nlschools
### Keywords: datasets

### ** Examples

## Don't show: 
op <- options(digits=5)
## End(Don't show)
nl1 <- within(nlschools, {
IQave <- tapply(IQ, class, mean)[as.character(class)]
IQ <- IQ - IQave
})
cen <- c("IQ", "IQave", "SES")
nl1[cen] <- scale(nl1[cen], center = TRUE, scale = FALSE)

nl.lme <- nlme::lme(lang ~ IQ*COMB + IQave + SES,
                    random = ~ IQ | class, data = nl1)
summary(nl.lme)
## Don't show: 
options(op)
## End(Don't show)



cleanEx()
nameEx("npk")
### * npk

flush(stderr()); flush(stdout())

### Name: npk
### Title: Classical N, P, K Factorial Experiment
### Aliases: npk
### Keywords: datasets

### ** Examples

options(contrasts = c("contr.sum", "contr.poly"))
npk.aov <- aov(yield ~ block + N*P*K, npk)
npk.aov
summary(npk.aov)
alias(npk.aov)
coef(npk.aov)
options(contrasts = c("contr.treatment", "contr.poly"))
npk.aov1 <- aov(yield ~ block + N + K, data = npk)
summary.lm(npk.aov1)
se.contrast(npk.aov1, list(N=="0", N=="1"), data = npk)
model.tables(npk.aov1, type = "means", se = TRUE)



base::options(contrasts = c(unordered = "contr.treatment",ordered = "contr.poly"))
cleanEx()
nameEx("oats")
### * oats

flush(stderr()); flush(stdout())

### Name: oats
### Title: Data from an Oats Field Trial
### Aliases: oats
### Keywords: datasets

### ** Examples

oats$Nf <- ordered(oats$N, levels = sort(levels(oats$N)))
oats.aov <- aov(Y ~ Nf*V + Error(B/V), data = oats, qr = TRUE)
summary(oats.aov)
summary(oats.aov, split = list(Nf=list(L=1, Dev=2:3)))
par(mfrow = c(1,2), pty = "s")
plot(fitted(oats.aov[[4]]), studres(oats.aov[[4]]))
abline(h = 0, lty = 2)
oats.pr <- proj(oats.aov)
qqnorm(oats.pr[[4]][,"Residuals"], ylab = "Stratum 4 residuals")
qqline(oats.pr[[4]][,"Residuals"])

par(mfrow = c(1,1), pty = "m")
oats.aov2 <- aov(Y ~ N + V + Error(B/V), data = oats, qr = TRUE)
model.tables(oats.aov2, type = "means", se = TRUE)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("parcoord")
### * parcoord

flush(stderr()); flush(stdout())

### Name: parcoord
### Title: Parallel Coordinates Plot
### Aliases: parcoord
### Keywords: hplot

### ** Examples

parcoord(state.x77[, c(7, 4, 6, 2, 5, 3)])

ir <- rbind(iris3[,,1], iris3[,,2], iris3[,,3])
parcoord(log(ir)[, c(3, 4, 2, 1)], col = 1 + (0:149)%/%50)



cleanEx()
nameEx("petrol")
### * petrol

flush(stderr()); flush(stdout())

### Name: petrol
### Title: N. L. Prater's Petrol Refinery Data
### Aliases: petrol
### Keywords: datasets

### ** Examples

library(nlme)
Petrol <- petrol
Petrol[, 2:5] <- scale(as.matrix(Petrol[, 2:5]), scale = FALSE)
pet3.lme <- lme(Y ~ SG + VP + V10 + EP,
                random = ~ 1 | No, data = Petrol)
pet3.lme <- update(pet3.lme, method = "ML")
pet4.lme <- update(pet3.lme, fixed = Y ~ V10 + EP)
anova(pet4.lme, pet3.lme)



cleanEx()
nameEx("plot.mca")
### * plot.mca

flush(stderr()); flush(stdout())

### Name: plot.mca
### Title: Plot Method for Objects of Class 'mca'
### Aliases: plot.mca
### Keywords: hplot multivariate

### ** Examples

plot(mca(farms, abbrev = TRUE))



cleanEx()
nameEx("plot.profile")
### * plot.profile

flush(stderr()); flush(stdout())

### Name: plot.profile
### Title: Plotting Functions for 'profile' Objects
### Aliases: plot.profile pairs.profile
### Keywords: models hplot

### ** Examples

## see ?profile.glm for an example using glm fits.

## a version of example(profile.nls) from R >= 2.8.0
fm1 <- nls(demand ~ SSasympOrig(Time, A, lrc), data = BOD)
pr1 <- profile(fm1, alpha = 0.1)
MASS:::plot.profile(pr1)
pairs(pr1) # a little odd since the parameters are highly correlated

## an example from ?nls
x <- -(1:100)/10
y <- 100 + 10 * exp(x / 2) + rnorm(x)/10
nlmod <- nls(y ~  Const + A * exp(B * x), start=list(Const=100, A=10, B=1))
pairs(profile(nlmod))



cleanEx()
nameEx("polr")
### * polr

flush(stderr()); flush(stdout())

### Name: polr
### Title: Ordered Logistic or Probit Regression
### Aliases: polr
### Keywords: models

### ** Examples

options(contrasts = c("contr.treatment", "contr.poly"))
house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
house.plr
summary(house.plr, digits = 3)
## slightly worse fit from
summary(update(house.plr, method = "probit", Hess = TRUE), digits = 3)
## although it is not really appropriate, can fit
summary(update(house.plr, method = "loglog", Hess = TRUE), digits = 3)
summary(update(house.plr, method = "cloglog", Hess = TRUE), digits = 3)

predict(house.plr, housing, type = "p")
addterm(house.plr, ~.^2, test = "Chisq")
house.plr2 <- stepAIC(house.plr, ~.^2)
house.plr2$anova
anova(house.plr, house.plr2)

house.plr <- update(house.plr, Hess=TRUE)
pr <- profile(house.plr)
confint(pr)
plot(pr)
pairs(pr)



base::options(contrasts = c(unordered = "contr.treatment",ordered = "contr.poly"))
cleanEx()
nameEx("predict.glmmPQL")
### * predict.glmmPQL

flush(stderr()); flush(stdout())

### Name: predict.glmmPQL
### Title: Predict Method for glmmPQL Fits
### Aliases: predict.glmmPQL
### Keywords: models

### ** Examples

fit <- glmmPQL(y ~ trt + I(week > 2), random = ~1 |  ID,
               family = binomial, data = bacteria)
predict(fit, bacteria, level = 0, type="response")
predict(fit, bacteria, level = 1, type="response")



cleanEx()
nameEx("predict.lda")
### * predict.lda

flush(stderr()); flush(stdout())

### Name: predict.lda
### Title: Classify Multivariate Observations by Linear Discrimination
### Aliases: predict.lda
### Keywords: multivariate

### ** Examples

tr <- sample(1:50, 25)
train <- rbind(iris3[tr,,1], iris3[tr,,2], iris3[tr,,3])
test <- rbind(iris3[-tr,,1], iris3[-tr,,2], iris3[-tr,,3])
cl <- factor(c(rep("s",25), rep("c",25), rep("v",25)))
z <- lda(train, cl)
predict(z, test)$class



cleanEx()
nameEx("predict.lqs")
### * predict.lqs

flush(stderr()); flush(stdout())

### Name: predict.lqs
### Title: Predict from an lqs Fit
### Aliases: predict.lqs
### Keywords: models

### ** Examples

set.seed(123)
fm <- lqs(stack.loss ~ ., data = stackloss, method = "S", nsamp = "exact")
predict(fm, stackloss)



cleanEx()
nameEx("predict.qda")
### * predict.qda

flush(stderr()); flush(stdout())

### Name: predict.qda
### Title: Classify from Quadratic Discriminant Analysis
### Aliases: predict.qda
### Keywords: multivariate

### ** Examples

tr <- sample(1:50, 25)
train <- rbind(iris3[tr,,1], iris3[tr,,2], iris3[tr,,3])
test <- rbind(iris3[-tr,,1], iris3[-tr,,2], iris3[-tr,,3])
cl <- factor(c(rep("s",25), rep("c",25), rep("v",25)))
zq <- qda(train, cl)
predict(zq, test)$class



cleanEx()
nameEx("profile.glm")
### * profile.glm

flush(stderr()); flush(stdout())

### Name: profile.glm
### Title: Method for Profiling glm Objects
### Aliases: profile.glm
### Keywords: regression models

### ** Examples

options(contrasts = c("contr.treatment", "contr.poly"))
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive = 20 - numdead)
budworm.lg <- glm(SF ~ sex*ldose, family = binomial)
pr1 <- profile(budworm.lg)
plot(pr1)
pairs(pr1)



base::options(contrasts = c(unordered = "contr.treatment",ordered = "contr.poly"))
cleanEx()
nameEx("qda")
### * qda

flush(stderr()); flush(stdout())

### Name: qda
### Title: Quadratic Discriminant Analysis
### Aliases: qda qda.data.frame qda.default qda.formula qda.matrix
###   model.frame.qda print.qda
### Keywords: multivariate

### ** Examples

tr <- sample(1:50, 25)
train <- rbind(iris3[tr,,1], iris3[tr,,2], iris3[tr,,3])
test <- rbind(iris3[-tr,,1], iris3[-tr,,2], iris3[-tr,,3])
cl <- factor(c(rep("s",25), rep("c",25), rep("v",25)))
z <- qda(train, cl)
predict(z,test)$class



cleanEx()
nameEx("rational")
### * rational

flush(stderr()); flush(stdout())

### Name: rational
### Title: Rational Approximation
### Aliases: rational .rat
### Keywords: math

### ** Examples

X <- matrix(runif(25), 5, 5)
zapsmall(solve(X, X/5)) # print near-zeroes as zero
rational(solve(X, X/5))



cleanEx()
nameEx("renumerate")
### * renumerate

flush(stderr()); flush(stdout())

### Name: renumerate
### Title: Convert a Formula Transformed by 'denumerate'
### Aliases: renumerate renumerate.formula
### Keywords: models

### ** Examples

denumerate(~(1+2+3)^3 + a/b)
## ~ (.v1 + .v2 + .v3)^3 + a/b
renumerate(.Last.value)
## ~ (1 + 2 + 3)^3 + a/b



cleanEx()
nameEx("rlm")
### * rlm

flush(stderr()); flush(stdout())

### Name: rlm
### Title: Robust Fitting of Linear Models
### Aliases: rlm rlm.default rlm.formula print.rlm predict.rlm psi.bisquare
###   psi.hampel psi.huber
### Keywords: models robust

### ** Examples

summary(rlm(stack.loss ~ ., stackloss))
rlm(stack.loss ~ ., stackloss, psi = psi.hampel, init = "lts")
rlm(stack.loss ~ ., stackloss, psi = psi.bisquare)



cleanEx()
nameEx("rms.curv")
### * rms.curv

flush(stderr()); flush(stdout())

### Name: rms.curv
### Title: Relative Curvature Measures for Non-Linear Regression
### Aliases: rms.curv print.rms.curv
### Keywords: nonlinear

### ** Examples

# The treated sample from the Puromycin data
mmcurve <- deriv3(~ Vm * conc/(K + conc), c("Vm", "K"),
                  function(Vm, K, conc) NULL)
Treated <- Puromycin[Puromycin$state == "treated", ]
(Purfit1 <- nls(rate ~ mmcurve(Vm, K, conc), data = Treated,
                start = list(Vm=200, K=0.1)))
rms.curv(Purfit1)
##Parameter effects: c^theta x sqrt(F) = 0.2121
##        Intrinsic: c^iota  x sqrt(F) = 0.092



cleanEx()
nameEx("rnegbin")
### * rnegbin

flush(stderr()); flush(stdout())

### Name: rnegbin
### Title: Simulate Negative Binomial Variates
### Aliases: rnegbin
### Keywords: distribution

### ** Examples

# Negative Binomials with means fitted(fm) and theta = 4.5
fm <- glm.nb(Days ~ ., data = quine)
dummy <- rnegbin(fitted(fm), theta = 4.5)



cleanEx()
nameEx("sammon")
### * sammon

flush(stderr()); flush(stdout())

### Name: sammon
### Title: Sammon's Non-Linear Mapping
### Aliases: sammon
### Keywords: multivariate

### ** Examples

swiss.x <- as.matrix(swiss[, -1])
swiss.sam <- sammon(dist(swiss.x))
plot(swiss.sam$points, type = "n")
text(swiss.sam$points, labels = as.character(1:nrow(swiss.x)))



cleanEx()
nameEx("stepAIC")
### * stepAIC

flush(stderr()); flush(stdout())

### Name: stepAIC
### Title: Choose a model by AIC in a Stepwise Algorithm
### Aliases: stepAIC extractAIC.gls terms.gls extractAIC.lme terms.lme
### Keywords: models

### ** Examples

quine.hi <- aov(log(Days + 2.5) ~ .^4, quine)
quine.nxt <- update(quine.hi, . ~ . - Eth:Sex:Age:Lrn)
quine.stp <- stepAIC(quine.nxt,
    scope = list(upper = ~Eth*Sex*Age*Lrn, lower = ~1),
    trace = FALSE)
quine.stp$anova

cpus1 <- cpus
for(v in names(cpus)[2:7])
  cpus1[[v]] <- cut(cpus[[v]], unique(quantile(cpus[[v]])),
                    include.lowest = TRUE)
cpus0 <- cpus1[, 2:8]  # excludes names, authors' predictions
cpus.samp <- sample(1:209, 100)
cpus.lm <- lm(log10(perf) ~ ., data = cpus1[cpus.samp,2:8])
cpus.lm2 <- stepAIC(cpus.lm, trace = FALSE)
cpus.lm2$anova

example(birthwt)
birthwt.glm <- glm(low ~ ., family = binomial, data = bwt)
birthwt.step <- stepAIC(birthwt.glm, trace = FALSE)
birthwt.step$anova
birthwt.step2 <- stepAIC(birthwt.glm, ~ .^2 + I(scale(age)^2)
    + I(scale(lwt)^2), trace = FALSE)
birthwt.step2$anova

quine.nb <- glm.nb(Days ~ .^4, data = quine)
quine.nb2 <- stepAIC(quine.nb)
quine.nb2$anova



cleanEx()
nameEx("summary.negbin")
### * summary.negbin

flush(stderr()); flush(stdout())

### Name: summary.negbin
### Title: Summary Method Function for Objects of Class 'negbin'
### Aliases: summary.negbin print.summary.negbin
### Keywords: models

### ** Examples

summary(glm.nb(Days ~ Eth*Age*Lrn*Sex, quine, link = log))



cleanEx()
nameEx("summary.rlm")
### * summary.rlm

flush(stderr()); flush(stdout())

### Name: summary.rlm
### Title: Summary Method for Robust Linear Models
### Aliases: summary.rlm print.summary.rlm
### Keywords: robust

### ** Examples

summary(rlm(calls ~ year, data = phones, maxit = 50))



cleanEx()
nameEx("theta.md")
### * theta.md

flush(stderr()); flush(stdout())

### Name: theta.md
### Title: Estimate theta of the Negative Binomial
### Aliases: theta.md theta.ml theta.mm
### Keywords: models

### ** Examples

quine.nb <- glm.nb(Days ~ .^2, data = quine)
theta.md(quine$Days, fitted(quine.nb), dfr = df.residual(quine.nb))
theta.ml(quine$Days, fitted(quine.nb))
theta.mm(quine$Days, fitted(quine.nb), dfr = df.residual(quine.nb))

## weighted example
yeast <- data.frame(cbind(numbers = 0:5, fr = c(213, 128, 37, 18, 3, 1)))
fit <- glm.nb(numbers ~ 1, weights = fr, data = yeast)
summary(fit)
mu <- fitted(fit)
theta.md(yeast$numbers, mu, dfr = 399, weights = yeast$fr)
theta.ml(yeast$numbers, mu, limit = 15, weights = yeast$fr)
theta.mm(yeast$numbers, mu, dfr = 399, weights = yeast$fr)



cleanEx()
nameEx("ucv")
### * ucv

flush(stderr()); flush(stdout())

### Name: ucv
### Title: Unbiased Cross-Validation for Bandwidth Selection
### Aliases: ucv
### Keywords: dplot

### ** Examples

ucv(geyser$duration)



cleanEx()
nameEx("waders")
### * waders

flush(stderr()); flush(stdout())

### Name: waders
### Title: Counts of Waders at 15 Sites in South Africa
### Aliases: waders
### Keywords: datasets

### ** Examples

plot(corresp(waders, nf=2))



cleanEx()
nameEx("whiteside")
### * whiteside

flush(stderr()); flush(stdout())

### Name: whiteside
### Title: House Insulation: Whiteside's Data
### Aliases: whiteside
### Keywords: datasets

### ** Examples

require(lattice)
xyplot(Gas ~ Temp | Insul, whiteside, panel =
  function(x, y, ...) {
    panel.xyplot(x, y, ...)
    panel.lmline(x, y, ...)
  }, xlab = "Average external temperature (deg. C)",
  ylab = "Gas consumption  (1000 cubic feet)", aspect = "xy",
  strip = function(...) strip.default(..., style = 1))

gasB <- lm(Gas ~ Temp, whiteside, subset = Insul=="Before")
gasA <- update(gasB, subset = Insul=="After")
summary(gasB)
summary(gasA)
gasBA <- lm(Gas ~ Insul/Temp - 1, whiteside)
summary(gasBA)

gasQ <- lm(Gas ~ Insul/(Temp + I(Temp^2)) - 1, whiteside)
coef(summary(gasQ))

gasPR <- lm(Gas ~ Insul + Temp, whiteside)
anova(gasPR, gasBA)
options(contrasts = c("contr.treatment", "contr.poly"))
gasBA1 <- lm(Gas ~ Insul*Temp, whiteside)
coef(summary(gasBA1))



base::options(contrasts = c(unordered = "contr.treatment",ordered = "contr.poly"))
cleanEx()
nameEx("width.SJ")
### * width.SJ

flush(stderr()); flush(stdout())

### Name: width.SJ
### Title: Bandwidth Selection by Pilot Estimation of Derivatives
### Aliases: width.SJ
### Keywords: dplot

### ** Examples

width.SJ(geyser$duration, method = "dpi")
width.SJ(geyser$duration)

width.SJ(galaxies, method = "dpi")
width.SJ(galaxies)



cleanEx()
nameEx("wtloss")
### * wtloss

flush(stderr()); flush(stdout())

### Name: wtloss
### Title: Weight Loss Data from an Obese Patient
### Aliases: wtloss
### Keywords: datasets

### ** Examples

wtloss.fm <- nls(Weight ~ b0 + b1*2^(-Days/th),
    data = wtloss, start = list(b0=90, b1=95, th=120))
wtloss.fm
plot(wtloss)
with(wtloss, lines(Days, fitted(wtloss.fm)))



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
