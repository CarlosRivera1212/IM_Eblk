options(digits = 22)

# ROTHAMSTED.BRUSSELS ----------------------------------------------------------------
library(agridat)
data(rothamsted.brussels)
dat <- rothamsted.brussels
head(dat)

d = as.matrix(dist(dat[,c('row','col')]))
di = 1/d
diag(di) = 0
We = di/rowSums(di)
dim(We)
# apply(We, 1, sum)

mod1 = aov(yield ~ trt, dat); summary(mod1)
mod2 = aov(yield ~ block + trt, dat); summary(mod2)

library(ape)
Moran.I(mod1$residuals, We)
Moran.I(mod2$residuals, We)


# MEAD.STRAWBERRY ---------------------------------------------------------------
library(agridat)
data(mead.strawberry)
dat <- mead.strawberry
head(dat)
dat$gen = as.factor(dat$gen)
dat$block = as.factor(dat$block)

table(dat$gen, dat$block)

d = as.matrix(dist(dat[,c('row','col')]))
di = 1/d
diag(di) = 0
We = di/rowSums(di)
dim(We)
# apply(We, 1, sum)

mod1 = aov(yield ~ gen, dat); summary(mod1)
mod2 = aov(yield ~ block + gen, dat); summary(mod2)

library(ape)
Moran.I(mod1$residuals, We)
Moran.I(mod2$residuals, We)


# FEDERER.TOBACCO ---------------------------------------------------------
library(agridat)
data(federer.tobacco)
dat <- federer.tobacco
head(dat)

d = as.matrix(dist(dat[,c('row','block')]))
di = 1/d
diag(di) = 0
We = di/rowSums(di)

mod1 = aov(height ~ dose, dat); summary(mod1)
mod2 = aov(height ~ factor(block) + dose, dat); summary(mod2)

library(ape)
Moran.I(mod1$residuals, We)
Moran.I(mod2$residuals, We)


# ROTHAMSTED.OATS --------------------------------------------------------------------
library(agridat)
data(rothamsted.oats)
dat <- rothamsted.oats
head(dat)

d = as.matrix(dist(dat[,c('row','col')]))
di = 1/d
diag(di) = 0
We = di/rowSums(di)

mod1 = aov(grain ~ trt, dat); summary(mod1)
mod2 = aov(grain ~ block + trt, dat); summary(mod2)

library(ape)
Moran.I(mod1$residuals, We)
Moran.I(mod2$residuals, We)


# DRUBAN.COMPETITION ------------------------------------------------------
library(agridat)
data(durban.competition)
# dat <- durban.competition[durban.competition$gen!='B',]
dat <- durban.competition
dat$block = as.factor(dat$block)
head(dat)

d = as.matrix(dist(dat[,c('block','col')]))
di = 1/d
diag(di) = 0
We = di/rowSums(di)

mod1 = aov(yield ~ gen, dat); summary(mod1)
mod2 = aov(yield ~ block + gen, dat); summary(mod2)

library(ape)
Moran.I(mod1$residuals, We)
Moran.I(mod2$residuals, We)

library(ggplot2)
ggplot(dat)+
  aes(block, col, fill = yield)+
  geom_tile()
