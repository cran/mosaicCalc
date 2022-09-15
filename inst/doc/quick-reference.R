## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(mosaic)
library(ggformula)
library(mosaicCalc)
library(palmerpenguins)
library(knitr)

## ----echo=FALSE, results="hide"-----------------------------------------------
g <- makeFun(2 + 3*x - 7*x^2 ~ .)

## ----echo=FALSE---------------------------------------------------------------
g

## ----echo=TRUE----------------------------------------------------------------
kable(penguins)

## ----echo=TRUE----------------------------------------------------------------
gf_point(flipper_length_mm ~ body_mass_g, data = palmerpenguins::penguins)

## ----echo=TRUE----------------------------------------------------------------
# linear axes: the default
gf_point(pressure ~ volume, data = Boyle) 

# Log-log axes
gf_point(pressure ~ volume, data = Boyle) %>%
  gf_refine(scale_x_log10(), scale_y_log10())

## ----echo=TRUE----------------------------------------------------------------
gf_point(flipper_length_mm ~ body_mass_g, 
         data = palmerpenguins::penguins,
         color = ~ sex)

## ----echo=TRUE----------------------------------------------------------------
gf_point(flipper_length_mm ~ body_mass_g, 
         data = palmerpenguins::penguins,
         color = ~ sex) %>%
  gf_lims(y=c(0,235), x=c(0,6500))


## ----echo=TRUE----------------------------------------------------------------
mod <- fitModel(pressure ~ a*volume^n, data = Boyle)

## ----echo=TRUE----------------------------------------------------------------
gf_point(pressure ~ volume, data = Boyle) %>%
  slice_plot(mod(volume) ~ volume, color="magenta")

## ----echo=TRUE----------------------------------------------------------------
coef(mod)

## ----echo=TRUE----------------------------------------------------------------
dist <- makeFun(v0*(t-t0) + g*(t-t0)^2 / 2 ~ ., g=-9.8, v0=0, t0=0)

## ----echo=TRUE----------------------------------------------------------------
on_Earth <- dist(2)
on_Mars  <- dist(2, g=-3.7)

## ----echo=TRUE----------------------------------------------------------------
df <- D(exp(t) * cos(t) ~ t)
df
slice_plot(df(t) ~ t, bounds(t=0:10))

## ----echo=TRUE----------------------------------------------------------------
antiD(sin(omega*t) ~ t)

## ----echo=TRUE----------------------------------------------------------------
solutions <- Zeros(sin(x) - 0.5 ~ x, bounds(x=0:10))
slice_plot(sin(x) ~ x, bounds(x=0:10)) %>%
  gf_point(0.5 ~ x, data = solutions)
solutions

## ----echo=TRUE----------------------------------------------------------------
f <- rfun(~ x, seed=943) # a random function
slice_plot(f(x) ~ x, bounds(x=-5:5))

## ----echo=TRUE----------------------------------------------------------------
df <- D(f(x) ~ x)

## ----echo=TRUE----------------------------------------------------------------
dzeros <- Zeros(df(x) ~ x, bounds(x=-5:5))
dzeros

## ----echo=TRUE----------------------------------------------------------------
ddf <- D(f(x) ~ x & x)

## ----echo=TRUE----------------------------------------------------------------
dzeros <- dzeros %>%
  mutate(val = f(x), convexity = sign(ddf(x)))
dzeros

## ----echo=TRUE----------------------------------------------------------------
slice_plot(f(x) ~ x, bounds(x=-5:5)) %>%
  gf_point(val ~ x, data = dzeros, color= ~ convexity)

## ----echo=TRUE----------------------------------------------------------------
argM(f(x) ~ x, bounds(x=-5:5))

## ----echo=TRUE----------------------------------------------------------------
f2 <- makeFun(x*sin(sqrt(3+x))*(1+cos(y))-y ~ .)
solns <- argM(f2(x, y) ~ x & y, bounds(x=c(-3,3), y=c(-3,3)))
solns

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  contour_plot(f2(x,y) ~ x & y, bounds(x=c(-3,3), y=c(-3,3))) %>%
#    gf_point(y ~ x, data = solns, color="red")

## ----create-functions, exercise=TRUE, exercise.cap="Basic modeling functions", exercise.lines=15----
identity_fun  <- makeFun(x ~ x)
constant_fun  <- makeFun(1 ~ x)
straight_line <- makeFun(m*x + b ~ x, b=0, m=1)
exponential   <- makeFun(exp(k * x) ~ x, k=1)
power_law     <- makeFun(x^p ~ x, p = 1/2)
sinusoid      <- makeFun(sin(2*pi*(t-t0)/P) ~ t, P=2, t0=0)
logarithm     <- makeFun(log(x, base=exp(1)) ~ x)
gaussian      <- makeFun(dnorm(x, mean, sd) ~ x, mean=0, sd=1)
sigmoid       <- makeFun(pnorm(x, mean, sd) ~ x, mean=0, sd=1)

identity_fun(3)
constant_fun(3)
power_law(3)

## ----assembling-functions, exercise=TRUE, exercise.cap="Assembling functions", exercise.lines=10----
# Linear combination (example)
f <- makeFun(a0 + a1*exp(k * x) ~ ., a0=30, a1=150, k=-0.5)
# Product (example)
g <- makeFun(dnorm(x, mean=0, sd=3) * sin(2*pi*t/P) ~ ., P=3)
# Composition (example)
h <- makeFun(exp(sin(2*pi*t/P) ~ x) ~ ., P = 3)

## ----graphing-functions, exercise=TRUE,  exercise.cap="Graphing functions", exercise.lines=10, eval=FALSE----
#  slice_plot(dnorm(x, mean=1, sd=2) ~ x,
#             bounds(x=-5:5))
#  contour_plot(dnorm(x, mean=1, sd=2) * pnorm(y, mean=-3, sd=1) ~ x + y,
#               bounds(x=-5:5, y=-5:5))
#  

## ----calculus-ops-der, exercise=TRUE, exercise.cap="Differentiation", exercise.lines=10----
f <- makeFun(exp(-0.5*x) * sin(2*pi*x/3) ~ .)
df <- D(f(x) ~ x)
slice_plot(df(x) ~ x, bounds(x=-5:5)) %>%
  slice_plot(f(x) ~ x, bounds(x=-5:5), color="orange3")

## ----calculus-ops-anti, exercise=TRUE, exercise.cap="Anti-differentiation", exercise.lines=10----
f <- makeFun(dnorm(x, mean=1, sd=2) ~ .)
F <- antiD(f(x) ~ x)
slice_plot(F(x) ~ x, bounds(x=-5:5)) %>%
  slice_plot(f(x) ~ x, color="orange3")
# Set "constant of integration"
slice_plot(F(x) ~ x, bounds(x=-5:5)) %>%
  slice_plot(f(x) ~ x, color="orange3") %>%
  slice_plot(F(x, C=0.25) ~ x, color="green")
# Definite integral
F(5) - F(-5)

## ----calculus-ops-solve, exercise=TRUE, exercise.cap="Zero finding", warning=FALSE, exercise.lines=10----
f <- makeFun(exp(sin(2*pi*x/3)) - 0.5 ~ .)
Zeros <- findZeros(f(x) ~ x, near=0, within=5)
Zeros
slice_plot(f(x) ~ x, bounds(x=-5:5)) %>%
  gf_hline(yintercept=0, color="orange3") %>%
  gf_vline(xintercept= ~ x, color="dodgerblue", data=Zeros)

## ----stans-data, exercise=TRUE, exercise.lines=10, exercise.cap="fitModel"----
gf_point(temp ~ time, data = CoolingWater)
# Eyeball half-life at 25
k0 <- -log(2)/25
mod <- fitModel(temp ~ A + B*exp(-k*time), data=CoolingWater,
                start=list(k=k0))
Plot <- gf_point(temp ~ time, data = CoolingWater) %>%
  slice_plot(mod(time) ~ time, color="dodgerblue", alpha=0.25, size=2) 

## ----linear-quad, exercise=TRUE, exercise.lines=10, warning=FALSE-------------
f <- makeFun(exp(-0.5*x)*sin(2*pi*x/3) ~ .)

df <- D(f(x) ~ x)
center_fun_on <- 0.9
ddf <- D(df(x) ~ x) # alternatively, D(f(x) ~ x + x)
lin_approx <- makeFun(f(x0) + df(x0)*(x-x0) ~ x, x0 = center_fun_on)
quad_approx <- makeFun(lin_approx(x) + 0.5*ddf(x0)*(x-x0)^2 ~ x, x0 = center_fun_on)
slice_plot(f(x) ~ x, bounds(x=0:1.5), size=2) %>%
  slice_plot(lin_approx(x) ~ x, color="blue") %>%
  slice_plot(quad_approx(x) ~ x, bounds(x=0:1.5), color="orange") %>%
  gf_vline(xintercept = center_fun_on, alpha=0.2, color="yellow", size=3)

