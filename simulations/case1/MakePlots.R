rm(list=ls())


library(latice)

load("sim1res.rda")

esc <- expand.grid(vars=c(50, 100), prop.miss=c(0.1,0.2,0.3), sigma.noise=c(0.125, 0.250))
