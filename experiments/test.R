# test that
library(HessianScreening)

d <- generateDesign(100, 10)

fit <- lassoPath(d$X, d$y)

write.csv(fit$beta, "results/test.csv")
