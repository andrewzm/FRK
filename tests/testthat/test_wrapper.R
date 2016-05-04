context("BAUs")

test_that("basic_wrapper",{
    library(sp)
    Z <- data.frame(x = runif(1000), y= runif(1000))
    Z$z <- sin(8*Z$x) + cos(8*Z$y) + 0.5*rnorm(100)
    coordinates(Z) = ~x+y
    S <- FRK(f = z~1,
             list(Z),
             cellsize = c(0.1,0.1),
             n_EM = 2)
    Pred <- SRE.predict(SRE_model = S,
                        obs_fs = TRUE)
})
