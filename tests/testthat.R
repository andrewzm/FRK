ops <- options("crayon.enabled" = FALSE)
library(testthat)
library("covr")

test_check("FRK")
options(ops)
