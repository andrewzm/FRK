ops <- options("crayon.enabled" = FALSE)
library(testthat)
require("covr")

test_check("FRK")
options(ops)
