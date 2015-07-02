context("linear algebra")
library(Matrix)
A <- sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
cholPermute(A)

test_that("cholsolve works as expected", {
    expect_identical(FRK:::cholPermute(A)$Qpermchol,t(chol(A)))
})
