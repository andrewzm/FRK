context("linear algebra")
library(Matrix)
A <- sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
X <- cholPermute(A)
y <- matrix(c(1,2),2,1)
b <- y %*% t(y)

test_that("cholsolve works as expected", {
    expect_identical(FRK:::cholPermute(A)$Qpermchol,t(chol(A)))
})

test_that("other linalg functions", {
  expect_equal(tr(A),sum(diag(A)))
  expect_equal(diag2(A,A),diag(A %*% A))
  expect_equal(Takahashi_Davis(A),solve(A))
  expect_equal(Takahashi_Davis(A),Takahashi_Davis(A,cholQp = X$Qpermchol, P = X$P))
  expect_equal(as.numeric(cholsolveAQinvAT(b,X$Qpermchol,X$P)),
               as.numeric(b  %*% solve(A) %*% t(b)))
  expect_is(Takahashi_Davis(A),"dgCMatrix")
  expect_is(amd_test(),"list")
  expect_equal(cholPermute(A),cholPermute(A,method="R"))
  expect_equal(cholsolve(A,b),solve(A) %*% b)
  expect_equal(cholsolve(A,b,perm = T),solve(A) %*% b)
  expect_equal(cholsolve(A,b,perm = T,cholQp = X$Qpermchol,P = X$P),solve(A) %*% b)

})
