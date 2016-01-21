print("Entering test_domains.R")
library(sp)

test_that("dimensions work", {
    S <- sphere(radius = 1)
    x <- matrix(c(-90,90,0,0),2,2)
    expect_equal(dimensions(sphere()), 2L)
    expect_equal(dimensions(plane()), 2L)
    expect_equal(dimensions(real_line()), 1L)
    expect_equal(dimensions(timeline()), 1L)
    expect_equal(dimensions(STplane()), 3L)
    expect_equal(dimensions(STsphere()), 3L)
})

test_that("sphere distance works", {
    S <- sphere(radius = 1)
    x <- matrix(c(-90,90,0,0),2,2)
    expect_identical(distance(S,x,x)[1,2], pi)
})


test_that("sphere space-time distance works",{
  ST <- STsphere(radius=1)
  x <- matrix(c(-90,90,0,0,0,2),2,3)
  expect_identical(distance(ST,x,x)[1,2],sqrt(pi^2 +2^2))
})
