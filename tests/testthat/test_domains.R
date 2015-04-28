S <- sphere(radius = 1)
x <- matrix(c(-90,90,0,0),2,2)

test_that("dimensions work", {
    expect_equal(dimensions(sphere()), 2L)
    expect_equal(dimensions(plane()), 2L)
    expect_equal(dimensions(real_line()), 1L)
    expect_equal(dimensions(domain(sphere())), 2L)
})

test_that("sphere distance works", {
    expect_identical(distance(S,x,x)[1,2], pi)
})
