m1 <- sphere(radius = 1)
m2 <- plane()
m3 <- real_line()

mu1 <- mu2 <- matrix(c(1,1),1,2)
mu3 <- matrix(1,1,1)
std <- 1

G1 <- GRBF_wrapper(m1,mu1,std)
G2 <- GRBF_wrapper(m2,mu2,std)
G3 <- GRBF_wrapper(m3,mu3,std)

test_that("basis work with all manifolds", {
    expect_identical(G2(s=matrix(c(0,1),1,2)), as.matrix(dnorm(mean=0,sd=1,x=1)*sqrt(2*pi)))
    expect_identical(G3(s=matrix(0,1,1)), as.matrix(dnorm(mean=0,sd=1,x=1)*sqrt(2*pi)))
    expect_is(G1(s=matrix(0,1,2)),"matrix")
    expect_is(G2(s=matrix(0,1,2)),"matrix")
    expect_is(G3(s=matrix(0,1,1)),"matrix")
})

mu <- matrix(-90 + 180*runif(100),50,2)
mu[,1] <- mu[,1]*2
std <- rep(500,50)
G_basis <- radial_basis(manifold = sphere(),loc = mu,scale=std,type="Gaussian")
test_that("we can show basis locations on sphere", {
    expect_true({show_basis(ggplot(),G_basis); TRUE})
})

mu <- matrix(rnorm(10),5,2)
std <- rep(0.3,5)
G_basis <- radial_basis(manifold = plane(),loc = mu,scale=std,type="Gaussian")
test_that("we can have multiple functions in a Basis object and plot them", {
    expect_equal(G_basis@n, 5)
    expect_equal(nrow(G_basis@pars$df), 5)
    expect_equal(G_basis@fn[[1]](mu)[1,1],1) # Value of basis at mean is 1
    expect_equal(diag(eval_basis(G_basis,s = mu)),rep(1,5))
    expect_true({show_basis(ggplot(),G_basis); TRUE})
})

mu <- matrix(runif(10),10,1)
std <- rep(0.3,10)
G_basis <- radial_basis(manifold = real_line(),loc = mu,scale=std,type="Gaussian")
test_that("we can have multiple functions in a Basis object and plot them", {
    expect_equal(G_basis@n, 10)
    expect_equal(nrow(G_basis@pars$df), 10)
    expect_equal(G_basis@fn[[1]](mu)[1,1],1) # Value of basis at mean is 1
    expect_equal(diag(eval_basis(G_basis,s = mu)),rep(1,10))
    expect_true({show_basis(ggplot(),G_basis,s1min=-1,s1max=2,ds1 = 0.01); TRUE})
})

