print("Entering test_basis.R")

m1 <- sphere(radius = 1)
m2 <- plane()
m3 <- real_line()

mu1 <- mu2 <- matrix(c(1,1),1,2)
mu3 <- matrix(1,1,1)
std <- 1

test_that("basis can be generated", {
    expect_true({ .GRBF_wrapper(m1,mu1,std); TRUE})
    expect_true({ .bisquare_wrapper(m1,mu1,std); TRUE})
    expect_true({ .exp_wrapper(m1,mu1,std); TRUE})
    expect_true({ .Matern32_wrapper(m1,mu1,std); TRUE})
})

G1 <- .GRBF_wrapper(m1,mu1,std)
G2 <- .GRBF_wrapper(m2,mu2,std)
G3 <- .GRBF_wrapper(m3,mu3,std)

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
G_basis <- local_basis(manifold = sphere(),loc = mu,scale=std,type="Gaussian")
test_that("we can show basis locations on sphere", {
    expect_true({show_basis(G_basis); TRUE})
})

mu <- matrix(rnorm(10),5,2)
std <- rep(0.3,5)
G_basis <- local_basis(manifold = plane(),loc = mu,scale=std,type="Gaussian")
test_that("we can have multiple functions in a Basis object on plane and plot them", {
    expect_equal(G_basis@n, 5)
    expect_equal(nrow(G_basis@df), 5)
    expect_equal(G_basis@fn[[1]](mu)[1,1],1) # Value of basis at mean is 1
    expect_equal(diag(eval_basis(G_basis,s = mu)),rep(1,5))
    expect_true({show_basis(G_basis); TRUE})
})

mu <- matrix(runif(10),10,1)
std <- rep(0.3,10)
G_basis <- local_basis(manifold = real_line(),loc = mu,scale=std,type="Gaussian")
test_that("we can have multiple functions in a Basis object on real line and plot them", {
    expect_equal(G_basis@n, 10)
    expect_equal(nrow(G_basis@df), 10)
    expect_equal(G_basis@fn[[1]](mu)[1,1],1) # Value of basis at mean is 1
    expect_equal(diag(eval_basis(G_basis,s = mu)),rep(1,10))
    expect_is(eval_basis(G_basis,s = mu),"Matrix")
    expect_true({show_basis(G_basis); TRUE})
})

test_that("can average basis over polygons in plane", {

    ## Get data
    library(sp)
    data(meuse)
    data(meuse.grid)
    coordinates(meuse) = ~x+y # change into an sp object
    gridded(meuse.grid) = ~x + y
    HexPts <- spsample(meuse.grid, type = "hexagonal", cellsize = 400)
    HexPols <- HexPoints2SpatialPolygons(HexPts)
    HexPols_df <- SpatialPolygonsDataFrame(HexPols,
                                           cbind(over(HexPols,meuse.grid),
                                                 coordinates(HexPts)))
    ## Generate basis functions
    G <- auto_basis(manifold = plane(),data=meuse,nres = 2,prune=10,type = "Gaussian")

        expect_true({eval_basis(G,coordinates(HexPts)); TRUE})
        expect_true({eval_basis(G,HexPols_df); TRUE})
        #plot(as.numeric(S1))
        #lines(as.numeric(S2),col='red')
})

test_that("can get ST basis using time repetition", {
    G_spatial <-  local_basis(manifold = sphere(),
                               loc=matrix(runif(20,min=-90,max=90),10,2),
                               scale=rep(20,10),
                               type="bisquare")
    G_space_time <- sp_to_ST_basis(G_spatial,1:10,manifold=STsphere())
    expect_is(G_space_time,"Basis")
    expect_is(manifold(G_space_time),"STsphere")
    expect_equal(nbasis(G_space_time),100)

})

test_that("can get ST basis using tensor product", {
    G_spatial <-  local_basis(manifold = sphere(),
                               loc=matrix(runif(20,min=-90,max=90),10,2),
                               scale=rep(20,10),
                               type="bisquare")

    G_temporal <- local_basis(manifold=real_line(),loc = matrix(c(2,7,12)),scale = rep(3,3))
    G_spacetime <- TensorP(G_spatial,G_temporal)
    expect_is(G_spacetime,"TensorP_Basis")
    expect_is(G_spacetime@Basis1,"Basis")
    expect_is(G_spacetime@Basis2,"Basis")
    expect_equal(nbasis(G_spacetime),30)

})
