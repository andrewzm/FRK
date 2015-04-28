#### GENERIC FUNCTIONS ######
#' @title Subset an MVST object
#' @description This function can be used to subset observations, meshes, and latent fields in the usual manner of \code{subset}.
#' @param x the MVST object
#' @param ... relations pertaining to the data frame within the MVST object
#' @return object of the same class as \code{x} unless \code{x} is a mesh, in which case a data frame is returned.
#' @export
#' @examples
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' icesat_sub <- subset(icesat_obs,t==3)
#' print(class(icesat_sub))
#'
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
#' Mesh_sub <- subset(Mesh,x < 0)
#' print(class(Mesh_sub))
setGeneric("subset")


#' @title Plot meshes and observations
#' @description This function is the general one used for plotting/visualising objects with MVST.
#' @param x an MVST object of class \code{Obs} or \code{Basis_GMRF}
#' @param y a character indicated which column to plot
#' @param ... other parameter used to configure the plot. These include \code{max} (upperbound on colour scale) and \code{min} (lowerbound on colour scale)
#' @return a ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
#' Mesh["z"] <- sin(Mesh["x"]/1000)*cos(Mesh["y"]/1000)
#' g <- plot(Mesh,"z")
#'
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' plot(subset(icesat_obs,t==2),"z",min=-0.5,max=0.5)}
setGeneric("plot")

#' @title Plots an interpolated field from a mesh
#' @docType methods
#' @description This function takes a Mesh and a column name in the mesh, to generate an interpolated
#' field of the indicated column
#' @param x a Mesh object
#' @param y a character indicated which column to plot
#' @param ds the resolution of the plot (defaulted to 40 x 40)
#' @param max upperbound on colour scale
#' @param min lowerbound on colour scale
#' @return a ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
#' Mesh["z"] <- sin(Mesh["x"]/1000)*cos(Mesh["y"]/1000)
#' g <- plot_interp(Mesh,"z",ds=200)}
setGeneric("plot_interp", function(x,y,ds,...) standardGeneric("plot_interp"))

setGeneric("show_basis", function(g,basis,...) standardGeneric("show_basis"))



#' @title Get data frame
#' @description Returns data frame from MVST block object. Can be used on objects of class \code{FEBasis,GMRF} and \code{Obs}. If used on an object of class \code{Graph}, the observation data frame is returned.
#' @param .Object object from which data frame will be extracted
#' @export
#' @examples
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' df <- getDf(icesat_obs)
setGeneric("getDf", function(.Object) standardGeneric("getDf"))

setGeneric("manifold", function(.Object) standardGeneric("manifold"))

#' @title Get precision matrix
#' @description Returns precision matrix from MVST object.
#' @param .Object object from which precision matrix will be extracted
#' @export
#' @examples
#' my_RW <- GMRF_RW(n=10, order=1, precinc =2, name="my_first_RW")
#' Q <- getPrecision(my_RW)
setGeneric("getPrecision", function(.Object) standardGeneric("getPrecision"))

#' @title Set precision matrix
#' @description Updates precision matrix from MVST object.
#' @param .Object object whose precision matrix will be changed
#' @param new precision matrix
#' @export
#' @examples
#' my_RW <- GMRF_RW(n=10, order=1, precinc =2, name="my_first_RW")
#' Q <- getPrecision(my_RW)
#' my_RW2 <- setPrecision(my_RW,Q + Imat(n=10))
#' Q2 <- getPrecision(my_RW2)
setGeneric("setPrecision", function(.Object,Q) standardGeneric("setPrecision"))

#' @title Get mean vector
#' @description Returns the mean from MVST object.
#' @param .Object object from which mean vector will be extracted
#' @export
#' @examples
#' my_RW <- GMRF_RW(n=10, order=1, precinc =2, name="my_first_RW")
#' mu <- getMean(my_RW)
setGeneric("getMean", function(.Object) standardGeneric("getMean"))

#' @title Get mass matrix
#' @description returns the mass matrix of a finite element basis
#' @param B an object of class \code{FEBasis}
#' @export
#' @examples
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
#' M <- mass_matrix(Mesh)
setGeneric("mass_matrix", function(B) standardGeneric("mass_matrix"))

#' @title Get stiffness matrix
#' @description returns the stiffness matrix of a finite element basis
#' @param B an object of class \code{FEBasis}
#' @export
#' @examples
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
#' M <- stiffness_matrix(Mesh)
setGeneric("stiffness_matrix", function(B) standardGeneric("stiffness_matrix"))

setGeneric("distance", function(d,x1,x2) standardGeneric("distance"))

setGeneric("eval_basis", function(basis,s,output="list") standardGeneric("eval_basis"))

#' @title Return size of object Get stiffness matrix
#' @description Returns the size of the object, where the 'size' is determined by the class of the object being passed. For a \code{GMRF} and \code{GMRF_basis}, this is the number of variables of the associated MRF, for an \code{FEBasis} this is the number of vertices in the mesh, and for an \code{Obs} this is the number of observations.
#' @param x an MVST object
#' @export
#' @examples
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' nrow(icesat_obs)
setGeneric("nrow", function(x) standardGeneric("nrow"))

#' @title Return head object data frame
#' @description Short for \code{head(getDf())}. Returns the top of the object data frame which contains most information on the object, see \code{getDf} for details.
#' @param x an MVST object
#' @param ... other parameters passed to \code{base::head}
#' @export
#' @examples
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' head(icesat_obs,n=10)
setGeneric("head", function(x,...) standardGeneric("head"))

#' @title Return tail object data frame
#' @description Short for \code{tail(getDf())}. Returns the tail of the object data frame which contains most information on the object, see \code{getDf} for details.
#' @param x an MVST object
#' @param ... other parameters passed to \code{base::tail}
#' @export
#' @examples
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' tail(icesat_obs,n=10)
setGeneric("tail", function(x,...) standardGeneric("tail"))

#' @title Concatenation
#' @description Concatenates MVST objects of the same class together. This is primarily used to join up \code{GMRF_basis} blocks and \code{Obs} blocks together.
#' @param ... a series of \code{MVST} objects
#' @export
#' @examples
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' icesat_2x <- concat(icesat_obs,icesat_obs)
setGeneric("concat", function(...) standardGeneric("concat"))

#' @title Compress graph
#' @description This function takes am object of class \code{Graph} and compresses it into a one-layer network of class \code{Graph_2nodes}.
#' The latter object can be then passed to \code{Infer} for a standard Gaussian update.
#' @param Graph object of class \code{Graph}.
#' @return Object of class \code{Graph_2nodes}.
#' @keywords Graph, compress
#' @export
#' @examples
#' \dontrun{
#' require(Matrix)
#' data(icesat)
#' data(surf_fe)
#'
#' ## First create observation object
#' icesat_obs <- Obs(df=icesat,
#'                  abs_lim = 5,
#'                  avr_method = "median",
#'                  box_size=100,
#'                  name="icesat")
#'
#' ## Now create GMRF defined over some FE basis
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
#'
#' mu <- matrix(0,nrow(Mesh),1)
#' Q <- sparseMatrix(i=1:nrow(surf_fe$p), j = 1:nrow(surf_fe$p), x = 1)
#'
#' my_GMRF <- GMRF(mu = mu, Q = Q,name="SURF",t_axis = 0:6)
#' SURF <-GMRF_basis(G = my_GMRF, Basis = Mesh)
#'
#' L1 <- link(SURF,icesat_obs)
#' e <- link_list(list(L1))
#' v <- block_list(list(O = icesat_obs, G = SURF))
#' G <- new("Graph",e=e,v=v)
#' G_reduced <- compress(G)
#' }
setGeneric("compress", function(Graph) standardGeneric("compress"))

#' @title Infer
#' @description Infer is a generic function which carries out inference for a given model. For now only a pure Gaussian model is considered; however
#' one has various options with which to carry out the Gaussian update in order to maximise use of resources, for example on linear combinations of
#' the state-space rather than the whole space/
#' @param Graph object of class \code{Graph_2nodes}.
#' @param SW if 1, the Shermany Woodbury is used for inference over linear combinations, see vignette for details. This option cannot be
#' set if linear combinations are not specified.
#' @param Comb an \code{m} \eqn{\times} \code{n} matrix where each row is a binary vector indicating which of the \code{n} states constitute the linear combination.
#' @return List with fields \code{Graph} (the original graph) and \code{Post_GMRF}. The latter is an object of class \code{GMRF} with mean and precision
#' given by the update. If a set of linear combinations is desired, then the list also contains a field \code{Comb_results}, a list with entries \code{mu}
#' and \code{cov}, the mean and covariance over the linear combinations respectively.
#' @keywords Graph, inference, Gaussian update
#' @export
#' @examples
#' \dontrun{
#' require(Matrix)
#' data(icesat)
#' data(surf_fe)
#'
#' ## First create observation object
#' icesat_obs <- Obs(df=icesat,
#'                  abs_lim = 5,
#'                  avr_method = "median",
#'                  box_size=100,
#'                  name="icesat")
#'
#' ## Now create GMRF defined over some FE basis
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
#'
#' mu <- matrix(0,nrow(Mesh),1)
#' Q <- sparseMatrix(i=1:nrow(surf_fe$p), j = 1:nrow(surf_fe$p), x = 1)
#'
#' my_GMRF <- GMRF(mu = mu, Q = Q,name="SURF",t_axis = 0:6)
#' SURF <-GMRF_basis(G = my_GMRF, Basis = Mesh)
#'
#' L1 <- link(SURF,icesat_obs)
#' e <- link_list(list(L1))
#' v <- block_list(list(O = icesat_obs, G = SURF))
#' G <- new("Graph",e=e,v=v)
#' G_reduced <- compress(G)
#' Results <- Infer(G_reduced)
#' }
setGeneric("Infer", function(Graph,...) standardGeneric("Infer"))

#' @title Set diffusion parameter on observations.
#' @description Accounts for averaging over observations, to account for signal leakage with nearby observations. See details.
#' @param .Object Object of class \code{Obs}.
#' @param alpha The diffusion coefficient. See details.
#' @param av_dist The distance within which observations are assumed to be coupled.
#' @return Object of class \code{Obj} with updated diffusion parameter \code{alpha0} and averaging matrix \code{P}.
#' @details  This function sets the diffusion parameter \eqn{\alpha} in the model \deqn{z = P(\alpha)y + e} where the matrix \eqn{P(\alpha)} is defined
#' as \deqn{P^{(i,j)} = \left\{ \begin{array}{ll}
#'                      \alpha, & i \sim j, (n_i)\theta > 0.9 \\ [2ex]
#'                                          1 - (n_i)\theta & i = j, (n_i)\theta > 0.9 \\ [2ex]
#'                                               \displaystyle \frac{1}{(n_i)+1} & \textrm{otherwise}
#'                                                                       \end{array} \right.
#'                                                                       }
#' where \eqn{n_i} denotes the number of neighbours of observation \eqn{i} and \eqn{\sim} denotes `neighbour of'.
#' The matrix describes the proportion of signal \eqn{n_i\alpha} which should be attributed to the spatial regions associated with the
#' neighbouring observations. If \eqn{n_i\theta} exceeds 0.9 (indicative of poor localisation), the observation
#' is assumed to be an equal average of itself and its neighbours. Two (observations are assigned as neighbours if their geometric
#' centres are distanced by less than \code{av_dist}.
#' @keywords Observation, average
#' @export
#' @examples
#' # Create three polygon 'footprints'
#' pol_df <- rbind(data.frame(id=1,x1=0,x2=0,x3=1,x4=1,y1=0,y2=1,y3=1,y4=0,t=0),
#'                 data.frame(id=2,x1=1,x2=1,x3=2,x4=2,y1=1,y2=2,y3=2,y4=1,t=0),
#'                 data.frame(id=3,x1=2,x2=2,x3=3,x4=3,y1=2,y2=3,y3=3,y4=2,t=0))
#' df <- rbind(data.frame(id=1,x=0.5,y=0.5,z=1,std=1,t=0),
#'             data.frame(id=2,x=1.5,y=1.5,z=2,std=1,t=0),
#'             data.frame(id=3,x=2.5,y=2.5,z=1.5,std=1,t=0))
#' O <- Obs_poly(df=df,pol_df=pol_df)
#' plot(O,"z")
#' Odiff <- setalpha(O,0.1,av_dist=2)
setGeneric("setalpha", function(.Object,alpha,av_dist) standardGeneric("setalpha"))

#' @title Set GMRF
#' @description This function takes an object of class \code{Graph_2nodes} and changes the GMRF (process) component accordingly. This can be used to
#' change the mean in of the process component, or the precision matrix when chaining Gaussian updates.
#' @param Graph object of class \code{Graph_2nodes}.
#' @param obj object of class \code{GMRF}.
#' @return Object of class \code{Graph_2nodes}.
#' @keywords Graph, two nodes
#' @export
#' @examples
#' \dontrun{
#' require(Matrix)
#' data(icesat)
#' data(surf_fe)
#'
#' ## First create observation object
#' icesat_obs <- Obs(df=icesat,
#'                  abs_lim = 5,
#'                  avr_method = "median",
#'                  box_size=100,
#'                  name="icesat")
#'
#' ## Now create GMRF defined over some FE basis
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
#'
#' mu <- matrix(0,nrow(Mesh),1)
#' Q <- sparseMatrix(i=1:nrow(surf_fe$p), j = 1:nrow(surf_fe$p), x = 1)
#'
#' my_GMRF <- GMRF(mu = mu, Q = Q,name="SURF",t_axis = 0:6)
#' SURF <-GMRF_basis(G = my_GMRF, Basis = Mesh)
#'
#' L1 <- link(SURF,icesat_obs)
#' e <- link_list(list(L1))
#' v <- block_list(list(O = icesat_obs, G = SURF))
#' G <- new("Graph",e=e,v=v)
#' G_reduced <- compress(G)
#'
#' mu2 <- matrix(1,nrow(Mesh),1)
#' my_GMRF2 <- GMRF(mu = mu2, Q = Q,name="SURF",t_axis = 0:6)
#'
#' G_reduced <- setGMRF(G_reduced, obj = my_GMRF2)
#'
#' }
setGeneric("setGMRF", function(Graph,obj) standardGeneric("setGMRF"))

#' @title Extract class
#' @description This function is used to extract either the observations or the processes from a list of vertices (i.e. a \code{block_list})
#' @param L an object of class \code{block_list}.
#' @param Cl A string identifying which object to extract; can either be \code{process} or \code{Obs}.
#' @return Object of class \code{Obs} or class \code{process}.
#' @keywords list extraction
#' @export
#' @examples
#' \dontrun{
#' require(Matrix)
#' data(icesat)
#' data(surf_fe)
#'
#' ## First create observation object
#' icesat_obs <- Obs(df=icesat,
#'                  abs_lim = 5,
#'                  avr_method = "median",
#'                  box_size=100,
#'                  name="icesat")
#'
#' ## Now create GMRF defined over some FE basis
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
#'
#' mu <- matrix(0,nrow(Mesh),1)
#' Q <- sparseMatrix(i=1:nrow(surf_fe$p), j = 1:nrow(surf_fe$p), x = 1)
#'
#' my_GMRF <- GMRF(mu = mu, Q = Q,name="SURF",t_axis = 0:6)
#' SURF <-GMRF_basis(G = my_GMRF, Basis = Mesh)
#'
#' L1 <- link(SURF,icesat_obs)
#' e <- link_list(list(L1))
#' v <- block_list(list(O = icesat_obs, G = SURF))
#' G <- new("Graph",e=e,v=v)
#' G_reduced <- compress(G)
#'
#' retrieved_proc_from_list <- extractClass(v,"process")
#' }
setGeneric("extractClass", function(L,Cl) standardGeneric("extractClass"))

#' @title Get incidence matrix
#' @description Retrieves the incidence matrix \code{C} from the observation model \eqn{z = Cx + e}.
#' @param .Object object of class \code{Graph_2nodes} or class \code{link}.
#' @return Object of class \code{Matrix}.
#' @keywords incidence matrix, observation model
#' @export
#' @examples
#' \dontrun{
#' require(Matrix)
#' data(icesat)
#' data(surf_fe)
#'
#' ## First create observation object
#' icesat_obs <- Obs(df=icesat,
#'                  abs_lim = 5,
#'                  avr_method = "median",
#'                  box_size=100,
#'                  name="icesat")
#'
#' ## Now create GMRF defined over some FE basis
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
#'
#' mu <- matrix(0,nrow(Mesh),1)
#' Q <- sparseMatrix(i=1:nrow(surf_fe$p), j = 1:nrow(surf_fe$p), x = 1)
#'
#' my_GMRF <- GMRF(mu = mu, Q = Q,name="SURF",t_axis = 0:6)
#' SURF <-GMRF_basis(G = my_GMRF, Basis = Mesh)
#'
#' L1 <- link(SURF,icesat_obs)
#' C_L1 <- getC(L1)
#' e <- link_list(list(L1))
#' v <- block_list(list(O = icesat_obs, G = SURF))
#' G <- new("Graph",e=e,v=v)
#' G_reduced <- compress(G)
#' C_Graph <- getC(G_reduced)
#' }
setGeneric("getC",function(.Object) standardGeneric("getC"))

#' @title Predictive variance for large data sets
#' @description This function takes inferential results from \code{Infer} and computes the predictive variance from the posterior GMRF at locations for a very large number of data points. This function only yields valid results if the non-zero locations in each row of the incidence matrix constructed using the validation data are a subset of those used to generate the posterior precision matrix, see references.
#' @param Results a list generated by the function \code{Infer}
#' @param G An object of class \code{Graph_2nodes} describing the relationship between the observations and the processes.
#' @return Object of class \code{Obs}.
#' @keywords predictive variance, validation
#' @references Jonathan C. Rougier, Andrew Zammit-Mangion and Nana Schoen (2014). Computation and visualisation for large-scale Gaussian updates. \url{http://arxiv.org/abs/1406.5005}.
#' @export
#' @examples
#' \dontrun{
#' require(Matrix)
#' data(icesat)
#' data(surf_fe)
#'
#' ## First create observation object
#' icesat_obs <- Obs(df=icesat,
#'                  abs_lim = 5,
#'                  avr_method = "median",
#'                  box_size=100,
#'                  name="icesat")
#'
#' ## Now create GMRF defined over some FE basis
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
#'
#' mu <- matrix(0,nrow(Mesh),1)
#' Q <- sparseMatrix(i=1:nrow(surf_fe$p), j = 1:nrow(surf_fe$p), x = 1)
#'
#' my_GMRF <- GMRF(mu = mu, Q = Q,name="SURF",t_axis = 0:6)
#' SURF <-GMRF_basis(G = my_GMRF, Basis = Mesh)
#'
#' L1 <- link(SURF,icesat_obs)
#' e <- link_list(list(L1))
#' v <- block_list(list(O = icesat_obs, G = SURF))
#' G <- new("Graph",e=e,v=v)
#' G_reduced <- compress(G)
#' Results <- Infer(G_reduced)
#'
#' Obs_test <- pred_variance_large(Results,G_reduced)
#' }
setGeneric("pred_variance_large", function(Results,G) standardGeneric("pred_variance_large"))

#' @title Split data into training/testing groups
#' @description \code{split_validation} takes an object of class \code{Obs} and returns a list of two objects of class \code{Obs}, where one object
#' contains the training data and the other contains the validation data. If the data is spatio-temporal, the validation data can be chosen to be
#' co-located in space.
#' @param .Object an object of class \code{Obs}.
#' @param samples an integer identifying the number of observations to use for validation.
#' @param common a flag which, if 1, indicates that validation data should be chosen to coincide spatially. This flag is not relevant if the
#' data is not spatio-temporal.
#' @param ... further arguments passed on to \code{subset}.
#' @return a list of two objects of class \code{Obs}.
#' @export
#' @examples
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' O2 <- split_validation(icesat_obs,100,common=0, t > 0)
setGeneric("split_validation", function(.Object,samples,common,...) standardGeneric("split_validation"))


#' @title Predictive statistics for small validation data
#' @description This function takes inferential results from \code{Infer} and computes predictive statistics on a validation data set. These include \code{RMS}, \code{PMCC}, \code{CR1} and \code{CR2} as defined in Sahu and Mardia (2005), and the Mahalanobis distance \code{Mahalanobis}, the eigenvalue errors \code{DE} and the pivoted Cholesky errors \code{DPC} as described in Bastos and O'Hagan (2009).
#' @param Results a list generated by the function \code{Infer}
#' @param G an object of class \code{Graph_2nodes} describing the relationship between the observations and the processes.
#' @param sim_obs if set to \code{T}, the observations are ignored and `ideal' simulated observations are used instead. This option can be used to effect a `Turing test', where the real-data case and the ideal case can be compared side-by-side.
#' @return A data frame with the predictive statistics as described above.
#' @keywords predictive variance, validation
#' @references Sahu, S. K., & Mardia, K. V. (2005). A Bayesian kriged Kalman model for short-term forecasting of air pollution levels. Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(1), 223-244. Bastos, L. S. and O'Hagan, A. (2008). Diagnostics for Gaussian process emulators. Technometrics 51, 425-438.
#' @export
#' @examples
#' require(Matrix)
#' data(icesat)
#' data(surf_fe)
#'
#' ## First create observation object
#' icesat_obs <- Obs(df=icesat,
#'                  abs_lim = 5,
#'                  avr_method = "median",
#'                  box_size=100,
#'                  name="icesat")
#'
#' ## Now split into a training/validation set
#' split_data <- split_validation(icesat_obs,sample=500,common=0, t==2)
#' icesat_validation <- split_data$O_val
#' icesat_training <- split_data$O_pruned
#'
#' ## Now create GMRF defined over some FE basis
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
#'
#' mu <- matrix(0,nrow(Mesh),1)
#' Q <- sparseMatrix(i=1:nrow(surf_fe$p), j = 1:nrow(surf_fe$p), x = 1)
#'
#' my_GMRF <- GMRF(mu = mu, Q = Q,name="SURF",t_axis = 0:6)
#' SURF <-GMRF_basis(G = my_GMRF, Basis = Mesh)
#'
#' L1 <- link(SURF,icesat_training)
#' e <- link_list(list(L1))
#' v <- block_list(list(O = icesat_training, G = SURF))
#' G <- new("Graph",e=e,v=v)
#' G_reduced <- compress(G)
#' Results <- Infer(G_reduced)
#'
#' ## Now we validate the results with icesat
#' L1 <- link(SURF,icesat_validation)
#' e <- link_list(list(L1))
#' v <- block_list(list(G1=SURF,O1=icesat_validation))
#' G <- Graph(e=e,v=v)
#' G_reduced <- compress(G)
#' val_results <- validate(Results,G_reduced,sim_obs=F)
setGeneric("validate", function(Results,G,sim_obs=F,...) standardGeneric("validate"))

#' @title Sample from a GMRF
#' @description Takes a GMRF object and, possibly, associated permutation matrix and Cholesky factor of permuted precision matrix, to generate
#' samples.
#' @param G object of class \code{GMRF}
#' @param L Cholesky factor of the precision matrix, if \code{P} is \code{NULL} then this is treated as the factor of an unpermuted matrix
#' @param reps number of samples to generate
#' @param P permutation matrix
#' @return a matrix of size \code{n} by \code{reps}
#' @export
#' @examples
#' G <- GMRF_RW()
#' G <- setPrecision(G, getPrecision(G) + 0.01*Imat(nrow(G)))
#' X <- sample_GMRF(G,reps=10)
setGeneric("sample_GMRF", function(G,L=NULL,reps=1,P=NULL) standardGeneric("sample_GMRF"))

setGeneric(".exist", function(L,to,from) standardGeneric(".exist"))
setGeneric("getData", function(.Object) standardGeneric("getData"))
setGeneric("basisinterp", function(G,s,weights) standardGeneric("basisinterp"))
setGeneric(".find_inc_matrix",  function(basis,obs,mulfun, mask, n_grid, muldata, md5_wrapper) standardGeneric(".find_inc_matrix"))
