#### CLASS DEFINITIONS ######

### Manifolds and measures
#'  @docType class
#'  @title measure
#'
#' @description Measure function
#'
#' @keywords Manifolds, spheres, planes
setClass("measure",representation(dist="function",dim="integer"),
         prototype(dist=fields::rdist,dim=2L))

#'  @docType class
#'  @title manifold
#'
#' @description Manifolds supported by the package.
#'
#' @keywords Manifolds, spheres, planes
setClass("manifold",representation(type="character", metric = "measure","VIRTUAL"))

#'  @docType class

#' @description Sphere manifold
#'
#' @keywords Manifolds, spheres, planes
setClass("sphere",representation(radius="numeric"),contains="manifold")

#'  @docType class
#'  @title plane
#'
#' @description Plane
#'
#' @keywords Manifolds, spheres, planes
setClass("plane",contains="manifold")

#'  @docType class
#'  @title real_line
#'
#' @description 1D real line
#'
#' @keywords Manifolds, spheres, planes, lines
setClass("real_line",contains="manifold")


#'  @docType class
#'  @title domain
#'
#' @description Domain
#'
#' @keywords Domain, bounding box
setClass("domain",representation(m = "manifold", bndary="matrix"))


####  Basis functions ####
#'  @docType class
#'  @title Parent class for basis functions
#'
#' @description A basis function contains three slots, the first is a list of parameters (pars), the second is the number (n) of basis functions and the third is a list
#' of functions (fn). The pars field usually contains parameters which appear in fn.
#'
#' @keywords Basis functions
#' @rdname Basisclass
setClass("Basis", representation(manifold="manifold",n = "numeric",fn="list",pars="list", df="data.frame"))

#'  @docType class
#'  @title Gaussian radial basis functions basis
#'
#' @description GRBFBasis inherits from the virtual class Basis. A GRBF basis functions is initialsied using the function \code{initGRBFbasis}.
#'
#' @keywords Basis functions, GRBF
setClass("GRBFBasis",contains="Basis")

#' @description bisquareBasis inherits from the virtual class Basis. A bisquare basis functions is initialsied using the function \code{initbisquarebasis}.
#'
#' @keywords Basis functions, GRBF
setClass("bisquareBasis",contains="Basis")


#'  @docType class
#'  @title Constant functions basis
#'
#' @description ConstBasis inherits from the virtual class Basis. A set of constant basis functions is initialsied using the function \code{initConstbasis}.
#'
#' @keywords Basis functions
#' #' @rdname ConstBasisclass
setClass("ConstBasis",contains="Basis")

#'  @docType class
#'  @title Finite element basis
#'
#' @description FEBasis inherits from the virtual class Basis. A fintie element basis is initialsied using the function \code{initFEbasis}.
#'
#' @keywords Basis functions
#' @rdname FEBasisclass
setClass("FEBasis",contains="Basis")

#'  @title SRE model
#' @description SRE model
#'
#' @keywords Spatial random effects, fixed rank kriging
setClass("SRE",representation(data="Spatial",
                              basis="Basis",
                              S = "Matrix",
                              V = "Matrix",
                              Z = "Matrix",
                              X = "Matrix"))

setClass("SRE.fit",contains="SRE",representation(
                              Khat = "Matrix",
                              alphahat = "numeric",
                              sigma2fshat = "numeric"))

####  Spatio-temporal blocks ####
#'  @docType class
#'  @title Parent (virtual) class for all blocks (latent fields and observations)
#'
#' @description A block is the most basic component of a multi-variate spatio-temporal model. It can either consist of a latent field or an observations. Each block
#' has a unique identifier (uid) and some cardinality n, the meaning of which depends on the specific class
#'
#' @keywords Block
#' @rdname blockclass
setClass("block", representation(uid="numeric",n="numeric","VIRTUAL"))


#'  @docType class
#'  @title Process block (virtual)
#'
#' @description A process block inherits from class \code{block}, however can only be used in certain ways in a graph. For example it can be linked to an observation
#' data block but not to another process block. It is the parent class of all latent processes.
#'
#' @keywords Block, process
#' @rdname processclass
setClass("process", contains="block",representation("VIRTUAL"))

#'  @docType class
#'  @title observation block
#'
#' @description An observation block inherits from class \code{block}, however can only be used in certain ways in a graph. For example it can be linked to an observation
#' process block but not to another observation block. To initialise use \code{initObs}.
#' @rdname Obsclass
setClass("Obs",contains=c("block"),representation(df="data.frame",args="list"))

#'  @docType class
#'  @title Link between latent process and observation
#'
#' @description The prime use of the \code{link} class is to generate an incidence matrix between observational data and the process.
#' This class is a virtual class, and is a parent for the specific \code{linkGO} class which specifies a link between a process block
#' and an observational block
#' @rdname linkclass
setClass("link", representation(from="block",to="block"))

#'  @docType class
#'  @title GMRF
#'
#' @description A basic GMRF object which can also take an intrinsic value for future definitions (intrinsic GMRFs not implemented yet). This
#' class inherits from \code{process} which in turn is a block.
#' @rdname GMRFclass
setClass("GMRF",contains="process",
         representation(mu="matrix", Q="dgCMatrix",intrinsic="numeric",rep="data.frame",
                        t_axis="numeric"),
         prototype(mu=matrix(0,2,1),
                   Q = sparseMatrix(i=c(1,2),j=c(1,2),x=1),
                   intrinsic=0,
                   n=2,rep=data.frame(),t_axis=0))

#'  @docType class
#'  @title VAR_Gauss
#'
#' @description A variable auto-regressive block which inherits from class \code{GMRF}. The primary difference is that this class is constructed using
#' temporal evolving functions. These functions are then used to construct a big precision matrix and mean vectors which are then passed on to the GMRF
#' constructor.
#' @rdname VAR_Gaussclass
setClass("VAR_Gauss",contains="GMRF",
         representation(mu_fun="function",A_fun = "function", B_fun = "function",
                        Qw_fun = "function",n="numeric",
                        Qb = "dgCMatrix"))

#'  @docType class
#'  @title GMRF_RW
#'
#' @description A random walk  which inherits from class \code{GMRF}. The primary difference is that this class is constructed using
#' a first-order auto-regressice structure. All random walks are intrinsic GMRFs.
#' @rdname GMRF_RWclass
setClass("GMRF_RW",contains="GMRF")

#'  @docType class
#'  @title GMRF_basis
#'
#' @description This class defines the amalgamation of a GMRF with a basis. It itself is a \code{process} block, i.e. it can be linked up with observations
#' using objects from the \code{link} class.
#' #' @rdname GMRF_basisclass
setClass("GMRF_basis",contains="process",  representation(G="GMRF",Basis="Basis"))

#'  @docType class
#'  @title Observations with non-trivial spatial footprints
#'
#' @description This class inherits from \code{Obs}. It extends \code{Obs} by providing the ability to encode spatial footprints with an additional input.
#' The use of this class requires installation of \code{gpclib}.
#' @rdname Obs_polyclass
setClass("Obs_poly",contains=c("Obs"),representation(pol="list",df2="data.frame"))

setClass("linkGO", contains="link",representation(Cmat="dgCMatrix"))
setClass("linkGG", contains="link",representation(cov_inter="matrix"))

#'  @docType class
#'  @title List of links
#'
#' @description The list of links is simply a collection of links, each with its own incidence matrix. The link list is combined with a block list in order to
#' compose a graph.
#' @rdname link_listclass
setClass("link_list", contains="list")

#'  @docType class
#'  @title List of blocks
#'
#' @description The list of blocks is simply a collection of blocks, which are either latent process or observations. The block list is combined with a link list in order to
#' compose a graph.
#' @rdname block_listclass
setClass("block_list", contains="list")

#' @docType class
#' @title Graph
#'
#' @description A collection of links (in a link list) and blocks (in a block list) composes a graph. The graph is the object over which we carry out inference.
#' @rdname graphclass
setClass("Graph", representation(e = "link_list", v = "block_list"),prototype(e = new("link_list"),v = new("block_list")))

#' @docType class
#' @title Graph_2nodes
#'
#' @description This class is a special case of \code{Graph}, where the link list and the block list only contain one node each. In general, a graph with many edges
#' can be compressed to a 2-node graph (one block latent process and one block observation set) using \code{compress}.
#' @rdname Graph_2nodesclass
setClass("Graph_2nodes",contains="Graph")
