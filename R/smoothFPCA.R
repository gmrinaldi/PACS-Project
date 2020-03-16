#' Smooth Functional Principal Component Analysis
#' 
#' @param datamatrix A matrix of dimensions #samples-by-#locations with the observed data values over the domain 
#' for each sample. The datamatrix needs to have zero mean. 
#' If the \code{locations} argument is left \code{NULL} the datamatrix has to be dimensions #samples-by-#nodes where #nodes
#' is the number of nodes of the mesh in the FEMbasis. In this case, each observation is associated to the corresponding 
#' node in the mesh. 
#' If the data are observed only on a subset of the mesh nodes, fill with \code{NA} the values of the 
#' \code{datamatrix} in correspondence of unobserved data. 
#' @param locations A #observations-by-2 matrix in the 2D case and #observations-by-3 matrix in the 2.5D and 3D case, where 
#' each row specifies the spatial coordinates \code{x} and \code{y} (and \code{z} in 2.5D and 3D) of the corresponding 
#' observation in the \code{datamatrix}. 
#' If the locations of the observations coincide with (or are a subset of) the nodes of the mesh in the \code{FEMbasis}, 
#' leave the parameter \code{locations = NULL} for a faster implementation. 
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param nPC An integer specifying the number of Principal Components to compute.
#' @param validation A string specifying the type of validation to perform. If \code{lambda} is a vector, it has to 
#' be specified as \code{"GCV"} or \code{"KFold"}. This parameter specify which method of cross-validation is used 
#' to select the best parameter \code{lambda} among those values of the smoothing parameter specified in \code{lambda} 
#' for each Principal Component.
#' @param NFolds This parameter is used only in case \code{validation = "KFold"}. It is an integer specifying 
#' the number of folds to use if the KFold cross-validation method for the 
#' selection of the best parameter \code{lambda} is chosen. Default value is 5. 
#' @param GCVmethod This parameter is considered only when \code{validation = "GCV"}. It can be either "Exact" or 
#' "Stochastic". If set to "Exact" the algoritm performs an exact (but possibly slow) computation 
#' of the GCV index. If set to "Stochastic" the GCV is approximated by a stochastic algorithm.
#' @param nrealizations The number of realizations to be used in the stochastic algorithm for the estimation of GCV.
#' @return A list with the following variables:
#' \itemize{
#' \item{\code{loadings.FEM}}{A \code{FEM} object that represents the L^2-normalized functional loadings for each 
#' Principal Component computed.}
#' \item{\code{scores}}{A #samples-by-#PrincipalComponents matrix that represents the unnormalized scores or PC vectors.}
#' \item{\code{lambda}}{A vector of length #PrincipalComponents with the values of the smoothing parameter \code{lambda} 
#' chosen for that Principal Component.}
#' \item{\code{variance_explained}}{A vector of length #PrincipalComponents where each value represent the variance explained by that component.}
#' \item{\code{cumsum_percentage}}{A vector of length #PrincipalComponents containing the cumulative percentage of the variance explained by the first components.}
#' }
#' @description This function implements a smooth functional principal component analysis over a planar mesh, 
#' a smooth manifold or a volume. 
#' @usage FPCA.FEM(locations = NULL, datamatrix, FEMbasis, lambda, nPC=1, 
#'          validation = NULL, NFolds = 5, GCVmethod = "Stochastic", nrealizations = 100)
#' @references Lila, E., Aston, J.A.D.,  Sangalli, L.M., 2016a. Smooth Principal Component Analysis over two-dimensional 
#' manifolds with an application to neuroimaging. Ann. Appl. Stat., 10(4), pp. 1854-1879. 
#' @export
#' @examples
#' library(fdaPDE)
#' 
#' ## Load the hub data
#' data(hub2.5D)
#' mesh = create.mesh.2.5D(nodes = hub2.5D.nodes, triangles = hub2.5D.triangles)
#' ## Create the Finite Element basis 
#' FEMbasis = create.FEM.basis(mesh)
#' ## Create a datamatrix
#' datamatrix = NULL
#' for(ii in 1:50){
#'   a1 = rnorm(1, mean = 1, sd = 1)
#'   a2 = rnorm(1, mean = 1, sd = 1)
#'   a3 = rnorm(1, mean = 1, sd = 1)
#'   
#'   func_evaluation = numeric(mesh$nnodes)
#'   for (i in 0:(mesh$nnodes-1)){
#'     func_evaluation[i+1] = a1* sin(2*pi*mesh$nodes[i+1,1]) +  a2* sin(2*pi*mesh$nodes[i+1,2]) +  a3*sin(2*pi*mesh$nodes[i+1,3]) +1
#'   }
#'   data = func_evaluation + rnorm(mesh$nnodes, mean = 0, sd = 0.5)
#'   datamatrix = rbind(datamatrix, data)
#' }
#' ## Compute the mean of the datamatrix and subtract it to the data
#' data_bar = colMeans(datamatrix)
#' data_demean = matrix(rep(data_bar,50), nrow=50, byrow=TRUE)
#' 
#' datamatrix_demeaned = datamatrix - data_demean
#' ## Set the smoothing parameter lambda
#' lambda = 0.00375
#' ## Estimate the first 2 Principal Components
#' FPCA_solution = FPCA.FEM(datamatrix = datamatrix_demeaned, 
#'                       FEMbasis = FEMbasis, lambda = lambda, nPC = 2)
#' 
#' ## Plot the functional loadings of the estimated Principal Components                           
#' plot(FPCA_solution$loadings.FEM)

FPCA.FEM<-function(locations = NULL, datamatrix, FEMbasis,lambda, nPC = 1, validation = NULL, NFolds = 5, 
                   GCVmethod = "Stochastic", nrealizations = 100)
{
  incidence_matrix=NULL # if areal fpca will be included in release, this should be put in the input
  
 if(class(FEMbasis$mesh) == "mesh.2D"){
 	ndim = 2
 	mydim = 2
 }else if(class(FEMbasis$mesh) == "mesh.2.5D"){
 	ndim = 3
 	mydim = 2
 }else if(class(FEMbasis$mesh) == "mesh.3D"){
 	ndim = 3
 	mydim = 3
 }else{
 	stop('Unknown mesh class')
 }
 
  if(GCVmethod=="Stochastic")
    GCVmethod=2
  else if(GCVmethod=="Exact")
    GCVmethod=1
  else{
    stop("GCVmethod must be either Stochastic or Exact")
  }
##################### Checking parameters, sizes and conversion ################################

  checkSmoothingParametersFPCA(locations, datamatrix, FEMbasis, incidence_matrix, lambda, nPC, validation, NFolds, GCVmethod ,nrealizations) 
  ## Coverting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  datamatrix = as.matrix(datamatrix)
  if(!is.null(incidence_matrix))
	incidence_matrix = as.matrix(incidence_matrix)
  lambda = as.matrix(lambda)
  
  checkSmoothingParametersSizeFPCA(locations, datamatrix, FEMbasis, incidence_matrix, lambda, ndim, mydim, validation, NFolds)
  
	  ################## End checking parameters, sizes and conversion #############################
  
  bigsol = NULL
  if(class(FEMbasis$mesh) == 'mesh.2D'){
	  print('C++ Code Execution')
	  bigsol = CPP_smooth.FEM.FPCA(locations, datamatrix, FEMbasis, incidence_matrix,
	  								lambda, ndim, mydim, nPC, validation, NFolds, 
									GCVmethod, nrealizations)
	  numnodes = nrow(FEMbasis$mesh$nodes)
  } else if(class(FEMbasis$mesh) == 'mesh.2.5D'){
	  print('C++ Code Execution')
	  bigsol = CPP_smooth.manifold.FEM.FPCA(locations, datamatrix, FEMbasis,
	  										incidence_matrix, lambda, ndim, mydim,
											nPC, validation, NFolds, GCVmethod, nrealizations)
	  numnodes = FEMbasis$mesh$nnodes
  } else if(class(FEMbasis$mesh) == 'mesh.3D'){
	  print('C++ Code Execution')
	  bigsol = CPP_smooth.volume.FEM.FPCA(locations, datamatrix, FEMbasis,
	  										incidence_matrix, lambda, ndim, mydim,
											nPC, validation, NFolds, GCVmethod, nrealizations)
	  numnodes = FEMbasis$mesh$nnodes
  }
  
  loadings=bigsol[[1]]
  loadings.FEM=FEM(loadings,FEMbasis)
  
  scores=bigsol[[2]]
  
  lambda=bigsol[[3]]
  
  variance_explained=bigsol[[4]]
  
  cumsum_percentage=bigsol[[5]]
  
  var=bigsol[[6]]
  
  reslist=list(loadings.FEM=loadings.FEM, scores=scores, lambda=lambda, variance_explained=variance_explained, cumsum_percentage=cumsum_percentage)
  return(reslist)
}
