#' Spatial regression with differential regularization
#' 
#' @param observations A vector of length #observations with the observed data values over the domain. 
#' If the \code{locations} argument is left NULL the vector of the observations have to be of length #nodes of the 
#' mesh in the FEMbasis. In this case, each observation is associated to the corresponding node in the mesh. 
#' If the observations are observed only on a subset of the mesh nodes, fill with \code{NA} the values of the vector 
#' \code{observations} in correspondence of unobserved data. 
#' @param locations A #observations-by-2 matrix in the 2D case and #observations-by-3 matrix in the 2.5D and 3D case, where 
#' each row specifies the spatial coordinates \code{x} and \code{y} (and \code{z} in 2.5D and 3D) of the corresponding 
#' observation in the vector \code{observations}.
#' If the locations of the observations coincide with (or are a subset of) the nodes of the mesh in the \code{FEMbasis}, 
#' leave the parameter \code{locations = NULL} for a faster implementation. 
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with 
#' the corresponding observed data value in \code{observations} and each column is a different covariate.
#' @param PDE_parameters A list specifying the parameters of the PDE in the regularizing term. Default is NULL, i.e. 
#' regularization is by means of the Laplacian (stationary, isotropic case). 
#' If the coefficients of the PDE are constant over the domain \code{PDE_parameters} must contain: 
#' \itemize{
#'    \item{\code{K}, a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic 
#' smoothing with a preferential direction that corresponds to the first eigenvector of the diffusion matrix K;}
#'    \item{\code{b}, a vector of length 2 of advection coefficients. This induces a 
#' smoothing only in the direction specified by the vector \code{b};} 
#'    \item{\code{c}, a scalar reaction coefficient. \code{c} induces a shrinkage of the surface to zero.}
#' }
#' If the coefficients of the PDE are space-varying \code{PDE_parameters} must contain: 
#' \itemize{
#' \item{\code{K}, a function that for each spatial location in the spatial domain (indicated by the vector of the 2 
#' spatial coordinates) returns a 2-by-2 matrix of diffusion coefficients. The function must support recycling for 
#' efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' an array with dimensions 2-by-2-by-#points.}
#' \item{\code{b}, a function that for each spatial location in the spatial domain returns 
#' a vector of length 2 of transport coefficients. The function must support recycling for efficiency reasons, thus 
#' if the input parameter is a #point-by-2 matrix, the output should be
#' a matrix with dimensions 2-by-#points;} 
#' \item{\code{c}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points;} 
#' \item{\code{u}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' \code{u} induces a reaction effect. The function must support recycling for efficiency reasons, thus if the input 
#' parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points.}
#' }
#' For 2.5D and 3D, only the Laplacian is available (\code{PDE_parameters=NULL}). 
#' @param incidence_matrix A #regions-by-#triangles/tetrahedrons matrix where the element (i,j) equals 1 if the j-th 
#' triangle/tetrahedron is in the i-th region and 0 otherwise.
#' This is needed only for areal data. In case of pointwise data, this parameter is set to \code{NULL}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param GCVmethod This parameter is considered only when \code{GCV=TRUE}. It can be either "Exact" or "Stochastic". 
#' If set to "Exact" the algoritm performs an exact (but possibly slow) computation 
#' of the GCV index. If set to "Stochastic" the GCV is approximated by a stochastic algorithm.
#' @param nrealizations This parameter is considered only when \code{GCV=TRUE} and \code{GCVmethod = "Stochastic"}. 
#' It is a positive integer that represents the number of uniform random variables used in stochastic GCV computation.      
#' @return A list with the following variables:
#' \itemize{
#'    \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#'    \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the Laplacian of the estimated spatial field.}
#'    \item{\code{beta}}{If covariates is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when 
#' the smoothing parameter is equal to \code{lambda[j]}.}
#'    \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'    \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#'    \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' }
#' @description This function implements a spatial regression model with differential regularization. 
#'  The regularizing term involves a Partial Differential Equation (PDE). In the simplest case the PDE involves only the 
#'  Laplacian of the spatial field, that induces an isotropic smoothing. When prior information about the anisotropy or 
#'  non-stationarity is available the PDE involves a general second order linear differential operator with possibly 
#'  space-varying coefficients. 
#'  The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions 
#'  can be imposed at the domain boundaries.
#' @usage smooth.FEM(locations = NULL, observations, FEMbasis, lambda, 
#'      covariates = NULL, PDE_parameters=NULL, incidence_matrix = NULL, 
#'      BC = NULL, GCV = FALSE, GCVmethod = "Stochastic", nrealizations = 100)
#' @export        

#' @references 
#' \itemize{
#'    \item{Sangalli, L. M., Ramsay, J. O., Ramsay, T. O. (2013). Spatial spline regression models. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(4), 681-703.}
#'    \item{Azzimonti, L., Sangalli, L. M., Secchi, P., Domanin, M., Nobile, F. (2015). Blood flow velocity field estimation 
#' via spatial regression with PDE penalization. Journal of the American Statistical Association, 110(511), 1057-1071.}
#' }
#' @examples
#' library(fdaPDE)
#' 
#' #### No prior information about anysotropy/non-stationarity (laplacian smoothing) ####
#' data(horseshoe2D)
#' FEMbasis = create.FEM.basis(mesh)
#' lambda = 10^-1
#' # no covariate
#' data = fs.test(mesh$nodes[,1], mesh$nodes[,2], exclude = FALSE) + rnorm(nrow(mesh$nodes), sd = 0.5)
#' 
#' solution = smooth.FEM.basis(observations = data, FEMbasis = FEMbasis, lambda = lambda)
#' plot(solution$fit.FEM)
#' 
#' # with covariates
#' covariate = covs.test(mesh$nodes[,1], mesh$nodes[,2])
#' data = fs.test(mesh$nodes[,1], mesh$nodes[,2], exclude = FALSE) + 2*covariate + rnorm(nrow(mesh$nodes), sd = 0.5)
#' 
#' solution = smooth.FEM.basis(observations = data, covariates = covariate, FEMbasis = FEMbasis, lambda = lambda)
#' # beta estimate:
#' solution$beta
#' # non-parametric estimate:
#' plot(solution$fit.FEM)
#' 
#' # Choose lambda with GCV:
#' lambda = 10^(-2:2)
#' solution = smooth.FEM.basis(observations = data, 
#'                             covariates = covariate, 
#'                             FEMbasis = FEMbasis, 
#'                             lambda = lambda, 
#'                             GCV = TRUE)
#' bestLambda = lambda[which.min(solution$GCV)]
#' 
#' 
#' #### Smoothing with prior information about anysotropy/non-stationarity and boundary conditions ####
#' # See Azzimonti et al. for reference to the current exemple
#' data(quasicircle2D)
#' mesh = create.MESH.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
#' FEMbasis = create.FEM.basis(mesh)
#' lambda = 10^-2
#' 
#' # Set the PDE parameters
#' R = 2.8
#' K1 = 0.1
#' K2 = 0.2
#' beta = 0.5
#' K_func<-function(points)
#' {
#'   output = array(0, c(2, 2, nrow(points)))
#'   for (i in 1:nrow(points))
#'     output[,,i] = 10*rbind(c(points[i,2]^2 + K1*points[i,1]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2),
#'                              (K1-1)*points[i,1]*points[i,2]),
#'                            c((K1-1)*points[i,1]*points[i,2],
#'                              points[i,1]^2 + K1*points[i,2]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2)))
#'   output
#' }
#' 
#' b_func<-function(points)
#' {
#'   output = array(0, c(2, nrow(points)))
#'   for (i in 1:nrow(points))
#'     output[,i] = 10*beta*c(points[i,1],points[i,2])
#'   output
#' }
#' 
#' c_func<-function(points)
#' {
#'   rep(c(0), nrow(points))
#' }
#' 
#' u_func<-function(points)
#' {
#'   rep(c(0), nrow(points))
#' }
#' PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
#' 
#' # Set the boundary conditions
#' BC = NULL
#' BC$BC_indices = which(mesh$nodesmarkers == 1) # b.c. on the complete boundary
#' BC$BC_values = rep(0,length(BC$BC_indices)) # homogeneus b.c.
#' 
#' # Since the data locations are a subset of the mesh nodes for a faster solution use:
#' dataNA = rep(NA, FEMbasis$nbasis)
#' dataNA[mesh$nodesmarkers == 0] = data
#' 
#' solution = smooth.FEM.basis(observations = dataNA, 
#'                             FEMbasis = FEMbasis, 
#'                             lambda = lambda, 
#'                             PDE_parameters = PDE_parameters, 
#'                             BC = BC)
#' plot(solution$fit.FEM)
#' image(solution$fit.FEM)
#' 
#' #### Smoothing with areal data ####
#' # See Azzimonti et al. for reference to the current exemple
#' data(quasicircle2Dareal)
#' 
#' FEMbasis = create.FEM.basis(mesh)
#' lambda = 10^-4
#' 
#' # Set the PDE parameters
#' R = 2.8
#' K1 = 0.1
#' K2 = 0.2
#' beta = 0.5
#' K_func<-function(points)
#' {
#'   output = array(0, c(2, 2, nrow(points)))
#'   for (i in 1:nrow(points))
#'     output[,,i] = 10*rbind(c(points[i,2]^2 + K1*points[i,1]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2),
#'                              (K1-1)*points[i,1]*points[i,2]),
#'                            c((K1-1)*points[i,1]*points[i,2],
#'                              points[i,1]^2 + K1*points[i,2]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2)))
#'   output
#' }
#' 
#' b_func<-function(points)
#' {
#'   output = array(0, c(2, nrow(points)))
#'   for (i in 1:nrow(points))
#'     output[,i] = 10*beta*c(points[i,1],points[i,2])
#'   output
#' }
#' 
#' c_func<-function(points)
#' {
#'   rep(c(0), nrow(points))
#' }
#' 
#' u_func<-function(points)
#' {
#'   rep(c(0), nrow(points))
#' }
#' PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
#' 
#' # Set the boundary conditions
#' BC = NULL
#' BC$BC_indices = which(mesh$nodesmarkers == 1) # b.c. on the complete boundary
#' BC$BC_values = rep(0,length(BC$BC_indices)) # homogeneus b.c.
#' 
#' solution = smooth.FEM.basis(observations = data, 
#'                             incidence_matrix = incidence_matrix,
#'                             FEMbasis = FEMbasis, 
#'                             lambda = lambda, 
#'                             PDE_parameters = PDE_parameters, 
#'                             BC = BC)
#' plot(solution$fit.FEM)
#' image(solution$fit.FEM)

smooth.FEM<-function(locations = NULL, observations, FEMbasis, lambda, 
                     covariates = NULL, PDE_parameters=NULL, incidence_matrix = NULL, 
                     BC = NULL, GCV = FALSE, GCVmethod = "Stochastic", nrealizations = 100)
{
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
  ##################### Checking parameters, sizes and conversion ################################
  
  if(GCVmethod=="Stochastic")
    GCVMETHOD=2
  else if(GCVmethod=="Exact")
    GCVMETHOD=1
  else{
    stop("GCVmethod must be either Stochastic or Exact")
  }
  
  space_varying=checkSmoothingParameters(locations, observations, FEMbasis, lambda, covariates, incidence_matrix, BC, GCV, PDE_parameters, GCVMETHOD , nrealizations)
  
  ## Coverting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  observations = as.matrix(observations)
  lambda = as.matrix(lambda)
  if(!is.null(covariates))
    covariates = as.matrix(covariates)
  if(!is.null(incidence_matrix))
    incidence_matrix = as.matrix(incidence_matrix)
  if(!is.null(BC))
  {
    BC$BC_indices = as.matrix(BC$BC_indices)
    BC$BC_values = as.matrix(BC$BC_values)
  }
  
  # if I have PDE non-sv case I need (constant) matrices as parameters
  
  if(!is.null(PDE_parameters) & space_varying==FALSE) 
  {
    PDE_parameters$K = as.matrix(PDE_parameters$K)
    PDE_parameters$b = as.matrix(PDE_parameters$b)
    PDE_parameters$c = as.matrix(PDE_parameters$c)
  }
  
  
  checkSmoothingParametersSize(locations, observations, FEMbasis, lambda, covariates, incidence_matrix, BC, GCV, space_varying, PDE_parameters,ndim, mydim)
  
  ################## End checking parameters, sizes and conversion #############################
  
  if(class(FEMbasis$mesh) == 'mesh.2D' & is.null(PDE_parameters)){	
    
    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
                                  covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                  BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations)
  
    numnodes = nrow(FEMbasis$mesh$nodes)
    
  } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==FALSE){
    
    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
                                      PDE_parameters = PDE_parameters,
                                      covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                      BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations)
    
    numnodes = nrow(FEMbasis$mesh$nodes)

  } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==TRUE){	
    
    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.sv.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
                                         PDE_parameters = PDE_parameters,
                                         covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                         BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations)
  
    numnodes = nrow(FEMbasis$mesh$nodes)
  
  }else if(class(FEMbasis$mesh) == 'mesh.2.5D'){
    
    bigsol = NULL  
    print('C++ Code Execution')
    if(!is.null(locations))
      stop("The option locations!=NULL for manifold domains is currently not implemented")
    bigsol = CPP_smooth.manifold.FEM.basis(locations, observations, FEMbasis, lambda, covariates, incidence_matrix, ndim, mydim, BC, GCV,GCVMETHOD, nrealizations)
    
    numnodes = FEMbasis$mesh$nnodes
    
  }else if(class(FEMbasis$mesh) == 'mesh.3D'){
    
    bigsol = NULL  
    print('C++ Code Execution')
    bigsol = CPP_smooth.volume.FEM.basis(locations, observations, FEMbasis, lambda, covariates, incidence_matrix, ndim, mydim, BC, GCV,GCVMETHOD, nrealizations)
    
    numnodes = FEMbasis$mesh$nnodes
  }
  
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)
  
  reslist = NULL
  beta = getBetaCoefficients(locations, observations, fit.FEM, covariates, incidence_matrix, ndim, mydim)
  if(GCV == TRUE)
  {
    seq=getGCV(locations = locations, observations = observations, fit.FEM = fit.FEM, covariates = covariates, incidence_matrix = incidence_matrix, edf = bigsol[[2]], ndim, mydim)
    reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, beta = beta, edf = bigsol[[2]], stderr = seq$stderr, GCV = seq$GCV)
  }else{
    reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, beta = beta)
  }
  
  return(reslist)
}

getBetaCoefficients<-function(locations, observations, fit.FEM, covariates, incidence_matrix = NULL, ndim, mydim)
{
  loc_nodes = NULL
  fnhat = NULL
  betahat = NULL
  
  if(!is.null(covariates))
  {
    if(is.null(locations))
    {
      loc_nodes = (1:length(observations))[!is.na(observations)]
      fnhat = as.matrix(fit.FEM$coeff[loc_nodes,])
    }else{
      loc_nodes = 1:length(observations)
      fnhat = eval.FEM(FEM = fit.FEM, locations = locations, incidence_matrix = incidence_matrix)
    }
    ## #row number of covariates, #col number of functions
    betahat = matrix(0, nrow = ncol(covariates), ncol = ncol(fnhat))
    for(i in 1:ncol(fnhat))
      betahat[,i] = as.vector(lm.fit(covariates,as.vector(observations-fnhat[,i]))$coefficients)
  }
  
  return(betahat)
}


getGCV<-function(locations, observations, fit.FEM, covariates = NULL, incidence_matrix = NULL, edf, ndim, mydim)
{
  loc_nodes = NULL
  fnhat = NULL
  
  edf = as.vector(edf)
  
  if(is.null(locations) && is.null(incidence_matrix))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    fnhat = as.matrix(fit.FEM$coeff[loc_nodes,])
  }else{
    loc_nodes = 1:length(observations)
    fnhat = eval.FEM(FEM = fit.FEM, locations = locations, incidence_matrix = incidence_matrix)
  }
  
  zhat = NULL
  zhat = matrix(nrow = length(loc_nodes), ncol = length(edf))
  if(!is.null(covariates))
  {
    desmatprod = ( solve( t(covariates) %*% covariates ) ) %*% t(covariates)
    for ( i in 1:length(edf))
    {
      betahat  = desmatprod %*% (observations-fnhat[,i])
      zhat[,i] = covariates %*% betahat + fnhat[,i]
    }
  }else{
    zhat = fnhat
  }
  
  np = length(loc_nodes)
  
  stderr2 = numeric(length(edf))
  GCV = numeric(length(edf))
  
  zhat <- as.matrix(zhat)
  
  if(any(np - edf <= 0))
  {
    warning("Some values of 'edf' are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.")  
  }
  
  for (i in 1:length(edf))
  {
    stderr2[i] = t(observations[loc_nodes] - zhat[,i]) %*% (observations[loc_nodes] - zhat[,i]) / ( np - edf[i] )
    GCV[i] = ( np / ( np - edf[i] )) * stderr2[i]
  }
  
  # NA if stderr2 is negative
  stderr = vector('numeric', length(stderr2));
  stderr[stderr2>=0] = sqrt(stderr2[stderr2>=0]);
  stderr[stderr2<0] = NaN;
  
  return(list(stderr = stderr, GCV = GCV))
}

