#' Evaluate a FEM object at a set of point locations
#' 
#' @param FEM A \code{FEM} object to be evaluated.
#' @param locations A 2-columns (in 2D) or 3-columns (in 2.5D and 3D) matrix with the spatial locations where the 
#' FEM object should be evaluated. 
#' @param incidence_matrix In case of areal evaluations, the #regions-by-#elements incidence matrix defining the regions 
#' where the FEM object should be evaluated.
#' @return 
#' A vector or a matrix of numeric evaluations of the \code{FEM} object. 
#' If the \code{FEM} object contains multiple finite element functions the output is a matrix, and 
#' each row corresponds to the location (or areal region) where the evaluation has been taken, while each column 
#' corresponds to the function evaluated. 
#' @description It evaluates a FEM object at the specified set of locations or areal regions. The locations are used for 
#' pointwise evaluations and incidence matrix for areal evaluations. 
#' The locations and the incidence matrix cannot be both NULL or both provided.  
#' @usage eval.FEM(FEM, locations = NULL, incidence_matrix = NULL)
#' @references 
#' \itemize{
#'    \item{Sangalli, L. M., Ramsay, J. O., & Ramsay, T. O. (2013). Spatial spline regression models. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(4), 681-703.}
#'    \item{Azzimonti, L., Sangalli, L. M., Secchi, P., Domanin, M., & Nobile, F. (2015). Blood flow velocity field estimation 
#' via spatial regression with PDE penalization. Journal of the American Statistical Association, 110(511), 1057-1071.}
#' }
#' @examples 
#' library(fdaPDE)
#' ## Upload the horseshoe2D data
#' data(horseshoe2D)
#' 
#' ## Create the 2D mesh
#' mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
#' ## Create the FEM basis
#' FEMbasis = create.FEM.basis(mesh)
#' ## Compute the coeff vector evaluating the desired function at the mesh nodes
#' ## In this case we consider the fs.test() function introduced by Wood et al. 2008
#' coeff = fs.test(mesh$nodes[,1], mesh$nodes[,2], exclude = FALSE)
#' ## Create the FEM object
#' FEMfunction = FEM(coeff, FEMbasis)
#' 
#' ## Evaluate the finite element function in the location (1,0.5)
#' eval.FEM(FEMfunction, locations = matrix(c(1, 0.5), ncol = 2))
#' 
#' ## Evaluate the mean of the finite element function over the fifth triangle of the mesh
#' incidence_matrix = matrix(0, ncol = nrow(mesh$triangles))
#' incidence_matrix[1,5] = 1
#' eval.FEM(FEMfunction, incidence_matrix = incidence_matrix)
#' @export

eval.FEM <- function(FEM, locations = NULL, incidence_matrix = NULL)
{
  if (is.null(FEM))
    stop("FEM required;  is NULL.")
  if(class(FEM) != "FEM")
    stop("'FEM' is not of class 'FEM'")
  if (is.null(locations) && is.null(incidence_matrix)) 
    stop("'locations' NOR 'incidence_matrix' required;  both are NULL.")
  if (!is.null(locations) && !is.null(incidence_matrix))
    stop("'locations' NOR 'incidence_matrix' required; both are given.")
  
  if(!is.null(locations))
   if(dim(locations)[1]==dim(FEM$FEMbasis$mesh$nodes)[1] & dim(locations)[2]==dim(FEM$FEMbasis$mesh$nodes)[2])
    warning("The locations matrix has the same dimensions as the mesh nodes. If you want to get the FEM object evaluation
            at the mesh nodes, use FEM$coeff instead")
  
  if (is.null(locations))
    locations <- matrix(nrow = 0, ncol = 2)
  else
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  
  res <- NULL
  
  if(class(FEM$FEMbasis$mesh)=='mesh.2D'){
    ndim = 2
    mydim = 2
    res = CPP_eval.FEM(FEM, locations, incidence_matrix, TRUE, ndim, mydim)
    
  }else if(class(FEM$FEMbasis$mesh)=='mesh.2.5D'){
    ndim = 3
    mydim = 2
    res = CPP_eval.manifold.FEM(FEM, locations, incidence_matrix, TRUE, ndim, mydim)
  }else if(class(FEM$FEMbasis$mesh)=='mesh.3D'){
    ndim = 3
    mydim = 3
    res = CPP_eval.volume.FEM(FEM, locations, incidence_matrix, TRUE, ndim, mydim)
  }
  
  return(as.matrix(res))
}


# R_eval.FEM <- function(FEM, locations)
# { 
#   if (is.vector(locations))
#   {
#     locations = t(as.matrix(locations))
#   }
#   
#   N = nrow(locations)
#   
#   # Augment Xvec and Yvec by ones for computing barycentric coordinates
#   Pgpts = cbind(matrix(1,N,1),locations[,1],locations[,2])
#   
#   # Get nodes and index
#   
#   FEMbasis = FEM$FEMbasis
#   mesh = FEMbasis$mesh
#   
#   nodes = mesh$nodes
#   triangles = mesh$triangles
#   coeff = FEM$coeff
#   nsurf = dim(coeff)[2]
#   
#   FEMbasis = FEM$FEMbasis
#   order = FEMbasis$order
#   #nodeindex = params$nodeindex
#   detJ = FEMbasis$detJ
#   
#   # 1st, 2nd, 3rd vertices of triangles
#   
#   v1 = nodes[triangles[,1],]
#   v2 = nodes[triangles[,2],]
#   v3 = nodes[triangles[,3],]
#   
#   if(order !=2 && order != 1)
#   {
#     stop('ORDER is neither 1 or 2.')
#   }
#   
#   # Denominator of change of coordinates chsange matrix
#   
#   modJ = FEMbasis$detJ
#   ones3 = matrix(1,3,1)
#   modJMat = modJ %*% t(ones3)
#   
#   M1 = cbind(v2[,1]*v3[,2] - v3[,1]*v2[,2], v2[,2] - v3[,2], v3[,1] - v2[,1])/(modJMat)
#   M2 = cbind(v3[,1]*v1[,2] - v1[,1]*v3[,2], v3[,2] - v1[,2], v1[,1] - v3[,1])/(modJMat)
#   M3 = cbind(v1[,1]*v2[,2] - v2[,1]*v1[,2], v1[,2] - v2[,2], v2[,1] - v1[,1])/(modJMat)
#   
#   ind = matrix(0,N,1)
#   for(i in 1:N)
#   {
#     ind[i] = R_insideIndex(mesh, as.numeric(locations[i,]))
#   }
#   
#   evalmat = matrix(NA, nrow=N, ncol=nsurf)
#   
#   for (isurf in 1:nsurf)
#   {
#     for(i in 1:N)
#     {
#       indi = ind[i]
#       
#       if(!is.nan(indi))
#       {
#         baryc1 = (M1[indi,]*Pgpts[i,]) %*% ones3
#         baryc2 = (M2[indi,]*Pgpts[i,]) %*% ones3
#         baryc3 = (M3[indi,]*Pgpts[i,]) %*% ones3
#         
#         if(order == 2)
#         {
#           c1 = coeff[triangles[indi,1],isurf]
#           c2 = coeff[triangles[indi,2],isurf]
#           c3 = coeff[triangles[indi,3],isurf]
#           c4 = coeff[triangles[indi,6],isurf]
#           c5 = coeff[triangles[indi,4],isurf]
#           c6 = coeff[triangles[indi,5],isurf]
#           
#           fval =  c1*(2* baryc1^2 - baryc1) +
#             c2*(2* baryc2^2 - baryc2) +
#             c3*(2* baryc3^2 - baryc3) +
#             c4*(4* baryc1 * baryc2) +
#             c5*(4* baryc2 * baryc3) +
#             c6*(4* baryc3 * baryc1)
#           evalmat[i,isurf] = fval
#         }
#         else
#         {
#           c1 = coeff[triangles[indi,1],isurf]
#           c2 = coeff[triangles[indi,2],isurf]
#           c3 = coeff[triangles[indi,3],isurf]
#           fval = c1*baryc1 + c2*baryc2 + c3*baryc3
#           evalmat[i,isurf] = fval
#         }
#       }
#     }
#   }
#   return(evalmat)
# }
# 
# R_insideIndex = function (mesh, location)
# {
#   #  insideIndex returns the index of the triangle containing the point
#   # (X,Y) if such a triangle exists, and NaN otherwise.
#   #  TRICOEF may have already been calculated for efficiency,
#   #  but if the function is called with four arguments, it is calculated.
#   
#   
#   eps=2.2204e-016
#   small = 10000*eps
#   
#   nodes = mesh$nodes
#   triangles = mesh$triangles
#   X = location[1]
#   Y = location[2]
#   
#   ntri   = dim(triangles)[[1]]
#   indtri   = matrix(1:ntri,ncol=1)
#   
#   #  compute coefficients for computing barycentric coordinates if needed
#   
#   tricoef = R_tricoefCal(mesh)
#   
#   #  compute barycentric coordinates
#   r3 = X - nodes[triangles[,3],1]
#   s3 = Y - nodes[triangles[,3],2]
#   lam1 = ( tricoef[,4]*r3 - tricoef[,2]*s3)
#   lam2 = (-tricoef[,3]*r3 + tricoef[,1]*s3)
#   lam3 = 1 - lam1 - lam2
#   
#   #  test these coordinates for a triple that are all between 0 and 1
#   int  = (-small <= lam1 & lam1 <= 1+small) & 
#     (-small <= lam2 & lam2 <= 1+small) & 
#     (-small <= lam3 & lam3 <= 1+small)
#   
#   #  return the index of this triple, or NaN if it doesn't exist
#   indi = indtri[int]
#   if (length(indi)<1)
#   {
#     ind = NA
#   }else{
#     ind = min(indi)
#   }
#   
#   ind
# }
