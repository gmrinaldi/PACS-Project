#' Create a FEM basis
#' 
#' @param mesh A \code{mesh.2D}, \code{mesh.2.5D} or \code{mesh.3D} object representing the domain triangulation. See \link{create.mesh.2D}, \link{create.mesh.2.5D}, \link{create.mesh.3D}.
#' @return A \code{FEMbasis} object. This contains the \code{mesh}, along with some additional quantities:
#' \itemize{
#' 	\item{\code{order}}{Either "1" or "2" for the 2D and 2.5D case, and "1" for the 3D case. 
#' 	Order of the Finite Element basis.} 
#' 	\item{\code{nbasis}}{Scalar. The number of basis.} 
#' 	\item{\code{transf_coord}}{It takes value only in the 2D case. It is a list of 4 vectors: diff1x, diff1y, diff2x and diff2y.
#' 	Each vector has length #triangles and encodes the information for the tranformation matrix that transforms the 
#' 	nodes of the reference triangle to the nodes of the i-th triangle. 
#' 	The tranformation matrix for the i-th triangle has the form [diff1x[i] diff2x[i]; diff1y[i] diff2y[i]].}
#' 	\item{\code{detJ}}{It takes value only in the 2D case. A vector of length #triangles. The ith element contains 
#' 	the determinant of the transformation from the reference triangle to the nodes of the i-th triangle. 
#' 	Its value is also the double of the area of each triangle of the basis.}
#' }
#' @description Sets up a Finite Element basis. It requires a \code{mesh.2D}, \code{mesh.2.5D} or \code{mesh.3D} object, 
#' as input. 
#' The basis' functions are globally continuos functions, that are polynomials once restricted to a triangle in the mesh. 
#' The current implementation includes linear finite elements (when \code{order = 1} in the input \code{mesh}) and 
#' quadratic finite elements (when \code{order = 2} in the input \code{mesh}).
#' @usage create.FEM.basis(mesh)
#' @seealso \code{\link{create.mesh.2D}}, \code{\link{create.mesh.2.5D}},\code{\link{create.mesh.3D}}
#' @examples 
#' ## Upload the quasicircle2D data
#' data(quasicircle2D)
#' 
#' ## Create the 2D mesh
#' mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
#' ## Plot it
#' plot(mesh)                   
#' ## Create the basis
#' FEMbasis = create.FEM.basis(mesh)
#' ## Upload the hub2.5D data
#' data(hub2.5D)
#' 
#' ## Create the 2.5D mesh
#' mesh = create.mesh.2.5D(nodes = hub2.5D.nodes, triangles = hub2.5D.triangles)
#' ## Plot it
#' plot(mesh)                   
#' ## Create the basis
#' FEMbasis = create.FEM.basis(mesh)
#' @export

create.FEM.basis = function(mesh)
{
  if(class(mesh)!='mesh.2D' & class(mesh)!='mesh.2.5D' & class(mesh)!='mesh.3D')
    stop("Unknown mesh class")
  
  if (class(mesh)=="mesh.2D"){

	  #  The number of basis functions corresponds to the number of vertices
	  #  for order = 1, and to vertices plus edge midpoints for order = 2

	  nbasis = dim(mesh$nodes)[[1]]
	  eleProp = R_elementProperties(mesh)

	  #eleProp = NULL
	  #if(CPP_CODE == FALSE)
	  #{
	  #  eleProp = R_elementProperties(mesh)
	  #}
  
  FEMbasis = list(mesh = mesh, order = as.integer(mesh$order), nbasis = nbasis, detJ=eleProp$detJ, transf_coord = eleProp$transf_coord)
  class(FEMbasis) = "FEMbasis"
  
  FEMbasis
  } else if (class(mesh) == "mesh.2.5D" || class(mesh) == "mesh.3D"){

  	  FEMbasis = list(mesh = mesh, order = as.integer(mesh$order),nbasis = mesh$nnodes)
  	  class(FEMbasis) = "FEMbasis"
  	  FEMbasis
  }
 }
#' Define a surface or spatial field by a Finite Element basis expansion
#' 
#' @param coeff A vector or a matrix containing the coefficients for the Finite Element basis expansion. The number of rows 
#' (or the vector's length) corresponds to the number of basis in \code{FEMbasis}. 
#' The number of columns corresponds to the number of functions. 
#' @param FEMbasis A \code{FEMbasis} object defining the Finite Element basis, created by \link{create.FEM.basis}.
#' @description This function defines a FEM object. 
#' @usage FEM(coeff,FEMbasis)
#' @return An \code{FEM} object. This contains a list with components \code{coeff} and \code{FEMbasis}.
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
#' ## Plot it
#' plot(FEMfunction)
#' @export

FEM<-function(coeff,FEMbasis)
{
  if (is.null(coeff)) 
    stop("coeff required;  is NULL.")
  if (is.null(FEMbasis)) 
    stop("FEMbasis required;  is NULL.")
  if(class(FEMbasis) != "FEMbasis")
    stop("FEMbasis not of class 'FEMbasis'")
  coeff = as.matrix(coeff)
  if(nrow(coeff) != FEMbasis$nbasis)
    stop("Number of row of 'coeff' different from number of basis")
  
  fclass = NULL
  fclass = list(coeff=coeff, FEMbasis=FEMbasis)
  class(fclass)<-"FEM"
  return(fclass)
}


R_elementProperties=function(mesh)
{
  nele = dim(mesh$triangles)[[1]]
  nodes = mesh$nodes
  triangles = mesh$triangles
  
  #detJ   = matrix(0,nele,1)      #  vector of determinant of transformations
  #metric = array(0,c(nele,2,2))  #  3-d array of metric matrices
  #transf = array(0,c(nele,2,2))
  
  transf_coord = NULL
  transf_coord$diff1x = nodes[triangles[,2],1] - nodes[triangles[,1],1]
  transf_coord$diff1y = nodes[triangles[,2],2] - nodes[triangles[,1],2]
  transf_coord$diff2x = nodes[triangles[,3],1] - nodes[triangles[,1],1]
  transf_coord$diff2y = nodes[triangles[,3],2] - nodes[triangles[,1],2]
  
  #  Jacobian or double of the area of triangle
  detJ = transf_coord$diff1x*transf_coord$diff2y - transf_coord$diff2x*transf_coord$diff1y
  
  #Too slow, computed only for stiff from diff1x,diff1y,..
  # for (i in 1:nele)
  # {
  #   #transf[i,,] = rbind(cbind(diff1x,diff2x),c(diff1y,diff2y))
  #   #  Compute controvariant transformation matrix OSS: This is (tranf)^(-T)
  #   Ael = matrix(c(diff2y, -diff1y, -diff2x,  diff1x),nrow=2,ncol=2,byrow=T)/detJ[i]
  #   
  #   #  Compute metric matrix
  #   metric[i,,] = t(Ael)%*%Ael
  # } 
  
  #FEStruct <- list(detJ=detJ, metric=metric, transf=transf)
  FEStruct <- list(detJ=detJ, transf_coord=transf_coord)
  return(FEStruct)
}
