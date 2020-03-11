#' Create a \code{MESH.2.5D} object from the connectivty matrix and nodes locations
#'
#' @param nodes A #nodes-by-3 matrix specifying the locations of each node
#' @param triangles A #triangles-by-3*order matrix specifying the indices of the nodes in each triangle
#' @param order Either "1" or "2". Order of the Finite Element basis default is order = 1
#' @return An object of the class \code{MESH.2.5D} with the following output:
#' \item{\code{nnodes}}{The #nodes contained in the mesh}
#' \item{\code{ntriangles}}{The #triangles contained in the mesh}
#' \item{\code{nodes}}{A #nodes-by-3 matrix containing the x,y and z coordinate for each point of the mesh}
#' \item{\code{triangles}}{A #triangles-by-3*order matrix specifying the indices of the nodes in each triangle of the mesh}
#' \item{\code{order}}{Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints). 
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements. Default is \code{order} = 1.}
#' @examples
#' #Load the matrix nodes and triangles
#'
#' library(fdaPDE)
#' data(sphereData)
#'
#' nodes=sphere$nodes
#' triangles=sphere$triangles
#'
#' #Create the triangulated mesh from the connectivity matrix and nodes locations
#' mesh=create.MESH.2.5D(nodes,triangles)
#' 

create.MESH.2.5D<- function(nodes, triangles, order = 1)
{
  nnodes = dim(nodes)[1]
  
  ntriangles = dim(triangles)[1]
  
  if(dim(triangles)[2]!= 3*order){
    if (order==1)
      stop("The matrix 'triangles' has the wrong number of columns. See second.order.mesh(...)")
    stop("The matrix 'triangles' has wrong number of columns. Should be 3*order \n")
  }
  out = list(nnodes=nnodes, ntriangles=ntriangles, nodes=as.matrix(nodes), triangles = as.matrix(triangles), order=as.integer(order))
  
  class(out)<-"MESH.2.5D"
  
  return(out)
}


#' Double the order of a fist order Finite Element mesh by adding middle points to each side of the triangles in the triangulation
#' @param mesh an object of class 'MESH.2.5D' that is the starting mesh of order 1
#' @param bc A vector specifying the indices of the nodes on which boundary conditions are applied.
#' @return If no boundaries conditions are passed:
#' \item{\code{mesh}}{An object of class 'MESH.2.5D' with the mesh of order 2.} Otherwhise a \code{list} with parameters:
#' \item{\code{mesh}}{An object of class 'MESH.2.5D' with the mesh of order 2.}
#' \item{\code{bc_index}}{An update of the vector specifying the indices of the nodes on which boundary conditions are applied.}
#' @usage second.order.MESH.2.5D(mesh,bc=NULL)
#' @seealso \code{\link{create.MESH.2.5D}}
#' @examples
#' #Loading mesh hub, a MESH.2.5D object of order=1 contained in the package
#' data(hub)
#' #Apply the function to generate a surface mesh of order 2
#' hub_order2 = second.order.MESH.2.5D(hub)

second.order.MESH.2.5D<-function(mesh,bc=NULL){
  if(class(mesh) != 'MESH.2.5D'){
    stop('This method is implemented only for a mesh of class MESH.2.5D')
  }else if(mesh$order != 1){
    stop('The object mesh must have order = 1')
  }else{
    toll=1e-5
    # T = matrix(mesh$triangles,nrow=mesh$ntriangles,ncol=3, byrow = TRUE)
    # V = matrix(mesh$nodes, nrow = mesh$nnodes, ncol= 3, byrow = TRUE)
    T=mesh$triangles
    V=mesh$nodes
    T <- cbind(T, matrix(0,nrow=nrow(T),ncol=3))
    nnodes=nrow(V)
    index=nrow(V)
    points = V[T[1,],]
    midpoints<-rbind((points[2,]+points[3,])/2,(points[1,]+points[3,])/2, (points[1,]+points[2,])/2);
    if(!is.null(bc)){
      isBC<-c( any(bc==T[1,2]) & any(bc==T[1,3]),
               any(bc==T[1,1]) & any(bc==T[1,3]),
               any(bc==T[1,2]) & any(bc==T[1,1]))
    }
    
    for (side in 1:3){
      point<-midpoints[side,]
      index<-index+1;
      V<-rbind(V,point)
      T[1,3+side]<-index;
      
      if(!is.null(bc)&&isBC[side]==1){
        bc<-c(bc,index)
      }
      
    }
    
    for (i in 2:nrow(T)){
      points = V[T[i,],]
      midpoints<-rbind((points[2,]+points[3,])/2,(points[1,]+points[3,])/2, (points[1,]+points[2,])/2);
      if(!is.null(bc)){
        isBC<-c( any(bc==T[i,2]) & any(bc==T[i,3]),
                 any(bc==T[i,1]) & any(bc==T[i,3]),
                 any(bc==T[i,2]) & any(bc==T[i,1]))
      }
      
      for (side in 1:3){
        point<-midpoints[side,]
        isthere<-apply(V[(nnodes+1):nrow(V),], 1, function(x) identical(as.vector(x), point))
        loc = which(isthere)
        if(length(loc)>0){
          loc = loc+nnodes
          T[i,3+side]<-loc[1]
        }else{
          index<-index+1;
          V<-rbind(V,point)
          T[i,3+side]<-index;
          
          if(!is.null(bc)&&isBC[side]==1){
            bc<-c(bc,index)
          }
        }
      }
    }
  }
  if(is.null(bc)){
    out = list(nnodes=nrow(V), ntriangles=nrow(T), nodes=V, triangles = T, order=2)
    class(out)<-"MESH.2.5D"
    return(out)
  }else{
    out = list(nnodes=nrow(V), ntriangles=nrow(T), nodes=V, triangles = T, order=2)
    class(out)<-"MESH.2.5D"
    retlist = list(mesh = out, bc_index=bc)
    return(retlist)
  }
}
