testN<-10^8
time<-0
nodes<-cbind(runif(testN),runif(testN),runif(testN))
triangles2_5D<-cbind(seq_len(testN-2),seq_len(testN-2)+1,seq_len(testN-2)+2)
# triangles2_5D<-rbind(triangles2_5D,cbind(seq_len(testN-2),seq_len(testN-2)+2,seq_len(testN-2)+3))
triangles3D<-cbind(seq_len(testN-3),seq_len(testN-3)+1,seq_len(testN-3)+2,seq_len(testN-3)+3)

library(fdaPDE)
library(Rcpp)

setwd("~/Desktop/PACS/Progetto/PACS-Project/R_mesh")
sourceCpp("R_mesh.cpp")

data(hub2.5D)
edges_CPP<-R_mesh_helper_2_5D(hub2.5D.triangles,hub2.5D.nodes,nrow(hub2.5D.nodes))

data(sphere3Ddata)

nodes=sphere3Ddata$nodes
tetrahedrons=sphere3Ddata$tetrahedrons

edges3D_CPP<-R_mesh_helper_3D(tetrahedrons,nrow(nodes),TRUE)

create.mesh.2.5D_CPP<- function(nodes, triangles = NULL, order = 1, nodesattributes = NULL, segments = NULL, holes = NULL)
{
  ##########################
  ###   Input checking   ###
  ##########################
  
  nodesmarkers <- vector(mode = "integer", 0)
  segmentsmarkers <- vector(mode = "integer", 0)
  edgesmarkers <- vector(mode = "integer", 0)
  
  nodes <- as.matrix(nodes)
  if (ncol(nodes) != 3)
    stop("Matrix of nodes should have 3 columns")
  if (anyDuplicated(nodes))
    stop("Duplicated nodes")
  
  ## If attributes not specified, set them to a matrix with zero columns
  if (any(is.null(nodesattributes))) {
    nodesattributes <- matrix(0, nrow(nodes), 0)
  }else{
    nodesattributes <- as.matrix(nodesattributes)
    if (nrow(nodesattributes) != nrow(nodes))
      stop("Point attribute matrix \'nodesattributes\' does not have same number of rows the point matrix \'nodes\'")
  }
  
  ## Deal with segments
  if (any(is.null(segments))) {
    segments <- matrix(0, 0, 2)
  } else {
    segments <- as.matrix(segments)
    if (ncol(segments) != 2) {
      stop("Matrix of segments should have 2 columns")
    }
  }
  
  if (any(is.null(holes)))
    holes <- matrix(0, 0, 2)
  holes <- as.matrix(holes)
  
  ## If triangles are not already specified
  if(any(is.null(triangles))){
    stop("Per il momento in questo caso serve triangles")
    # triangles = matrix(0,nrow = 0, ncol = 3)
  } else {
    triangles = as.matrix(triangles)
  }
  
  outCPP<-R_mesh_helper_2_5D(triangles[,1:3],nrow(nodes), order==2 && ncol(triangles) == 3)
  
  out<-NULL
  
  if(order==1 && ncol(triangles) == 3){
    out <- list(nodes=nodes, nodesmarkers=outCPP$nodesmarkers, nodesattributes=nodesattributes,
                triangles=triangles, segments=segments, segmentsmarkers=segmentsmarkers,
                edges=matrix(outCPP$edges,ncol=2), edgesmarkers=outCPP$edgesmarkers, neighbors=matrix(outCPP$neighbors,ncol=3), holes=holes, order=order)
  }
  else if(order==2 && ncol(triangles) == 6){ # triangles matrix contains both the true triangles and the midpoints ones
    out <- list(nodes=nodes, nodesmarkers=outCPP$nodesmarkers, nodesattributes=nodesattributes,
                triangles=triangles, segments=segments, segmentsmarkers=segmentsmarkers,
                edges=matrix(outCPP$edges,ncol=2), edgesmarkers=outCPP$edgesmarkers, neighbors=matrix(outCPP$neighbors,ncol=3), holes=holes, order=order)
  }
  else if(order==2 && ncol(triangles) == 3){
    print("You set order=2 but passed a matrix of triangles with just 3 columns. The midpoints for each edge will be computed.")
    
    edges<-matrix(outCPP$edges,ncol=2)
    midpoints <- nodes[edges[,1],]/2+nodes[edges[,2],]/2

    nodes<-rbind(nodes,midpoints)
    nodesmarkers<-c(nodesmarkers,rep(0,nrow(midpoints)))
    
    out <- list(nodes=nodes, nodesmarkers=outCPP$nodesmarkers, nodesattributes=nodesattributes,
                triangles=cbind(triangles,matrix(outCPP$extra,ncol=3)), segments=segments, segmentsmarkers=segmentsmarkers,
                edges=edges, edgesmarkers=outCPP$edgesmarkers, neighbors=matrix(outCPP$neighbors,ncol=3), holes=holes, order=order)
  }
  else{
    stop("The number of columns of triangles matrix is not consistent with the order parameter")
  }
  
  class(out)<-"mesh.2.5D"
  
  return(out)
}

create.mesh.3D_CPP<- function(nodes, tetrahedrons, order = 1, nodesattributes = NULL, segments = NULL, holes = NULL)
{
  ##########################
  ###   Input checking   ###
  ##########################
  
  nodesmarkers <- vector(mode = "integer", 0)
  segmentsmarkers <- vector(mode = "integer", 0)
  edgesmarkers <- vector(mode = "integer", 0)
  facesmarkers <- vector(mode = "integer", 0)
  
  nodes = as.matrix(nodes)
  if (ncol(nodes) != 3)
    stop("Matrix of nodes should have 3 columns")
  if (anyDuplicated(nodes))
    stop("Duplicated nodes")
  
  ## If attributes not specified, set them to a matrix with zero columns
  if (any(is.null(nodesattributes))) {
    nodesattributes <- matrix(0, nrow(nodes), 0)
  }else{
    nodesattributes <- as.matrix(nodesattributes)
    if (nrow(nodesattributes) != nrow(nodes))
      stop("Point attribute matrix \'nodesattributes\' does not have same number of rows the point matrix \'nodes\'")
  }
  
  ## Deal with segments
  if (any(is.null(segments))) {
    segments <- matrix(0, 0, 2)
  } else {
    segments <- as.matrix(segments)
    if (ncol(segments) != 2) {
      stop("Matrix of segments should have 2 columns")
    }
  }
  
  if (any(is.null(holes)))
    holes <- matrix(0, 0, 2)
  holes <- as.matrix(holes)
  
  ## If tetrahedrons are not already specified
  if(any(is.null(tetrahedrons))){
    stop("Per il momento in questo caso serve tetrahedrons")
    # triangles = matrix(0,nrow = 0, ncol = 3)
  } else {
    tetrahedrons <- as.matrix(tetrahedrons)
  }

  outCPP<-R_mesh_helper_3D(tetrahedrons[,1:4],nrow(nodes),order==2 && ncol(tetrahedrons) == 4)  
  out<-NULL
  
  if(order==1 && ncol(tetrahedrons) == 4){
    out <- list(nodes=nodes, nodesmarkers=outCPP$nodesmarkers, nodesattributes=nodesattributes,
                tetrahedrons=tetrahedrons, segments=segments, segmentsmarkers=segmentsmarkers,
                faces=matrix(outCPP$faces,ncol=3), facesmarkers=outCPP$facesmarkers, neighbors=matrix(outCPP$neighbors,ncol=4),
                holes=holes, order=order)
  }
  else if(order==2 && ncol(tetrahedrons) == 10){ # tetrahedrons matrix contains both the true tetrahedrons and the midpoints ones
    out <- list(nodes=nodes, nodesmarkers=outCPP$nodesmarkers, nodesattributes=nodesattributes,
                tetrahedrons=tetrahedrons, segments=segments, segmentsmarkers=segmentsmarkers,
                faces=matrix(outCPP$faces,ncol=3), facesmarkers=facesmarkers, neighbors=matrix(outCPP$neighbors,ncol=4),
                holes=holes, order=order)
  }
  else if(order==2 && ncol(tetrahedrons) == 4){
    print("You set order=2 but passed a matrix of tetrahedrons with just 4 columns. The midpoints for each edge will be computed.")
    
    edges<-matrix(outCPP$edges,ncol=2)
    midpoints <- nodes[edges[,1],]/2+nodes[edges[,2],]/2

    nodes<-rbind(nodes,midpoints)
    nodesmarkers<-c(outCPP$nodesmarkers,rep(FALSE,nrow(midpoints)))

    out <- list(nodes=nodes, nodesmarkers=nodesmarkers, nodesattributes=nodesattributes,
                tetrahedrons=cbind(tetrahedrons,matrix(outCPP$extra,ncol=6)), segments=segments, segmentsmarkers=segmentsmarkers,
                faces=matrix(outCPP$faces,ncol=3), facesmarkers=outCPP$facesmarkers, neighbors=matrix(outCPP$neighbors,ncol=4),
                holes=holes, order=order)
  }
  else{
    stop("The number of columns of tetrahedrons matrix is not consistent with the order parameter")
  }
  
  class(out)<-"mesh.3D"
  
  return(out)
}


