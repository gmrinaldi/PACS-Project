triangulate_native <- function(P, PB, PA, S, SB,H, TR, flags) {
  ## It is necessary to check for NAs and NaNs, as the triangulate C
  ## code crashes if fed with them

  P  <- as.matrix(P)
  PB <- as.integer(PB)
  PA <- as.matrix(PA)
  S  <- as.matrix(S)
  SB <- as.integer(SB)
  H  <- as.matrix(H)
  TR  <- as.matrix(TR)

  storage.mode(P)  <- "double"
  storage.mode(PA) <- "double"
  #storage.mode(PB) <- "integer"
  storage.mode(S)  <- "integer"
  #storage.mode(SB) <- "integer"
  storage.mode(H)  <- "double"
  storage.mode(TR) <- "integer"
  storage.mode(flags) <- 'character'
  ## Call the main routine
  out <- .Call("R_triangulate_native",
               t(P),
               PB,
               PA,
               t(S),
               SB,
               H,
               t(TR),
               flags,
               PACKAGE="fdaPDE")
  names(out) <- c("P", "PB", "PA", "T", "S", "SB", "E", "EB","TN", "VP", "VE", "VN", "VA")
  class(out) <- "triangulation"
  return(out)
}

#' Create a 2D triangular mesh
#'
#' @param nodes A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.
#' @param nodesattributes A matrix with #nodes rows containing nodes' attributes.
#' These are passed unchanged to the output. If a node is added during the triangulation process or mesh refinement, its attributes are computed
#' by linear interpolation using the attributes of neighboring nodes. This functionality is for instance used to compute the value
#' of a Dirichlet boundary condition at boundary nodes added during the triangulation process.
#' @param segments A #segments-by-2 matrix. Each row contains the row's indices in \code{nodes} of the vertices where the segment starts from and ends to.
#' Segments are edges that are not splitted during the triangulation process. These are for instance used to define the boundaries
#' of the domain. If this is input is NULL, it generates a triangulation over the
#' convex hull of the points specified in \code{nodes}.
#' @param holes A #holes-by-2 matrix containing the x and y coordinates of a point internal to each hole of the mesh. These points are used to carve holes
#' in the triangulation, when the domain has holes.
#' @param triangles A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix.
#' This option is used when a triangulation is already available. It specifies the triangles giving the row's indices in \code{nodes} of the triangles' vertices and (when \code{nodes} = 2) also if the triangles' edges midpoints. The triangles' vertices and midpoints are ordered as described
#' at \cr https://www.cs.cmu.edu/~quake/triangle.highorder.html.
#' In this case the function \code{create.mesh.2D} is used to produce a complete mesh.2D object.
#' @param order Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints).
#' These are
#' respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements. Default is \code{order} = 1.
#' @param verbosity This can be '0', '1' or '2'. It indicates the level of verbosity in the triangulation process. When \code{verbosity} = 0 no message is returned
#' during the triangulation. When \code{verbosity} = 2 the triangulation process is described step by step by displayed messages.
#' Default is \code{verbosity} = 0.
#' @description This function is a wrapper of the Triangle library (http://www.cs.cmu.edu/~quake/triangle.html). It can be used
#' to create a triangulation of the domain of interest starting from a list of points, to be used as triangles' vertices, and a list of segments, that define the domain boundary. The resulting
#' mesh is a Constrained Delaunay triangulation. This is constructed in a way to preserve segments provided in the input \code{segments} without splitting them. This imput can be used to define the boundaries
#' of the domain. If this imput is NULL, it generates a triangulation over the
#' convex hull of the points.
#' It is also possible to create a mesh.2D from the nodes locations and the connectivity matrix.
#' @usage create.mesh.2D(nodes, nodesattributes = NULL, segments = NULL, holes = NULL,
#'                      triangles = NULL, order = 1, verbosity = 0)
#' @seealso \code{\link{refine.mesh.2D}}, \code{\link{create.FEM.basis}}
#' @return An object of the class mesh.2D with the following output:
#' \itemize{
#' \item{\code{nodes}}{A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.}
#' \item{\code{nodesmarkers}}{A vector of length #nodes, with entries either '1' or '0'. An entry '1' indicates that the corresponding node is a boundary node; an entry '0' indicates that the corresponding node is not a boundary node.}
#' \item{\code{nodesattributes}}{A matrix with #nodes rows containing nodes' attributes.
#' These are passed unchanged from the input.}
#' \item{\code{triangles}}{A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix.
#' This option is used when a triangulation is already available. It specifies the triangles giving the indices in \code{nodes} of the triangles' vertices and (when \code{nodes} = 2) also if the triangles' edges midpoints. The triangles' vertices and midpoints are ordered as described
#' at  \cr https://www.cs.cmu.edu/~quake/triangle.highorder.html.}
#' \item{\code{segmentsmarker}}{A vector of length #segments with entries either '1' or '0'. An entry '1' indicates that the corresponding element in \code{segments} is a boundary segment;
#' an entry '0' indicates that the corresponding segment is not a boundary segment.}
#' \item{\code{edges}}{A #edges-by-2 matrix containing all the edges of the triangles in the output triangulation. Each row contains the row's indices in \code{nodes}, indicating the nodes where the edge starts from and ends to.}
#' \item{\code{edgesmarkers}}{A vector of lenght #edges with entries either '1' or '0'. An entry '1' indicates that the corresponding element in \code{edge} is a boundary edge;
#' an entry '0' indicates that the corresponding edge is not a boundary edge.}
#' \item{\code{neighbors}}{A #triangles-by-3 matrix. Each row contains the indices of the three neighbouring triangles. An entry '-1' indicates that
#' one edge of the triangle is a boundary edge.}
#' \item{\code{holes}}{A #holes-by-2 matrix containing the x and y coordinates of a point internal to each hole of the mesh. These points are used to carve holes
#' in the triangulation, when the domain has holes.}
#' \item{\code{order}}{Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints).
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements.}
#' }
#' @export
#' @examples
#' library(fdaPDE)
#'
#' ## Upload the quasicirle2D data
#' data(quasicircle2D)
#'
#' ## Create mesh from boundary
#' ## if the domain is convex it is sufficient to call:
#' mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations))
#' plot(mesh)
#'
#' ## if the domain is not convex, pass in addition the segments the compose the boundary:
#' mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
#'
#' ## Create mesh from data locations (without knowing the boundary)
#' mesh = create.mesh.2D(nodes = locations)
#' plot(mesh)
#' ## In this case the domain is the convex hull of the data locations.
#' ## Do this only if you do not have any information about the shape of the domain of interest.

create.mesh.2D <- function(nodes, nodesattributes = NULL, segments = NULL, holes = NULL, triangles = NULL, order = 1, verbosity = 0)
{
  ##########################
  ###   Input checking   ###
  ##########################

  # Triangle finds out which are on the border (see https://www.cs.cmu.edu/~quake/triangle.help.html)
  nodesmarkers = vector(mode = "integer", 0)
  segmentsmarkers = vector(mode = "integer", 0)

  nodes = as.matrix(nodes)
  if (ncol(nodes) != 2)
    stop("Matrix of nodes should have 2 columns")
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

  ## If hole not specified, set it to empty matrix
  if (any(is.null(holes)))
    holes <- matrix(0, 0, 2)
  holes = as.matrix(holes)

  ## If triangles are not already specified
  if(any(is.null(triangles)))
    triangles = matrix(0,nrow = 0, ncol = 3)
  triangles = as.matrix(triangles)

  ## Set meshing parameters ##
  flags="ven"
  if(nrow(segments) == 0){
    flags = paste(flags,"c",sep = '')
  }

  if(nrow(segments)>0){
    flags = paste(flags,"p",sep = '')
  }

  #If order=2 add flag for second order nodes
  if(order == 2){
    flags = paste(flags,"o2",sep = '')
  }
  if(order < 1 || order >2){
    print('Order must be 1 or 2')
  }

  if(nrow(triangles) > 0){
    flags = paste(flags,"r",sep = '')
  }

  if (verbosity == 0) {
    flags = paste(flags,"Q",sep = '')
  }
  if (verbosity == 1) {
    flags = paste(flags,"V",sep = '')
  }
  if (verbosity == 2) {
    flags = paste(flags,"VV",sep = '')
  }

  out<-NULL
  #If triangles is null it makes the trianglulation
  #If triangle is not null it makes a refinement with no parameter, to compose the mesh object
  out <- triangulate_native(
    nodes,
    nodesmarkers,
    nodesattributes,
    segments,
    segmentsmarkers,
    t(holes),
    triangles,
    flags
  )

  names(out)[1]<-"nodes"
  names(out)[2]<-"nodesmarkers"
  names(out)[3]<-"nodesattributes"
  names(out)[4]<-"triangles"
  names(out)[5]<-"segments"
  names(out)[6]<-"segmentsmarkers"
  names(out)[7]<-"edges"
  names(out)[8]<-"edgesmarkers"
  names(out)[9]<-"neighbors"

  out[13]<-NULL
  out[12]<-NULL
  out[11]<-NULL
  out[10]<-NULL

  out[[10]] = holes
  names(out)[10]<-"holes"
  out[[11]] = order
  names(out)[11]<-"order"


  class(out)<-"mesh.2D"

  return(out)
}

#' Refine a 2D triangular mesh
#'
#' @param mesh A mesh.2D object representing the triangular mesh, created by \link{create.mesh.2D}.
#' @param minimum_angle A scalar specifying a minimun value for the triangles angles.
#' @param maximum_area A scalar specifying a maximum value for the triangles areas.
#' @param delaunay A boolean parameter indicating whether or not the output mesh should satisfy the Delaunay condition.
#' @param verbosity This can be '0', '1' or '2'. It indicates the level of verbosity in the triangulation process.
#' @description This function refines a Constrained Delaunay triangulation into a Conforming Delaunay triangulation. This is a wrapper of the Triangle library (http://www.cs.cmu.edu/~quake/triangle.html). It can be used to
#' refine a mesh previously created with \link{create.mesh.2D}. The algorithm can add Steiner points (points through which the \code{segments} are splitted)
#' in order to meet the imposed refinement conditions.
#' @usage refine.mesh.2D(mesh, minimum_angle, maximum_area, delaunay, verbosity)
#' @seealso \code{\link{create.mesh.2D}}, \code{\link{create.FEM.basis}}
#' @return A mesh.2D object representing the refined triangular mesh,  with the following output:
#' \itemize{
#' \item{\code{nodes}}{A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.}
#' \item{\code{nodesmarkers}}{A vector of length #nodes, with entries either '1' or '0'. An entry '1' indicates that the corresponding node is a boundary node; an entry '0' indicates that the corresponding node is not a boundary node.}
#' \item{\code{nodesattributes}}{nodesattributes A matrix with #nodes rows containing nodes' attributes.
#' These are passed unchanged to the output. If a node is added during the triangulation process or mesh refinement, its attributes are computed
#' by linear interpolation using the attributes of neighboring nodes. This functionality is for instance used to compute the value
#' of a Dirichlet boundary condition at boundary nodes added during the triangulation process.}
#' \item{\code{triangles}}{A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix.}
#' \item{\code{edges}}{A #edges-by-2 matrix. Each row contains the row's indices of the nodes where the edge starts from and ends to.}
#' \item{\code{edgesmarkers}}{A vector of lenght #edges with entries either '1' or '0'. An entry '1' indicates that the corresponding element in \code{edge} is a boundary edge;
#' an entry '0' indicates that the corresponding edge is not a boundary edge.}
#' \item{\code{neighbors}}{A #triangles-by-3 matrix. Each row contains the indices of the three neighbouring triangles. An entry '-1' indicates that
#' one edge of the triangle is a boundary edge.}
#' \item{\code{holes}}{A #holes-by-2 matrix containing the x and y coordinates of a point internal to each hole of the mesh. These points are used to carve holes
#' in the triangulation, when the domain has holes.}
#' \item{\code{order}}{Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints).
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements.}
#' }
#' @examples
#' library(fdaPDE)
#'
#' ## Upload the quasicircle2D data
#' data(quasicircle2D)
#'
#' ## Create mesh from boundary:
#' mesh = create.mesh.2D(nodes = boundary_nodes, segments = boundary_segments)
#' plot(mesh)
#' ## Refine the mesh with the maximum area criterion:
#' finemesh = refine.mesh.2D(mesh = mesh, maximum_area = 0.1)
#' plot(finemesh)
#' ## Refine the mesh with the minimum angle criterion:
#' finemesh2 = refine.mesh.2D(mesh = mesh, minimum_angle = 30)
#' plot(finemesh2)
#' @export

refine.mesh.2D<-function(mesh, minimum_angle = NULL, maximum_area = NULL, delaunay = FALSE, verbosity = 0)
{
  if(class(mesh) !="mesh.2D")
    stop("Sorry, this function is implemented just for mesh.2D class ")

  flags="rpven"

  if(!is.null(minimum_angle)){
    flags <- paste(flags, "q", sprintf("%.12f", minimum_angle), sep='')
  }

  if(!is.null(maximum_area)){
    flags <- paste(flags, "a", sprintf("%.12f", maximum_area), sep='')
  }

  if(delaunay){
    flags <- paste(flags, "D", sep='')
  }

  if(mesh$order==2){
    flags <- paste(flags, "o2", sep='')
  }

  if (verbosity == 0) {
    flags = paste(flags,"Q",sep = '')
  }
  if (verbosity == 1) {
    flags = paste(flags,"V",sep = '')
  }
  if (verbosity == 2) {
    flags = paste(flags,"VV",sep = '')
  }

  # Triangle finds out which are on the border (see https://www.cs.cmu.edu/~quake/triangle.help.html)
  mesh$nodesmarkers = vector(mode = "integer", 0)
  mesh$segmentsmarkers = vector(mode = "integer", 0)

  out<-NULL
  #If triangles is null it makes the trianglulation
  #If triangle is not null it makes a refinement with no parameter, to compose the mesh object
  out <- triangulate_native(
    mesh$nodes,
    mesh$nodesmarkers,
    mesh$nodesattributes,
    mesh$segments,
    mesh$segmentmarkers,
    t(mesh$holes),
    mesh$triangles,
    flags
  )

  names(out)[1]<-"nodes"
  names(out)[2]<-"nodesmarkers"
  names(out)[3]<-"nodesattributes"
  names(out)[4]<-"triangles"
  names(out)[5]<-"segments"
  names(out)[6]<-"segmentsmarkers"
  names(out)[7]<-"edges"
  names(out)[8]<-"edgesmarkers"
  names(out)[9]<-"neighbors"

  out[13]<-NULL
  out[12]<-NULL
  out[11]<-NULL
  out[10]<-NULL

  out[[10]] = mesh$holes
  names(out)[10]<-"holes"
  out[[11]] = mesh$order
  names(out)[11]<-"order"

  class(out)<-"mesh.2D"

  return(out)
}

#' Create a \code{mesh.2.5D} object from the nodes locations and the connectivty matrix
#'
#' @param nodes A #nodes-by-3 matrix specifying the locations of each node.
#' @param triangles A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix,
#' specifying the indices of the nodes in each triangle.
#' @param order Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes
#' (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints).
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements.
#' Default is \code{order} = 1.
#' @return An object of the class \code{mesh.2.5D} with the following output:
#' \itemize{
#' \item{\code{nnodes}}{The #nodes in the mesh.}
#' \item{\code{ntriangles}}{The #triangles in the mesh.}
#' \item{\code{nodes}}{A #nodes-by-3 matrix containing the x,y and z coordinate for each point of the mesh.}
#' \item{\code{triangles}}{A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix,
#' specifying the indices of the nodes in each triangle.}
#' \item{\code{order}}{Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes
#' (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints).
#' It is passed unchanged from the input.}
#' }
#' @export
#' @examples
#' library(fdaPDE)
#'
#' ## Upload the hub2.5D the data
#' data(hub2.5D)
#'
#' ## Create mesh from nodes and connectivity matrix:
#' mesh = create.mesh.2.5D(nodes = hub2.5D.nodes, triangles = hub2.5D.triangles)
#' plot(mesh)

create.mesh.2.5D<- function(nodes, triangles = NULL, order = 1, nodesattributes = NULL, segments = NULL, holes = NULL)
{
  ##########################
  ###   Input checking   ###
  ##########################

  nodesmarkers = vector(mode = "integer", 0)
  segmentsmarkers = vector(mode = "integer", 0)
  edgesmarkers = vector(mode = "integer", 0)

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
  holes = as.matrix(holes)

  ## If triangles are not already specified
  if(any(is.null(triangles))){
    stop("Per il momento in questo caso serve triangles")
    # triangles = matrix(0,nrow = 0, ncol = 3)
  } else {
    triangles = as.matrix(triangles)
  }

  #Compute neighbors matrix (see Ueng, Sikorski 1996)
  neighbors <- matrix(-1, nrow = nrow(triangles), ncol = 3)

  #Note: keep only first 3 cols of triangles to avoid issues in case order=2
  ordered_triangles<-cbind(t(apply(triangles[,1:3], 1, sort)),1:nrow(triangles))

  edge_list<-apply(ordered_triangles, 1, function(x){
    E1<-c(x[c(1,2)],x[4],1)
    E2<-c(x[c(1,3)],x[4],2)
    E3<-c(x[c(2,3)],x[4],3)
    list(E1,E2,E3)
  })

  edge_list<-unlist(edge_list, recursive = FALSE)

  for (level in 2:1){
      bin_list <- vector(mode = "list", length = nrow(nodes))
      for (i in 1:length(edge_list))
          bin_list[[edge_list[[i]][level]]]<-c(bin_list[[edge_list[[i]][level]]],edge_list[i])
      edge_list <- unlist(bin_list, recursive = FALSE)
  }
  
  #Repeated[i] is true if the i-th edge is a duplicate of the i-1-th edge
  #Caution: <<- is used here to fill out the neighbors matrix (outside of the function scope)
  # (no problem here because there is no variable called neighbors in the scope)
  repeated<-mapply(function(X,Y){
    out<-identical(X[1:2],Y[1:2])
    if (out){
      adjacent<-rbind(X[3:4],Y[3:4])
      neighbors[adjacent]<<-adjacent[2:1,1]
    }
    return(out)
  },edge_list[-length(edge_list)],edge_list[-1])
  #First edge is not considered (since there is no "0-th edge" set repeated[1]=FALSE)
  repeated<-c(FALSE,repeated)

  #Remove duplicates and set edges
  edges <- matrix(unlist(edge_list),ncol=4, byrow = TRUE)[,1:2]
  edges <- edges[!repeated,]

  #Set edgesmarkers and nodemarkers
  #Boundary edges are not shared among triangles hence they are never repeated
  #Meaning that both repeated[i] and repeated[i+1] have to be FALSE
  #FALSE at the end accounts for the last edge
  edgesmarkers <- !c(repeated[-1],FALSE)[!repeated]
  nodesmarkers <- 1:nrow(nodes) %in% edges[edgesmarkers]
  edgesmarkers <- as.numeric(edgesmarkers)
  nodesmarkers <- as.numeric(nodesmarkers)

  out<-NULL

  if(order==1 && ncol(triangles) == 3){
    out = list(nodes=nodes, nodesmarkers=nodesmarkers, nodesattributes=nodesattributes,
        triangles=triangles, segments=segments, segmentsmarkers=segmentsmarkers,
          edges=edges, edgesmarkers=edgesmarkers, neighbors=neighbors, holes=holes, order=order)
  }
  else if(order==2 && ncol(triangles) == 6){ # triangles matrix contains both the true triangles and the midpoints ones
    out = list(nodes=nodes, nodesmarkers=nodesmarkers, nodesattributes=nodesattributes,
        triangles=triangles, segments=segments, segmentsmarkers=segmentsmarkers,
          edges=edges, edgesmarkers=edgesmarkers, neighbors=neighbors, holes=holes, order=order)
  }
  else if(order==2 && ncol(triangles) == 3){
    print("You set order=2 but passed a matrix of triangles with just 3 columns. The midpoints for each edge will be computed.")

    midpoints <- nodes[edges[,1],]/2+nodes[edges[,2],]/2
    triangle_labels <- matrix(unlist(edge_list),ncol=4, byrow = TRUE)[,3:4]
    triangle_labels[,2]<-triangle_labels[,2]+3
    indices<-nrow(nodes)+cumsum(!repeated)

    triangles<-cbind(triangles,matrix(0,nrow=nrow(triangles),ncol=3))
    triangles[triangle_labels]<-indices
    
    nodes<-rbind(nodes,midpoints)
    nodesmarkers<-c(nodesmarkers,rep(0,nrow(midpoints)))
    
    out = list(nodes=nodes, nodesmarkers=nodesmarkers, nodesattributes=nodesattributes,
        triangles=triangles, segments=segments, segmentsmarkers=segmentsmarkers,
               edges=edges, edgesmarkers=edgesmarkers, neighbors=neighbors, holes=holes, order=order)
  }
  else{
    stop("The number of columns of triangles matrix is not consistent with the order parameter")
  }

  class(out)<-"mesh.2.5D"

  return(out)
}

#' Create a \code{mesh.3D} object from the connectivity matrix and nodes locations
#'
#' @param nodes A #nodes-by-3 matrix specifying the locations of each node.
#' @param tetrahedrons A #tetrahedrons-by-4 matrix specifying the indices of the nodes in each tetrahedrons.
#' @param order Order of the Finite Element basis. Only order = 1 is currently implemented.
#' @return An object of the class \code{mesh.3D} with the following output:
#' \itemize{
#' \item{\code{nnodes}}{The #nodes in the mesh.}
#' \item{\code{ntetrahedrons}}{The #tetrahedrons in the mesh.}
#' \item{\code{nodes}}{A #nodes-by-3 matrix containing the x,y and z coordinate for each point of the mesh.}
#' \item{\code{tetrahedrons}}{A #tetrahedrons-by-4 matrix specifying the indices of the nodes in each tetrahedron of the mesh.}
#' \item{\code{order}}{It specifies the order of the associated Finite Element basis. When order = 1, each
#' mesh tetrahedron is represented by 4 nodes (the tetrahedron vertices).}
#' }
#' @export
#' @examples
#' library(fdaPDE)
#'
#' ##Load the matrix nodes and tetrahedrons
#' data(sphere3Ddata)
#'
#' nodes=sphere3Ddata$nodes
#' tetrahedrons=sphere3Ddata$tetrahedrons
#'
#' ##Create the triangulated mesh from the connectivity matrix and nodes locations
#' mesh=create.mesh.3D(nodes,tetrahedrons)

create.mesh.3D<- function(nodes, tetrahedrons, order = 1, nodesattributes = NULL, segments = NULL, holes = NULL)
{
  ##########################
  ###   Input checking   ###
  ##########################

  nodesmarkers = vector(mode = "integer", 0)
  segmentsmarkers = vector(mode = "integer", 0)
  edgesmarkers = vector(mode = "integer", 0)
  facesmarkers = vector(mode = "integer", 0)

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
  holes = as.matrix(holes)

  ## If tetrahedrons are not already specified
  if(any(is.null(tetrahedrons))){
    stop("Per il momento in questo caso serve tetrahedrons")
    # triangles = matrix(0,nrow = 0, ncol = 3)
  } else {
    tetrahedrons = as.matrix(tetrahedrons)
  }

  #Compute neighbors (see Ueng, Sikorski 1996)
  neighbors <- matrix(-1, nrow = nrow(tetrahedrons), ncol = 4)

  #Note: select first 4 cols of tetrahedrons to avoid issues in case order=2
  ordered_tetrahedrons<-cbind(t(apply(tetrahedrons[,1:4], 1, sort)),1:nrow(tetrahedrons))

  faces_list<-apply(ordered_tetrahedrons, 1, function(x){
    F1<-c(x[c(1,2,3)],x[5],1)
    F2<-c(x[c(1,2,4)],x[5],2)
    F3<-c(x[c(1,3,4)],x[5],3)
    F4<-c(x[c(2,3,4)],x[5],4)

    list(F1,F2,F3,F4)
  })

  faces_list<-unlist(faces_list, recursive = FALSE)

  for (level in 3:1){
      bin_list <- vector(mode = "list", length = nrow(nodes))
      for (i in 1:length(faces_list))
          bin_list[[faces_list[[i]][level]]]<-c(bin_list[[faces_list[[i]][level]]],faces_list[i])
    faces_list <- unlist(bin_list, recursive = FALSE)
  }

  #Repeated[i] is true if the i-th face is a duplicate of the i-1-th face
  #Caution: <<- is used here to fill out the neighbors matrix (outside of the function scope)
  # (no problem here because there is no variable called neighbors in the scope)
  repeated_faces<-mapply(function(X,Y){
    out<-identical(X[1:3],Y[1:3])
    if (out){
      adjacent<-rbind(X[4:5],Y[4:5])
      neighbors[adjacent]<<-adjacent[2:1,1]
    }
    return(out)
  },faces_list[-length(faces_list)],faces_list[-1])
  #First face is not considered (since there is no "0-th face" set repeated[1]=FALSE)
  repeated_faces<-c(FALSE,repeated_faces)

  #Remove duplicates and set faces
  faces <- matrix(unlist(faces_list),ncol=5, byrow = TRUE)[,1:3]
  faces <- faces[!repeated_faces,]

  #Set facesmarkers and nodemarkers
  facesmarkers <- !c(repeated_faces[-1],FALSE)[!repeated_faces]
  nodesmarkers <- 1:nrow(nodes) %in% faces[facesmarkers]

  #Set edges
  edge_list <- apply(ordered_tetrahedrons, 1, function(x){
    E1<-c(x[c(1,2)],x[5],1)
    E2<-c(x[c(1,3)],x[5],2)
    E3<-c(x[c(1,4)],x[5],3)
    E4<-c(x[c(2,3)],x[5],4)
    E5<-c(x[c(2,4)],x[5],5)
    E6<-c(x[c(3,4)],x[5],6)
    list(E1,E2,E3,E4,E5,E6)
  })
  
  edge_list<-unlist(edge_list, recursive = FALSE)

  for (level in 2:1){
      bin_list <- vector(mode = "list", length = nrow(nodes))
      for (i in 1:length(edge_list))
          bin_list[[edge_list[[i]][level]]]<-c(bin_list[[edge_list[[i]][level]]],edge_list[i])
      edge_list <- unlist(bin_list, recursive = FALSE)
  }

  repeated_edges<-mapply(function(X,Y){
      identical(X[1:2],Y[1:2])
    },edge_list[-length(edge_list)],edge_list[-1])
  repeated_edges<-c(FALSE,repeated_edges)

  #Remove duplicates and set edges
  edges <- matrix(unlist(edge_list),ncol=4, byrow = TRUE)[,1:2]
  edges <- edges[!repeated_edges,]

  #Look for boundary edges (check which edges connect boundary nodes) 
  edgesmarkers <- matrix(edges %in% (1:nrow(nodes))[nodesmarkers], ncol=2)
  edgesmarkers <- edgesmarkers[,1] & edgesmarkers[,2]
    
  edgesmarkers <- as.numeric(edgesmarkers)
  facesmarkers <- as.numeric(facesmarkers)
  nodesmarkers <- as.numeric(nodesmarkers)

  out<-NULL

  if(order==1 && ncol(tetrahedrons) == 4){
    out = list(nodes=nodes, nodesmarkers=nodesmarkers, nodesattributes=nodesattributes,
        tetrahedrons=tetrahedrons, segments=segments, segmentsmarkers=segmentsmarkers,
          edges=edges, edgesmarkers=edgesmarkers, faces=faces, facesmarkers=facesmarkers,
            neighbors=neighbors, holes=holes, order=order)
  }
  else if(order==2 && ncol(tetrahedrons) == 10){ # tetrahedrons matrix contains both the true tetrahedrons and the midpoints ones
    out = list(nodes=nodes, nodesmarkers=nodesmarkers, nodesattributes=nodesattributes,
        tetrahedrons=tetrahedrons, segments=segments, segmentsmarkers=segmentsmarkers,
          edges=edges, edgesmarkers=edgesmarkers, faces=faces, facesmarkers=facesmarkers,
            neighbors=neighbors, holes=holes, order=order)
  }
  else if(order==2 && ncol(tetrahedrons) == 4){
    print("You set order=2 but passed a matrix of tetrahedrons with just 4 columns. The midpoints for each edge will be computed.")

    midpoints <- nodes[edges[,1],]/2+nodes[edges[,2],]/2
    
    tetrahedron_labels <- matrix(unlist(edge_list),ncol=4, byrow = TRUE)[,3:4]
    tetrahedron_labels[,2] <- tetrahedron_labels[,2]+4
    indices<-nrow(nodes)+cumsum(!repeated_edges)
    
    tetrahedrons<-cbind(tetrahedrons,matrix(0,nrow=nrow(tetrahedrons),ncol=6))
    tetrahedrons[tetrahedron_labels]<-indices
    
    nodes<-rbind(nodes,midpoints)
    nodesmarkers<-c(nodesmarkers,rep(0,nrow(midpoints)))
    
    out = list(nodes=nodes, nodesmarkers=nodesmarkers, nodesattributes=nodesattributes,
        tetrahedrons=tetrahedrons, segments=segments, segmentsmarkers=segmentsmarkers,
               edges=edges, edgesmarkers=edgesmarkers, faces=faces, facesmarkers=facesmarkers,
                  neighbors=neighbors, holes=holes, order=order)
  }
  else{
    stop("The number of columns of tetrahedrons matrix is not consistent with the order parameter")
  }

  class(out)<-"mesh.3D"

  return(out)
}
