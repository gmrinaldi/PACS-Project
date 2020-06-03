#ifndef MESH_IMP_H_
#define MESH_IMP_H_

#ifdef R_VERSION_
template <UInt ORDER, UInt mydim, UInt ndim>
MeshHandler<ORDER,mydim,ndim>::MeshHandler(SEXP mesh)
{
	mesh_ 		= mesh;
	points_ 	= REAL(VECTOR_ELT(mesh_, 0));
	sides_ 		= INTEGER(VECTOR_ELT(mesh_, 6));
	elements_  = INTEGER(VECTOR_ELT(mesh_, 3));
	neighbors_  = INTEGER(VECTOR_ELT(mesh_, 8));

	num_nodes_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 0), R_DimSymbol))[0];
	num_sides_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 6), R_DimSymbol))[0];
	num_elements_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 3), R_DimSymbol))[0];

}
#endif

template <UInt ORDER, UInt mydim, UInt ndim>
inline Point<ndim> MeshHandler<ORDER,mydim,ndim>::getPoint(UInt id) const
{
	return Point<ndim>(id, points_, num_nodes_);
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::getElement(UInt id) const
{
	typename meshElement::elementPoints elPoints;
	for (int i=0; i<how_many_nodes(ORDER,mydim); ++i)
		elPoints[i] = getPoint(elements_[id + i*num_elements_]);
	return meshElement(id,elPoints);
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::getNeighbors(UInt id_element, UInt number) const
{
	UInt id_neighbor{neighbors_[id_element + number * num_elements_]};
	//return empty element if "neighbor" not present (out of boundary!)
	return (id_neighbor==-1) ? meshElement() : getElement(id_neighbor);
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::findLocationNaive(const Point<ndim>& point) const
{
	for(UInt id=0; id < num_elements_; ++id){
		meshElement current_element{getElement(id)};
		if(current_element.isPointInside(point))
			return current_element;
	}
	return meshElement(); //default element with NVAL ID
}

// Visibility walk algorithm which uses barycentric coordinate [Sundareswara et al]
//Starting triangles usually n^(1/3) points
template <UInt ORDER, UInt mydim, UInt ndim>
template <UInt m, UInt n>																																		//vvvvvvvvv actual return type if enabled
typename std::enable_if<n==m && n==ndim && m==mydim, typename MeshHandler<ORDER,mydim,ndim>::meshElement>::type
MeshHandler<ORDER,mydim,ndim>::findLocationWalking(const Point<ndim>& point, const Element<how_many_nodes(ORDER,mydim),mydim,ndim>& starting_element) const
{
	meshElement current_element{starting_element};
	//Test for found Element, or out of border
	while(current_element.hasValidId() && !current_element.isPointInside(point))
		current_element = getNeighbors(current_element.getId(), current_element.getPointDirection(point));

	return current_element;
}

// This function finds the closest mesh nodes to a given set of 3D points.
template <UInt ORDER, UInt mydim, UInt ndim>
template <UInt m, UInt n>																//vvvvvvvvv actual return type
typename std::enable_if<(m<n) && n==ndim && m==mydim, std::vector<UInt> >::type
MeshHandler<ORDER,mydim,ndim>::find_closest(const std::vector<Point<3> > &points) const{

	std::vector<UInt> closest_ID;
	closest_ID.reserve(points.size());

	//exclude midpoints in order 2
	const UInt num_actual_nodes = (ORDER==1) ? num_nodes_ : num_nodes_ - num_sides_;

	for(auto const &point : points){
		UInt min_pos;
		Real min_dist{std::numeric_limits<Real>::max()};
		for(int i=0; i<num_actual_nodes; ++i){
			Real distance{point.dist2(getPoint(i))};
			if(distance<min_dist){
				min_dist=distance;
				min_pos=i;
			}
		}
		closest_ID.push_back(min_pos);
	}
	return closest_ID;
}

// Naive function for projection onto surface
template <UInt ORDER, UInt mydim, UInt ndim>
template <UInt m, UInt n>																//vvvvvvvvv actual return type
typename std::enable_if<m!=n && n==ndim && m==mydim, std::vector<Point<3> > >::type
MeshHandler<ORDER,mydim,ndim>::project(const std::vector<Point<3> > &points) const{

	std::vector<Point<3> > projections;
	projections.reserve(points.size());

	for(auto const &point : points){
		// First find the closest node for the current point
		UInt closest_node{find_closest(point)};

		// Second build up a patch of elements onto which to project
		std::vector<UInt> patch;
		// Conservative estimate to avoid reallocation in most cases
		patch.reserve(18);
		for(int i=0; i<3*num_elements_; ++i)
			if(closest_node==elements_[i])
				patch.push_back(i%num_elements_);

		// Third compute the projections on the elements in the patch and keep the closest
		Real min_dist{std::numeric_limits<Real>::max()};
		Point<3> proj_point;
		for(auto const &i : patch){
			Element<how_many_nodes(ORDER,2),2,3> current_element{getElement(i)};
			Point<3> curr_proj{current_element.computeProjection(point)};
			Real distance{point.dist2(curr_proj)};
			if(distance<min_dist){
				proj_point = curr_proj;
				min_dist = distance;
			}
		}
		projections.push_back(proj_point);
	}
	return projections;
}



template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printPoints(std::ostream& os)
{
	os<<"# Nodes: "<<num_nodes_<<std::endl;
	for(UInt i=0; i<num_nodes_; ++i)
		os<<getPoint(i);
}


template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printElements(std::ostream& os)
{
	os << "# Triangles: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i )
		os<<getElement(i);
}

template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printNeighbors(std::ostream& os)
{
	os << "# Neighbors list: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i ){
		for( UInt k = 0; k < mydim+1; ++k)
			os<<neighbors_[i+k*num_elements_]<<" ";
		os<<std::endl;
	}
}

#endif
