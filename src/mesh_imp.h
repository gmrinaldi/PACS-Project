#ifndef MESH_IMP_H_
#define MESH_IMP_H_

#ifdef R_VERSION_
template <UInt ORDER, UInt mydim, UInt ndim>
MeshHandlerCore<ORDER,mydim,ndim>::MeshHandlerCore(SEXP mesh)
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
MeshHandlerCore<ORDER,mydim,ndim>::~MeshHandlerCore() {}

template <UInt ORDER, UInt mydim, UInt ndim>
inline Point<ndim> MeshHandlerCore<ORDER,mydim,ndim>::getPoint(Id id) const
{
	return Point<ndim>(id, points_, this->num_nodes());
}

template <UInt ORDER, UInt mydim, UInt ndim>
Element<how_many_nodes(ORDER,mydim),mydim,ndim> MeshHandlerCore<ORDER,mydim,ndim>::getElement(Id id) const
{
	std::array<Point<ndim>, how_many_nodes(ORDER,mydim)> element_points;
	for (int i=0; i<how_many_nodes(ORDER,mydim); ++i)
	{
		Id curr{elements_[id + i*num_elements_]};
		element_points[i]= this->getPoint(curr);
	}
	return meshElement(id, element_points);
}

template <UInt ORDER, UInt mydim, UInt ndim>
Element<how_many_nodes(ORDER,mydim),mydim,ndim> MeshHandlerCore<ORDER,mydim,ndim>::getNeighbors(Id id_element, UInt number) const
{
	Id id_neighbor{neighbors_[id_element + number * num_elements_]};
	//return empty element if "neighbor" not present (out of boundary!)
	return (id_neighbor==-1) ? meshElement() : this->getElement(id_neighbor);
}

template <UInt ORDER, UInt mydim, UInt ndim>
Element<how_many_nodes(ORDER,mydim),mydim,ndim> MeshHandlerCore<ORDER,mydim,ndim>::findLocationNaive(const Point<ndim>& point) const
{
	for(Id id=0; id < num_elements_; ++id){
		meshElement current_element{this->getElement(id)};
		if(current_element.isPointInside(point))
			return current_element;
	}
	return meshElement(); //default element with NVAL ID
}

// Visibility walk algorithm which uses barycentric coordinate [Sundareswara et al]
//Starting triangles usually n^(1/3) points
template <UInt ORDER, UInt mydim, UInt ndim>
Element<how_many_nodes(ORDER,mydim),mydim,ndim> MeshHandler<ORDER,mydim,ndim>::findLocationWalking(const Point<ndim>& point, const Element<how_many_nodes(ORDER,mydim),mydim,ndim>& starting_element) const
{
	static_assert(ndim==mydim,
								"ERROR: Walking algorithm does not work for manifold data! see mesh_imp.h");
	//Walking algorithm to the point
	meshElement current_element{starting_element};

	//Test for found Element, or out of border
	while(current_element.hasValidId() && !current_element.isPointInside(point))
		current_element = this->getNeighbors(current_element.getId(), current_element.getPointDirection(point));

	return current_element;
}

// This function finds the closest mesh nodes to a given set of 3D points.
template <UInt ORDER>
std::vector<UInt> MeshHandler<ORDER,2,3>::find_closest(const std::vector<Point<3> > &points) const{

	std::vector<UInt> closest_ID;
	closest_ID.reserve(points.size());

	//exclude midpoints in order 2
	const UInt num_actual_nodes = (ORDER==1) ? this->num_nodes() : this->num_nodes() - this->num_edges();

	for(auto const &point : points){
		UInt min_pos;
		Real min_dist{std::numeric_limits<Real>::max()};
		for(int i=0; i<num_actual_nodes; ++i){
			Real distance{point.dist2(this->getPoint(i))};
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
template <UInt ORDER>
std::vector<Point<3> > MeshHandler<ORDER,2,3>::project(const std::vector<Point<3> > &points) const{

	std::vector<Point<3> > projections;
	projections.reserve(points.size());

	for(auto const &point : points){
		// First find the closest node for the current point
		UInt closest_node{this->find_closest(point)};

		// Second build up a patch of elements onto which to project
		std::vector<UInt> patch;
		// Conservative estimate to avoid reallocation in most cases
		patch.reserve(18);
		for(int i=0; i<3*this->num_elements(); ++i)
			if(closest_node==this->elements_[i])
				patch.push_back(i%this->num_elements());

		// Third compute the projections on the elements in the patch and keep the closest
		Real min_dist{std::numeric_limits<Real>::max()};
		Point<3> proj_point;
		for(auto const &i : patch){
			Element<how_many_nodes(ORDER,2),2,3> current_element{this->getElement(i)};
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
void MeshHandlerCore<ORDER,mydim,ndim>::printPoints(std::ostream& os)
{
	os<<"# Nodes: "<<num_nodes_<<std::endl;
	for(UInt i=0; i<num_nodes_; ++i)
		os<<getPoint(i);
}


template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandlerCore<ORDER,mydim,ndim>::printElements(std::ostream& os)
{
	os << "# Triangles: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i )
		os<<getElement(i);
}

template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandlerCore<ORDER,mydim,ndim>::printNeighbors(std::ostream& os)
{
	os << "# Neighbors list: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i ){
		for( UInt k = 0; k < mydim+1; ++k)
			os<<neighbors_[i+k*num_elements_]<<" ";
		os<<std::endl;
	}
}

#endif
