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
Point<ndim> MeshHandlerCore<ORDER,mydim,ndim>::getPoint(Id id) const
{
	std::array<Real,ndim> coord;
	for(int i=0; i<ndim; ++i)
		coord[i]=points_[id+i*num_nodes_];
	return Point<ndim>(id, Identifier::NVAL, coord);
}

template <UInt ORDER, UInt mydim, UInt ndim>
Element<how_many_nodes(ORDER,mydim),mydim,ndim> MeshHandlerCore<ORDER,mydim,ndim>::getElement(Id id) const
{
	std::array<Point<ndim>, how_many_nodes(ORDER,mydim)> element_points;
	Id curr_point;
	for (int i=0; i<how_many_nodes(ORDER,mydim); ++i)
	{
		curr_point = elements_[i*num_elements_ + id];
		element_points[i]= this->getPoint(curr_point);
	}
	return meshElement(id, element_points);
}

template <UInt ORDER, UInt mydim, UInt ndim>
Element<how_many_nodes(ORDER,mydim),mydim,ndim> MeshHandlerCore<ORDER,mydim,ndim>::getNeighbors(Id id_element, UInt number) const
{
	Id id_neighbor = neighbors_[number * num_elements_ + id_element];
	//return empty element if "neighbor" not present (out of boundary!)
	return (id_neighbor==-1) ? meshElement() : this->getElement(id_neighbor);
}

template <UInt ORDER, UInt mydim, UInt ndim>
Element<how_many_nodes(ORDER,mydim),mydim,ndim> MeshHandlerCore<ORDER,mydim,ndim>::findLocationNaive(const Point<ndim> point) const
{
	meshElement current_element;
	for(Id id=0; id < num_elements_; ++id)
	{
		current_element = this->getElement(id);
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
	meshElement current_element = starting_element;
	int direction=0;

	//Test for found Element, or out of border
	while(current_element.getId() != Identifier::NVAL && !current_element.isPointInside(point))
	{
		direction = current_element.getPointDirection(point);
		current_element = this->getNeighbors(current_element.getId(), direction);
	}

	return current_element;
}

// template <UInt ORDER>
// std::vector<UInt> MeshHandler<ORDER,2,3>::find_closest(const std::vector<Point<3> > &inpoints) const{
// 	// Container for ID and distances of closest points
// 	std::vector<std::pair<UInt,Real> > closest_to;
// 	closest_to.reserve(inpoints.size());
//
// 	Point<3> curr_node{this->getPoint(0)};
// 	for(auto const point : inpoints)
// 		closest_to.emplace_back(0,distance(curr_node,point));
//
// 	for(int i=1; i<this->num_nodes(); ++i){
// 		curr_node=this->getPoint(i);
// 		for(int j=0; j<inpoints.size(); ++j){
// 			Real dist=distance(curr_node,inpoints[j]);
// 			if (dist<closest_to[j].second)
// 				closest_to[j]=std::make_pair(i,dist);
// 		}
// 	}
//
// 	std::vector<UInt> closest_ID;
// 	closest_ID.reserve(inpoints.size());
// 	for (auto const c : closest_to)
// 		closest_ID.push_back(c.first);
//
// 	return closest_ID;
// }
//

template <UInt ORDER>
std::vector<UInt> MeshHandler<ORDER,2,3>::find_closest(const std::vector<Point<3> > &inpoints) const{

	std::vector<UInt> closest_ID(inpoints.size());
	std::vector<Real> distances;
	distances.reserve(inpoints.size());

	Point<3> curr_node{this->getPoint(0)};
	for (auto const point : inpoints)
		distances.push_back(distance(curr_node,point));

	for(int i=1; i<this->num_nodes(); ++i){
		curr_node=this->getPoint(i);
		for(int j=0; j<inpoints.size(); ++j){
			Real dist=distance(curr_node,inpoints[j]);
			if (dist<distances[j]){
				distances[j]=dist;
				closest_ID[j]=i;
			}
		}
	}

	return closest_ID;
}



template <UInt ORDER>
std::vector<Point<3> > MeshHandler<ORDER,2,3>::project(const std::vector<Point<3> > &inpoints) const{

	const UInt num_elements = this->num_elements();
	// First find the closest node for each point
	std::vector<UInt> closest_nodes{this->find_closest(inpoints)};

	// Second build up a patch of elements on which to project
	std::vector<std::vector<UInt> > patches(inpoints.size());

	// Loop over the elements (first consider only true nodes)
	for(int i=0; i<3*num_elements; ++i){
		for(int j=0; j<closest_nodes.size(); ++j){
			if(closest_nodes[j]==this->elements_[i]){
				patches[j].push_back(i%num_elements);
				// Add to the patch the elements that lie opposite the closest node
				patches[j].push_back(this->neighbors_[i]);
			}
		}
	}
	// In the 2nd order case midpoints are also considered
	// Note: in the 1st order case this loop is skipped altogether
	for(int i=3*num_elements; i<this->num_elements_; ++i){
		for(int j=0; j<closest_nodes.size(); ++j){
			if(closest_nodes[j]==this->elements_[i]){
				patches[j].push_back(i%num_elements);
				// Add to the patch the neighboring elements
				patches[j].push_back(this->neighbors_[(i-num_elements)%(3*num_elements)]);
				patches[j].push_back(this->neighbors_[(i-2*num_elements)%(3*num_elements)]);
			}
		}
	}

	// Third compute the projections on the elements in each patch and keep the closest
	std::vector<Point<3> > projections(inpoints.size());
	Element<how_many_nodes(ORDER,2),2,3> current_element;
	Point<3> proj_point;

	for(int i=0; i<inpoints.size(); ++i){
		current_element = this->getElement(patches[i][0]);
		projections[i] = current_element.computeProjection(inpoints[i]);
		Real dist = distance(projections[i], inpoints[i]);

		for(int j=1; j<patches[i].size(); ++j){
			current_element = this->getElement(patches[i][j]);
			proj_point = current_element.computeProjection(inpoints[i]);
			if(distance(proj_point, inpoints[i])<dist){
				projections[i] = proj_point;
				dist = distance(projections[i], inpoints[i]);
			}
		}
	}

	return projections;

}


// Note: this print functions need to be checked!
template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandlerCore<ORDER,mydim,ndim>::printPoints(std::ostream & out)
{
	for(UInt i = 0; i < num_nodes_; ++i)
	{
		out<<"-"<< i <<"-"<<"("<<points_[i]<<","<<points_[num_nodes_+i]<<")"<<std::endl<<"------"<<std::endl;
	}
}


template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandlerCore<ORDER,mydim,ndim>::printElements(std::ostream & out)
{

	out << "# Triangles: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i )
	{
		out<<"-"<< i <<"- ";
		for( UInt k = 0; k < 3*ORDER; ++k)
			out<<elements_[k*num_elements_ + i]<<"   ";
		out<<std::endl;
	}

}

template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandlerCore<ORDER,mydim,ndim>::printNeighbors(std::ostream & out)
{

	out << "# Neighbors list: "<< num_elements_ <<std::endl;
	for (UInt i = 0; i < num_elements_; ++i )
	{
		out<<"-"<< i <<"- ";
		for( UInt k = 0; k < 3; ++k)
			out<<neighbors_[k*num_elements_ + i]<<"   ";
		out<<std::endl;
	}

}

#endif
