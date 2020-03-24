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
Point<ndim> MeshHandler<ORDER,mydim,ndim>::getPoint(Id id) const
{
	std::array<Real,ndim> coord;
	for(int i=0; i<ndim; ++i)
		coord[i]=points_[id+i*num_nodes_];
	return Point<ndim>(id, Identifier::NVAL, coord);
}

template <UInt ORDER, UInt mydim, UInt ndim>
Element<how_many_nodes(ORDER,mydim),mydim,ndim> MeshHandler<ORDER,mydim,ndim>::getElement(Id id) const
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
Element<how_many_nodes(ORDER,mydim),mydim,ndim> MeshHandler<ORDER,mydim,ndim>::getNeighbors(Id id_element, UInt number) const
{
	Id id_neighbor = neighbors_[number * num_elements_ + id_element];
	//return empty element if "neighbor" not present (out of boundary!)
	return (id_neighbor==-1) ? meshElement() : this->getElement(id_neighbor);
}

template <UInt ORDER, UInt mydim, UInt ndim>
Element<how_many_nodes(ORDER,mydim),mydim,ndim> MeshHandler<ORDER,mydim,ndim>::findLocationNaive(const Point<ndim> point) const
{
	meshElement current_element;
	for(Id id=0; id < num_elements_; ++id)
	{
		current_element = getElement(id);
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
		current_element = getNeighbors(current_element.getId(), direction);
	}

	return current_element;
}


// Note: this print functions need to be checked!
template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printPoints(std::ostream & out)
{
	for(UInt i = 0; i < num_nodes_; ++i)
	{
		out<<"-"<< i <<"-"<<"("<<points_[i]<<","<<points_[num_nodes_+i]<<")"<<std::endl<<"------"<<std::endl;
	}
}


template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printElements(std::ostream & out)
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
void MeshHandler<ORDER,mydim,ndim>::printNeighbors(std::ostream & out)
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
