//#include "mesh_objects.hpp"
#ifndef __MESH_OBJECTS_IMP_HPP__
#define __MESH_OBJECTS_IMP_HPP__

template <UInt NNODES, UInt mydim, UInt ndim>
Eigen::Matrix<Real,mydim+1,1> ElementCore<NNODES,mydim,ndim>::getBaryCoordinates(const Point<ndim>& point) const
{
	ElementCore<NNODES,mydim,ndim> t = *this;
	Eigen::Matrix<Real,mydim+1,1> lambda;

	std::array<Real,ndim> diff = point - t[0];

	lambda.tail(mydim) = M_invJ_ * Eigen::Map<Eigen::Matrix<Real,ndim,1> >(diff.data());

  lambda(0) = 1 - lambda.tail(mydim).sum();

	return lambda;

}

template <UInt NNODES, UInt mydim, UInt ndim>
void ElementCore<NNODES,mydim,ndim>::print(std::ostream & out) const
{
	out<<"ElementCore -"<< id_ <<"- ";
	for (UInt i=0; i<NNODES; ++i)
		out<<points_[i].getId()<<"  ";
	out<<std::endl;
}


// This function is called to construct elements in 2D and 3D
template <UInt NNODES, UInt mydim, UInt ndim>
void Element<NNODES,mydim,ndim>::computeProperties()
{
	ElementCore<NNODES,mydim,ndim> &t = *this;
	std::array<Real,ndim> diff;

	for (int i=0; i<mydim; ++i){
			diff=t[i+1]-t[0];
			this->M_J_.col(i) = Eigen::Map<Eigen::Matrix<Real,ndim,1> >(diff.data());
	}
	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!
	this->detJ_ = this->M_J_.determinant();
	this->M_invJ_ = this->M_J_.inverse();
	this->metric_ = this->M_invJ_ * this->M_invJ_.transpose();
}


//Primary template member definition
template <UInt NNODES, UInt mydim, UInt ndim>
bool Element<NNODES,mydim,ndim>::isPointInside(const Point<ndim>& point) const
{
	static constexpr Real eps = 	std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,mydim+1,1> lambda = this->getBaryCoordinates(point);

	return (-tolerance<=lambda.array()).all();

}



template <UInt NNODES, UInt mydim, UInt ndim>
int Element<NNODES,mydim,ndim>::getPointDirection(const Point<ndim>& point) const
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,mydim+1,1> lambda = this->getBaryCoordinates(point);

	//Find the minimum coordinate (if negative stronger straight to the point searched)
	int min_index;
	lambda.minCoeff(&min_index);

	return (lambda[min_index] < -tolerance)	? min_index : -1;
}

// Template member specializations for 2.5D

template <UInt NNODES>
void Element<NNODES,2,3>::computeProperties()
{
	ElementCore<NNODES,2,3> &t = *this;
	std::array<Real,3> diff;

	for (int i=0; i<2; ++i){
			diff=t[i+1]-t[0];
			this->M_J_.col(i) = Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());
	}

	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!

	Eigen::Matrix<Real,2,2> G_J_=this->M_J_.transpose()*this->M_J_;
	this->detJ_ = G_J_.determinant();
	this->metric_ = G_J_.inverse();
	this->M_invJ_ = this->metric_ * this->M_J_.transpose();
}


template <UInt NNODES>
bool Element<NNODES,2,3>::isPointInside(const Point<3>& point) const
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Element<NNODES,2,3> t=*this;

 	Eigen::Matrix<Real,3,3> A;

	std::array<Real,3> diff = point - t[0];

 	A << M_J_, Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());

 	// NOTE: this method is as fast as ColPivHouseholderQR for such small matrices
 	// but this is optimized for rank computations (see eigen documentation)
 	// Check if point belongs to plane!
 	if(A.fullPivHouseholderQr().isInvertible())
			return false;

	Eigen::Matrix<Real,3,1> lambda = this->getBaryCoordinates(point);

	return (-tolerance<=lambda.array()).all();

}
template <UInt NNODES>
Point<3> Element<NNODES,2,3>::computeProjection(const Point<3>& point) const
{
	// Note: no need for tolerances here

	if(isPointInside(point))
		return point;

	Eigen::Matrix<Real,3,1> lambda = this->getBaryCoordinates(point);
	// Convention: (+,+,+) means that all lambda are positive and so on
  // For visual reference: (remember that edges are numbered wrt the node in front)
  //
	//\ (-,-,+)|
  //  \      |
  //    \    |
  //      \  |
  //        \|
  //         |3
  //         | \
  //         |   \
  // (+,-,+) |     \       (-,+,+)
  //         |       \
  //         |         \
  //         |           \
  //         |  (+,+,+)    \
  //  _____1 |_______________\ 2 _____________
  //         |                 \
  // (+,-,-) |    (+,+,-)        \  (-,+,-)
  //         |                     \

	// If (+,-,-) the projection lies beyond node 1
	// Simply return node 1 (same for the others)
	if(lambda(0)>0 && lambda(1)<0 && lambda(2)<0)
		return this->points_[0];
	else if (lambda(0)<0 && lambda(1)>0 && lambda(2)<0)
		return this->points_[1];
	else if (lambda(0)<0 && lambda(1)<0 && lambda(2)>0)
		return this->points_[2];

	Eigen::Matrix<Real,3,1> coords3D;

	// If (+,+,-) the projection lies beyond edge 3
	// Simply scale it back on the edge
	if(lambda(0)>0 && lambda(1)>0)
		coords3D = this->M_J_.col(0)/lambda.head(2).sum();
	// If (+,-,+) the projection lies beyond edge 2
	else if (lambda(0)>0 && lambda(2)>0)
		coords3D = this->M_J_.col(1)/(lambda(0)+lambda(2));
	// If (+,+,+) the projection lies inside the element
	// So just convert back to 3D coords
	if(lambda(0)>0 && lambda(1)>0 && lambda(2)>0)
		coords3D = this->M_J_ * lambda.tail(2);
	// If (-,+,+) the projection lies beyond edge 1
	// Scale back and convert
	else
		coords3D = this->M_J_ * lambda.tail(2)/lambda.tail(2).sum();

	//Add point1 for translation and return
	std::array<Real,3> proj;

	for(int i=0; i<3; ++i)
		proj[i] = coords3D[i] + this->points_[0][i];

	return Point<3>(proj);
}


#endif
