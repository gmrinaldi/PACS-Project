//#include "mesh_objects.hpp"
#ifndef __MESH_OBJECTS_IMP_HPP__
#define __MESH_OBJECTS_IMP_HPP__

// Member distance funtion for class Point
template <UInt ndim>
Real Point<ndim>::distance(const Point<ndim> &other) const {
	Real dist{0};
	pointCoords diff{*this - other};
	for (auto const &i : diff)
		dist+=i*i;
	return std::sqrt(dist);
}

// From here on implementation of Element/ElementCore
template <UInt NNODES, UInt mydim, UInt ndim>
Eigen::Matrix<Real,mydim+1,1> ElementCore<NNODES,mydim,ndim>::getBaryCoordinates(const Point<ndim>& point) const
{
	Eigen::Matrix<Real,mydim+1,1> lambda;

	std::array<Real,ndim> diff = point - points_[0];

	lambda.tail(mydim).noalias() = M_invJ_ * Eigen::Map<Eigen::Matrix<Real,ndim,1> >(diff.data());

  lambda(0) = 1 - lambda.tail(mydim).sum();

	return lambda;

}

//Primary template member definition
template <UInt NNODES, UInt mydim, UInt ndim>
bool ElementCore<NNODES,mydim,ndim>::isPointInside(const Point<ndim>& point) const
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,mydim+1,1> lambda = getBaryCoordinates(point);

	return (-tolerance<=lambda.array()).all();

}

// This function is called to construct elements in 2D and 3D
template <UInt NNODES, UInt mydim, UInt ndim>
void Element<NNODES,mydim,ndim>::computeProperties()
{
	for (int i=0; i<mydim; ++i){
			std::array<Real,ndim> diff(this->points_[i+1] - this->points_[0]);
			this->M_J_.col(i) = Eigen::Map<Eigen::Matrix<Real,ndim,1> >(diff.data());
	}
	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!
	this->detJ_ = this->M_J_.determinant();
	this->M_invJ_ = this->M_J_.inverse();
	this->metric_.noalias() = this->M_invJ_ * this->M_invJ_.transpose();
}


template <UInt NNODES, UInt mydim, UInt ndim>
int Element<NNODES,mydim,ndim>::getPointDirection(const Point<ndim>& point) const
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,mydim+1,1> lambda = this->getBaryCoordinates(point);

	//Find the minimum coordinate (if negative stronger straight to the point searched)
	UInt min_index;
	lambda.minCoeff(&min_index);

	return (lambda[min_index] < -tolerance)	? min_index : -1;
}

// Template member specializations for 2.5D

template <UInt NNODES>
void Element<NNODES,2,3>::computeProperties()
{
	for (int i=0; i<2; ++i){
		  std::array<Real,3> diff(this->points_[i+1] - this->points_[0]);
			this->M_J_.col(i) = Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());
	}

	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!

	Eigen::Matrix<Real,2,2> G_J_=this->M_J_.transpose()*this->M_J_;
	this->detJ_ = G_J_.determinant();
	this->metric_ = G_J_.inverse();
	this->M_invJ_.noalias() = this->metric_ * this->M_J_.transpose();
}


template <UInt NNODES>
bool Element<NNODES,2,3>::isPointInside(const Point<3>& point) const
{
	std::array<Real,3> diff(point - this->points_[0]);

 	Eigen::Matrix<Real,3,3> A;
 	A << this->M_J_, Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());

 	// NOTE: this method is as fast as ColPivHouseholderQR for such small matrices
 	// but this is optimized for rank computations (see eigen documentation)
	return !A.fullPivHouseholderQr().isInvertible() && ElementCore<NNODES,2,3>::isPointInside(point);

}
template <UInt NNODES>
Point<3> Element<NNODES,2,3>::computeProjection(const Point<3>& point) const
{
	// Note: no need for tolerances here

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
	// If (+,+,+) the projection lies inside the element
	// So just convert back to 3D coords
	if(lambda(0)>0 && lambda(1)>0 && lambda(2)>0)
		coords3D = this->M_J_ * lambda.tail<2>();
	// If (+,+,-) the projection lies beyond edge 3
	// Simply scale it back on the edge and convert
	else if(lambda(0)>0 && lambda(1)>0)
		coords3D = this->M_J_.col(0)/lambda.head<2>().sum();
	// If (+,-,+) the projection lies beyond edge 2
	else if (lambda(0)>0 && lambda(2)>0)
		coords3D = this->M_J_.col(1)/(lambda(0)+lambda(2));
	// If (-,+,+) the projection lies beyond edge 1
	else
		coords3D = this->M_J_ * lambda.tail<2>()/lambda.tail<2>().sum();

	// Translate back by the first point and return
	return Point<3>({coords3D[0],coords3D[1],coords3D[2]}) + (this->points_[0]);
}


// This covers all order 1 cases
template <UInt NNODES, UInt mydim, UInt ndim>
inline Real ElementCore<NNODES,mydim,ndim>::evaluate_point(const Point<ndim>& point, const Eigen::Matrix<Real,NNODES,1>& coefficients) const
{
  return(coefficients.dot(this->getBaryCoordinates(point)));
}

// Full specialization for order 2 in 2D
template <>
inline Real ElementCore<6,2,2>::evaluate_point(const Point<2>& point, const Eigen::Matrix<Real,6,1>& coefficients) const
{
	Eigen::Matrix<Real,3,1> lambda = this->getBaryCoordinates(point);
  return( coefficients[0]*lambda[0]*(2*lambda[0] - 1) +
            coefficients[1]*lambda[1]*(2*lambda[1] - 1) +
            coefficients[2]*lambda[2]*(2*lambda[2] - 1) +
            coefficients[3]*4*lambda[1]*lambda[2] +
            coefficients[4]*4*lambda[2]*lambda[0] +
            coefficients[5]*4*lambda[0]*lambda[1] );
}

 // Full specialization for order 2 in 2.5D
template <>
inline Real ElementCore<6,2,3>::evaluate_point(const Point<3>& point, const Eigen::Matrix<Real,6,1>& coefficients) const
{
	Eigen::Matrix<Real,3,1> lambda = this->getBaryCoordinates(point);
	return( coefficients[0]*lambda[0]*(2*lambda[0] - 1) +
            coefficients[1]*lambda[1]*(2*lambda[1] - 1) +
            coefficients[2]*lambda[2]*(2*lambda[2] - 1) +
            coefficients[3]*4*lambda[1]*lambda[2] +
            coefficients[4]*4*lambda[2]*lambda[0] +
            coefficients[5]*4*lambda[0]*lambda[1] );
}

// Full specialization for order 2 in 3D
// MEMO: this works assuming edges are ordered like so: (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
template <>
inline Real ElementCore<10,3,3>::evaluate_point(const Point<3>& point, const Eigen::Matrix<Real,10,1>& coefficients) const
{
 Eigen::Matrix<Real,4,1> lambda = this->getBaryCoordinates(point);
 return( coefficients[0]*lambda[0]*(2*lambda[0] - 1) +
           coefficients[1]*lambda[1]*(2*lambda[1] - 1) +
           coefficients[2]*lambda[2]*(2*lambda[2] - 1) +
           coefficients[3]*lambda[3]*(2*lambda[3] - 1) +
           coefficients[4]*(4*lambda[1]*lambda[0]) +
           coefficients[5]*(4*lambda[2]*lambda[0]) +
           coefficients[6]*(4*lambda[3]*lambda[0]) +
           coefficients[7]*(4*lambda[1]*lambda[2]) +
           coefficients[8]*(4*lambda[2]*lambda[3]) +
           coefficients[9]*(4*lambda[3]*lambda[1]) );
}


// This covers all order 1 cases
template <UInt NNODES, UInt mydim, UInt ndim>
inline Eigen::Matrix<Real,ndim,1> ElementCore<NNODES,mydim,ndim>::evaluate_der_point(const Point<ndim>& point, const Eigen::Matrix<Real,NNODES,1>& coefficients) const
{
  Eigen::Matrix<Real,mydim,mydim+1> B1;
  B1.col(0).setConstant(-1);
  B1.rightCols(mydim).setIdentity();

  return(this->M_invJ_.transpose()*B1*coefficients);
}

// Full specialization for order 2 in 2D
template <>
inline Eigen::Matrix<Real,2,1> ElementCore<6,2,2>::evaluate_der_point(const Point<2>& point, const Eigen::Matrix<Real,6,1>& coefficients) const
{
  // Multiply by 4 for convenience in assemblying B2
	Eigen::Matrix<Real,3,1> bar4 = 4*this->getBaryCoordinates(point);
	Eigen::Matrix<Real,2,6> B2;
  B2 << 1-bar4(0), bar4(1)-1,           0, bar4(2),        -bar4(2), bar4(0)-bar4(1),
        1-bar4(0),           0, bar4(2)-1, bar4(1), bar4(0)-bar4(2),        -bar4(1);
	return(this->M_invJ_.transpose()*B2*coefficients);
}

// Full specialization for order 2 in 2.5D
template <>
inline Eigen::Matrix<Real,3,1> ElementCore<6,2,3>::evaluate_der_point(const Point<3>& point, const Eigen::Matrix<Real,6,1>& coefficients) const
{
  // Multiply by 4 for convenience in assemblying B2
	Eigen::Matrix<Real,3,1> bar4 = 4*this->getBaryCoordinates(point);
	Eigen::Matrix<Real,2,6> B2;
  B2 << 1-bar4(0), bar4(1)-1,           0, bar4(2),        -bar4(2), bar4(0)-bar4(1),
        1-bar4(0),           0, bar4(2)-1, bar4(1), bar4(0)-bar4(2),        -bar4(1);
	return(this->M_invJ_.transpose()*B2*coefficients);
}

// Full specialization for order 2 in 3D
// MEMO: this works assuming edges are ordered like so: (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
template <>
inline Eigen::Matrix<Real,3,1> ElementCore<10,3,3>::evaluate_der_point(const Point<3>& point, const Eigen::Matrix<Real,10,1>& coefficients) const
{
  // Multiply by 4 for convenience in assemblying B2
	Eigen::Matrix<Real,4,1> bar4 = 4*this->getBaryCoordinates(point);
	Eigen::Matrix<Real,3,10> B2;
  B2 << 1-bar4(0),  bar4(1)-1,          0,         0, bar4(0)-bar4(1),        -bar4(2),        -bar4(3),  bar4(2),        0, bar4(3),
        1-bar4(0),           0, bar4(2)-1,         0,        -bar4(1), bar4(0)-bar4(2),        -bar4(3),  bar4(1),  bar4(3),       0,
        1-bar4(0),           0,         0, bar4(3)-1,        -bar4(1),        -bar4(2), bar4(0)-bar4(3),        0,  bar4(2), bar4(1);

	return(this->M_invJ_.transpose()*B2*coefficients);
}


// This covers all order 1 cases
template <UInt NNODES, UInt mydim, UInt ndim>
inline Real ElementCore<NNODES,mydim,ndim>::integrate(const Eigen::Matrix<Real,NNODES,1>& coefficients) const
{
	// This works because of linearity!
	return this->getMeasure() * coefficients.mean();
}

// Full specialization for order 2 in 2D
template <>
inline Real ElementCore<6,2,2>::integrate(const Eigen::Matrix<Real,6,1>& coefficients) const
{
	// Evaluate the function on midpoints and weight each term equally
	return this->getMeasure() * coefficients.tail<3>().mean();
}

// Full specialization for order 2 in 2.5D
template <>
inline Real ElementCore<6,2,3>::integrate(const Eigen::Matrix<Real,6,1>& coefficients) const
{
	// Evaluate the function on midpoints and weight each term equally
	return this->getMeasure() * coefficients.tail<3>().mean();
}

// Full specialization for order 2 in 3D
// MEMO: this works assuming edges are ordered like so: (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
template <>
inline Real ElementCore<10,3,3>::integrate(const Eigen::Matrix<Real,10,1>& coefficients) const
{
	// In this case a more complicated integration scheme is needed!
	static constexpr Real shape_fun[]={0, -1./9, -1./9, -1./9, 1./3, 1./3, 1./3, 1./9, 1./9, 1./9,
																		-1./9, 0, -1./9, -1./9, 1./3, 1./9, 1./9, 1./3, 1./9, 1./3,
																		-1./9, -1./9, 0, -1./9, 1./9, 1./3, 1./9, 1./3, 1./3, 1./9,
																		-1./9, -1./9, -1./9, 0, 1./9, 1./9, 1./3, 1./9, 1./3, 1./3};

	return this->getMeasure() * (-.8*(-.125*coefficients.head<4>().sum()+.25*coefficients.tail<6>().sum())+
			.45*coefficients.replicate<4,1>().dot(Eigen::Map<const Eigen::Matrix<Real,40,1> >(shape_fun)));
}





#endif
