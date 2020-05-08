//#include "mesh_objects.hpp"
#ifndef __MESH_OBJECTS_IMP_HPP__
#define __MESH_OBJECTS_IMP_HPP__

// Member functions for class Point
// Constructor specializations
template <>
inline Point<2>::Point(Id id, const Real* const points, const UInt num_points) :
		Point(id, {points[id], points[id+num_points]}) {}

template <>
inline Point<3>::Point(Id id, const Real* const points, const UInt num_points) :
		Point(id, {points[id], points[id+num_points], points[id+2*num_points]}) {}

template <>
constexpr Point<2>::Point(const Real(&coord)[2]) :
		coord_({coord[0], coord[1]}) {}

template <>
constexpr Point<3>::Point(const Real(&coord)[3]) :
		coord_({coord[0], coord[1], coord[2]}) {}

template <>
inline Point<2>::Point(const pointEigen &coord) :
		coord_({coord[0], coord[1]}) {}

template <>
inline Point<3>::Point(const pointEigen &coord) :
		coord_({coord[0], coord[1], coord[2]}) {}

// This function returns the squared euclidean distance between "this" point and other
template <UInt ndim>
inline Real Point<ndim>::dist2(const Point<ndim> &other) const {
	Real dist2{0.};
	for (UInt i=0; i<ndim; ++i)
		dist2+=(coord_[i]-other[i])*(coord_[i]-other[i]);
	return dist2;
}

// This function returns the euclidean distance between "this" point and other
template <UInt ndim>
inline Real Point<ndim>::dist(const Point<ndim> &other) const {
	return std::sqrt(this->dist2(other));
}

// Overloaded +=/-= operators
template <UInt ndim>
inline Point<ndim>& Point<ndim>::operator+=(const Point &other){
	for (UInt i=0; i<ndim; ++i)
		coord_[i]+=other[i];
	return *this;
}

template <UInt ndim>
inline Point<ndim>& Point<ndim>::operator-=(const Point &other){
	for (UInt i=0; i<ndim; ++i)
		coord_[i]-=other[i];
	return *this;
}

// From here on implementation of Element/ElementCore
template <UInt NNODES, UInt mydim, UInt ndim>
Eigen::Matrix<Real,mydim+1,1> ElementCore<NNODES,mydim,ndim>::getBaryCoordinates(const Point<ndim> &point) const
{
	Eigen::Matrix<Real,mydim+1,1> lambda;

	// .template is needed! See Eigen documentation regarding
	// the template and typename keywords in C++
	lambda.template tail<mydim>().noalias() = M_invJ_ * (EigenMap2Const_t(&point[0])-EigenMap2Const_t(&points_[0][0]));

  lambda(0) = 1 - lambda.template tail<mydim>().sum();

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
	{
		EigenMap2Const_t basePoint(&(this->points_[0][0]));
		for (int i=0; i<mydim; ++i)
			this->M_J_.col(i) = EigenMap2Const_t(&(this->points_[i+1][0]))-basePoint;
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
	{
		EigenMap2Const_t basePoint(&(this->points_[0][0]));
		for (int i=0; i<2; ++i)
			this->M_J_.col(i) = EigenMap2Const_t(&(this->points_[i+1][0]))-basePoint;
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
 	Eigen::Matrix<Real,3,3> A;
 	A << this->M_J_, (EigenMap2Const_t(&point[0])-EigenMap2Const_t(&(this->points_[0][0])));

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
	return Point<3>(coords3D)+=(this->points_[0]);
}


// This covers all order 1 cases
template <UInt NNODES, UInt mydim, UInt ndim>
inline Real ElementCore<NNODES,mydim,ndim>::evaluate_point(const Point<ndim>& point, const Eigen::Matrix<Real,NNODES,1>& coefficients) const
{
  return coefficients.dot(this->getBaryCoordinates(point));
}

// Full specialization for order 2 in 2D
template <>
inline Real ElementCore<6,2,2>::evaluate_point(const Point<2>& point, const Eigen::Matrix<Real,6,1>& coefficients) const
{
	Eigen::Matrix<Real,3,1> lambda = this->getBaryCoordinates(point);
  return coefficients[0] * lambda[0] * (2*lambda[0]-1) +
         coefficients[1] * lambda[1] * (2*lambda[1]-1) +
         coefficients[2] * lambda[2] * (2*lambda[2]-1) +
         coefficients[3] * 4 * lambda[1] * lambda[2] +
         coefficients[4] * 4 * lambda[2] * lambda[0] +
         coefficients[5] * 4 * lambda[0] * lambda[1];
}

 // Full specialization for order 2 in 2.5D
template <>
inline Real ElementCore<6,2,3>::evaluate_point(const Point<3>& point, const Eigen::Matrix<Real,6,1>& coefficients) const
{
	Eigen::Matrix<Real,3,1> lambda = this->getBaryCoordinates(point);
	return coefficients[0] * lambda[0] * (2*lambda[0]-1) +
         coefficients[1] * lambda[1] * (2*lambda[1]-1) +
         coefficients[2] * lambda[2] * (2*lambda[2]-1) +
         coefficients[3] * 4 * lambda[1]*lambda[2] +
         coefficients[4] * 4 * lambda[2]*lambda[0] +
         coefficients[5] * 4 * lambda[0]*lambda[1];
}

// Full specialization for order 2 in 3D
// MEMO: this works assuming edges are ordered like so: (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
template <>
inline Real ElementCore<10,3,3>::evaluate_point(const Point<3>& point, const Eigen::Matrix<Real,10,1>& coefficients) const
{
 Eigen::Matrix<Real,4,1> lambda = this->getBaryCoordinates(point);
 return coefficients[0] * lambda[0] * (2*lambda[0]-1) +
        coefficients[1] * lambda[1] * (2*lambda[1]-1) +
        coefficients[2] * lambda[2] * (2*lambda[2]-1) +
        coefficients[3] * lambda[3] * (2*lambda[3]-1) +
        coefficients[4] * 4 * lambda[1] * lambda[0] +
        coefficients[5] * 4 * lambda[2] * lambda[0] +
        coefficients[6] * 4 * lambda[3] * lambda[0] +
        coefficients[7] * 4 * lambda[1] * lambda[2] +
        coefficients[8] * 4 * lambda[2] * lambda[3] +
        coefficients[9] * 4 * lambda[3] * lambda[1];
}


// This covers all order 1 cases
template <UInt NNODES, UInt mydim, UInt ndim>
inline Eigen::Matrix<Real,ndim,1> ElementCore<NNODES,mydim,ndim>::evaluate_der_point(const Point<ndim>& point, const Eigen::Matrix<Real,NNODES,1>& coefficients) const
{
  Eigen::Matrix<Real,mydim,mydim+1> B1;
  B1.col(0).setConstant(-1);
  B1.rightCols(mydim).setIdentity();

  return this->M_invJ_.transpose()*B1*coefficients;
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
	return this->M_invJ_.transpose()*B2*coefficients;
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
	return this->M_invJ_.transpose()*B2*coefficients;
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

	return this->M_invJ_.transpose()*B2*coefficients;
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
	static constexpr Real basis_fun[]={2, -1, -1, 1, 4, 4,
																		 -1, 2, -1, 4, 1, 4,
																		 -1, -1, 2, 4, 4, 1
																		};

	return this->getMeasure()/27 * coefficients.replicate<3,1>().dot(Eigen::Map<const Eigen::Matrix<Real,18,1> >(basis_fun));
}

// Full specialization for order 2 in 2.5D
template <>
inline Real ElementCore<6,2,3>::integrate(const Eigen::Matrix<Real,6,1>& coefficients) const
{
	static constexpr Real basis_fun[]={2, -1, -1, 1, 4, 4,
																		 -1, 2, -1, 4, 1, 4,
																		 -1, -1, 2, 4, 4, 1
																		};

	return this->getMeasure()/27 * coefficients.replicate<3,1>().dot(Eigen::Map<const Eigen::Matrix<Real,18,1> >(basis_fun));
}

// Full specialization for order 2 in 3D
// MEMO: this works assuming edges are ordered like so: (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
template <>
inline Real ElementCore<10,3,3>::integrate(const Eigen::Matrix<Real,10,1>& coefficients) const
{
	// In this case a more complicated integration scheme is needed!
	static constexpr Real basis_fun[]={0.1, -0.1, -0.1, -0.1, 0.323606797749979, 0.323606797749979, 0.323606797749979, 0.076393202250021, 0.076393202250021, 0.076393202250021,
																		 -0.1, 0.1, -0.1, -0.1, 0.323606797749979, 0.076393202250021, 0.076393202250021, 0.323606797749979, 0.076393202250021, 0.323606797749979,
																		 -0.1, -0.1, 0.1, -0.1, 0.076393202250021, 0.323606797749979, 0.076393202250021, 0.323606797749979, 0.323606797749979, 0.076393202250021,
																		 -0.1, -0.1, -0.1, 0.1, 0.076393202250021, 0.076393202250021, 0.323606797749979, 0.076393202250021, 0.323606797749979, 0.323606797749979
																		};

	return this->getMeasure() * 0.25*coefficients.replicate<4,1>().dot(Eigen::Map<const Eigen::Matrix<Real,40,1> >(basis_fun));
}





#endif
