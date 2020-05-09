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

// From here on implementation of Element
template <UInt NNODES, UInt mydim, UInt ndim>
Eigen::Matrix<Real,mydim+1,1> Element<NNODES,mydim,ndim>::getBaryCoordinates(const Point<ndim> &point) const
{
	Eigen::Matrix<Real,mydim+1,1> lambda;

	// .template is needed! See Eigen documentation regarding
	// the template and typename keywords in C++
	lambda.template tail<mydim>().noalias() = M_invJ_ * (EigenMap2Const_t(&point[0])-EigenMap2Const_t(&points_[0][0]));

  lambda(0) = 1 - lambda.template tail<mydim>().sum();

	return lambda;

}

template <UInt NNODES, UInt mydim, UInt ndim>
bool Element<NNODES,mydim,ndim>::isPointInsideImpl(const Point<ndim>& point, std::false_type) const
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,mydim+1,1> lambda = this->getBaryCoordinates(point);

	return (-tolerance<=lambda.array()).all();

}

// Implementation for manifold data
template <UInt NNODES, UInt mydim, UInt ndim>
bool Element<NNODES,mydim,ndim>::isPointInsideImpl(const Point<ndim>& point, std::true_type) const
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

 	Eigen::Matrix<Real,ndim,ndim> A;
 	A << M_J_, (EigenMap2Const_t(&point[0])-EigenMap2Const_t(&points_[0][0]));

	Eigen::Matrix<Real,mydim+1,1> lambda = this->getBaryCoordinates(point);

 	// NOTE: this method is as fast as ColPivHouseholderQR for such small matrices
 	// but this is optimized for rank computations (see eigen documentation)
	return !A.fullPivHouseholderQr().isInvertible() && (-tolerance<=lambda.array()).all();

}

// This function is called to construct elements in 2D and 3D
template <UInt NNODES, UInt mydim, UInt ndim>
void Element<NNODES,mydim,ndim>::computeProperties(std::false_type)
{
	{
		EigenMap2Const_t basePoint(&points_[0][0]);
		for (int i=0; i<mydim; ++i)
			M_J_.col(i) = EigenMap2Const_t(&points_[i+1][0])-basePoint;
	}
	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!
	M_invJ_ = M_J_.inverse();
	element_measure = std::abs(M_J_.determinant())/factorial(ndim);
}

// implementation for manifold data
template <UInt NNODES, UInt mydim, UInt ndim>
void Element<NNODES,mydim,ndim>::computeProperties(std::true_type)
{
	std::cout<<"ciao"<<std::endl;
	{
		EigenMap2Const_t basePoint(&points_[0][0]);
		for (int i=0; i<2; ++i)
			M_J_.col(i) = EigenMap2Const_t(&points_[i+1][0])-basePoint;
	}

	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!
	M_invJ_.noalias() = (M_J_.transpose()*M_J_).inverse() * M_J_.transpose();
	// Area of 3D triangle is half the norm of cross product of two sides!
	element_measure = M_J_.col(0).cross(M_J_.col(1)).norm()/2;
}


template <UInt NNODES, UInt mydim, UInt ndim>
template <UInt m, UInt n>
typename std::enable_if<n==m && n==ndim && m==mydim, int>::type Element<NNODES,mydim,ndim>::getPointDirection(const Point<ndim>& point) const
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,mydim+1,1> lambda = this->getBaryCoordinates(point);

	//Find the minimum coordinate (if negative stronger straight to the point searched)
	UInt min_index;
	lambda.minCoeff(&min_index);

	return (lambda[min_index] < -tolerance)	? min_index : -1;
}

template <UInt NNODES, UInt mydim, UInt ndim>
template <UInt m, UInt n>
typename std::enable_if<(m<n) && n==ndim && m==mydim, Point<3> >::type Element<NNODES,mydim,ndim>::computeProjection(const Point<3>& point) const
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
		return points_[0];
	else if (lambda(0)<0 && lambda(1)>0 && lambda(2)<0)
		return points_[1];
	else if (lambda(0)<0 && lambda(1)<0 && lambda(2)>0)
		return points_[2];

	Eigen::Matrix<Real,3,1> coords3D;
	// If (+,+,+) the projection lies inside the element
	// So just convert back to 3D coords
	if(lambda(0)>0 && lambda(1)>0 && lambda(2)>0)
		coords3D = M_J_ * lambda.tail<2>();
	// If (+,+,-) the projection lies beyond edge 3
	// Simply scale it back on the edge and convert
	else if(lambda(0)>0 && lambda(1)>0)
		coords3D = M_J_.col(0)/lambda.head<2>().sum();
	// If (+,-,+) the projection lies beyond edge 2
	else if (lambda(0)>0 && lambda(2)>0)
		coords3D = M_J_.col(1)/(lambda(0)+lambda(2));
	// If (-,+,+) the projection lies beyond edge 1
	else
		coords3D = M_J_ * lambda.tail<2>()/lambda.tail<2>().sum();

	// Translate back by the first point and return
	return Point<3>(coords3D)+=points_[0];
}

// This covers all order 1 cases
template <UInt NNODES, UInt mydim, UInt ndim>
inline Real Element<NNODES,mydim,ndim>::evaluate_point(const Point<ndim>& point, const Eigen::Matrix<Real,NNODES,1>& coefficients) const
{
  return coefficients.dot(this->getBaryCoordinates(point));
}

// Full specialization for order 2 in 2D
template <>
inline Real Element<6,2,2>::evaluate_point(const Point<2>& point, const Eigen::Matrix<Real,6,1>& coefficients) const
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
inline Real Element<6,2,3>::evaluate_point(const Point<3>& point, const Eigen::Matrix<Real,6,1>& coefficients) const
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
inline Real Element<10,3,3>::evaluate_point(const Point<3>& point, const Eigen::Matrix<Real,10,1>& coefficients) const
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
inline Real Element<NNODES,mydim,ndim>::integrate(const Eigen::Matrix<Real,NNODES,1>& coefficients) const
{
	// This works because of linearity!
	return this->getMeasure() * coefficients.mean();
}

// Full specialization for order 2 in 2D
template <>
inline Real Element<6,2,2>::integrate(const Eigen::Matrix<Real,6,1>& coefficients) const
{
	static constexpr Real basis_fun[]={2, -1, -1, 1, 4, 4,
																		 -1, 2, -1, 4, 1, 4,
																		 -1, -1, 2, 4, 4, 1
																		};

	return this->getMeasure()/27 * coefficients.replicate<3,1>().dot(Eigen::Map<const Eigen::Matrix<Real,18,1> >(basis_fun));
}

// Full specialization for order 2 in 2.5D
template <>
inline Real Element<6,2,3>::integrate(const Eigen::Matrix<Real,6,1>& coefficients) const
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
inline Real Element<10,3,3>::integrate(const Eigen::Matrix<Real,10,1>& coefficients) const
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
