//#include "mesh_objects.hpp"
#ifndef __MESH_OBJECTS_IMP_HPP__
#define __MESH_OBJECTS_IMP_HPP__

//Primary template member definition
template <UInt NNODES, UInt mydim, UInt ndim>
void Element<NNODES,mydim,ndim>::computeProperties()
{
	Element<NNODES,mydim,ndim> &t = *this;
	std::array<Real,ndim> diff;

	for (int i=0; i<mydim; ++i){
			diff=point_diff<ndim>(t[i+1],t[0]);
			M_J_.col(i) = Eigen::Map<Eigen::Matrix<Real,ndim,1> >(diff.data());
	}
	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!
	detJ_ = M_J_.determinant();
	M_invJ_ = M_J_.inverse();
	metric_ = M_invJ_ * M_invJ_.transpose();
}

template <UInt NNODES, UInt mydim, UInt ndim>
Eigen::Matrix<Real,mydim+1,1> Element<NNODES,mydim,ndim>::getBaryCoordinates(const Point<ndim>& point) const
{
	Element<NNODES,mydim,ndim> t = *this;
	Eigen::Matrix<Real,mydim+1,1> lambda;

	std::array<Real,ndim> diff = point_diff<ndim>(point, t[0]);

	lambda.tail(mydim) = M_invJ_ * Eigen::Map<Eigen::Matrix<Real,ndim,1> >(diff.data());

  lambda(0) = 1 - lambda.tail(mydim).sum();

	return lambda;

}

//Primary template member definition
template <UInt NNODES, UInt mydim, UInt ndim>
bool Element<NNODES,mydim,ndim>::isPointInside(const Point<ndim>& point) const
{
	Real eps = 	std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;
	Eigen::Matrix<Real,mydim+1,1> lambda = getBaryCoordinates(point);

	return (-tolerance<=lambda.array()).all();

}


template <UInt NNODES, UInt mydim, UInt ndim>
void Element<NNODES,mydim,ndim>::print(std::ostream & out) const
{
	out<<"Element -"<< id_ <<"- ";
	for (UInt i=0; i<NNODES; ++i)
		out<<points_[i].getId()<<"  ";
	out<<std::endl;
}

template <UInt NNODES, UInt mydim, UInt ndim>
int Element<NNODES,mydim,ndim>::getPointDirection(const Point<ndim>& point) const
{
	Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,mydim+1,1> lambda = getBaryCoordinates(point);

	//Find the minimum coordinate (if negative stronger straight to the point searched)
	int min_index;
	lambda.minCoeff(&min_index);

	if(lambda[min_index] < -tolerance)
			return min_index;
	else
			return -1;
}

// Template member specializations for 2.5D

template <UInt NNODES>
void Element<NNODES,2,3>::computeProperties()
{
	Element<NNODES,2,3> &t = *this;
	std::array<Real,3> diff;

	for (int i=0; i<2; ++i){
			diff=point_diff<3>(t[i+1],t[0]);
			M_J_.col(i) = Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());
	}

	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!

	Eigen::Matrix<Real,2,2> G_J_=M_J_.transpose()*M_J_;
	detJ_ = G_J_.determinant();
	metric_ = G_J_.inverse();
	M_invJ_ = metric_ * M_J_.transpose();
}

template <UInt NNODES>
Eigen::Matrix<Real,3,1> Element<NNODES,2,3>::getBaryCoordinates(const Point<3>& point) const
{
	Element<NNODES,2,3> t = *this;
	Eigen::Matrix<Real,3,1> lambda;

	std::array<Real,3> diff = point_diff<3>(point, t[0]);

	lambda.tail<2>() = M_invJ_ * Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());

  lambda(0) = 1 - lambda.tail<2>().sum();

	return lambda;

}

template <UInt NNODES>
bool Element<NNODES,2,3>::isPointInside(const Point<3>& point) const
{
	Real eps = 	std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Element<NNODES,2,3> t=*this;

 	Eigen::Matrix<Real,3,3> A;

 	std::array<Real,3> diff = point_diff<3>(point, t[0]);

 	A << M_J_, Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());

 	// NOTE: this method is as fast as ColPivHouseholderQR for such small matrices
 	// but this is optimized for rank computations (see eigen documentation)
 	// Check if point belongs to plane!
 	if(A.fullPivHouseholderQr().isInvertible())
			return false;

	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);

	return (-tolerance<=lambda.array()).all();

}

template <UInt NNODES>
void Element<NNODES,2,3>::print(std::ostream & out) const
{
	out<<"Element -"<< id_ <<"- ";
	for (UInt i=0; i<NNODES; ++i)
		out<<points_[i].getId()<<"  ";
	out<<std::endl;
}


#endif
