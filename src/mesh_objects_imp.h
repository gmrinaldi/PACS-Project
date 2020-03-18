//#include "mesh_objects.hpp"
#ifndef __MESH_OBJECTS_IMP_HPP__
#define __MESH_OBJECTS_IMP_HPP__

// Template member specialization declaration for 2.5D
template <>
void Element<2,3>::computeProperties();

//Primary template member definition
template <UInt mydim, UInt ndim>
void Element<mydim,ndim>::computeProperties()
{
	Element<mydim,ndim> &t = *this;
	std::vector<Real> diff;

	for (int i=0; i<mydim; ++i){
			diff=point_diff(t[i+1],t[0]);
			M_J_.col(i) = Eigen::Map<Eigen::Matrix<Real,ndim,1> >(diff.data());
	}
	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!
	detJ_ = M_J_.determinant();
	M_pseudo_invJ_ = M_J_.inverse();
	metric_ = M_pseudo_invJ_ * M_pseudo_invJ_.transpose();
}

template <UInt mydim, UInt ndim>
Eigen::Matrix<Real,mydim+1,1> Element<mydim,ndim>::getBaryCoordinates(const Point& point) const
{
	Element<mydim,ndim> t = *this;
	Eigen::Matrix<Real,mydim+1,1> lambda;
	Eigen::Matrix<Real,ndim,1> b;

	std::vector<Real> diff = point_diff(point, t[0]);
	b = Eigen::Map<Eigen::Matrix<Real,ndim,1> >(diff.data());
	lambda.tail(mydim) = M_pseudo_invJ_*b;

  lambda(0) = 1 - lambda.tail(mydim).sum();

	return lambda;

}

// Template member specialization declaration for 2.5D
template <>
bool Element<2,3>::isPointInside(const Point& point) const;

//Primary template member definition
template <UInt mydim, UInt ndim>
bool Element<mydim,ndim>::isPointInside(const Point& point) const
{
	Real eps = 	std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;
	Eigen::Matrix<Real,mydim+1,1> lambda = getBaryCoordinates(point);

	return (-tolerance<=lambda.array()).all();

}


template <UInt mydim, UInt ndim>
void Element<mydim,ndim>::print(std::ostream & out) const
{
	out<<"Element -"<< id_ <<"- ";
	for (UInt i=0; i<NNODES; ++i)
		out<<points_[i].getId()<<"  ";
	out<<std::endl;
}

// This works in general! But is it the right thing?
template <UInt mydim, UInt ndim>
int Element<mydim,ndim>::getPointDirection(const Point& point) const
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

#endif
