/*
 * mesh_object.cpp
 *
 *  Created on: Aug 16, 2015
 *      Author: eardi
 */

#include "mesh_objects.h"

void Point::print(std::ostream & out) const
{
	out<<"Point -"<< id_ <<"- "<<"("<<coord_[0]<<","<<coord_[1]<<","<<coord_[2]<<")"<<std::endl<<"------"<<std::endl;
}

void Edge::print(std::ostream & out) const
{
	out<<"Edge -"<< id_ <<"- "<<"("<<points_[0].getId()<<","<<points_[1].getId()<<")"<<std::endl;
}

std::vector<Real> point_diff(const Point &lhs, const Point &rhs){
		std::vector<Real> diff;
		for (int i=0; i<lhs.ndim; ++i)
				diff.push_back(lhs[i]-rhs[i]);
		return diff;
};

// Template member specialization for 2.5D
template <>
void Element<2,3>::computeProperties()
{
	Element<2,3> &t = *this;
	std::vector<Real> diff;

	for (int i=0; i<2; ++i){
			diff=point_diff(t[i+1],t[0]);
			M_J_.col(i) = Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());
	}

	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!

	Eigen::Matrix<Real,2,2> G_J_=M_J_.transpose()*M_J_;
	metric_ = G_J_.inverse();
	M_pseudo_invJ_ = metric_ * M_J_.transpose();
}


// Template member specialization for 2.5D
template <>
bool Element<2,3>::isPointInside(const Point& point) const
{
	Real eps = 	std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;

	Element<2,3> t=*this;

 	Eigen::Matrix<Real,3,3> A;
	Eigen::Matrix<Real,3,1> lambda;

 	std::vector<Real> diff = point_diff(point, t[0]);

 	A << M_J_, Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());

 	// NOTE: this method is as fast as ColPivHouseholderQR for such small matrices
 	// but this is optimized for rank computations (see eigen documentation)
 	// Check if point belongs to plane!
 	if(!A.fullPivHouseholderQr().isInvertible())
			return false;

	lambda.tail<2>() = M_pseudo_invJ_ * Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());

	lambda(0) = 1 - lambda.tail<2>().sum();

	return (-tolerance<=lambda.array()).all();

}
