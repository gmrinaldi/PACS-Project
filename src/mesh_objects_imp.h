//#include "mesh_objects.hpp"
#ifndef __MESH_OBJECTS_IMP_HPP__
#define __MESH_OBJECTS_IMP_HPP__

template <UInt NNODES>
void Element<NNODES,2,2>::computeProperties()
{

	Element<NNODES,2,2> &t = *this;

	std::vector<Real> d1=point_diff(t[1],t[0]);
	std::vector<Real> d2=point_diff(t[2],t[0]);

	// M_J_(0,0) = d1[0];			// (x2-x1)
	// M_J_(1,0) = d1[1];			// (y2-y1)
	// M_J_(0,1) = d2[0];			// (x3-x1)
	// M_J_(1,1) = d2[1];			// (y3-y1)

	M_J_.col(0) = Eigen::Map<Eigen::Matrix<Real,2,1> >(d1.data());
	M_J_.col(1) = Eigen::Map<Eigen::Matrix<Real,2,1> >(d2.data());

	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!

	detJ_ = M_J_.determinant();

	M_invJ_ = M_J_.inverse();

	metric_ = M_invJ_ * M_invJ_.transpose();
}

template <UInt NNODES>
Eigen::Matrix<Real,3,1> Element<NNODES,2,2>::getBaryCoordinates(const Point& point) const{

	Element<NNODES,2,2> t = *this;
	Eigen::Matrix<Real,3,1> lambda;
	Eigen::Matrix<Real,2,1> b;

	//Compute barycentric coordinates for the point
	std::vector<Real> diff = point_diff(point, t[0]);
	b = Eigen::Map<Eigen::Matrix<Real,2,1> >(diff.data());
	lambda.tail<2>() = M_invJ_*b;
  lambda(0) = 1 - lambda.tail<2>().sum();

	return lambda;

}


template <UInt NNODES>
bool Element<NNODES,2,2>::isPointInside(const Point& point) const
{
	Real eps = 	std::numeric_limits<Real>::epsilon,
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);

	return (-tolerance<=lambda.array()).all();

}


// TO BE FIXED: if one dir -1, try with others
template <UInt NNODES>
int Element<NNODES,2,2>::getPointDirection(const Point& point) const
{
	Real eps = std::numeric_limits<Real>::epsilon,
		 tolerance = 10 * eps;

	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);

	//Find the minimum coordinate (if negative stronger straight to the point searched)
	int min_index;
	lambda.minCoeff(&min_index);

	if(lambda[min_index] < -tolerance)
			return min_index;
	else
			return -1;
}


template <UInt NNODES>
void Element<NNODES,2,2>::print(std::ostream & out) const
{
	out<<"Triangle -"<< id_ <<"- ";
	for (UInt i=0; i<NNODES; ++i)
		out<<points_[i].getId()<<"  ";
	out<<std::endl;
}


//IMPLEMENTAZIONE myDim=2, nDim=3

template <UInt NNODES>
void Element<NNODES,2,3>::computeProperties()
{

	Element<NNODES,2,3> &t = *this;

	std::vector<Real> d1=point_diff(t[1],t[0]);
	std::vector<Real> d2=point_diff(t[2],t[0]);

	// M_J_(0,0) = d1[0];			// (x2-x1)
	// M_J_(1,0) = d1[1];			// (y2-y1)
	// M_J_(2,0) = d1[2];			// (z2-z1)
	// M_J_(0,1) = d2[0];			// (x3-x1)
	// M_J_(1,1) = d2[1];			// (y3-y1)
	// M_J_(2,1) = d2[2];			// (z3-z1)

	M_J_.col(0) = Eigen::Map<Eigen::Matrix<Real,3,1> >(d1.data());
	M_J_.col(1) = Eigen::Map<Eigen::Matrix<Real,3,1> >(d2.data());

	G_J_=M_J_.transpose()*M_J_;

	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!

	detJ_ = G_J_.determinant();

	metric_ = G_J.inverse();

}

template <UInt NNODES>
Eigen::Matrix<Real,3,1> Element<NNODES,2,3>::getBaryCoordinates(const Point& point) const{

  Element<NNODES,2,3> t=*this;
  Eigen::Matrix<Real,3,1> lambda;
	Eigen::Matrix<Real,3,1> b;

	std::vector<Real> diff = point_diff(point, t[0]);
	b = Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());

	lambda.tail<2>() = metric_ * M_J_.transpose() * b;

	// lambda.tail<2>() = M_J_.colPivHouseholderQr().solve(b);

	lambda(0) = 1 - lambda.tail<2>().sum();

	return lambda;

}


// THIS COMMENT IS FROM BERAHA, COSMO: We solve 3 scalar equation in 2 unknowns(u,v)
// u*(P1-P0)+v*(P2-P0)=P-P0
// if the system is solveable, P is in the plane (P1,P2,P0), if in addition
// u,v>=0 and u+v<=1 then P is inside the triangle

template <UInt NNODES>
bool Element<NNODES,2,3>::isPointInside(const Point& point) const
{
	Real eps = std::numeric_limits<Real>::epsilon;
	Real tolerance = 10 * eps;

//THIS COMMENT IS FROM BERAHA, COSMO First: check consistency trough Rouchè-Capelli theorem

	Element<NNODES,2,3> t=*this;

	Eigen::Matrix<Real,3,3> A;

	std::vector<Real> diff = point_diff(point, t[0]);

	A = M_J_, Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());

 // NOTE: this method is as fast as ColPivHouseholderQR for such small matrices
 // but this is optimized for rank computations (see eigen documentation)
	if (!A.FullPivHouseholderQR().isInvertible())
			return FALSE;

	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);

	return (-tolerance<=lambda.array()).all();

}


/*
template <UInt NNODES>
int Triangle<NNODES,2,3>::getPointDirection(const Point& point) const
{

	//da implementare
	std::cerr<<"ancora da implementare";
}
*/

template <UInt NNODES>
void Element<NNODES,2,3>::print(std::ostream & out) const
{
	out<<"Triangle -"<< id_ <<"- ";
	for (UInt i=0; i<NNODES; ++i)
		out<<points_[i].getId()<<"  ";
	out<<std::endl;
}


//IMPLEMENTAZIONE myDim=3, nDim=3

template <UInt NNODES>
void Element<NNODES,3,3>::computeProperties()
{

	Element<NNODES,3,3> &t = *this;
	std::vector<Real> d1=point_diff(t[1],t[0]);
	std::vector<Real> d2=point_diff(t[2],t[0]);
	std::vector<Real> d3=point_diff(t[3],t[0]);


	// M_J_(0,0) = d1[0];			// (x2-x1)
	// M_J_(1,0) = d1[1];			// (y2-y1)
	// M_J_(2,0) = d1[2];			// (z2-z1)
	// M_J_(0,1) = d2[0];			// (x3-x1)
	// M_J_(1,1) = d2[1];			// (y3-y1)
	// M_J_(2,1) = d2[2];			// (z3-z1)
	// M_J_(0,2) = d3[0];			// (x4-x1)
	// M_J_(1,2) = d3[1];			// (y4-y1)
	// M_J_(2,2) = d3[2];			// (z4-z1)

	M_J_.col(0) = Eigen::Map<Eigen::Matrix<Real,3,1> >(d1.data());
	M_J_.col(1) = Eigen::Map<Eigen::Matrix<Real,3,1> >(d2.data());
	M_J_.col(2) = Eigen::Map<Eigen::Matrix<Real,3,1> >(d3.data());

	// NOTE: for small (not bigger than 4x4) matrices eigen directly calculates
	// determinants and inverses, it is very efficient!

	Real detMJ_ = M_J_.determinant();

	M_invJ_=M_J_.inverse();

	G_J_=M_J_.transpose()*M_J_;

	detJ_ = G_J_.determinant();

	metric_ = G_J_.inverse();

}


template <UInt NNODES>
Eigen::Matrix<Real,4,1> Element<NNODES,3,3>::getBaryCoordinates(const Point& point) const{


	Element<NNODES,3,3> t=*this;
	Eigen::Matrix<Real,4,1> lambda;
	Eigen::Matrix<Real,3,1> b;

	std::vector<Real> diff = point_diff(point, t[0]);
	b = Eigen::Map<Eigen::Matrix<Real,3,1> >(diff.data());

	lambda.tail<3>() = M_invJ_ * b;

	lambda[0]=1-lambda.tail<3>.sum();

	return lambda;

}


template <UInt NNODES>
bool Element<NNODES,3,3>::isPointInside(const Point& point) const
{
	Real eps = std::numeric_limits<Real>::epsilon;
	Real tolerance = 10 * eps;

	Eigen::Matrix<Real,4,1> lambda = getBaryCoordinates(point);

	return (-tolerance<=lambda.array()).all();
}


template <UInt NNODES>
void Element<NNODES,3,3>::print(std::ostream & out) const
{
	out<<"Tetrahedron -"<< id_ <<"- ";
	for (UInt i=0; i<NNODES; ++i)
		out<<points_[i].getId()<<"  ";
	out<<std::endl;
}




#endif
