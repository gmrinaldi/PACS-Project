#ifndef __MESH_OBJECTS_HPP__
#define __MESH_OBJECTS_HPP__

#include "fdaPDE.h"

typedef UInt Id;
typedef UInt BcId;

//Simple function to evaluate factorial at compile time
//Needed for getVolume member of Element
constexpr int factorial(int n) {
    return n ? (n * factorial(n - 1)) : 1;
}


// //Define the NotValid meaning value
// constexpr UInt Identifier::NVAL=std::numeric_limits<UInt>::max();


//!  This class gives some common methods to all mesh objects.
class Identifier{
public:

	//! An static const Unisgned Integer.
    /*! Needed to identify the Not Valid Id. */
	static constexpr UInt NVAL=std::numeric_limits<UInt>::max();

	Identifier(UInt id):id_(id),bcId_(NVAL){}
	Identifier(UInt id, UInt bcId):id_(id),bcId_(bcId){}

	bool unassignedId()const {return id_==NVAL;}
	bool unassignedBc()const {return bcId_==NVAL;}

	Id id() const {return id_;}
	BcId bcId() const {return bcId_;}
	Id getId() const {return id_;}


	protected:
	Id id_;
	BcId bcId_;
};


//!  This class implements a 3D point, the default is z=0 => 2D point
class Point : public Identifier{
public:

	friend std::vector<Real> point_diff(const Point &, const Point &);

  Point(Real x, Real y) :
					Identifier(NVAL, NVAL), coord_({x,y}), ndim(2) {}
	Point(Real x, Real y, Real z) :
					Identifier(NVAL, NVAL), coord_({x,y,z}), ndim(3) {}
	Point(Id id, BcId bcId, Real x, Real y) :
					Identifier(id, bcId), coord_({x,y}), ndim(2) {}
	Point(Id id, BcId bcId, Real x, Real y, Real z) :
					Identifier(id, bcId), coord_({x,y,z}), ndim(3) {}

	void print(std::ostream & out) const;
	Real operator[](UInt i) const {return coord_[i];}
private:
	const std::vector<Real> coord_;
	const UInt ndim;
};


//!  This class implements an Edge, as an objects composed by two 2D points.
class Edge : public Identifier{
  public:
    static const UInt NNODES=2;
    static const UInt numSides=1;
    static const UInt myDim=1;

    Edge(Id id, BcId bcId, const Point& start, const Point& end) :
					Identifier(id, bcId), points_({start,end}) {}

    void print(std::ostream & out) const;

    Point operator[](UInt i) const {return points_[i];}

 private:
    const std::vector<Point> points_;
  };

//! This is an abstract template class called Element
/*!
 * mydim is the dimension of the object: e.g. a triangle has mydim=2, a tethraedron
 *       has mydim = 3
 *
 * ndim is the dimension of the space in which the object is embedded
 *
*/


template <UInt mydim, UInt ndim>
class Element : public Identifier {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	static const UInt numVertices=ndim+1;
	static const UInt numSides=ndim*(ndim+1)/2;
	static const UInt myDim=mydim;

	//! This constructor creates an Element, given its Id and an std array with the three object Point the will define the Element
		Element(Id id, const std::vector<Point>& points) :
					Identifier(id), points_(points), NNODES(points_.size()) {this->computeProperties();}

	//! Overloading of the operator [],  taking the Node number and returning a node as Point object.
		/*!
		 * For node numbering convention see:
			\param i an integer argument.
			\return the Point object
		*/
	Point operator[](UInt i) const {return points_[i];}

	//! A member that computes the barycentric coordinates.
		/*!
			\param point a Point object
			\return The three baricentric coordinates of the point
		*/

	Real getDetJ() const {return detJ_;}
	const Eigen::Matrix<Real,ndim,mydim>& getM_J() const {return M_J_;}
	const Eigen::Matrix<Real,ndim,mydim>& getM_pseudo_invJ() const {return M_pseudo_invJ_;}
	const Eigen::Matrix<Real,ndim,mydim>& getMetric() const {return metric_;}

	//! A member returning the area of the finite element
			/*!
				\return a Real value representing the area of the triangle from which we updated the element
				\sa  updateElement(Element<Integrator::NNODES> t)
			*/
	Real getVolume() const {return std::abs(detJ_)/factorial(ndim);}

	Eigen::Matrix<Real,mydim+1,1> getBaryCoordinates(const Point& point) const;

	//! A member that tests if a Point is located inside an Element.
		/*!
			\param point a Point object.
			\return True if the point is inside the triangle
		*/
	bool isPointInside(const Point& point) const;

	//! A memeber that verifies which edge separates the Triangle from a Point.
		/*!
			\param point a Point object.
			\return The number of the Edge that separates the point
			from the triangle and -1 if the point is inside the triangle.
		*/
	int getPointDirection(const Point& point) const;

	//! A member that prints the main properties of the triangle
		/*!
			\param out a std::outstream.
		*/
	void print(std::ostream & out) const;

private:
	// IMPORTANT DO NOT CHANGE THE ORDER OF points_ and NNODES!
	const std::vector<Point> points_;
	const UInt NNODES;
	Eigen::Matrix<Real,ndim,mydim> M_J_;
	Eigen::Matrix<Real,mydim,ndim> M_pseudo_invJ_; //this is invM_J_ in 2D and 3D
	Eigen::Matrix<Real,mydim,mydim> metric_;
	Real detJ_;
	void computeProperties();
};

//!  This class implements a Triangle as an objects composed by three or six nodes.
/*!
 *  The first three nodes represent the vertices, the others the internal nodes,
 *  following this enumeration: !IMPORTANT! different from Sangalli code!
 *
 * 		        3
 * 			    *
 * 		     /    \
 * 		  5 *	   * 4
 * 		  /	        \
 * 		 *_____*_____*
 * 		1	   6	  2
*/

//!  This class implements a Triangle as an objects composed by three or six nodes, embedded in a 3-dimensional space
/*!
 *  The first three nodes represent the vertices, the others the internal nodes,
 *  following this enumeration: !IMPORTANT! different from Sangalli code!
 *
 * 			    3
 * 			    *
 * 		     /    \
 * 		  5 *	   * 4
 * 		  /	        \
 * 		 *_____*_____*
 * 		1	   6	  2
*/



//!  This class implements a Tetrahedron as an objects composed by four or ten nodes, embedded in a 3-dimensional space. Currently, only the 4 nodes version is implemented. The tetrahedron is an Element with mydim=3 and ndim=3

//
// template <UInt Nodes, UInt mydim, UInt ndim>
// inline Real evaluate_point(const Element<Nodes,mydim,ndim>& t, const Point& point, const Eigen::Matrix<Real,Nodes,1>& coefficients)
// {
// 	//std::cerr<< "TRYING TO EVALUATE ORDER NOT IMPLEMENTED" << std::endl;
// 	return 0;
// }
//
// template <>
// inline Real evaluate_point<3,2,2>(const Element<3,2,2>& t, const Point& point, const Eigen::Matrix<Real,3,1>& coefficients)
// {
// 	Eigen::Matrix<Real,3,1> bary_coeff = t.getBaryCoordinates(point);
// 	//std::cout<< "B-coord: "<<bary_coeff<<std::endl;
// 	return(coefficients.dot(bary_coeff));
// }
//
// template <>
// inline Real evaluate_point<6,2,2>(const Element<6,2,2>& t, const Point& point, const Eigen::Matrix<Real,6,1>& coefficients)
// {
// 	Eigen::Matrix<Real,3,1> bary_coeff = t.getBaryCoordinates(point);
// 	return( coefficients[0]*(2*bary_coeff[0]*bary_coeff[0] - bary_coeff[0]) +
//             coefficients[1]*(2*bary_coeff[1]*bary_coeff[1] - bary_coeff[1]) +
//             coefficients[2]*(2*bary_coeff[2]*bary_coeff[2] - bary_coeff[2]) +
//             coefficients[3]*(4*bary_coeff[1]* bary_coeff[2]) +
//             coefficients[4]*(4*bary_coeff[2]* bary_coeff[0]) +
//             coefficients[5]*(4*bary_coeff[0]* bary_coeff[1]) );
// }
//
//
//
// /*! THIS COMMENT COMES FROM BERAHA, COSMO: in this case, the implementation is not as trivial
//  first solve the linear system (p-p0) = (p1-p0)*alpha + (p2-p0)*beta + N*gamma
//  where p0,p1,p2 are the vertices of the triangle, p is the point
//  (observe that, if the point is inside the triangle, gamma=0)
//  then the solution u(p) = u(p0) + alpa*(u(p1) - u(p0) + beta*(u(p2)-u(p0))
//  */
// template <>
// inline Real evaluate_point<3,2,3>(const Element<3,2,3>& t, const Point& point, const Eigen::Matrix<Real,3,1>& coefficients)
// {
// 	Eigen::Matrix<Real,3,1> bary_coeff = t.getBaryCoordinates(point);
// 	return(coefficients.dot(bary_coeff));
// }
//
// template <>
// inline Real evaluate_point<6,2,3>(const Element<6,2,3>& t, const Point& point, const Eigen::Matrix<Real,6,1>& coefficients)
// {
// 	Eigen::Matrix<Real,3,1> bary_coeff = t.getBaryCoordinates(point);
// 	return( coefficients[0]*(2*bary_coeff[0]*bary_coeff[0] - bary_coeff[0]) +
//             coefficients[1]*(2*bary_coeff[1]*bary_coeff[1] - bary_coeff[1]) +
//             coefficients[2]*(2*bary_coeff[2]*bary_coeff[2] - bary_coeff[2]) +
//             coefficients[3]*(4*bary_coeff[1]*bary_coeff[2]) +
//             coefficients[4]*(4*bary_coeff[2]*bary_coeff[0]) +
//             coefficients[5]*(4*bary_coeff[0]*bary_coeff[1]) );
// }
//
// //! Implementation for tetrahedrons
// template <>
// inline Real evaluate_point<4,3,3>(const Element<4,3,3>& t, const Point& point, const Eigen::Matrix<Real,4,1>& coefficients)
// {
// 	Eigen::Matrix<Real,4,1> bary_coeff=t.getBaryCoordinates(point);
// 	return(coefficients.dot(bary_coeff));
// }
//
// template <UInt Nodes,UInt mydim, UInt ndim>
// inline Eigen::Matrix<Real,ndim,1> evaluate_der_point(const Element<Nodes,mydim,ndim>& t, const Point& point, const Eigen::Matrix<Real,Nodes,1>& coefficients)
// {
// 	//std::cerr<< "TRYING TO EVALUATE ORDER NOT IMPLEMENTED" << std::endl;
// 	Eigen::Matrix<Real,ndim,1> null;
// 	return(null);
// }
//
// template <>
// inline Eigen::Matrix<Real,2,1> evaluate_der_point<3,2,2>(const Element<3,2,2>& t, const Point& point, const Eigen::Matrix<Real,3,1>& coefficients)
// {
// 	Eigen::Matrix<Real,2,3> B1;
// 	B1 << t[1][1] - t[2][1], t[2][1] - t[0][1], t[0][1] - t[1][1],
// 		t[2][0] - t[1][0], t[0][0] - t[2][0], t[1][0] - t[0][0];
//
// 	B1 = B1 / (2 * t.getVolume());
//
// 	return(B1*coefficients);
//
// }
//
// template <>
// inline Eigen::Matrix<Real,2,1> evaluate_der_point<6,2,2>(const Element<6,2,2>& t, const Point& point, const Eigen::Matrix<Real,6,1>& coefficients)
// {
// 	Eigen::Matrix<Real,3,1> L = t.getBaryCoordinates(point);
// 	Eigen::Matrix<Real,2,3> B1;
// 	B1 << t[1][1] - t[2][1], t[2][1] - t[0][1], t[0][1] - t[1][1],
// 		t[2][0] - t[1][0], t[0][0] - t[2][0], t[1][0] - t[0][0];
// 	B1 = B1 / (2 * t.getVolume());
// 	Eigen::Matrix<Real,3,6> B2;
// 	B2 << 4*L[0]-1, 0       , 0       , 0        , 4*L[2], 4*L[1],
// 		  0       , 4*L[1]-1, 0       , 4*L[2]   , 0     , 4*L[0],
// 		  0       , 0       , 4*L[2]-1, 4*L[1]   , 4*L[0], 0     ;
// 	return(B1*B2*coefficients);
// }
//
//
//
// template <>
// inline Eigen::Matrix<Real,3,1> evaluate_der_point<4,3,3>(const Element<4,3,3>& t, const Point& point, const Eigen::Matrix<Real,4,1>& coefficients)
// {
//
// 	Eigen::Matrix<Real,3,4> B1;
// 	B1 << -1,1,0,0,
// 	      -1,0,1,0,
// 	      -1,0,0,1;
// 	B1 = B1 / (6 * t.getVolume());
// 	/*Eigen::Matrix<Real,3,3> B1;
// 	B1(0,0)=-t.getM_J()(1,2)*t.getM_J()(2,1) + t.getM_J()(1,1)*t.getM_J()(2,2);
// 	B1(0,1)= t.getM_J()(0,2)*t.getM_J()(2,1) - t.getM_J()(0,1)*t.getM_J()(2,2);
// 	B1(0,2)=-t.getM_J()(0,2)*t.getM_J()(1,1) + t.getM_J()(0,1)*t.getM_J()(1,2);
// 	B1(1,0)= t.getM_J()(1,2)*t.getM_J()(2,0) - t.getM_J()(1,0)*t.getM_J()(2,2);
// 	B1(1,1)=-t.getM_J()(0,2)*t.getM_J()(2,0) + t.getM_J()(0,0)*t.getM_J()(2,2);
// 	B1(1,2)= t.getM_J()(0,2)*t.getM_J()(1,0) - t.getM_J()(0,0)*t.getM_J()(1,2);
// 	B1(2,0)=-t.getM_J()(1,1)*t.getM_J()(2,0) + t.getM_J()(1,0)*t.getM_J()(2,1);
// 	B1(2,1)= t.getM_J()(0,1)*t.getM_J()(2,0) - t.getM_J()(0,0)*t.getM_J()(1,2);
// 	B1(2,2)=-t.getM_J()(0,1)*t.getM_J()(1,0) + t.getM_J()(0,0)*t.getM_J()(1,1);
//
// 	B1 = B1 / (6*std::sqrt(t.getDetJ()));
// 	*/
//
// 	return(B1*coefficients);
//
// }



#include "mesh_objects_imp.h"
#endif
