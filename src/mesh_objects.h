#ifndef __MESH_OBJECTS_HPP__
#define __MESH_OBJECTS_HPP__


#include "fdaPDE.h"

typedef UInt Id;
typedef UInt BcId;

//Simple function to evaluate factorial at compile time
//Needed for getVolume/getArea member of Element
constexpr UInt factorial(UInt n) {
    return n ? (n * factorial(n - 1)) : 1;
}


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

// This class implements a n-dimensional point
// Doesn't allow creation of points in dimension other than 2 or 3
template<UInt ndim>
class Point : public Identifier{
  static_assert(ndim==2 || ndim==3,
    "ERROR: Trying to create a Point object in dimension different than 2D or 3D; see mesh_objects.h");

  public:
    using pointCoords=std::array<Real,ndim>;

    Point(): Identifier(NVAL, NVAL) {};
    Point(pointCoords coord) :
              Identifier(NVAL, NVAL), coord_(coord) {}
    Point(Id id, BcId bcId, pointCoords coord) :
              Identifier(id, bcId), coord_(coord) {}

    Real operator[](UInt i) const {return coord_[i];}

    // Overload the "-" operator to take 2 points of the same dimension and compute
    // the coordinate difference.
    // It returns an array for convenience, keep it in mind!
    friend pointCoords operator -(const Point& lhs, const Point& rhs){
      pointCoords diff(lhs.coord_);
      for (int i=0; i<ndim; ++i)
          diff[i]-=rhs[i];
      return diff;
    };

    friend Real distance(const Point& lhs, const Point& rhs){
      pointCoords diff=lhs-rhs;
      Real dist=0;
      for (auto const &i : diff)
        dist+=diff*diff;
      return std::sqrt(dist);
    };

  private:
    pointCoords coord_;
  };

//!  This class implements an Edge, as an objects composed by two points.
template <UInt ndim>
class Edge : public Identifier{
  public:
    static constexpr UInt NNODES=2;
    static constexpr UInt numSides=1;
    static constexpr UInt myDim=1;

    Edge(Id id, BcId bcId, const Point<ndim>& start, const Point<ndim>& end) :
					Identifier(id, bcId), points_({start,end}) {}

    void print(std::ostream & out) const;

    Point<ndim> operator[](UInt i) const {return points_[i];}

 private:
    const std::array<Point<ndim>,2> points_;
  };


//! This is an abstract template class called Element
/*!
 * mydim is the dimension of the object: e.g. a triangle has mydim=2, a tethraedron
 *       has mydim = 3
 *
 * ndim is the dimension of the space in which the object is embedded
 *
*/

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

//!  This class implements a Tetrahedron as an objects composed by four or ten nodes, embedded in a 3-dimensional space.
// For NNODES=10 the edges are ordered like so: (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
// The midpoints are also expected to follow this convention!

template <UInt NNODES, UInt mydim, UInt ndim>
class ElementCore : public Identifier {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	static constexpr UInt numVertices=mydim+1;
	static constexpr UInt numSides=mydim*(mydim+1)/2;
	static constexpr UInt myDim=mydim;

  using elementPoints=std::array<Point<ndim>,NNODES>;

  //! This constructor creates an "empty" Element, with an Id Not Valid
  ElementCore() :
          Identifier(NVAL) {}

	//! This constructor creates an Element, given its Id and an std array with the three object Point the will define the Element
	ElementCore(Id id, const elementPoints& points) :
					Identifier(id), points_(points) {}

	//! Overloading of the operator [],  taking the Node number and returning a node as Point object.
	Point<ndim> operator[](UInt i) const {return points_[i];}

	Real getDetJ() const {return detJ_;}
	const Eigen::Matrix<Real,ndim,mydim>& getM_J() const {return M_J_;}
	const Eigen::Matrix<Real,mydim,ndim>& getM_invJ() const {return M_invJ_;} //in 2.5D this is actually the pseudoinverse!
	const Eigen::Matrix<Real,mydim,mydim>& getMetric() const {return metric_;}

  //! A member that computes the barycentric coordinates.
	Eigen::Matrix<Real,mydim+1,1> getBaryCoordinates(const Point<ndim>& point) const;

	//! A member that prints the main properties of the triangle
	void print(std::ostream & out) const;

protected:
	elementPoints points_;
	Eigen::Matrix<Real,ndim,mydim> M_J_;
	Eigen::Matrix<Real,mydim,ndim> M_invJ_; //in 2.5D this is actually the pseudoinverse!
	Eigen::Matrix<Real,mydim,mydim> metric_;
	Real detJ_;
};



// This class adds some useful methods for 2D and 3D elements
template <UInt NNODES, UInt mydim, UInt ndim>
class Element : public ElementCore<NNODES,mydim,ndim> {
public:
  //! This constructor creates an "empty" Element, with an Id Not Valid
  Element() :
        ElementCore<NNODES,mydim,ndim>() {}

  //! This constructor creates an Element, given its Id and an std array with the three object Point the will define the Element
  Element(Id id, const std::array<Point<ndim>,NNODES>& points) :
        ElementCore<NNODES,mydim,ndim>(id,points) {this->computeProperties();}

  //! A member returning the area/volume of the element
  // Note: both are needed for compatibility issues with previous implementation
  Real getArea() const {return std::abs(detJ_)/factorial(ndim);}
  Real getVolume() const {return std::abs(detJ_)/factorial(ndim);}

  //! A member that tests if a Point is located inside an Element.
  bool isPointInside(const Point<ndim>& point) const;

  //! A member that verifies which edge/face separates the Triangle/Tetrahedron from a Point.
  int getPointDirection(const Point<ndim>& point) const;

private:
  void computeProperties();

};


// This partial specialization deals with the 2.5D case
template <UInt NNODES>
class Element<NNODES, 2,3> : public ElementCore<NNODES,2,3> {
public:
  //! This constructor creates an "empty" Element, with an Id Not Valid
  Element() :
        ElementCore<NNODES,2,3>() {}

  //! This constructor creates an Element, given its Id and an std array with the three object Point the will define the Element
  Element(Id id, const std::array<Point<3>,NNODES>& points) :
        ElementCore<NNODES,2,3>(id,points) {this->computeProperties();}

  //! A member returning the area of the element
  // Mind the sqrt!
  Real getArea() const {return std::sqrt(detJ_)/2;}

  //! A member that tests if a Point is located inside an Element.
  // Note: this is implemented (slightly) differently in this case
  // In particular one has to check that the point lies on the same plane as the triangle!
  bool isPointInside(const Point<3>& point) const;

  // This function projects a 3D point (XYZ coordinates!) onto the element
  // Note: if the projection lies outside the element the function returns
  // the closest point on the boundary of the element instead
  Point<3> computeProjection(const Point<3>& point) const;

private:
  void computeProperties();

};


// This covers all the order 1 cases!
template <UInt Nodes, UInt mydim, UInt ndim>
inline Real evaluate_point(const Element<Nodes,mydim,ndim>& t, const Point<ndim>& point, const Eigen::Matrix<Real,Nodes,1>& coefficients)
{
  return(coefficients.dot(t.getBaryCoordinates(point)));
	return 0;
}

// Full specialization for order 2 in 2D
template <>
inline Real evaluate_point<6,2,2>(const Element<6,2,2>& t, const Point<2>& point, const Eigen::Matrix<Real,6,1>& coefficients)
{
	Eigen::Matrix<Real,3,1> lambda = t.getBaryCoordinates(point);
  return( coefficients[0]*lambda[0]*(2*lambda[0] - 1) +
            coefficients[1]*lambda[1]*(2*lambda[1] - 1) +
            coefficients[2]*lambda[2]*(2*lambda[2] - 1) +
            coefficients[3]*4*lambda[1]*lambda[2] +
            coefficients[4]*4*lambda[2]*lambda[0] +
            coefficients[5]*4*lambda[0]*lambda[1] );
}



 // Full specialization for order 2 in 2.5D
template <>
inline Real evaluate_point<6,2,3>(const Element<6,2,3>& t, const Point<3>& point, const Eigen::Matrix<Real,6,1>& coefficients)
{
	Eigen::Matrix<Real,3,1> lambda = t.getBaryCoordinates(point);
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
inline Real evaluate_point<10,3,3>(const Element<10,3,3>& t, const Point<3>& point, const Eigen::Matrix<Real,10,1>& coefficients)
{
 Eigen::Matrix<Real,4,1> lambda = t.getBaryCoordinates(point);
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


// This covers the order 1 case
template <UInt Nodes,UInt mydim, UInt ndim>
inline Eigen::Matrix<Real,ndim,1> evaluate_der_point(const Element<Nodes,mydim,ndim>& t, const Point<ndim>& point, const Eigen::Matrix<Real,Nodes,1>& coefficients)
{
  Eigen::Matrix<Real,mydim,mydim+1> B1;
  B1.col(0).setConstant(-1);
  B1.rightCols(mydim).setIdentity();

  return(t.getM_invJ().transpose()*B1*coefficients);
}


template <>
inline Eigen::Matrix<Real,2,1> evaluate_der_point<6,2,2>(const Element<6,2,2>& t, const Point<2>& point, const Eigen::Matrix<Real,6,1>& coefficients)
{
  // Multiply by 4 for convenience in assemblying B2
	Eigen::Matrix<Real,3,1> bar4 = 4*t.getBaryCoordinates(point);
	Eigen::Matrix<Real,2,6> B2;
  B2 << 1-bar4(0), bar4(1)-1,           0, bar4(2),        -bar4(2), bar4(0)-bar4(1),
        1-bar4(0),           0, bar4(2)-1, bar4(1), bar4(0)-bar4(2),        -bar4(1);
	return(t.getM_invJ().transpose()*B2*coefficients);
}

template <>
inline Eigen::Matrix<Real,3,1> evaluate_der_point<6,2,3>(const Element<6,2,3>& t, const Point<3>& point, const Eigen::Matrix<Real,6,1>& coefficients)
{
  // Multiply by 4 for convenience in assemblying B2
	Eigen::Matrix<Real,3,1> bar4 = 4*t.getBaryCoordinates(point);
	Eigen::Matrix<Real,2,6> B2;
  B2 << 1-bar4(0), bar4(1)-1,           0, bar4(2),        -bar4(2), bar4(0)-bar4(1),
        1-bar4(0),           0, bar4(2)-1, bar4(1), bar4(0)-bar4(2),        -bar4(1);
	return(t.getM_invJ().transpose()*B2*coefficients);
}

// Full specialization for order 2 in 3D
// MEMO: this works assuming edges are ordered like so: (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
template <>
inline Eigen::Matrix<Real,3,1> evaluate_der_point<10,3,3>(const Element<10,3,3>& t, const Point<3>& point, const Eigen::Matrix<Real,10,1>& coefficients)
{
  // Multiply by 4 for convenience in assemblying B2
	Eigen::Matrix<Real,4,1> bar4 = 4*t.getBaryCoordinates(point);
	Eigen::Matrix<Real,3,10> B2;
  B2 << 1-bar4(0),  bar4(1)-1,          0,         0, bar4(0)-bar4(1),        -bar4(2),        -bar4(3),  bar4(2),        0, bar4(3),
        1-bar4(0),           0, bar4(2)-1,         0,        -bar4(1), bar4(0)-bar4(2),        -bar4(3),  bar4(1),  bar4(3),       0,
        1-bar4(0),           0,         0, bar4(3)-1,        -bar4(1),        -bar4(2), bar4(0)-bar4(3),        0,  bar4(2), bar4(1);

	return(t.getM_invJ().transpose()*B2*coefficients);
}




#include "mesh_objects_imp.h"
#endif
