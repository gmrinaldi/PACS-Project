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
    Point(Id id, pointCoords coord) :
              Identifier(id), coord_(coord) {}
    Point(Id id, BcId bcId, pointCoords coord) :
              Identifier(id, bcId), coord_(coord) {}

    Real operator[](UInt i) const {return coord_[i];}
    Real distance(const Point&) const;

    // Overload the "+" ("-") operator to take 2 points of the same dimension and compute
    // the coordinate sum (difference).
    // Note: they return an array for convenience, but since an array is convertible to Point
    // something like this is fine Point<ndim> a, b; // Point<ndim> c=a-b;
    friend pointCoords operator +(const Point& lhs, const Point& rhs){
      pointCoords diff{lhs.coord_};
      for (int i=0; i<ndim; ++i)
          diff[i]+=rhs[i];
      return diff;
    };

    friend pointCoords operator -(const Point& lhs, const Point& rhs){
      pointCoords diff{lhs.coord_};
      for (int i=0; i<ndim; ++i)
          diff[i]-=rhs[i];
      return diff;
    };


    friend Real distance(const Point& lhs, const Point& rhs){
      return lhs.distance(rhs);
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

  // Default destructor
  virtual ~ElementCore()=default;

	//! Overloading of the operator [],  taking the Node number and returning a node as Point object.
	Point<ndim> operator[](UInt i) const {return points_[i];}

	Real getDetJ() const {return detJ_;}
	Eigen::Matrix<Real,ndim,mydim>& getM_J() const {return M_J_;}
	Eigen::Matrix<Real,mydim,ndim>& getM_invJ() const {return M_invJ_;} //in 2.5D this is actually the pseudoinverse!
	Eigen::Matrix<Real,mydim,mydim>& getMetric() const {return metric_;}

  //! A member that computes the barycentric coordinates.
	Eigen::Matrix<Real,mydim+1,1> getBaryCoordinates(const Point<ndim>& point) const;

  //! A member that tests if a Point is located inside an Element.
  virtual bool isPointInside(const Point<ndim>& point) const;

  //! A member returning the area/volume of the element
  virtual Real getMeasure() const =0;

  //! A member to evaluate functions in a point inside the element
  Real evaluate_point(const Point<ndim>& point, const Eigen::Matrix<Real,NNODES,1>&) const;
  //! A member to evaluate derivatives at a point inside the element
  Eigen::Matrix<Real,ndim,1> evaluate_der_point(const Point<ndim>& point, const Eigen::Matrix<Real,NNODES,1>&) const;
  //! A member to evaluate integrals on the element
  Real integrate(const Eigen::Matrix<Real,NNODES,1>&) const;

	//! Overload the << operator to easily print element info (note: define it in class
  // to avoid a forward declaration)
  friend std::ostream& operator<<(std::ostream& os, const ElementCore& el){
    os<<"Element "<< el.getId() <<": ";
    for (UInt i=0; i<NNODES; ++i)
      os<<el[i].getId()<<" ";
    os<<std::endl;
    return os;
  }

protected:
	elementPoints points_;
	Eigen::Matrix<Real,ndim,mydim> M_J_;
	Eigen::Matrix<Real,mydim,ndim> M_invJ_; //in 2.5D this is actually the pseudoinverse!
	Eigen::Matrix<Real,mydim,mydim> metric_;
	Real detJ_;

  // pure virtual member to make ElementCore an abstract base class
  virtual void computeProperties()=0;
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
  Real getMeasure() const {return std::abs(this->detJ_)/factorial(ndim);}
  //! Additional members returning the area/volume of the element
  // Note: both are needed for ease of use
  Real getArea() const {return this->getMeasure();}
  Real getVolume() const {return this->getMeasure();}

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
  Real getMeasure() const {return std::sqrt(this->detJ_)/2;}
  //! Additional member returning the area added for ease of use
  Real getArea() const {return this->getMeasure();}

  //! A member that tests if a Point is located inside an Element.
  // Note: this is implemented (slightly) differently in this case
  // because one has to check that the point lies on the same plane as the triangle!
  bool isPointInside(const Point<3>& point) const;

  // This function projects a 3D point (XYZ coordinates!) onto the element
  // Note: if the projection lies outside the element the function returns
  // the closest point on the boundary of the element instead
  Point<3> computeProjection(const Point<3>& point) const;

private:
  void computeProperties();

};



#include "mesh_objects_imp.h"
#endif
