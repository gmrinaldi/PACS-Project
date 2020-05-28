#ifndef __MESH_OBJECTS_HPP__
#define __MESH_OBJECTS_HPP__


#include "fdaPDE.h"

//Simple function to evaluate factorial at compile time
//Needed for getVolume/getArea member of Element
constexpr UInt factorial(UInt n) {
    return n ? (n * factorial(n - 1)) : 1;
}

//Simple function to convert from order to NNODES at compile time
//Computes how many nodes in an element of given order and dimension
constexpr UInt how_many_nodes(UInt ORDER, UInt mydim) {
    return factorial(ORDER+mydim)/(factorial(ORDER)*factorial(mydim));
}

//!  This class gives some common methods to all mesh objects.
struct Identifier{

	//! An static const Unsigned Integer.
    /*! Needed to identify the Not Valid Id. */
	static constexpr UInt NVAL = std::numeric_limits<UInt>::max();

	bool unassignedId() const {return id_==NVAL;}
  bool hasValidId() const {return !unassignedId();}
	bool unassignedBc() const {return bcId_==NVAL;}

	UInt id() const {return id_;}
	UInt bcId() const {return bcId_;}
	UInt getId() const {return id_;}

protected:
  UInt id_=NVAL;
	UInt bcId_=NVAL;

  constexpr Identifier()=default;
  constexpr Identifier(UInt id) : id_(id) {}
  constexpr Identifier(UInt id, UInt bcId) : id_(id), bcId_(bcId) {}

  ~Identifier()=default;
};


// This class implements a n-dimensional point
// Doesn't allow creation of points in dimension other than 2 or 3
template<UInt ndim>
class Point : public Identifier{
  static_assert(ndim==2 || ndim==3,
								 "ERROR! TRYING TO INSTANTIATE POINT IN UNIMPLEMENTED DIMENSION! See mesh_objects.h");
  public:
    using pointCoords = std::array<Real,ndim>;
    using EigenCoords = Eigen::Matrix<Real,ndim,1>;
    // Default constructor initializing the origin
    constexpr Point()=default;
    // Full constructor initializing every member explicitly
    constexpr Point(UInt id, UInt bcId, const pointCoords& coord) :
              Identifier(id, bcId), coord_(coord) {}
    // Additional constructors initializing only some members
    constexpr Point(const pointCoords& coord) :
              coord_(coord) {}
    constexpr Point(UInt id, const pointCoords& coord) :
              Identifier(id), coord_(coord) {}
    // Additional constructor (needed for disambiguation when constructing from
    // bracketed initializer lists)
    constexpr Point(const Real(&coord)[ndim]);

    // Additional constructor allowing conversion from eigen objects
    Point(const EigenCoords &coord);

    // Additional constructor for convenience in dealing with meshes
    Point(UInt id, const Real* const points, const UInt num_points);

    Real& operator[](UInt i) {return coord_[i];}
    const Real& operator[](UInt i) const {return coord_[i];}

    // Member/friend functions returning the distance/squared distance between two points
    Real dist(const Point&) const;
    Real dist2(const Point&) const;

    friend Real dist(const Point& lhs, const Point& rhs) {return lhs.dist(rhs);};
    friend Real dist2(const Point& lhs, const Point& rhs) {return lhs.dist2(rhs);};

    Point& operator+=(const Point&);
    Point& operator-=(const Point&);
    // Overload the "+" ("-") operator to take 2 points of the same dimension and compute
    // the coordinate sum (difference).
    friend Point operator+(Point lhs, const Point& rhs) {return lhs+=rhs;};
    friend Point operator-(Point lhs, const Point& rhs) {return lhs-=rhs;};


    //! Overload the << operator to easily print Point info (note: define it in class
    // to avoid a forward declaration)
    friend std::ostream& operator<<(std::ostream& os, const Point& p){
      if(p.hasValidId())
        os<<p.getId()<<":";
      for (const auto &c : p.coord_)
        os<<" "<<c;
      return os<<std::endl;
    }

  private:
    pointCoords coord_;
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
class Element : public Identifier {

  static_assert((mydim==2 || mydim==3) &&
                 mydim <= ndim &&
                 (NNODES==how_many_nodes(1,mydim) || NNODES==how_many_nodes(2,mydim)),
                 "ERROR! TRYING TO INSTANTIATE ELEMENT WITH WRONG NUMBER OF NODES AND/OR DIMENSIONS! See mesh_objects.h");

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  using elementPoints = std::array<Point<ndim>,NNODES>;
  using iterator = typename elementPoints::iterator;
  using const_iterator = typename elementPoints::const_iterator;

  //! This constructor creates an "empty" Element, with an Id Not Valid
  Element()=default;

	//! This constructor creates an Element, given its Id and an std array with the Points
	Element(UInt id, const elementPoints& points) :
					Identifier(id), points_(points) {computeProperties();}

	//! Overloading of the operator [],  taking the Node number and returning a node as a reference to a Point object.
  // Point<ndim>& operator[](UInt i) {return points_[i];}
	const Point<ndim>& operator[](UInt i) const {return points_[i];}

  //! Define begin and end iterators (this also gives "ranged for" for free)
  // iterator begin() {return points_.begin();}
  // iterator end() {return points_.end();}

  //! Define const begin and end iterators
  const_iterator begin() const {return points_.begin();}
  const_iterator end() const {return points_.end();}

	const Eigen::Matrix<Real,ndim,mydim>& getM_J() const {return M_J_;}
	const Eigen::Matrix<Real,mydim,ndim>& getM_invJ() const {return M_invJ_;} //in 2.5D this is actually the pseudoinverse!

  //! A member that computes the barycentric coordinates.
	Eigen::Matrix<Real,mydim+1,1> getBaryCoordinates(const Point<ndim>& point) const;

  //! A member that tests if a Point is located inside an Element.
  bool isPointInside(const Point<ndim>& point) const;

  //! A member that verifies which edge/face separates the Triangle/Tetrahedron from a Point.
  // This function is not implemented for manifold data (we make sure of this at compile time using enable_if)
  int getPointDirection(const Point<ndim>&) const;

  //! Some members returning the area/volume of the element
  Real getMeasure() const {return element_measure;}
  Real getArea() const {return this->getMeasure();}
  Real getVolume() const {return this->getMeasure();}

  //! A member to evaluate functions at a point inside the element
  Real evaluate_point(const Point<ndim>&, const Eigen::Matrix<Real,NNODES,1>&) const;
  //! A member to evaluate integrals on the element
  Real integrate(const Eigen::Matrix<Real,NNODES,1>&) const;

	//! Overload the << operator to easily print element info (note: define it in class
  // to avoid a forward declaration)
  friend std::ostream& operator<<(std::ostream& os, const Element& el){
    os<< el.getId() <<":";
    for (const auto &p : el)
      os<<" "<<p.getId();
    return os<<std::endl;
  }

private:
  // Some useful type aliases

  // Data members
	elementPoints points_;
	Eigen::Matrix<Real,ndim,mydim> M_J_;
	Eigen::Matrix<Real,mydim,ndim> M_invJ_; //in 2.5D this is actually the pseudoinverse!
  Real element_measure;

  void computeProperties();
};

template <UInt NNODES>
class Element<NNODES, 2, 3> : public Identifier {

  static_assert(NNODES==3 || NNODES==6,
                 "ERROR! TRYING TO INSTANTIATE SURFACE ELEMENT WITH WRONG NUMBER OF NODES! See mesh_objects.h");

  // Some useful type aliases

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // Some more useful type aliases (these could be useful outside the class
  // so they are public)
  using elementPoints = std::array<Point<3>, NNODES>;
  using iterator = typename elementPoints::iterator;
  using const_iterator = typename elementPoints::const_iterator;

  //! This constructor creates an "empty" Element, with an Id Not Valid
  Element()=default;

	//! This constructor creates an Element, given its Id and an std array with the Points
	Element(UInt id, const elementPoints& points) :
					Identifier(id), points_(points) {computeProperties();}

	//! Overloading of the operator [],  taking the Node number and returning a node as a reference to a Point object.
  // Point<ndim>& operator[](UInt i) {return points_[i];}
	const Point<3>& operator[](UInt i) const {return points_[i];}

  //! Define begin and end iterators (this also gives "ranged for" for free)
  // iterator begin() {return points_.begin();}
  // iterator end() {return points_.end();}

  //! Define const begin and end iterators
  const_iterator begin() const {return points_.begin();}
  const_iterator end() const {return points_.end();}

	const Eigen::Matrix<Real,3,2>& getM_J() const {return M_J_;}
	const Eigen::Matrix<Real,2,3>& getM_invJ() const {return M_invJ_;} //this is actually the pseudoinverse!

  //! A member that computes the barycentric coordinates.
	Eigen::Matrix<Real,3,1> getBaryCoordinates(const Point<3>& point) const;

  //! A member that tests if a Point is located inside an Element.
  bool isPointInside(const Point<3>& point) const;

  //! Some members returning the area/volume of the element
  Real getMeasure() const {return element_measure;}
  Real getArea() const {return this->getMeasure();}
  Real getVolume() const {return this->getMeasure();}

  //! A member to evaluate functions at a point inside the element
  Real evaluate_point(const Point<3>&, const Eigen::Matrix<Real,NNODES,1>&) const;
  //! A member to evaluate integrals on the element
  Real integrate(const Eigen::Matrix<Real,NNODES,1>&) const;

  // This funtion is implemented (and needed) only for manifold elements!
  // This function projects a 3D point (XYZ coordinates!) onto the element
  // Note: if the projection lies outside the element the function returns
  // the closest point on the boundary of the element instead
  Point<3> computeProjection(const Point<3>&) const;

	//! Overload the << operator to easily print element info (note: define it in class
  // to avoid a forward declaration)
  friend std::ostream& operator<<(std::ostream& os, const Element& el){
    os<< el.getId() <<":";
    for (const auto &p : el)
      os<<" "<<p.getId();
    return os<<std::endl;
  }

private:
  // Data members
	elementPoints points_;
	Eigen::Matrix<Real,3,2> M_J_;
	Eigen::Matrix<Real,2,3> M_invJ_; //this is actually the pseudoinverse!
  Real element_measure;

  void computeProperties();

};


#include "mesh_objects_imp.h"
#endif
