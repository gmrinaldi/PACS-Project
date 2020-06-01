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

//! This class contains the ID of the objects and gives some common methods
// to get the ID and check if it is valid
struct Identifier{

  // The biggest UInt number represents a not valid ID (NVAL)
	static constexpr UInt NVAL = std::numeric_limits<UInt>::max();

  // Some methods to check if the object has a valid ID
	bool unassignedId() const {return id_==NVAL;}
  bool hasValidId() const {return !unassignedId();}
	bool unassignedBc() const {return bcId_==NVAL;}

  // Some methods to get the object ID
	UInt id() const {return id_;}
	UInt bcId() const {return bcId_;}
	UInt getId() const {return id_;}

protected:
  // Default object ID is NVAL if not explicitly set
  UInt id_=NVAL;
  // Note: bcId was probably intended to mark boundary objects but is never actually
  // used, so it could be removed
	UInt bcId_=NVAL;

  // Protected constructors make the class act as an abstract base class
  // without any runtime cost due to virtual functions
  // Note: these constructors are constexpr so that any derived class with a
  // constexpr constructor can work as a literal type (this is particularly relevant for Point)
  constexpr Identifier()=default;
  constexpr Identifier(UInt id) : id_(id) {}
  constexpr Identifier(UInt id, UInt bcId) : id_(id), bcId_(bcId) {}

  ~Identifier()=default;
};



//! This class implements a Point in dimension 2 or 3
template<UInt ndim>
class Point : public Identifier{
  static_assert(ndim==2 || ndim==3,
								 "ERROR! TRYING TO INSTANTIATE POINT IN UNIMPLEMENTED DIMENSION! See mesh_objects.h");
  public:
    using pointCoords = std::array<Real,ndim>;
    using EigenCoords = Eigen::Matrix<Real,ndim,1>;

    // Note: some of the constructors are declared constexpr so that Point can be
    // used as a literal type (see integration.h)
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

    // Additional constructor allowing construction from an eigen vector
    Point(const EigenCoords &coord);

    // Additional constructor for convenience in dealing with R data (e.g. meshes)
    Point(UInt id, const Real* const points, const UInt num_points);

    // Overloaded subscript operator
    Real& operator[](UInt i) {return coord_[i];}
    const Real& operator[](UInt i) const {return coord_[i];}

    // Member functions returning the distance/squared distance between two points
    // Note: both are defined because the squared distance function is more efficient
    // and should be preferred if one is only interested in comparing distances
    Real dist(const Point&) const;
    Real dist2(const Point&) const;
    // Friend functions returning the distance/squared distance between two points
    friend Real dist(const Point& lhs, const Point& rhs) {return lhs.dist(rhs);};
    friend Real dist2(const Point& lhs, const Point& rhs) {return lhs.dist2(rhs);};

    // Overloaded "+="/"-=" operators
    // These operators add/subtract the coordinates elementwise
    Point& operator+=(const Point&);
    Point& operator-=(const Point&);
    // Overloaded "+"/"-" operator
    // These operators add/subtract the coordinates elementwise and return a new point
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


//! This class implements an Element (i.e. a triangle or tetrahedron)
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


//!  This class implements a Tetrahedron as an objects composed by four or ten nodes, embedded in a 3-dimensional space.
// For NNODES=10 the edges are ordered like so: (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
// The midpoints are also expected to follow this convention!

// Note: only order 1 and 2 are implemented at the moment

template <UInt NNODES, UInt mydim, UInt ndim>
class Element : public Identifier {

  static_assert((mydim==2 || mydim==3) &&
                 mydim <= ndim &&
                 (NNODES==how_many_nodes(1,mydim) || NNODES==how_many_nodes(2,mydim)),
                 "ERROR! TRYING TO INSTANTIATE ELEMENT WITH WRONG NUMBER OF NODES AND/OR DIMENSIONS! See mesh_objects.h");

public:
  // Note: this macro is needed to avoid alignment issues!
  // See Eigen documentation "Structures Having Eigen Members"
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  using elementPoints = std::array<Point<ndim>,NNODES>;
  using iterator = typename elementPoints::iterator;
  using const_iterator = typename elementPoints::const_iterator;

  //! This constructor creates an "empty" Element, with an Id Not Valid
  Element()=default;

	//! This constructor creates an Element, given its Id and an array containing the Points
	Element(UInt id, const elementPoints& points) :
					Identifier(id), points_(points) {computeProperties();}

	//! Overloaded subscript operator
  // Note: only the const version is available because any change in points_
  // would require a call to computeProperties() to keep the element in a valid state
	const Point<ndim>& operator[](UInt i) const {return points_[i];}

  //! Define begin and end iterators (this also gives "ranged for" for free)
  // Note: only the const version is available because any change in points_
  // would require a call to computeProperties() to keep the element in a valid state
  const_iterator begin() const {return points_.begin();}
  const_iterator end() const {return points_.end();}

	const Eigen::Matrix<Real,ndim,mydim>& getM_J() const {return M_J_;}
	const Eigen::Matrix<Real,mydim,ndim>& getM_invJ() const {return M_invJ_;}

  //! A member function that computes the barycentric coordinates of a given Point
	Eigen::Matrix<Real,mydim+1,1> getBaryCoordinates(const Point<ndim>& point) const;

  //! A member function that tests if a given Point is located inside the Element
  bool isPointInside(const Point<ndim>& point) const;

  //! A member function that returns beyond which edge/face a given Point lies
  // This function returns -1 if the point is inside the element
  int getPointDirection(const Point<ndim>&) const;

  //! Some members returning the area/volume of the element
  // Note: all three are available for convenience and compatibility reasons with existing code
  Real getMeasure() const {return element_measure;}
  Real getArea() const {return getMeasure();}
  Real getVolume() const {return getMeasure();}

  //! A member to evaluate a function at a point given the function's coefficients
  // on the element's basis functions
  // Note: this function assumes that the point is inside the element!
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
  // Data members
	elementPoints points_;
  // A matrix encoding a linear transformation from barycentric coordinates to
  // standard coordinates (modulo a translation)
	Eigen::Matrix<Real,ndim,mydim> M_J_;
  // A matrix which is the inverse of M_J_
	Eigen::Matrix<Real,mydim,ndim> M_invJ_;
  // A Real storing the area/volume of the element
  Real element_measure;

  // A member initializing M_J_, M_invJ_ and element_measure at construction
  void computeProperties();
};


//! Partial template specialization for surface elements
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

template <UInt NNODES>
class Element<NNODES, 2, 3> : public Identifier {

  static_assert(NNODES==3 || NNODES==6,
                 "ERROR! TRYING TO INSTANTIATE SURFACE ELEMENT WITH WRONG NUMBER OF NODES! See mesh_objects.h");

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  using elementPoints = std::array<Point<3>, NNODES>;
  using iterator = typename elementPoints::iterator;
  using const_iterator = typename elementPoints::const_iterator;

  //! This constructor creates an "empty" Element, with an Id Not Valid
  Element()=default;

	//! This constructor creates an Element, given its Id and an std array with the Points
	Element(UInt id, const elementPoints& points) :
					Identifier(id), points_(points) {computeProperties();}

  //! Overloaded subscript operator
  // Note: only the const version is available because any change in points_
  // would require a call to computeProperties() to keep the element in a valid state
	const Point<3>& operator[](UInt i) const {return points_[i];}

  //! Define begin and end iterators (this also gives "ranged for" for free)
  // Note: only the const version is available because any change in points_
  // would require a call to computeProperties() to keep the element in a valid state
  const_iterator begin() const {return points_.begin();}
  const_iterator end() const {return points_.end();}

	const Eigen::Matrix<Real,3,2>& getM_J() const {return M_J_;}
	const Eigen::Matrix<Real,2,3>& getM_invJ() const {return M_invJ_;} //this is actually the pseudoinverse!

  //! A member function that computes the barycentric coordinates of a given point
	Eigen::Matrix<Real,3,1> getBaryCoordinates(const Point<3>& point) const;

  //! A member function that tests if a Point is located inside the Element
  bool isPointInside(const Point<3>& point) const;

  //! Some members returning the area/volume of the element
  // Note: all three are available for convenience and compatibility reasons with existing code
  Real getMeasure() const {return element_measure;}
  Real getArea() const {return getMeasure();}
  Real getVolume() const {return getMeasure();}

  //! A member to evaluate functions at a point inside the element
  Real evaluate_point(const Point<3>&, const Eigen::Matrix<Real,NNODES,1>&) const;
  //! A member to evaluate integrals on the element
  Real integrate(const Eigen::Matrix<Real,NNODES,1>&) const;

  // A member function that projects a 3D point (XYZ coordinates!) onto the element
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
  // A matrix encoding a linear transformation from barycentric coordinates to
  // standard coordinates (modulo a translation)
	Eigen::Matrix<Real,3,2> M_J_;
  // A matrix which is the pseudoinverse of M_J_
	Eigen::Matrix<Real,2,3> M_invJ_;
  // A Real storing the area/volume of the element
  Real element_measure;

  // A member initializing M_J_, M_invJ_ and element_measure at construction
  void computeProperties();

};


#include "mesh_objects_imp.h"
#endif
