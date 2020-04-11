#ifndef MESH_H_
#define MESH_H_

#include "fdaPDE.h"
#include "mesh_objects.h"


//Simple function to convert from order to NNODES at compile time
//Needed to uniform implementation of meshes
//Note: factorial function is defined in mesh_objects.h
constexpr UInt how_many_nodes(UInt order, UInt mydim) {
    return factorial(order+mydim)/(factorial(order)*factorial(mydim));
}


template <UInt ORDER, UInt mydim, UInt ndim>
class MeshHandlerCore{
public:
	using meshElement=Element<how_many_nodes(ORDER,mydim),mydim,ndim>;
	//! A constructor.
		/*!
			* The constructor permits the initialization of the mesh from an R object
			* constructed with the TriLibrary (our R wrapper for the Triangle library)
			* in 2D (in 2.5D and 3D R functions can produce a compatible object if the
			* triangulation is already available)
		*/

	MeshHandlerCore(Real* points, UInt* sides, UInt* elements, UInt* neighbors, UInt num_nodes, UInt num_sides, UInt num_elements):
			points_(points), sides_(sides), elements_(elements), neighbors_(neighbors), num_nodes_(num_nodes), num_sides_(num_sides), num_elements_(num_elements) {};

	#ifdef R_VERSION_
	MeshHandlerCore(SEXP Rmesh);
	#endif

  // No copy, no assignmente for this class!
  // They are not needed and since this class will contain pointers better safe than sorry!
  MeshHandlerCore(const MeshHandlerCore&) = delete;
  MeshHandlerCore &operator=(const MeshHandlerCore&) = delete;

	virtual ~MeshHandlerCore()=0;

	//! A normal member returning an unsigned integer value.
		/*!
			\return The number of nodes in the mesh
		*/
	UInt num_nodes() const {return num_nodes_/ndim;}

	//! A normal member returning an unsigned integer value.
		/*!
			\return The number of elements in the mesh
		*/
	UInt num_elements() const {return num_elements_/how_many_nodes(ORDER,mydim);}

	//! A normal member returning an unsigned integer value.
		/*!
			\return The number of distinct sides (edges for mydim=2,
			faces for mydim=3) in the mesh
		*/
	UInt num_sides() const {return num_sides_/mydim;}

	//! A normal member returning a n-dimensional Point
		/*!
		 * \param id an Id argument
			\return The point with the specified id
		*/
	Point<ndim> getPoint(Id id) const;

	//! A normal member returning an Element
		/*!
		 * \param id an Id argument
			\return The element with order coerent to that of the mesh with the specified id
		*/
	meshElement getElement(Id id) const;

	//The "number" neighbor of element i is opposite the "number" vertex of element i
		//! A normal member returning the Neighbors of a element
		/*!
		 * \param id the id of the element
		 * \param number the number of the vertex
			\return The element that has as a side the one opposite to the specified
			vertex
		*/
	meshElement getNeighbors(Id id_element, UInt number) const;

	void printPoints(std::ostream & out);
	void printElements(std::ostream & out);
	void printNeighbors(std::ostream & out);

	//! A normal member returning the element on which a point is located
		/*!
		 * This method implements a simply research between all the elements of the mesh
		 * \param point the point we want to locate
			\return The element that contains the point
		*/
	meshElement findLocationNaive(const Point<ndim> point) const;


protected:
	#ifdef R_VERSION_
	SEXP mesh_;
	#endif
	Real *points_;
	UInt *sides_;
	UInt *elements_;
	UInt *neighbors_;

	UInt num_nodes_, num_sides_, num_elements_;

};

// Additional methods for 2D and 3D
template <UInt ORDER, UInt mydim, UInt ndim>
class MeshHandler : public MeshHandlerCore<ORDER, mydim, ndim>{
public:
  using meshElement=Element<how_many_nodes(ORDER,mydim),mydim,ndim>;

  MeshHandler(Real* points, UInt* sides, UInt* elements, UInt* neighbors, UInt num_nodes, UInt num_sides, UInt num_elements) :
      MeshHandlerCore<ORDER,mydim,ndim>(points, sides, elements, neighbors, num_nodes, num_sides, num_elements) {}

  //! A normal member returning the element on which a point is located
    /*!
    * This method implements a Visibility Walk Algorithm (further details in: Walking in a triangulation, Devillers et al)
    * \param point the point we want to locate
    * \param starting_elements a vector of points that specifies the poposed starting
    * points for the walking algorithm
    \return The element that contains the point
    */
   meshElement findLocationWalking(const Point<ndim>& point, const meshElement& starting_element) const;

};

// Useful to add some methods peculiar to surface meshes
template <UInt ORDER>
class MeshHandler<ORDER,2,3> : public MeshHandlerCore<ORDER,2,3>{
public:
  MeshHandler(Real* points, UInt* sides, UInt* elements, UInt* neighbors, UInt num_nodes, UInt num_sides, UInt num_elements) :
      MeshHandlerCore<ORDER,2,3>(points, sides, elements, neighbors, num_nodes, num_sides, num_elements) {}

  // This function projects points onto the mesh
  std::vector<Point<3> > project(const std::vector<Point<3> >&) const;
private:
  // This function computes the closest nodes to the given points and returns their index
  std::vector<UInt> find_closest(const std::vector<Point<3> >&) const;
};



#include "mesh_imp.h"

#endif
