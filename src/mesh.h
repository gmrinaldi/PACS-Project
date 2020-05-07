#ifndef MESH_H_
#define MESH_H_

#include "fdaPDE.h"
#include "mesh_objects.h"


template <UInt ORDER, UInt mydim, UInt ndim>
class MeshHandlerCore{
static_assert(ORDER==1 || ORDER==2,
  "ERROR: Only first and second order case implemented for now; see mesh.h");
public:
  // Note: how_many_nodes function is defined in mesh_objects.h
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

	virtual ~MeshHandlerCore()=0;

	//! A normal member returning an unsigned integer value.
		/*!
			\return The number of nodes in the mesh
		*/
	UInt num_nodes() const {return num_nodes_;}

	//! A normal member returning an unsigned integer value.
		/*!
			\return The number of elements in the mesh
		*/
	UInt num_elements() const {return num_elements_;}

	//! A normal member returning an unsigned integer value.
		/*!
			\return The number of distinct sides (edges for mydim=2,
			faces for mydim=3) in the mesh
		*/
	UInt num_sides() const {return num_sides_;}

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

	void printPoints(std::ostream &);
	void printElements(std::ostream &);
	void printNeighbors(std::ostream &);

	//! A normal member returning the element on which a point is located
		/*!
		 * This method implements a simply research between all the elements of the mesh
		 * \param point the point we want to locate
			\return The element that contains the point
		*/
	meshElement findLocationNaive(const Point<ndim>&) const;


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

  //! Additional members returning the number of edges/faces of the element for convenience
  UInt num_edges() const {return this->num_sides();}
  UInt num_faces() const {return this->num_sides();}

  //! A normal member returning the element on which a point is located
    /*!
    * This method implements a Visibility Walk Algorithm (further details in: Walking in a triangulation, Devillers et al)
    * \param point the point we want to locate
    * \param starting_elements a vector of points that specifies the poposed starting
    * points for the walking algorithm
    \return The element that contains the point
    */
   meshElement findLocationWalking(const Point<ndim>&, const meshElement&) const;

};

// Useful to add some methods peculiar to surface meshes
template <UInt ORDER>
class MeshHandler<ORDER,2,3> : public MeshHandlerCore<ORDER,2,3>{
public:
  MeshHandler(Real* points, UInt* sides, UInt* elements, UInt* neighbors, UInt num_nodes, UInt num_sides, UInt num_elements) :
      MeshHandlerCore<ORDER,2,3>(points, sides, elements, neighbors, num_nodes, num_sides, num_elements) {}

  // Additional member returning the number of edges added for convenience
  UInt num_edges() const {return this->num_sides();}

  // This function projects points onto the mesh
  std::vector<Point<3> > project(const std::vector<Point<3> >&) const;
private:
  // This function computes the closest nodes to the given points and returns their index
  std::vector<UInt> find_closest(const std::vector<Point<3> >&) const;
};



#include "mesh_imp.h"

#endif
