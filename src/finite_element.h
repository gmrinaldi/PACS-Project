#ifndef __FINITE_ELEMENT_HPP__
#define __FINITE_ELEMENT_HPP__

// Needed for IntegratorHelper
#include "integration.h"

// Needed for element, point and how_many_nodes
#include "mesh_objects.h"

// This is an abstract base class that wraps Element objects
// It stores the data needed by a triangular or tetrahedral finite element
template <UInt ORDER, UInt mydim, UInt ndim>
class FiniteElementData{
	static_assert((ORDER==1 || ORDER==2) &&
								(mydim==2 || mydim==3) &&
								 mydim <= ndim,
								 "ERROR! TRYING TO INSTANTIATE FINITE ELEMENT WITH WRONG NUMBER OF NODES AND/OR DIMENSIONS! See finite_element.h");
public:
  using meshSide = std::array<UInt, mydim>;

	// This type encodes an appropriate quadrature rule depending on order and dimension
	using Integrator = typename IntegratorHelper::Integrator<ORDER,mydim>;

	// Number of basis function on the element
	static constexpr UInt NBASES = how_many_nodes(ORDER,mydim);

	// A default constructor
	FiniteElementData();

	// No move/copy constructors/operations
	// Since the class acts as a wrapper they are not needed!
	// Moreover this guarantees that the class can't be passed by value
	// (it is less error prone)
	FiniteElementData(const FiniteElementData&) = delete;
	FiniteElementData(FiniteElementData&&) = delete;
	FiniteElementData& operator=(const FiniteElementData&) = delete;
	FiniteElementData& operator=(FiniteElementData&&) = delete;

	// Pure virtual destructor making the class an ABC
	// Note: this has negligible runtime cost for this class
	virtual ~FiniteElementData()=0;

	// A member that accepts a new element to wrap
	void updateElement(const Element<NBASES,mydim,ndim>& t);

	// Overloaded subscript operator
	// It returns the i-th point of the underlying element
	// Note: read-only access because changing a point in an element
	// currently invalidates its state (see: "mesh_objects.h")
	const Point<ndim>& operator[] (UInt i) const {return t_[i];}

	// Members returning the area/volume of the underlying element
	Real getMeasure() const {return t_.getMeasure();}
	Real getArea() const {return t_.getMeasure();}
	Real getVolume() const {return t_.getMeasure();}

	// A member returning the ID of the underlying element
	Real getId() const {return t_.getId();}

	// A member returning the iq-th quadrature point as a point in ndim-ensional space
	Point<ndim> coorQuadPt(UInt iq){
		Point<ndim> quadPt{t_.getM_J() * Eigen::Map<const Eigen::Matrix<Real,mydim,1> >(&Integrator::NODES[iq][0])};
		return quadPt += t_[0];
	}

	// A member returning ^phi(i,iq)
	Real getPhi(UInt i, UInt iq) const {return referencePhi(i, iq);}
	// A member returning the global index of a quadrature node
	UInt getGlobalIndex(UInt iq) const {return Integrator::NNODES * t_.getId() + iq;}

protected:
	// The underlying element
	Element<NBASES,mydim,ndim> t_;
	// A matrix Phi
	// Phi(i,iq) is ^phi_i(node_iq), i.e. the i-th basis function evaluated at the
	// iq-th quadrature node on the reference element
	Eigen::Matrix<Real, NBASES, Integrator::NNODES> referencePhi;
	// A block matrix PhiDer
	// Each block iq is made of nabla ^phi(node_iq), i.e. it stores the gradients
	// of the basis functions evaluated at the quadrature nodes on the reference element
	Eigen::Matrix<Real, mydim, NBASES*Integrator::NNODES> referencePhiDer;
	// A block matrix elementPhiDer
	// Each block iq is made of nabla phi(node_iq), i.e. it stores the gradients
	// of the basis functions evaluated at the quadrature nodes on the underlying element
	Eigen::Matrix<Real, ndim, NBASES*Integrator::NNODES> elementPhiDer;

	// Members initializing the corresponding matrices at construction
	void setPhi();
	void setPhiDer();
	// A member updating elementPhiDer after each element update
	void setElementPhiDer();

};


// This class implements all the needed methods to assemble the FE matrix

template <UInt ORDER, UInt mydim, UInt ndim>
struct FiniteElement : public FiniteElementData<ORDER, mydim, ndim> {

	using Integrator = typename FiniteElementData<ORDER, mydim, ndim>::Integrator;

	static constexpr UInt NBASES = FiniteElementData<ORDER, mydim, ndim>::NBASES;

	FiniteElement()=default;

	Real stiff_impl(UInt iq, UInt i, UInt j) const {
		return this->elementPhiDer.col(iq*NBASES+i).dot(this->elementPhiDer.col(iq*NBASES+j));
	}

	Real stiff_impl(UInt iq, UInt i, UInt j, const Eigen::Matrix<Real, ndim, ndim>& K) const {
		return this->elementPhiDer.col(iq*NBASES+i).dot(K * this->elementPhiDer.col(iq*NBASES+j));
	}

	Real mass_impl(UInt iq, UInt i, UInt j) const {
		return this->referencePhi(i,iq) * this->referencePhi(j,iq);
	}

	Real grad_impl(UInt iq, UInt i, UInt j) const {
		return this->referencePhi(i,iq) * this->elementPhiDer(0,iq*NBASES+j);
	}

	Real grad_impl(UInt iq, UInt i, UInt j, const Eigen::Matrix<Real, ndim, 1>& b) const {
		return this->referencePhi(i,iq) * b.dot(this->elementPhiDer.col(iq*NBASES+j));
	}

};

#include "finite_element_imp.h"

#endif
