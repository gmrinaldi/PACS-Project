#ifndef __FINITE_ELEMENT_HPP__
#define __FINITE_ELEMENT_HPP__


#include "integration.h"
#include "mesh_objects.h"
#include "param_functors.h"
#include <type_traits>


//!  This class implements all properties of a Triangular or Tetrahedral Finite Element
/*!
 * This class is the most important one of the entire code
 * and implements everything needed by a triangular or tetrahedral finite elemnt
 *
 * It takes as a template parameter a class that implements the method used
 * for determining the mass, stiff and grad matrices
*/

template <UInt ORDER, UInt mydim, UInt ndim>
class FiniteElement{
public:
	using Integrator = typename std::conditional<mydim==2,
												typename std::conditional<ORDER==1, IntegratorTriangleP2, IntegratorTriangleP4>::type,
												typename std::conditional<ORDER==1, IntegratorTetrahedronP2, IntegratorTetrahedronP4>::type>::type;

	static constexpr UInt NBASES=how_many_nodes(ORDER,mydim);

	using return_t = Eigen::Matrix<Real,NBASES,NBASES>;


private:

	Element<NBASES,mydim,ndim> t_;
	Eigen::Matrix<Real, NBASES, Integrator::NNODES> phiMapMaster_;
	//Num local bases x num coords x num quad nodes
	Eigen::Matrix<Real, NBASES, Integrator::NNODES*mydim> phiDerMapMaster_;
	Eigen::Matrix<Real, NBASES, Integrator::NNODES*ndim> invTrJPhiDerMapMaster_;

	void setPhiMaster();
	void setPhiDerMaster();
	void setInvTrJPhiDerMaster();

public:
	//! This is an empty constructor
    /*!
        For efficiency and Expression Templates organization of the
        code, the use of this class is based on the updateElement class
    */

	FiniteElement();

	//! A member updating the Finite Element properties
    /*!
      \param t an element from which to update the finite element properties
    */
	void updateElement(const Element<NBASES,mydim,ndim> &t);

	const Point<ndim>& operator[](UInt i) const {return t_[i];}
	// Both Area and Volume getter added for convenience
	Real getMeasure() const {return t_.getMeasure();}
	Real getArea() const {return t_.getMeasure();}
	Real getVolume() const {return t_.getMeasure();}

	Real getId() const {return t_.getId();}

	Point<ndim> coorQuadPt(UInt iq){
		return Point<ndim>(t_.getM_J()*Eigen::Map<const Eigen::Matrix<Real,mydim,1> >(&Integrator::NODES[iq][0]))+=t_[0];
	}

	//Returns \hat{phi}
	Real phiMaster(UInt i, UInt iq) const {return phiMapMaster_(i, iq);}
	UInt getGlobalIndex(UInt iq) {return Integrator::NNODES * t_.getId() + iq;}

	auto stiff_impl(UInt iq) const -> decltype(invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq) * invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq).transpose()) {
		// compute nabla(phi) x nabla(phi) at node iq
		return invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq) * invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq).transpose();
	}

	auto stiff_impl(UInt iq, const Diffusion<ndim>& K) const -> decltype(invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq) * K() * invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq).transpose()) {
		// compute nabla(phi) x K nabmyla(phi) at node iq
		// Memo: K is symmetric!
		return invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq) * K() * invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq).transpose();
	}

	auto stiff_impl(UInt iq, const DiffusionSV<ndim>& K) const -> decltype(invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq) * K(0) * invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq).transpose()) {
		// compute nabla(phi) x K nabmyla(phi) at node iq
		// Memo: K is symmetric!
		UInt globalIndex = this->getGlobalIndex(iq);
		return invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq) * K(globalIndex) * invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq).transpose();
	}


	auto mass_impl(UInt iq) const -> decltype(phiMapMaster_.col(iq)*phiMapMaster_.col(iq).transpose()) {
		// compute phi x phi at node iq
		return phiMapMaster_.col(iq)*phiMapMaster_.col(iq).transpose();
	}

	auto grad_impl(UInt iq) const -> decltype(phiMapMaster_.col(iq) * invTrJPhiDerMapMaster_.col(ndim*iq).transpose()) {
		// compute phi x nabla(phi) x b at node iq
		// default is b=u_x
		return phiMapMaster_.col(iq) * invTrJPhiDerMapMaster_.col(ndim*iq).transpose();
	}

	auto grad_impl(UInt iq, const Advection<ndim>& b) const -> decltype(phiMapMaster_.col(iq) * b().transpose() * invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq).transpose()) {
		// compute phi x nabla(phi) x b at node iq
		return phiMapMaster_.col(iq) * b(0).transpose() * invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq).transpose();
	}

	auto grad_impl(UInt iq, const AdvectionSV<ndim>& b) const -> decltype(phiMapMaster_.col(iq) * b(0).transpose() * invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq).transpose()) {
		// compute phi x nabla(phi) x b at node iq
		UInt globalIndex = this->getGlobalIndex(iq);
		return phiMapMaster_.col(iq) * b(iq).transpose() * invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq).transpose();
	}




};

#include "finite_element_imp.h"

#endif
