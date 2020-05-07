#ifndef __FINITE_ELEMENT_HPP__
#define __FINITE_ELEMENT_HPP__

#include "fdaPDE.h"
#include "integration.h"
#include "mesh_objects.h"

//!  This class implements all properties of a Triangular or Tetrahedral Finite Element
/*!
 * This class is the most important one of the entire code
 * and implements everything needed by a triangular or tetrahedral finite elemnt
 *
 * It takes as a template parameter a class that implements the method used
 * for determining the mass, stiff and grad matrices
*/

template <class Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FiniteElement{
private:
	static constexpr UInt NBASES=how_many_nodes(ORDER,mydim);
	Element<NBASES,mydim,ndim> t_;
	Eigen::Matrix<Real, NBASES, Integrator::NNODES> phiMapMaster_;
	//Numero basi locali x Num coordinate x numero nodi integrazione
	Eigen::Matrix<Real, NBASES, Integrator::NNODES*mydim> phiDerMapMaster_;
	Eigen::Matrix<Real, NBASES, Integrator::NNODES*mydim> invTrJPhiDerMapMaster_;

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

	// Both Area and Volume getter added for convenience
	constexpr Real getAreaReference()
	{
		return 1./factorial(mydim);
	}

	constexpr Real getVolumeReference()
	{
		return 1./factorial(mydim);
	}

	Real getDet()
	{
		return t_.getDetJ();
	}

	Point<ndim> coorQuadPt(UInt iq)
	{
		return Point<ndim>(t_.getM_J()*Eigen::Map<const Eigen::Matrix<Real,mydim,1> >(&Integrator::NODES[iq][0])) + t_[0];
	}

	UInt getGlobalIndex(UInt iq)
	{
		return Integrator::NNODES * t_.getId() + iq;
	}

	//Returns \hat{phi}
	Real phiMaster(UInt i, UInt iq) const;

	//Returns \nabla \hat{phi}
	Real phiDerMaster(UInt i, UInt ic, UInt iq) const;

	//Returns J^{-1} \nabla \hat{phi}
	Real invTrJPhiDerMaster(UInt i, UInt ic, UInt iq) const;

};


// Partial specialization for 2.5D
// Could be removed adjusting the code
template <class Integrator ,UInt ORDER>
class FiniteElement<Integrator, ORDER, 2,3>{
private:
	static constexpr UInt NBASES=how_many_nodes(ORDER,2);
	Element<NBASES,2,3> t_;
	Eigen::Matrix<Real, NBASES, Integrator::NNODES> phiMapMaster_;
	//Numero basi locali x Num coordinate x numero nodi integrazione
	Eigen::Matrix<Real, NBASES, Integrator::NNODES*2> phiDerMapMaster_;
	//Eigen::Matrix<Real,3*ORDER, Integrator::NNODES*2> invTrJPhiDerMapMaster_;
	Eigen::Matrix<Real,2,2> metric_;

	void setPhiMaster();
	void setPhiDerMaster();
	//void setInvTrJPhiDerMaster();

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
	void updateElement(const Element<NBASES,2,3> &t);

	constexpr Real getAreaReference()
	{
		return 1./2;
	}

	Real getDet()
	{
		return t_.getDetJ();
	}

	Point<3> coorQuadPt(UInt iq)
	{
		return Point<3>(t_.getM_J()*Eigen::Map<const Eigen::Matrix<Real,3,1> >(&Integrator::NODES[iq][0])) + t_[0];
	}

	UInt getGlobalIndex(UInt iq)
	{
		return Integrator::NNODES * t_.getId() + iq;
	}

	//Returns \hat{phi}
	Real phiMaster(UInt i, UInt iq) const;

	//Returns \nabla \hat{phi}
	Real phiDerMaster(UInt i, UInt ic, UInt iq) const;

	//Returns J^{-1} \nabla \hat{phi}
	//Real invTrJPhiDerMaster(UInt i, UInt ic, UInt iq) const;

	Eigen::Matrix<Real,2,2> metric()const {return metric_;};

};

#include "finite_element_imp.h"

#endif
