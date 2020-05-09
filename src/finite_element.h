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
	//Num local bases x num coords x num quad nodes
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

	const Point<ndim>& operator[](UInt i) const {return t_[i];}
	// Both Area and Volume getter added for convenience
	Real getMeasure() const {t_.getMeasure();}
	Real getArea() const {t_.getMeasure();}
	Real getVolume() const {t_.getMeasure();}

	Real getId() const {return t_.getId();}

	Point<ndim> coorQuadPt(UInt iq){
		return Point<ndim>(t_.getM_J()*Eigen::Map<const Eigen::Matrix<Real,mydim,1> >(&Integrator::NODES[iq][0]))+=t_[0];
	}

	//Returns \hat{phi}
	Real phiMaster(UInt i, UInt iq) const {return phiMapMaster_(i, iq);}
	UInt getGlobalIndex(UInt iq) {return Integrator::NNODES * t_.getId() + iq;}

	// These function are used in matrix_assembler.h
	Real stiff_impl(UInt i, UInt j, UInt iq);
	Real stiff_anys_impl(UInt i, UInt j, UInt iq);
	Real mass_impl(UInt i, UInt j, UInt iq);
	Real grad_impl(UInt i, UInt j, UInt iq);

};

#include "finite_element_imp.h"

#endif
