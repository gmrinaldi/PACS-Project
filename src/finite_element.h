#ifndef __FINITE_ELEMENT_HPP__
#define __FINITE_ELEMENT_HPP__

#include "integration.h"

#include "mesh_objects.h"

template <UInt ORDER, UInt mydim, UInt ndim>
class FiniteElementData{
	static_assert((ORDER==1 || ORDER==2) &&
								(mydim==2 || mydim==3) &&
								 mydim <= ndim,
								 "ERROR! TRYING TO INSTANTIATE FINITE ELEMENT WITH WRONG NUMBER OF NODES AND/OR DIMENSIONS! See finite_element.h");
public:

	using Integrator = typename IntegratorHelper::Integrator<ORDER,mydim>;

	static constexpr UInt NBASES = how_many_nodes(ORDER,mydim);

	FiniteElementData();

	FiniteElementData(const FiniteElementData&) = delete;
	FiniteElementData(FiniteElementData&&) = delete;
	FiniteElementData& operator=(const FiniteElementData&) = delete;
	FiniteElementData& operator=(FiniteElementData&&) = delete;

	virtual ~FiniteElementData()=0;

	void updateElement(const Element<NBASES,mydim,ndim>& t);

	const Point<ndim>& operator[] (UInt i) const {return t_[i];}
	// Both Area and Volume getter added for convenience
	Real getMeasure() const {return t_.getMeasure();}
	Real getArea() const {return t_.getMeasure();}
	Real getVolume() const {return t_.getMeasure();}

	Real getId() const {return t_.getId();}

	Point<ndim> coorQuadPt(UInt iq){
		Point<ndim> quadPt{t_.getM_J() * Eigen::Map<const Eigen::Matrix<Real,mydim,1> >(&Integrator::NODES[iq][0])};
		return quadPt += t_[0];
	}

	//Returns \hat{phi}
	Real phiMaster(UInt i, UInt iq) const {return referencePhi(i, iq);}
	UInt getGlobalIndex(UInt iq) const {return Integrator::NNODES * t_.getId() + iq;}

protected:

	Element<NBASES,mydim,ndim> t_;
	Eigen::Matrix<Real, NBASES, Integrator::NNODES> referencePhi;
	//Num local bases x num coords x num quad nodes
	Eigen::Matrix<Real, NBASES, mydim*Integrator::NNODES> referencePhiDer;
	Eigen::Matrix<Real, NBASES, ndim*Integrator::NNODES> elementPhiDer;

	void setPhi();
	void setPhiDer();
	void setElementPhiDer();

};

template <UInt ORDER, UInt mydim, UInt ndim>
struct FiniteElement : public FiniteElementData<ORDER, mydim, ndim> {

	using Integrator = typename FiniteElementData<ORDER, mydim, ndim>::Integrator;

	static constexpr UInt NBASES = FiniteElementData<ORDER, mydim, ndim>::NBASES;

	FiniteElement()=default;

	auto stiff_impl(UInt iq) const -> decltype(this->elementPhiDer.template block<NBASES,ndim>(0,0) * this->elementPhiDer.template block<NBASES,ndim>(0,0).transpose()) {
		// compute nabla(phi) x nabla(phi) at node iq
		return this->elementPhiDer.template block<NBASES,ndim>(0,ndim*iq) * this->elementPhiDer.template block<NBASES,ndim>(0,ndim*iq).transpose();
	}

	auto stiff_impl(UInt iq, const Eigen::Matrix<Real, ndim, ndim>& K) const -> decltype(this->elementPhiDer.template block<NBASES,ndim>(0,0) * K * this->elementPhiDer.template block<NBASES,ndim>(0,0).transpose()) {
		// compute nabla(phi) x K nabmyla(phi) at node iq
		// Memo: K is symmetric!
		return this->elementPhiDer.template block<NBASES,ndim>(0,ndim*iq) * K * this->elementPhiDer.template block<NBASES,ndim>(0,ndim*iq).transpose();
	}

	auto mass_impl(UInt iq) const -> decltype(this->referencePhi.col(0) * this->referencePhi.col(0).transpose()) {
		// compute phi x phi at node iq
		return this->referencePhi.col(iq) * this->referencePhi.col(iq).transpose();
	}

	auto grad_impl(UInt iq) const -> decltype(this->referencePhi.col(0) * this->elementPhiDer.col(0).transpose()) {
		// compute phi x nabla(phi) x b at node iq
		// default is b=u_x
		return this->referencePhi.col(iq) * this->elementPhiDer.col(ndim*iq).transpose();
	}

	auto grad_impl(UInt iq, const Eigen::Matrix<Real, ndim, 1>& b) const -> decltype(this->referencePhi.col(0) * b.transpose() * this->elementPhiDer.template block<NBASES,ndim>(0,0).transpose()) {
		// compute phi x nabla(phi) x b at node iq
		return this->referencePhi.col(iq) * b.transpose() * this->elementPhiDer.template block<NBASES,ndim>(0,ndim*iq).transpose();
	}

};

#include "finite_element_imp.h"

#endif
