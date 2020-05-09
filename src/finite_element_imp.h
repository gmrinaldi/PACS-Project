#ifndef __FINITE_ELEMENT_IMP_HPP__
#define __FINITE_ELEMENT_IMP_HPP__

#include<iostream>

// Template auxiliary function declaration
template<UInt NBASES, UInt mydim>
Eigen::Matrix<Real, NBASES, 1> reference_eval_point(const Point<mydim> &node);

template<UInt NBASES, UInt mydim>
Eigen::Matrix<Real, NBASES,mydim> reference_eval_der_point(const Point<mydim> &node);

//! FE implementation
template <class Integrator, UInt ORDER, UInt mydim, UInt ndim>
FiniteElement<Integrator, ORDER, mydim, ndim>::FiniteElement()
{
	//How it will be used, it does not depend on J^-1 -> set one time
	setPhiMaster();
	setPhiDerMaster();
}

template <class Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FiniteElement<Integrator, ORDER, mydim, ndim>::updateElement(const Element<NBASES,mydim,ndim> &t)
{
	t_ = t;

	//it does depend on J^-1 -> set for each element
	setInvTrJPhiDerMaster();

}

template <class Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FiniteElement<Integrator, ORDER, mydim, ndim>::setPhiMaster()
{
	for(UInt iq=0; iq<Integrator::NNODES; ++iq)
		phiMapMaster_.col(iq)=reference_eval_point<NBASES,mydim>(Integrator::NODES[iq]);
}

template <class Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FiniteElement<Integrator, ORDER, mydim, ndim>::setPhiDerMaster()
{
	// .template is needed! See Eigen documentation regarding
	// the template and typename keywords in C++
	for(UInt iq=0; iq<Integrator::NNODES; ++iq)
		phiDerMapMaster_.template block<NBASES, mydim>(0, mydim*iq)=reference_eval_der_point<NBASES,mydim>(Integrator::NODES[iq]);
}

template <class Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FiniteElement<Integrator, ORDER, mydim, ndim>::setInvTrJPhiDerMaster()
{
	// we need J^(-T) nabla( phi)
	for (auto iq=0; iq < Integrator::NNODES; ++iq)
			invTrJPhiDerMapMaster_.template block<NBASES,ndim>(0,ndim*iq).noalias() = phiDerMapMaster_.template block<NBASES,mydim>(0,mydim*iq)*t_.getM_invJ();
}

template <class Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real FiniteElement<Integrator, ORDER, mydim, ndim>::stiff_impl(UInt i, UInt j, UInt iq){
	// compute nabla(phi_i) dot nabla(phi_j) at node iq
	return invTrJPhiDerMapMaster_.template block<1,ndim>(i,ndim*iq).dot(invTrJPhiDerMapMaster_.template block<1,ndim>(j,ndim*iq));
}

template <class Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real FiniteElement<Integrator, ORDER, mydim, ndim>::stiff_anys_impl(UInt i, UInt j, UInt iq, const Eigen::Matrix<Real,2,2>& K){
	// compute nabla(phi_i) dot K nabla(phi_j) at node iq
	// Memo: K is symmetric!
	return invTrJPhiDerMapMaster_.template block<1,ndim>(i,ndim*iq).transpose()*K*invTrJPhiDerMapMaster_.template block<1,ndim>(j,ndim*iq);
}

template <class Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real FiniteElement<Integrator, ORDER, mydim, ndim>::mass_impl(UInt i, UInt j, UInt iq){
	// compute phi_i x phi_j at node iq
	return phiMapMaster_(i,iq)*phiMapMaster_(j,iq);
}

template <class Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real FiniteElement<Integrator, ORDER, mydim, ndim>::grad_impl(UInt i, UInt j, UInt iq){
	// compute phi_i x nabla(phi_j) at node iq
	return phiMapMaster_(i,iq)*invTrJPhiDerMapMaster_.template block<1,ndim>(j,ndim*iq);
}


// Templates for auxiliary functions
// This function evaluates the basis function on the reference element
// at the quadrature nodes
// This function covers all order 1 cases
template<UInt NBASES, UInt mydim>
Eigen::Matrix<Real, NBASES, 1> reference_eval_point(const Point<mydim> &node){
	Eigen::Matrix<Real, NBASES, 1> phi;
	phi.template tail<mydim>()=Eigen::Map<const Eigen::Matrix<Real,mydim,1> >(&node[0]);
	phi(0) = 1 - phi.template tail<mydim>().sum();
	return phi;
}

// Full specialization for order 2 in 2D and 2.5D (same reference element!)
template<>
Eigen::Matrix<Real, 6, 1> reference_eval_point<6,2>(const Point<2> &node){
	Eigen::Matrix<Real, 6, 1> phi;
	phi << (1-node[0]-node[1]) * (1-2*node[0]-2*node[1]),
										 node[0] * (2*node[0]-1),
										 node[1] * (2*node[1]-1),
									 4*node[0] * node[1],
								   4*node[1] * (1-node[0]-node[1]),
									 4*node[0] * (1-node[0]-node[1]);
	return phi;
}

// Full specialization for order 2 in 3D
template<>
Eigen::Matrix<Real, 10, 1> reference_eval_point<10,3>(const Point<3> &node){
	Eigen::Matrix<Real, 10, 1> phi;
	phi << (1-node[0]-node[1]-node[2]) * (1-2*node[0]-2*node[1]-2*node[2]),
														 node[0] * (2*node[0]-1),
														 node[1] * (2*node[1]-1),
														 node[2] * (2*node[2]-1),
													 4*node[0] * (1-node[0]-node[1]-node[2]),
													 4*node[1] * (1-node[0]-node[1]-node[2]),
													 4*node[2] * (1-node[0]-node[1]-node[2]),
													 4*node[0] * node[1],
													 4*node[1] * node[2],
													 4*node[2] * node[0];
	return phi;
}

// This function evaluates the ndim-gradient of basis function on the reference element
// at the quadrature nodes
// This function covers all order 1 cases
template<UInt NBASES, UInt mydim>
Eigen::Matrix<Real, NBASES,mydim> reference_eval_der_point(const Point<mydim> &node){
	Eigen::Matrix<Real,NBASES,mydim> B1;
	B1.row(0).setConstant(-1);
	B1.bottomRows(mydim).setIdentity();
	return B1;
}

// Full specialization for order 2 in 2D and 2.5D
template<>
Eigen::Matrix<Real, 6,2> reference_eval_der_point<6,2>(const Point<2> &node){
	Eigen::Matrix<Real,6,2> B2;
	B2 << 1-4*(1-node[0]-node[1]), 1-4*(1-node[0]-node[1]),
										4*node[0]-1, 											 0,
															0, 						 4*node[1]-1,
											4*node[1],  						 4*node[0],
										 -4*node[1], 4*(1-node[0]-2*node[1]),
				4*(1-2*node[0]-node[1]), 						  -4*node[0];
	return B2;
}

// Full specialization for order 2 in 3D
template<>
Eigen::Matrix<Real, 10,3> reference_eval_der_point<10,3>(const Point<3> &node){
	Eigen::Matrix<Real,10,3> B2;
	B2 << 1-4*(1-node[0]-node[1]-node[2]), 	1-4*(1-node[0]-node[1]-node[2]), 1-4*(1-node[0]-node[1]-node[2]),
														4*node[0]-1, 		 			 						          0,															 0,
																			0, 										  4*node[1]-1, 															 0,
				4*(1-2*node[0]-node[1]-node[2]),  						  			 -4*node[0], 										  -4*node[0],
														 -4*node[1],  4*(1-node[0]-2*node[1]-node[2]), 											-4*node[1],
														 -4*node[2], 											 -4*node[2], 4*(1-node[0]-node[1]-2*node[2]),
															4*node[1],												4*node[0], 															 0,
																			0, 												4*node[2], 											 4*node[1],
															4*node[2],																0,											 4*node[0];
	return B2;
}



#endif
