#ifndef MATRIX_ASSEMBLER_IMP_H_
#define MATRIX_ASSEMBLER_IMP_H_


template<UInt ORDER, typename Integrator, UInt mydim, UInt ndim, typename A>
void Assembler::operKernel(EOExpr<A> oper, const MeshHandler<ORDER,mydim,ndim>& mesh,
	                     FiniteElement<Integrator,ORDER,mydim,ndim>& fe, SpMat& OpMat)
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;
	static constexpr UInt NBASES = FiniteElement<Integrator,ORDER,mydim,ndim>::NBASES;
	using local_matr_t = typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t;

	std::vector<coeff> triplets;
	triplets.reserve(NBASES*NBASES*mesh.num_elements());

	std::vector<UInt> identifiers;
	identifiers.reserve(NBASES);

	local_matr_t loc_matr = local_matr_t::Zero();

  for(int t=0; t<mesh.num_elements(); ++t){

		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		for(int i=0; i<NBASES; ++i)
			identifiers.push_back(fe[i].id());

		for(int iq = 0; iq < Integrator::NNODES; ++iq)
				loc_matr += oper(fe, iq) * Integrator::WEIGHTS[iq];

		loc_matr *= fe.getMeasure();

		for (int j=0; j<NBASES; ++j)
			for (int i=0; i<NBASES; ++i)
				triplets.push_back(coeff(identifiers[i],identifiers[j],loc_matr(i,j)));

		identifiers.clear();
		loc_matr.setZero();
	}

  UInt nnodes = mesh.num_nodes();
  OpMat.resize(nnodes, nnodes);
	OpMat.setFromTriplets(triplets.begin(),triplets.end());
	OpMat.prune(tolerance);
}

template<UInt ORDER, typename Integrator, UInt mydim, UInt ndim>
void Assembler::forcingTerm(const MeshHandler<ORDER,mydim,ndim>& mesh,
	                     FiniteElement<Integrator, ORDER,mydim,ndim>& fe, const ForcingTerm& u, VectorXr& forcingTerm)
{
	static constexpr UInt NBASES = FiniteElement<Integrator,ORDER,mydim,ndim>::NBASES;

	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  for(int t=0; t<mesh.num_elements(); ++t){

		fe.updateElement(mesh.getElement(t));

		for(int i=0; i<NBASES; ++i){
			Real s=0;
			for(int iq = 0; iq < Integrator::NNODES; ++iq){
				UInt globalIndex = fe.getGlobalIndex(iq);
				s +=  fe.phiMaster(i,iq)* u(globalIndex) * Integrator::WEIGHTS[iq];//(*)
			}
			forcingTerm[fe[i].id()] += s * fe.getMeasure();
		}
	}
}

#endif
