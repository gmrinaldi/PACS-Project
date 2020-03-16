#ifndef __MIXEDFEREGRESSION_IMP_HPP__
#define __MIXEDFEREGRESSION_IMP_HPP__

#include <iostream>
#include <chrono>
#include <random>
#include <fstream>

#include "R_ext/Print.h"

//#include <libseq/mpi.h>
#include "../inst/include/dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::addDirichletBC()
{
	UInt id1,id3;

	UInt nnodes = mesh_.num_nodes();

	const std::vector<UInt>& bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real>& bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices.size();

	Real pen=10e20;

	for( auto i=0; i<nbc_indices; i++)
	 {
			id1=bc_indices[i];
			id3=id1+nnodes;

			matrixNoCov_.coeffRef(id1,id1)=pen;
			matrixNoCov_.coeffRef(id3,id3)=pen;


			_rightHandSide(id1)=bc_values[i]*pen;
			_rightHandSide(id3)=0;
	 }

	matrixNoCov_.makeCompressed();
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::setPsi(){

	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	psi_.resize(nlocations, nnodes);
	if (regressionData_.isLocationsByNodes()) //pointwise data
	{
		std::vector<coeff> tripletAll;
		auto k = regressionData_.getObservationsIndices();
		tripletAll.reserve(k.size());
		for (int i = 0; i< k.size(); ++i){
			tripletAll.push_back(coeff(i,k[i],1.0));
		}
		psi_.setFromTriplets(tripletAll.begin(),tripletAll.end());
		psi_.makeCompressed();
	}
	else if (regressionData_.getNumberOfRegions() == 0)
	{
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> tri_activated;
		Eigen::Matrix<Real,Nodes,1> coefficients;

		Real evaluator;

		for(UInt i=0; i<nlocations;i++)
		{
			tri_activated = mesh_.findLocationNaive(regressionData_.getLocations()[i]);
			if(tri_activated.getId() == Identifier::NVAL)
			{
				#ifdef R_VERSION_
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
				#else
				std::cout << "ERROR: Point " << i+1 <<" is not in the domain\n";
				#endif
			}else
			{
				for(UInt node = 0; node < Nodes ; ++node)
				{
					coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
					coefficients(node) = 1; //Activates only current base
					evaluator = evaluate_point<Nodes,mydim,ndim>(tri_activated, regressionData_.getLocations()[i], coefficients);
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		}
		psi_.makeCompressed();
	}
	else //areal data
	{
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;

		Real *tab; //Psi_i
		tab = (Real*) malloc(sizeof(Real)*nnodes);
		for(UInt i=0; i<nlocations;i++) //nlocations = number of regions
		{
			for (UInt k=0; k<nnodes; k++) {tab[k]=0;}
			for (UInt j=0; j<mesh_.num_elements(); j++)
			{
				if (regressionData_.getIncidenceMatrix()(i,j) == 1) //element j is in region i
				{
					Element<Nodes, mydim, ndim> tri = mesh_.getElement(j);
					for (UInt k=0; k<Nodes; k++)
					{
						tab[tri[k].getId()] += integratePsi(tri,k); // integral over tri of psi_k
					}
				}
			}
			for (int k=0; k<nnodes; k++)
			{
				if (tab[k] != 0)
				{
					psi_.insert(i,k) = tab[k]/A_(i); //divide by |D_i|
				}
			}
		}
		free(tab);
		psi_.makeCompressed();
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
MatrixXr MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::LeftMultiplybyQ(const MatrixXr& u)
{	
	if (regressionData_.getCovariates().rows() == 0){
		return u;
	}
	else{
		MatrixXr W(this->regressionData_.getCovariates());
		if (isWTWfactorized_ == false ){
			WTW_.compute(W.transpose()*W);
			isWTWfactorized_=true;
		}
		MatrixXr Pu= W*WTW_.solve(W.transpose()*u);
		return u-Pu;
	}

}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::buildMatrixNoCov(const SpMat& Psi,  const SpMat& R1,  const SpMat& R0) {

	UInt nnodes = mesh_.num_nodes();
	
	SpMat DMat;
	
	if(regressionData_.getNumberOfRegions()==0) // pointwise data
	    DMat=Psi.transpose()*Psi;
	else                                        // areal data: need to add the diag(|D_1|,...,|D_N|)
	    DMat=Psi.transpose()*A_.asDiagonal()*Psi;    
	

	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*R1.nonZeros() + R0.nonZeros());

	for (int k=0; k<DMat.outerSize(); ++k)
		for (SpMat::InnerIterator it(DMat,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
		}
	for (int k=0; k<R0.outerSize(); ++k)
		for (SpMat::InnerIterator it(R0,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
		}
	for (int k=0; k<R1.outerSize(); ++k)
	  for (SpMat::InnerIterator it(R1,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
	  }
	for (int k=0; k<R1.outerSize(); ++k)
	  for (SpMat::InnerIterator it(R1,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
	  }

	matrixNoCov_.setZero();
	matrixNoCov_.resize(2*nnodes,2*nnodes);
	matrixNoCov_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	matrixNoCov_.makeCompressed();
	//std::cout<<"Coefficients' Matrix Set Correctly"<<std::endl;
}

//template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
//void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::buildMatrixOnlyCov(const SpMat& Psi,  const MatrixXr& H) {
//
//	UInt nnodes=mesh_.num_nodes();
//	
//	MatrixXr NWblock= MatrixXr::Zero(2*nnodes,2*nnodes);
//	
//	if(regressionData_.getNumberOfRegions()==0) 
//    	NWblock.topLeftCorner(nnodes,nnodes)=Psi.transpose()*(-H)*Psi;
//    else
//	    NWblock.topLeftCorner(nnodes,nnodes)=Psi.transpose()*A_.asDiagonal()*(-H)*Psi;	
//
////    matrixOnlyCov_.setZero();
////    matrixOnlyCov_.resize(2*nnodes,2*nnodes);
//    
//    matrixOnlyCov_=NWblock.sparseView();
////    for(int i=0; i<nnodes; ++i)
////     for(int j=0; j<nnodes; ++j)
////       matrixOnlyCov_.coeffRef(i,j)=NWblock(i,j);
//
//	matrixOnlyCov_.makeCompressed();
//}
//

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::system_factorize() {

	UInt nnodes = mesh_.num_nodes();

	// First phase: Factorization of matrixNoCov
	matrixNoCovdec_.compute(matrixNoCov_);
	
	if (regressionData_.getCovariates().rows() != 0) {
		// Second phase: factorization of matrix  G =  C + [V * matrixNoCov^-1 * U]= C + D

		// Definition of matrix U = [ psi^T * A * W | 0 ]^T and V= [ W^T*psi| 0]
		
		MatrixXr W(this->regressionData_.getCovariates());
		
		U_ = MatrixXr::Zero(2*nnodes, W.cols());
		
		V_=MatrixXr::Zero(W.cols(),2*nnodes);
		V_.leftCols(nnodes)=W.transpose()*psi_;
			
		if(regressionData_.getNumberOfRegions()==0){ // pointwise data
		  U_.topRows(nnodes) = psi_.transpose()*W;
		} 
		else{                                          //areal data
		  U_.topRows(nnodes) = psi_.transpose()*A_.asDiagonal()*W;
        }
		
		MatrixXr D = V_*matrixNoCovdec_.solve(U_);
		
		// G = C + D
		MatrixXr G = -W.transpose()*W + D;
		Gdec_.compute(G);
		
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename Derived>
MatrixXr MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::system_solve(const Eigen::MatrixBase<Derived> &b) {

	// Resolution of the system matrixNoCov * x1 = b
	MatrixXr x1 = matrixNoCovdec_.solve(b);
	
	if (regressionData_.getCovariates().rows() != 0) {
		// Resolution of G * x2 = V * x1
		
		MatrixXr x2 = Gdec_.solve(V_*x1); 
	
		// Resolution of the system matrixNoCov * x3 = U * x2
		x1 -= matrixNoCovdec_.solve(U_*x2);
	}
	
	return x1;
}

 template<typename InputHandler, typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::setQ()
 {
 	//std::cout<<"Computing Orthogonal Space Projection Matrix"<<std::endl;
 	Q_.resize(H_.rows(),H_.cols());
 	Q_ = -H_;
 	for (int i=0; i<H_.rows();++i)
 	{
 		Q_(i,i) += 1;
 	}
 }
 
 template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim, ndim>::setH()
 {

 	UInt nlocations = regressionData_.getNumberofObservations();
 

 	MatrixXr W(this->regressionData_.getCovariates());

 	if(regressionData_.isLocationsByNodes())
 	{
 		MatrixXr W_reduced(regressionData_.getNumberofObservations(), W.cols());
 		for (auto i=0; i<nlocations;++i)
 		{
 			auto index_i = regressionData_.getObservationsIndices()[i];
 			for (auto j=0; j<W.cols();++j)
 			{
 				W_reduced(i,j) = W(index_i,j);
 			}
 		}
 		W = W_reduced;
 	}
 
 
 	MatrixXr WTW(W.transpose()*W);
 
 	H_=W*WTW.ldlt().solve(W.transpose()); // using cholesky LDLT decomposition for computing hat matrix
 }

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim, ndim>::setA()
{
	UInt nRegions = regressionData_.getNumberOfRegions();
	A_.resize(nRegions,1);
	for (int i=0; i<nRegions; i++)
	{
		A_(i)=0;
		for (int j=0; j<regressionData_.getIncidenceMatrix().cols(); j++)
		{
			if (regressionData_.getIncidenceMatrix()(i,j) == 1)
			{
				A_(i)+=mesh_.elementMeasure(j);
			}
		}
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::getRightHandData(VectorXr& rightHandData)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	rightHandData = VectorXr::Zero(nnodes);

	if (regressionData_.getCovariates().rows() == 0) //no covariate
	{
		if (regressionData_.isLocationsByNodes())
		{
			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				rightHandData(index_i) = regressionData_.getObservations()[i];
			}
		}
		else if (regressionData_.getNumberOfRegions() == 0) //pointwise data
		{
			rightHandData=psi_.transpose()*regressionData_.getObservations();
		}
		else //areal data
		{
			rightHandData=psi_.transpose()*A_.asDiagonal()*regressionData_.getObservations();
		}
	}
	else if (regressionData_.getNumberOfRegions() == 0) //with covariates, pointwise data
	{
		rightHandData=psi_.transpose()*LeftMultiplybyQ(regressionData_.getObservations());
	}
	else //with covariates, areal data
	{
		rightHandData=psi_.transpose()*A_.asDiagonal()*LeftMultiplybyQ(regressionData_.getObservations());
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::computeDegreesOfFreedom(UInt output_index, Real lambda)
{
	int GCVmethod = regressionData_.getGCVmethod();
	switch (GCVmethod) {
		case 1:
			computeDegreesOfFreedomExact(output_index, lambda);
			break;
		case 2:
			computeDegreesOfFreedomStochastic(output_index, lambda);
			break;
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::computeDegreesOfFreedomExact(UInt output_index, Real lambda)
{
    
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	Real degrees=0;

	// Case 1: MUMPS
	if (regressionData_.isLocationsByNodes() && regressionData_.getCovariates().rows() == 0 )
	{
		auto k = regressionData_.getObservationsIndices();
		DMUMPS_STRUC_C id;
		//int myid, ierr;
        //int argc=0;
        //char ** argv= NULL;
        //MPI_Init(&argc,&argv);
		//ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

		id.sym=0;
		id.par=1;
		id.job=JOB_INIT;
		id.comm_fortran=USE_COMM_WORLD;
		dmumps_c(&id);

		std::vector<int> irn;
		std::vector<int> jcn;
		std::vector<double> a;
		std::vector<int> irhs_ptr;
		std::vector<int> irhs_sparse;
		double* rhs_sparse= (double*)malloc(nlocations*sizeof(double));
		
		//if( myid==0){
			id.n=2*nnodes;
			for (int j=0; j<matrixNoCov_.outerSize(); ++j){
				for (SpMat::InnerIterator it(matrixNoCov_,j); it; ++it){
					irn.push_back(it.row()+1);
					jcn.push_back(it.col()+1);
					a.push_back(it.value());
				}
			}
		//}
		id.nz=irn.size();
		id.irn=irn.data();
		id.jcn=jcn.data();
		id.a=a.data();
		id.nz_rhs=nlocations;
		id.nrhs=2*nnodes;
		int j = 1;
		irhs_ptr.push_back(j);
		for (int l=0; l<k[0]-1; ++l) {
			irhs_ptr.push_back(j);
		}
		for (int i=0; i<k.size()-1; ++i) {
			++j;
			for (int l=0; l<k[i+1]-k[i]; ++l) {
				irhs_ptr.push_back(j);
			}
			
		}
		++j;
		for (int i=k[k.size()-1]; i < id.nrhs; ++i) {
			irhs_ptr.push_back(j);
		}
		for (int i=0; i<nlocations; ++i){
			irhs_sparse.push_back(k[i]+1);
		}
		id.irhs_sparse=irhs_sparse.data();
		id.irhs_ptr=irhs_ptr.data();
		id.rhs_sparse=rhs_sparse;

		#define ICNTL(I) icntl[(I)-1]
		//Output messages suppressed
		id.ICNTL(1)=-1;
		id.ICNTL(2)=-1;
		id.ICNTL(3)=-1;
		id.ICNTL(4)=0;
		id.ICNTL(20)=1;
		id.ICNTL(30)=1;
		id.ICNTL(14)=200;

		id.job=6;
		dmumps_c(&id);
		id.job=JOB_END;
		dmumps_c(&id);

		//if (myid==0){
			for (int i=0; i< nlocations; ++i){
				//std::cout << "rhs_sparse" << rhs_sparse[i] << std::endl;
				degrees+=rhs_sparse[i];
			}
		//}
		free(rhs_sparse);

		//MPI_Finalize();
	}
	// Case 2: Eigen
	else{
		MatrixXr X1;
		if (regressionData_.getNumberOfRegions() == 0){ //pointwise data
			X1 = psi_.transpose() * LeftMultiplybyQ(psi_);
		}else{ //areal data
			X1 = psi_.transpose() * A_.asDiagonal() * LeftMultiplybyQ(psi_);
		}

		if (isRcomputed_ == false){
			isRcomputed_ = true;
			Eigen::SparseLU<SpMat> solver;
			solver.compute(R0_);
			auto X2 = solver.solve(R1_);
			R_ = R1_.transpose() * X2;
		}

		MatrixXr X3 = X1 + lambda * R_;
		Eigen::LDLT<MatrixXr> Dsolver(X3);

		auto k = regressionData_.getObservationsIndices();

		if(regressionData_.isLocationsByNodes() && regressionData_.getCovariates().rows() != 0) {
			degrees += regressionData_.getCovariates().cols();

			// Setup rhs B
			MatrixXr B;
			B = MatrixXr::Zero(nnodes,nlocations);
			// B = I(:,k) * Q
			for (auto i=0; i<nlocations;++i) {
				VectorXr ei = VectorXr::Zero(nlocations);
				ei(i) = 1;
				VectorXr Qi = LeftMultiplybyQ(ei);
				for (int j=0; j<nlocations; ++j) {
					B(k[i], j) = Qi(j);
				}
			}
			// Solve the system TX = B
			MatrixXr X;
			X = Dsolver.solve(B);
			// Compute trace(X(k,:))
			for (int i = 0; i < k.size(); ++i) {
				degrees += X(k[i], i);
			}
		}

		if (!regressionData_.isLocationsByNodes()){
			MatrixXr X;
			X = Dsolver.solve(MatrixXr(X1));

			if (regressionData_.getCovariates().rows() != 0) {
				degrees += regressionData_.getCovariates().cols();
			}
			for (int i = 0; i<nnodes; ++i) {
				degrees += X(i,i);
			}
		}
	}
	_dof[output_index] = degrees;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::computeDegreesOfFreedomStochastic(UInt output_index, Real lambda)
{	
	
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	std::default_random_engine generator;
	// Creation of the random matrix
	std::bernoulli_distribution distribution(0.5);
	UInt nrealizations = regressionData_.getNrealizations();
	MatrixXr u(nlocations, nrealizations);
	for (int j=0; j<nrealizations; ++j) {
		for (int i=0; i<nlocations; ++i) {
			if (distribution(generator)) {
				u(i,j) = 1.0;
			}
			else {
				u(i,j) = -1.0;
			}
		}
	}

	// Define the first right hand side : | I  0 |^T * psi^T * A * Q * u
	MatrixXr b = MatrixXr::Zero(2*nnodes,u.cols());
	if (regressionData_.getNumberOfRegions() == 0){
		b.topRows(nnodes) = psi_.transpose() * LeftMultiplybyQ(u);
	}else{
		b.topRows(nnodes) = psi_.transpose() * A_.asDiagonal() * LeftMultiplybyQ(u);
	}

	// Resolution of the system
	MatrixXr x = system_solve(b);
	
//	Eigen::SparseLU<SpMat> solver;
//	
//	if(!regressionData_.getCovariates().rows()==0){
//	
//		this->buildMatrixOnlyCov(psi_, H_);
//		
//		SpMat coeffMatrix_= matrixNoCov_ + matrixOnlyCov_;
//
//	    solver.compute(coeffMatrix_); //matrixNoCov_+matrixOnlyCov_ = full system matrix when there are covariates
//
//	}
//	else{
//		solver.compute(matrixNoCov_);
//	}
//	
//	auto x = solver.solve(b);
//	if(solver.info()!=Eigen::Success)
//		{
//			#ifdef R_VERSION_
//	        Rprintf("Solving system for stoch gcv failed!!!\n");
//	        #endif
//		}

	MatrixXr uTpsi = u.transpose()*psi_;
	VectorXr edf_vect(nrealizations);
	Real q = 0;

	// Degrees of freedom = q + E[ u^T * psi * | I  0 |* x ]
	if (regressionData_.getCovariates().rows() != 0) {
		q = regressionData_.getCovariates().cols();
	}
	// For any realization we compute the degrees of freedom
	for (int i=0; i<nrealizations; ++i) {
		
		edf_vect(i) = uTpsi.row(i).dot(x.col(i).head(nnodes)) + q;
	}

	// Estimates: sample mean, sample variance
	Real mean = edf_vect.sum()/nrealizations;
	_dof[output_index] = mean;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename A>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::apply(EOExpr<A> oper, const ForcingTerm & u)
{
	UInt nnodes=mesh_.num_nodes();
	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	setA();
	setPsi();
	
	if(!regressionData_.getCovariates().rows() == 0)
	{
		setH();
		setQ();
	}

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, mesh_, fe, R1_);
	Assembler::operKernel(mass, mesh_, fe, R0_);
	VectorXr rightHandData;
	getRightHandData(rightHandData); //updated
	this->_rightHandSide = VectorXr::Zero(2*nnodes);
	this->_rightHandSide.topRows(nnodes)=rightHandData;
	
	VectorXr forcingTerm;
	
	if(this->isSpaceVarying)
	{
        //u=ForcingTerm(VectorXr::Zero(nnodes));
		Assembler::forcingTerm(mesh_, fe, u, forcingTerm);
	
	}
	
	this->_solution.resize(regressionData_.getLambda().size());
	this->_dof.resize(regressionData_.getLambda().size());

	for(UInt i = 0; i<regressionData_.getLambda().size(); ++i)
	{
		Real lambda = regressionData_.getLambda()[i];
		SpMat R1_lambda = (-lambda)*R1_;
		SpMat R0_lambda = (-lambda)*R0_;
		
		this->buildMatrixNoCov(psi_, R1_lambda, R0_lambda);
		//this->buildCoeffMatrix(DMat_, R1_lambda, R0_lambda);
				
		if(this->isSpaceVarying)
		{
		    _rightHandSide.bottomRows(nnodes)=lambda*forcingTerm;	
		}
		
		//Applying boundary conditions if necessary
		if(regressionData_.getDirichletIndices().size() != 0)  // if areal data NO BOUNDARY CONDITIONS
			addDirichletBC();
	
		system_factorize();	

	    _solution[i] = this->template system_solve(this->_rightHandSide);

		if(regressionData_.computeDOF())
		{
			computeDegreesOfFreedom(i,lambda);
		}
			
		else
			_dof[i] = -1;
	}
}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionData, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionData& regressionData):MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	    MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>::apply(stiff, ForcingTerm(std::vector<Real>(1)));
	}
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataElliptic, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionDataElliptic& regressionData):MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
	if(mydim!=2 || ndim !=2){

	#ifdef R_VERSION_
		Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
	#else
		std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
	#endif

	}else{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	    const Real& c = this->regressionData_.getC();
	    const Eigen::Matrix<Real,2,2>& K = this->regressionData_.getK();
	    const Eigen::Matrix<Real,2,1>& b = this->regressionData_.getBeta();

	    MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>::apply(c*mass+stiff[K]+dot(b,grad), ForcingTerm(std::vector<Real>(1)));
	}
	}
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionDataEllipticSpaceVarying& regressionData):MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
	if(mydim!=2 || ndim !=2){

	#ifdef R_VERSION_
		Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
	#else
		std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
	#endif

	}else{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

		const Reaction& c = this->regressionData_.getC();
		const Diffusivity& K = this->regressionData_.getK();
		const Advection& b = this->regressionData_.getBeta();
		const ForcingTerm& u= this->regressionData_.getU();
		
		this->isSpaceVarying=TRUE;

		MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>::apply(c*mass+stiff[K]+dot(b,grad), u);
	}
	}
};

#endif

