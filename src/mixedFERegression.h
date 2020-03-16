#ifndef __MIXEDFEREGRESSION_HPP__
#define __MIXEDFEREGRESSION_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "regressionData.h"
#include "solver.h"
#include "integratePsi.h"
#include <memory>

/*! A base class for the smooth regression. 
*/
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegressionBase
{
	protected:
		
	const MeshHandler<ORDER, mydim, ndim> &mesh_;
	const InputHandler& regressionData_;
	
	//  system matrix= 	|psi^T * A *psi | lambda R1^T  |   +  |psi^T * A * (-H) * psi |  O |   =  matrixNoCov + matrixOnlyCov
	//	                |     R1        | R0	      |      |         O             |  O |
	
	SpMat matrixNoCov_;	//! System matrix with psi^T*psi or psi^T*A*psi in north-west block  (is the full system matrix if no covariates)
	//SpMat matrixOnlyCov_; //! coeffmatrix=matrixNoCov+matrixOnlyCov
	SpMat R1_;	//! North-east block of system matrix matrixNoCov_
	SpMat R0_;	//! South-east block of system matrix matrixNoCov_
	SpMat psi_; //! Psi matrix of the model
	MatrixXr U_;	//! psi^T * W or psi^T * A * W padded with zeros, needed for Woodbury decomposition
	MatrixXr V_;   //! W^T*psi, if pointwise data is U^T, needed for Woodbury decomposition
	VectorXr z_; //! Observations
		
		
	Eigen::SparseLU<SpMat> matrixNoCovdec_; // Stores the factorization of matrixNoCov_
	Eigen::PartialPivLU<MatrixXr> Gdec_;	// Stores factorization of G =  C + [V * matrixNoCov^-1 * U]
	Eigen::PartialPivLU<MatrixXr> WTW_;	// Stores the factorization of W^T * W
	bool isWTWfactorized_ = false;
	bool isRcomputed_ = false;
	MatrixXr R_; //R1 ^T * R0^-1 * R1
	
	MatrixXr Q_;  //! Identity - H, projects onto the orthogonal subspace
 	MatrixXr H_; //! The hat matrix of the regression
	
	VectorXr A_; //A_.asDiagonal() = diag(|D_1|,...,|D_N|) areal matrix, = identity nnodesxnnodes if pointwise data

	VectorXr _rightHandSide;                     //!A Eigen::VectorXr: Stores the system right hand side.
	std::vector<VectorXr> _solution; //!A Eigen::VectorXr: Stores the system solution.
	std::vector<Real> _dof; //! A vector storing the computed dofs
	
	bool isSpaceVarying=false; // used to distinguish whether to use the forcing term u in apply() or not
	
	//! A member function computing the Psi matrix
	void setPsi();
	//! A member function computing the no-covariates version of the system matrix 
	void buildMatrixNoCov(const SpMat& Psi,  const SpMat& R1,  const SpMat& R0);
	// A member function computing the matrix to be added to matrixNoCov_ to obtain the full system matrix
	//void buildMatrixOnlyCov(const SpMat& Psi, const MatrixXr& H);
	//! A function that given a vector u, performs Q*u efficiently
	MatrixXr LeftMultiplybyQ(const MatrixXr& u);
	//! A function which adds Dirichlet boundary conditions before solving the system ( Remark: BC for areal data are not implemented!)
	void addDirichletBC();
	//! A member function which builds the Q matrix
	void setQ();
	//! A member function which builds the H matrix
 	void setH();
 	//! A member function which builds the A vector containing the areas of the regions in case of areal data
	void setA(); 
	
	//! A member function returning the system right hand data
	void getRightHandData(VectorXr& rightHandData);
	
	//! A member function computing the dofs
	void computeDegreesOfFreedom(UInt output_index, Real lambda);
	//! A function computing dofs in case of exact GCV, it is called by computeDegreesOfFreedom
	void computeDegreesOfFreedomExact(UInt output_index, Real lambda);
	//! A function computing dofs in case of stochastic GCV, it is called by computeDegreesOfFreedom
	void computeDegreesOfFreedomStochastic(UInt output_index, Real lambda);
    
    //! A function to factorize the system, using Woodbury decomposition when there are covariates
	void system_factorize();
	//! A function which solves the factorized system 
	template<typename Derived>
	MatrixXr system_solve(const Eigen::MatrixBase<Derived>&);
	
	public:
	//!A Constructor.
	MixedFERegressionBase(const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& regressionData): mesh_(mesh), regressionData_(regressionData) {};

	//! The function solving the system, used by the children classes. Saves the result in _solution
	/*!
	    \param oper an operator, which is the Stiffness operator in case of Laplacian regularization
	    \param u the forcing term, will be used only in case of anysotropic nonstationary regression
	*/
	template<typename A>
	void apply(EOExpr<A> oper,const ForcingTerm & u);
	
	//! A inline member that returns a VectorXr, returns the whole solution_.
	inline std::vector<VectorXr> const & getSolution() const{return _solution;};
	//! A function returning the computed dofs of the model
	inline std::vector<Real> const & getDOF() const{return _dof;};
};

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression : public MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, ndim, mydim>& mesh, const InputHandler& regressionData):MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
		std::cout << "Option not implemented! \n";
	}
};

#include "mixedFERegression_imp.h"

#endif
