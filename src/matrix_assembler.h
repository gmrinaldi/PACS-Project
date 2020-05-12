#ifndef MATRIX_ASSEMBLER_H_
#define MATRIX_ASSEMBLER_H_


#include "fdaPDE.h"
#include "finite_element.h"
#include "mesh_objects.h"
#include "param_functors.h"
#include "mesh.h"

//! A Stiff class: a class for the stiffness operator.
struct Stiff{
	//! A definition of operator () taking three arguments.
    /*!
     * Evaluates the stiffness operator (i,j) of the current planar finite element.
     * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,2,2>, current planar finite element
     * returns a matrix!.
     */
	template<class Integrator, UInt ORDER, UInt mydim, UInt ndim>
	typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
  operator() (FiniteElement<Integrator,ORDER,mydim,ndim>& currentfe_, UInt iq){
	  return currentfe_.stiff_impl(iq);
	}
};

struct Grad{
	//! A definition of operator () taking three arguments.
	/*!
	* Evaluates the component ic of the vGrad operator (i,j) on the current finite elemente.
	* \param i is an unsigned int, current finite element local index
	* \param j is an unsigned int, current finite element local index
	* \param ic is an unsigned int, vGrad component to be evaluated
	* returns a double.
	*/

	template<class Integrator, UInt ORDER, UInt mydim, UInt ndim>
	typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
	operator() (FiniteElement<Integrator,ORDER,mydim,ndim>& currentfe_, UInt iq){
		return currentfe_.grad_impl(iq, Eigen::Matrix<Real,ndim,1>::Unit(0));
	}
};

//! A Mass class: a class for the mass operator.
struct Mass{
  //! A definition of operator () taking three arguments.
  /*!
  * Evaluates the mass operator (i,j) of the current finite element.
  * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,2,2>, current planar finite element
  * \param i is an unsigned int, current finite element local index
  * \param j is an unsigned int, current finite element local index
  * returns a double.
  */

	template<class Integrator, UInt ORDER, UInt mydim, UInt ndim>
	typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
  operator() (FiniteElement<Integrator,ORDER,mydim,ndim>& currentfe_, UInt iq){
	  return currentfe_.mass_impl(iq);
	}
};

template <class T>
class StiffAnys{

	template<UInt ndim>
	using is_T_SV = typename std::enable_if<std::is_same<T, Diffusion<ndim, true> >::value, bool>::type;
	template<UInt ndim>
	using is_T_not_SV = typename std::enable_if<std::is_same<T, Diffusion<ndim, false> >::value, bool>::type;

	template <typename A>
	friend class EOExpr;

	//Note: the constructor is private because only EOExpr<Stiff> will need to use this type
	StiffAnys(const T& K): K_(K){};

	const T& K_;

public:
	//! A definition of operator () taking four arguments.
	/*!
	* Evaluates the product of: the derivative of basis(i) with respect to coordinate ic1 and the derivative of basis(j) with respect
	* to coordinate ic2 ,on current finite elemente.
	* \param i is an unsigned int, current finite element local index
	* \param j is an unsigned int, current finite element local index
	* \param ic1 is an unsigned int, the variable respect whom the derivative is take: ic1=0 abscissa, ic1=1 ordinata
	* \param ic1 is an unsigned int, the variable respect whom the derivative is take: ic1=0 abscissa, ic1=1 ordinata
	* returns a double.
	*/

	template<class Integrator, UInt ORDER, UInt mydim, UInt ndim, is_T_not_SV<ndim> = true>
	typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
	operator() (FiniteElement<Integrator,ORDER,mydim,ndim>& currentfe_, UInt iq){
		return currentfe_.stiff_anys_impl(iq, K_());
	}

	template<class Integrator, UInt ORDER, UInt mydim, UInt ndim, is_T_SV<ndim> = true>
	typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
	operator() (FiniteElement<Integrator,ORDER,mydim,ndim>& currentfe_, UInt iq){
		UInt globalIndex = currentfe_.getGlobalIndex(iq);
		return currentfe_.stiff_anys_impl(iq, K_(globalIndex));
	}

};

//! A bGrad class: a class for Gradient operator.with given b
template <class T>
class bGrad{
	template<UInt ndim>
	using is_T_SV = typename std::enable_if<std::is_same<T, Advection<ndim, true> >::value, bool>::type;
	template<UInt ndim>
	using is_T_not_SV = typename std::enable_if<std::is_same<T, Advection<ndim, false> >::value, bool>::type;

	template <typename A>
	friend class EOExpr;

	//Note: the constructor is private because only EOExpr<Grad> will need to use this type

	bGrad(const T& b): b_(b){};

	const T& b_;

public:
	//! A definition of operator () taking three arguments.
	/*!
	* Evaluates the component ic of the vGrad operator (i,j) on the current finite elemente.
	* \param i is an unsigned int, current finite element local index
	* \param j is an unsigned int, current finite element local index
	* \param ic is an unsigned int, vGrad component to be evaluated
	* returns a double.
	*/

	template<class Integrator, UInt ORDER, UInt mydim, UInt ndim, is_T_not_SV<ndim> = true>
	typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
	operator() (FiniteElement<Integrator,ORDER,mydim,ndim>& currentfe_, UInt iq){
		return currentfe_.grad_impl(iq,b_());
	}

	template<class Integrator, UInt ORDER, UInt mydim, UInt ndim, is_T_SV<ndim> = true>
	typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
	operator() (FiniteElement<Integrator,ORDER,mydim,ndim>& currentfe_, UInt iq){
		UInt globalIndex = currentfe_.getGlobalIndex(iq);
		return currentfe_.grad_impl(iq,b_(globalIndex));
	}

};

//generic template class wrapper
//! A ETWrapper class: Expression Template Wrapper.
/*!
 * Class that mimic the behaviour of a generic operator defined above: following
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */



template<typename A>
class EOExpr{
	  //! "A" is a generic type
	  A a_;
public:
	//! A constructor.
	/*!
	 * \param object is a constant reference to a generic operator.
	 */
	 EOExpr(const A& a):a_(a){};

	 template<typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
	 typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
 	 operator() (FiniteElement<Integrator, ORDER,mydim,ndim>& currentfe_, UInt iq){
		 return a_(currentfe_, iq);
	 }
};

template<>
class EOExpr<Stiff>{

	Stiff a_;

public:

	EOExpr(const Stiff& a):a_(a){};

	// The subscript operator is used to convert from Stiff to StiffAnys
	template<UInt ndim, bool is_SV>
	EOExpr<StiffAnys<Diffusion<ndim,is_SV> > >  operator[] (const Diffusion<ndim, is_SV>& K){
		typedef EOExpr<StiffAnys<Diffusion<ndim,is_SV> > > ExprT;
		StiffAnys<Diffusion<ndim,is_SV> > anys(K);
		return ExprT(anys);
	}

	template<typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
	typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
	operator() (FiniteElement<Integrator, ORDER,mydim,ndim>& currentfe_, UInt iq){
		return a_(currentfe_, iq);
	}
};

template<>
class EOExpr<Grad>{
private:

	Grad a_;

public:
	//! A constructor.
	/*!
	* \param object is a constant reference to a generic operator.
	*/
	EOExpr(const Grad& a):a_(a){};

	// The subscript operator is used to convert from Grad to bGrad
	template<UInt ndim, bool is_SV>
	EOExpr<bGrad<Advection<ndim, is_SV> > > operator[] (const Advection<ndim, is_SV>& b){
		typedef EOExpr<bGrad<Advection<ndim, is_SV> > > ExprT;
		bGrad<Advection<ndim, is_SV> > bgrad(b);
		return ExprT(bgrad);
	}

	template<typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
	typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
 	operator() (FiniteElement<Integrator, ORDER,mydim,ndim>& currentfe_, UInt iq){
		return a_(currentfe_, iq);
	}
};


//composition of two wrappers (operator)
template<typename A, typename B, typename Op>
class EOBinOp{
		//! "A" is a generic type.
		/*!
		 * Stores the first operand.
		 */
		A a_;
		//! "B" is a generic type.
		/*!
		 * Stores the second operand.
		 */
		B b_;
public:
	//! A constructor.
	/*!
	 * \param a is a constant reference to a generic type.
	 * \param b is a constant reference to a generic type.
	 */
		EOBinOp(const A& a ,const B& b): a_(a), b_(b){};
	 //! A definition of operator () taking two arguments.
	 /*!
     * \param i is an unsigned int
     * \param j is an unsigned int
     * applies the generic operation defined by the type Op to the two generic objects a_, b_;
     * returns a type P variable
	 */
		template<typename Integrator, UInt ORDER,UInt mydim,UInt ndim>
		typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
 		operator () (FiniteElement<Integrator, ORDER,mydim,ndim>& currentfe_, UInt iq){
		  return Op::apply(a_(currentfe_, iq), b_(currentfe_, iq));
	  }
	};

template<class B, class Op>
class EOBinOp<Real, B, Op>{
	Real M_a;
	B M_b;
public:
	EOBinOp(Real a, const B& b):M_a(a),M_b(b) {};

	template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
	typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
 	operator()(FiniteElement<Integrator, ORDER,mydim,ndim>& currentfe_, int iq){
		return Op::apply(M_a, M_b(currentfe_, iq));
	}
};

template<class B, class Op>
class EOBinOp<Reaction, B, Op>{
	const Reaction& M_a;
	B M_b;
public:
	EOBinOp(const Reaction& a, const B& b):M_a(a), M_b(b) {};

	template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
	typename FiniteElement<Integrator, ORDER, mydim, ndim>::return_t
 	operator()(FiniteElement<Integrator, ORDER, mydim, ndim>& currentfe_, int iq){
		UInt globalIndex = currentfe_.getGlobalIndex(iq);
		return Op::apply(M_a(globalIndex),M_b(currentfe_, iq));
	}
};

//wrappers addition
//! A ETWAdd class: Expression Template Wrapper Addition
/*!
 * Class that defines Addition operation, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */

struct EOAdd{
	//! A static inline method taking two arguments.
	/*!
	 *The actual addition operation
	 * \param a is of P type, first addend
	 * \param b is of P type, second addend
	 */
	template<UInt N>
 	static Eigen::Matrix<Real,N,N> apply(const Eigen::Matrix<Real,N,N> &a, const Eigen::Matrix<Real,N,N> &b){
 		return a+b;
 	}
};

//multiplication by real scalar
//! A ETWMult class: Expression Template Wrapper Multiplication.
/*!
 * Class that defines Multiplication operation, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */

struct EOMult{
	 //! A static inline method taking two arguments.
	/*!
	 * The actual multiplication operation.
	 * \param a is of P type, first operand
	 * \param b is a Real, second operand
	 */
	 template<UInt N>
   static Eigen::Matrix<Real,N,N> apply(Real a, const Eigen::Matrix<Real,N,N> &b){
  		return a*b;
   }
};

//operator +
//! Overloading of operator +.
/*!
 * Following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 * Takes two arguments:
 * \param a is const reference ETWrapper<P, A>
 * \param b is const reference ETWrapper<P, A>
 * \return a ETWrapper<P,ETWBinOp<P, ETWrapper<P,A>, ETWrapper<P, B>, ETWAdd<P> > which is resolved at compile time.
 */
template<typename A, typename B>
EOExpr<EOBinOp<EOExpr<A>, EOExpr<B>, EOAdd > >
operator + (const EOExpr<A>&  a, const EOExpr<B>&  b){

	  typedef EOBinOp<EOExpr<A>, EOExpr<B>, EOAdd > ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}

template<typename B>
EOExpr<EOBinOp<Reaction, EOExpr<B>, EOMult > >
operator * (const Reaction&  a, const EOExpr<B>&  b){

	  typedef EOBinOp<Reaction, EOExpr<B>, EOMult> ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}

template<typename B>
EOExpr<EOBinOp<Real, EOExpr<B>, EOMult > >
operator * (Real a, const EOExpr<B>&  b){

	  typedef EOBinOp<Real, EOExpr<B>, EOMult > ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}


//!A Assmbler class: discretize a generic differential operator in a sparse matrix
//template<UInt mydim, UInt ndim>
struct Assembler{
  //! A template member taking three arguments: discretize differential operator
  /*!
   * \param oper is a template expression : the differential operator to be discretized.
   * \param mesh is const reference to a MeshHandler<ORDER,2,2>: the mesh where we want to discretize the operator.
   * \param fe is a const reference to a FiniteElement
   * stores the discretization in SPoper_mat_
   */

  //Return triplets vector
  template<UInt ORDER, typename Integrator, UInt mydim, UInt ndim, typename A>
  static void operKernel(EOExpr<A> oper,const MeshHandler<ORDER,mydim,ndim>& mesh,
  	                     FiniteElement<Integrator, ORDER,mydim,ndim>& fe, SpMat& OpMat);

  template<UInt ORDER, typename Integrator, UInt mydim, UInt ndim>
  static void forcingTerm(const MeshHandler<ORDER,mydim,ndim>& mesh, FiniteElement<Integrator, ORDER,mydim,ndim>& fe, const ForcingTerm& u, VectorXr& forcingTerm);

};


#include "matrix_assembler_imp.h"

#endif
