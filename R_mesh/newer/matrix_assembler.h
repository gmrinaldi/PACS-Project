#ifndef MATRIX_ASSEMBLER_H_
#define MATRIX_ASSEMBLER_H_

#include "finite_element.h"
#include "mesh_objects.h"
#include "param_functors.h"
#include "mesh.h"


struct Stiff {
  template<UInt ORDER, UInt mydim, UInt ndim>
  auto operator() (FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(fe_.stiff_impl(iq)) {
    return fe_.stiff_impl(iq);
  }
};

struct Grad {
  template<UInt ORDER, UInt mydim, UInt ndim>
  auto operator() (FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(fe_.grad_impl(iq)) {
    return fe_.grad_impl(iq);
  }
};

struct Mass {
  template<UInt ORDER, UInt mydim, UInt ndim>
  auto operator() (FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(fe_.mass_impl(iq)) {
    return fe_.mass_impl(iq);
  }
};

template <class T>
class StiffAnys{
  template <typename A>
	friend class EOExpr;

  StiffAnys(const T& K) : K_(K) {}
  const T& K_;

  template<UInt ORDER, UInt mydim, UInt ndim>
  auto operator() (FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(fe_.stiff_impl(iq, K_)) {
    return fe_.stiff_impl(iq, K_);
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
	 EOExpr(const A& a) : a_(a){};

   template<class T>
   typename std::enable_if<std::is_same<A,Stiff>::value, EOExpr<StiffAnys<T> > >::type operator[] (const T& K){
     typedef EOExpr<StiffAnys<T> > ExprT;
     StiffAnys<T> anys(K);
     return ExprT(anys);
   }

	 template<UInt ORDER, UInt mydim, UInt ndim>
 	 auto operator() (FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(a_(fe_, iq)) {
     // std::cout<<type_name<decltype(a_(helper_, iq))>()<<std::endl;
		 return a_(fe_, iq);
	 }
};


class cMass{

  friend EOExpr<cMass> operator * (const Reaction&, const EOExpr<Mass>&);

  cMass(const Reaction& c) : c_(c) {}
  const Reaction& c_;
public:
  template<UInt ORDER, UInt mydim, UInt ndim>
  auto operator() (FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(c_(0)*fe_.mass_impl(iq)){
    UInt globalIndex=fe_.getGlobalIndex(iq);
    return c_(globalIndex)*fe_.mass_impl(iq);
  }
};

EOExpr<cMass> operator * (const Reaction&  c, const EOExpr<Mass>&  mass){
	  typedef EOExpr<cMass> ExprT;
    cMass cmass(c);
	  return ExprT(cmass);
}

template <class T>
class bGrad{

  template<class A>
  friend EOExpr<bGrad<A> > dot(const A&, const EOExpr<Grad>&);

  bGrad(const T& b) : b_(b) {}
  const T& b_;
public:
  template<UInt ORDER, UInt mydim, UInt ndim>
  auto operator() (FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(fe_.grad_impl(iq, b_)) {
    return fe_.grad_impl(iq, b_);
  }
};


template<class T>
EOExpr<bGrad<T> > dot(const T& b, const EOExpr<Grad>& grad){
  typedef EOExpr<bGrad<T> > ExprT;
  bGrad<T> bgrad(b);
  return ExprT(bgrad);
}




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
		EOBinOp(const A& a ,const B& b) : a_(a), b_(b) {};
	 //! A definition of operator () taking two arguments.
	 /*!
     * \param i is an unsigned int
     * \param j is an unsigned int
     * applies the generic operation defined by the type Op to the two generic objects a_, b_;
     * returns a type P variable
	 */
		template<UInt ORDER,UInt mydim,UInt ndim>
 		auto operator () (FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(Op::apply(a_(fe_, iq), b_(fe_, iq))) {
		  return Op::apply(a_(fe_, iq), b_(fe_, iq));
	  }
	};

template<class B, class Op>
class EOBinOp<Real, B, Op>{
	Real M_a;
	B M_b;
public:
	EOBinOp(Real a, const B& b) : M_a(a), M_b(b) {};

	template<UInt ORDER, UInt mydim, UInt ndim>
  auto operator () (FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(Op::apply(M_a, M_b(fe_, iq))) {
		return Op::apply(M_a, M_b(fe_, iq));
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
	template<class T, class V>
	static auto apply(const T& a, const V& b) -> decltype(a+b) {
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
	template<class T>
 	static auto apply(Real a, const T &b) -> decltype(a*b) {
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
EOExpr<EOBinOp<Real, EOExpr<B>, EOMult > >
operator * (Real a, const EOExpr<B>&  b){
	  typedef EOBinOp<Real, EOExpr<B>, EOMult > ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}

//!A Assembler class: discretize a generic differential operator in a sparse matrix
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
  template<UInt ORDER, UInt mydim, UInt ndim, typename A>
  static void operKernel(EOExpr<A> oper,const MeshHandler<ORDER,mydim,ndim>& mesh,
  	                     FiniteElement<ORDER,mydim,ndim>& fe, SpMat& OpMat);

  template<UInt ORDER, UInt mydim, UInt ndim>
  static void forcingTerm(const MeshHandler<ORDER,mydim,ndim>& mesh, FiniteElement<ORDER,mydim,ndim>& fe, const ForcingTerm& u, VectorXr& forcingTerm);

};




#include "matrix_assembler_imp.h"

#endif
