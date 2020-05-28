#ifndef EXPRESSION_HH
#define EXPRESSION_HH

template <UInt ORDER, UInt mydim, UInt ndim>
class FiniteElement;

template <UInt ndim, bool is_SV>
struct Diffusion;

struct Reaction;


struct Stiff {
  template<UInt ORDER, UInt mydim, UInt ndim>
  auto operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(fe_.stiff_impl(iq)) {
    return fe_.stiff_impl(iq);
  }
};

struct Grad {
  template<UInt ORDER, UInt mydim, UInt ndim>
  auto operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(fe_.grad_impl(iq)) {
    return fe_.grad_impl(iq);
  }
};

struct Mass {
  template<UInt ORDER, UInt mydim, UInt ndim>
  auto operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(fe_.mass_impl(iq)) {
    return fe_.mass_impl(iq);
  }
};

template<typename A>
class EOExpr{
	  //! "A" is a generic type
	  A a_;
public:
	//! A constructor.
	/*!
	 * \param object is a constant reference to a generic operator.
	 */
	 EOExpr(const A& a) : a_(a) {};

	 template<UInt ORDER, UInt mydim, UInt ndim>
 	 auto operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(a_(fe_, iq)) {
		 return a_(fe_, iq);
	 }
};

template<>
class EOExpr<Stiff>{
	  //! "A" is a generic type
	  Stiff a_;
public:
	//! A constructor.
	/*!
	 * \param object is a constant reference to a generic operator.
	 */
	 EOExpr(const Stiff& a) : a_(a) {};

   template<UInt ndim, bool is_SV>
   EOExpr<const Diffusion<ndim, is_SV>&> operator[] (const Diffusion<ndim, is_SV>& K){
     typedef EOExpr<const Diffusion<ndim, is_SV>&> ExprT;
     return ExprT(K);
   }

	 template<UInt ORDER, UInt mydim, UInt ndim>
 	 auto operator() (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(a_(fe_, iq)) {
		 return a_(fe_, iq);
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
		EOBinOp(const A& a ,const B& b) : a_(a), b_(b) {};
	 //! A definition of operator () taking two arguments.
	 /*!
     * \param i is an unsigned int
     * \param j is an unsigned int
     * applies the generic operation defined by the type Op to the two generic objects a_, b_;
     * returns a type P variable
	 */
		template<UInt ORDER,UInt mydim,UInt ndim>
 		auto operator () (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(Op::apply(a_(fe_, iq), b_(fe_, iq))) {
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
  auto operator () (const FiniteElement<ORDER,mydim,ndim>& fe_, UInt iq) -> decltype(Op::apply(M_a, M_b(fe_, iq))) {
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
EOExpr<EOBinOp<EOExpr<A>, EOExpr<B>, EOAdd> >
operator + (const EOExpr<A>&  a, const EOExpr<B>&  b){

	  typedef EOBinOp<EOExpr<A>, EOExpr<B>, EOAdd > ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}

template<typename B>
EOExpr<EOBinOp<Real, EOExpr<B>, EOMult> >
operator * (Real a, const EOExpr<B>& b){
	  typedef EOBinOp<Real, EOExpr<B>, EOMult> ExprT;
	  return EOExpr<ExprT> (ExprT(a,b));
}

EOExpr<const Reaction&> operator * (const Reaction&  c, const EOExpr<Mass>&  mass){
    typedef EOExpr<const Reaction&> ExprT;
    return ExprT(c);
}



#endif
