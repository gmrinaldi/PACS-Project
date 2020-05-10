/*
 * param_functors.hpp
 *
 *  Created on: Jun 2, 2015
 *      Author: eardi
 */
#ifndef PARAM_FUNCTORS_H_
#define PARAM_FUNCTORS_H_

//#include "matrix_assembler.hpp"
class Diffusivity{
  using diffusion_matr = Eigen::Matrix<Real,2,2>;
	using diff_matr_container = std::vector<diffusion_matr, Eigen::aligned_allocator<diffusion_matr> >;
	diff_matr_container K_;
public:
	Diffusivity(const diffusion_matr& K) :
			K_({K}) {}
	Diffusivity(const diff_matr_container& K):
			K_(K){}
	#ifdef R_VERSION_
	Diffusivity(SEXP RGlobalVector){
		UInt num_int_nodes = Rf_length(RGlobalVector)/4;
		K_.reserve(num_int_nodes);
		for(UInt i=0; i<num_int_nodes; ++i)
				K_.push_back(Eigen::Map<const diffusion_matr>(&REAL(RGlobalVector)[4*i]));
	}
	#endif

	const diffusion_matr& operator()() const {return K_[0];}
	const diffusion_matr& operator()(UInt globalNodeIndex) const {return K_[globalNodeIndex];}

	diff_matr_container::const_iterator begin() const {return K_.begin();}
	diff_matr_container::const_iterator end() const {return K_.end();}
};

class Advection{
	using advection_vec = Eigen::Matrix<Real,2,1>;
	using adv_vec_container = std::vector<advection_vec, Eigen::aligned_allocator<advection_vec> >;
	adv_vec_container beta_;
public:
	Advection(const advection_vec& beta):
		beta_({beta}) {}
	Advection(const adv_vec_container& beta):
		beta_(beta){}
	#ifdef R_VERSION_
	Advection(SEXP RGlobalVector){
		UInt num_int_nodes = Rf_length(RGlobalVector)/2;
		beta_.reserve(num_int_nodes);
		for(UInt i=0; i<num_int_nodes; ++i)
			beta_.push_back(Eigen::Map<const advection_vec>(&REAL(RGlobalVector)[2*i]);

	}
	#endif
	const advection_vec& operator()() const {return beta_[0];}
	const advection_vec& operator()(UInt globalNodeIndex) const {return beta_[globalNodeIndex];}

	adv_vec_container::const_iterator begin() const {return beta_.begin();}
	adv_vec_container::const_iterator end() const {return beta_.end();}

};

class Reaction{
	std::vector<Real> c_;
public:
	Reaction(const Real &c)
		c_({c}) {}
	Reaction(const std::vector<Real>& c):
		c_(c) {}
#ifdef R_VERSION_
	Reaction(SEXP RGlobalVector){
		UInt num_int_nodes = Rf_length(RGlobalVector);
		c_.reserve(num_int_nodes);
		for(UInt i=0; i<num_int_nodes; ++i)
			c_.push_back(REAL(RGlobalVector)[i]);
	}
	#endif
	const Real& operator()() const {return c_[0];}
	const Real& operator()(UInt globalNodeIndex) const {return c_[globalNodeIndex];}

	std::vector<Real>::const_iterator begin() const {return c_.begin();}
	std::vector<Real>::const_iterator end() const {return c_.end();}

};

class ForcingTerm{
	std::vector<Real> u_;
public:
	ForcingTerm(const std::vector<Real>& u):
		u_(u) {}
	#ifdef R_VERSION_
	ForcingTerm(SEXP RGlobalVector){
		UInt num_int_nodes = Rf_length(RGlobalVector);
		u_.reserve(num_int_nodes);
		for(UInt i=0; i<num_int_nodes; ++i)
			u_.push_back(REAL(RGlobalVector)[i]);
	}
	#endif
	Real operator()(UInt globalNodeIndex) const {return u_[globalNodeIndex];}
};


#endif /* PARAM_FUNCTORS_H_ */
