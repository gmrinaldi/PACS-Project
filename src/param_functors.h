/*
 * param_functors.hpp
 *
 *  Created on: Jun 2, 2015
 *      Author: eardi
 */
#ifndef PARAM_FUNCTORS_H_
#define PARAM_FUNCTORS_H_

template <UInt ndim, bool is_space_varying = false>
struct Diffusion{
  using diffusion_matr = Eigen::Matrix<Real,ndim,ndim>;

	Diffusion(const diffusion_matr& K):
			K_(K){}
	const diffusion_matr& operator()() const {return K_;}

private:
  diffusion_matr K_;

};

template <UInt ndim>
struct Diffusion<ndim, true>{
  using diffusion_matr = Eigen::Matrix<Real,ndim,ndim>;
	using diff_matr_container = std::vector<diffusion_matr, Eigen::aligned_allocator<diffusion_matr> >;

	Diffusion(const diff_matr_container& K):
			K_(K){}
  #ifdef R_VERSION_
	Diffusion(SEXP RGlobalVector){
		UInt num_int_nodes = Rf_length(RGlobalVector)/(ndim*ndim);
		K_.reserve(num_int_nodes);
		for(UInt i=0; i<num_int_nodes; ++i)
				K_.push_back(Eigen::Map<const diffusion_matr>(&REAL(RGlobalVector)[ndim*ndim*i]));
	}
	#endif
	const diffusion_matr& operator()(UInt globalNodeIndex) const {return K_[globalNodeIndex];}

private:
  diff_matr_container K_;

};

template <UInt ndim, bool is_space_varying = false>
struct Advection{
	using advection_vec = Eigen::Matrix<Real,ndim,1>;

	Advection(const advection_vec& beta):
		beta_(beta){}

	const advection_vec& operator()() const {return beta_;}

private:
  advection_vec beta_;

};

template<UInt ndim>
struct Advection<ndim, true>{
	using advection_vec = Eigen::Matrix<Real,ndim,1>;
	using adv_vec_container = std::vector<advection_vec, Eigen::aligned_allocator<advection_vec> >;

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
	const advection_vec& operator()(UInt globalNodeIndex) const {return beta_[globalNodeIndex];}

private:
  adv_vec_container beta_;

};

class Reaction{
	std::vector<Real> c_;
public:
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
	const Real& operator()(UInt globalNodeIndex) const {return c_[globalNodeIndex];}

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
