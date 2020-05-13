#ifndef __REGRESSIONDATA_IMP_HPP__
#define __REGRESSIONDATA_IMP_HPP__

template<UInt ndim>
RegressionData<ndim>::RegressionData(std::vector<Point<ndim> >& locations, VectorXr& observations, UInt order, std::vector<Real> lambda, MatrixXr& covariates, MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF):
					locations_(locations), observations_(observations), covariates_(covariates), incidenceMatrix_(incidenceMatrix),
					order_(order), lambda_(lambda), bc_values_(bc_values), bc_indices_(bc_indices), DOF_(DOF)
{
	nRegions_ = incidenceMatrix_.rows();
	if(locations_.size()==0 && nRegions_==0)
	{
		locations_by_nodes_ = true;
		for(int i = 0; i<observations_.size();++i)
			observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}
}

template<UInt ndim>
RegressionDataElliptic<ndim>::RegressionDataElliptic(std::vector<Point<ndim> >& locations, VectorXr& observations, UInt order,
												std::vector<Real> lambda, Eigen::Matrix<Real,ndim,ndim>& K,
												Eigen::Matrix<Real,ndim,1>& beta, Real c, MatrixXr& covariates,
												MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,
												std::vector<Real>& bc_values, bool DOF) :
		 RegressionData<ndim>(locations, observations, order, lambda, covariates, incidenceMatrix, bc_indices, bc_values, DOF), K_(K), beta_(beta), c_(c) {}

template<UInt ndim>
RegressionDataEllipticSpaceVarying<ndim>::RegressionDataEllipticSpaceVarying(std::vector<Point<ndim> >& locations,
									VectorXr& observations, UInt order, std::vector<Real> lambda,
									const std::vector<Eigen::Matrix<Real,ndim,ndim>, Eigen::aligned_allocator<Eigen::Matrix<Real,ndim,ndim> > >& K,
									const std::vector<Eigen::Matrix<Real,ndim,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,ndim,1> > >& beta,
									const std::vector<Real>& c, const std::vector<Real>& u,
									MatrixXr& covariates, MatrixXi& incidenceMatrix,
									std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF):
		RegressionData<ndim>(locations, observations, order, lambda, covariates, incidenceMatrix, bc_indices, bc_values, DOF), K_(K), beta_(beta), c_(c), u_(u) {}


#ifdef R_VERSION_
template<UInt ndim>
RegressionData<ndim>::RegressionData(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP Rcovariates,
							SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF, SEXP RGCVmethod,
							SEXP Rnrealizations)
{
	setLocations(Rlocations);
	setIncidenceMatrix(RincidenceMatrix);
	setObservations(Robservations);
	setCovariates(Rcovariates);
	setNrealizations(Rnrealizations);

	GCVmethod_ = INTEGER(RGCVmethod)[0];

	order_ =  INTEGER(Rorder)[0];
	DOF_ = INTEGER(DOF)[0];
	UInt length_indexes = Rf_length(RBCIndices);

	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) +  length_indexes);

	bc_values_.assign(REAL(RBCValues),REAL(RBCValues) + Rf_length(RBCIndices));

  UInt length_lambda = Rf_length(Rlambda);
  for (UInt i = 0; i<length_lambda; ++i)
		lambda_.push_back(REAL(Rlambda)[i]);

}

template<UInt ndim>
RegressionDataElliptic<ndim>::RegressionDataElliptic(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta,
				 SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations) :
				 		RegressionData<ndim>(Rlocations, Robservations, Rorder, Rlambda, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, DOF,RGCVmethod, Rnrealizations),
						K_(RK),	beta_(Rbeta), c_(REAL(Rc)[0]) {}

template<UInt ndim>
RegressionDataEllipticSpaceVarying<ndim>::RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta,
				 SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations):
					 RegressionData<ndim>(Rlocations, Robservations, Rorder, Rlambda, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, DOF,RGCVmethod, Rnrealizations),
					 K_(RK), beta_(Rbeta), c_(Rc), u_(Ru) {}



// template<UInt ndim>
// void RegressionDataEllipticSpaceVarying<ndim>::print(std::ostream & out) const
// {
// }

template<UInt ndim>
void RegressionData<ndim>::setObservations(SEXP Robservations)
{
	UInt n_obs_ = Rf_length(Robservations);
	observations_.resize(n_obs_);
	observations_indices_.reserve(n_obs_);

	UInt count = 0;
	if(locations_.size() == 0 && nRegions_ == 0)
	{
		locations_by_nodes_ = true;
		for(auto i=0;i<n_obs_;++i)
		{
			if(!ISNA(REAL(Robservations)[i]))
			{
				observations_[count] = REAL(Robservations)[i];
				count++;
				observations_indices_.push_back(i);
			}
		}
		observations_.conservativeResize(count, Eigen::NoChange);
	}
	else // locations_.size() > 0 NOR nRegions_ > 0
	{
		locations_by_nodes_ = false;
		for(auto i=0;i<n_obs_;++i)
		{
			observations_[i] = REAL(Robservations)[i];
		}
	}

	//std::cout<<"Observations #"<<observations_.size()<<std::endl<<observations_<<std::endl;
	//for(auto i=0;i<observations_indices_.size();++i)	std::cout<<observations_indices_[i]<<std::endl;
}

template<UInt ndim>
void RegressionData<ndim>::setCovariates(SEXP Rcovariates)
{
	n_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[0];
	p_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[1];

	covariates_.resize(n_, p_);

	for(auto i=0; i<n_; ++i)
	{
		for(auto j=0; j<p_ ; ++j)
		{
			covariates_(i,j)=REAL(Rcovariates)[i+ n_*j];
		}
	}
}

template<UInt ndim>
void RegressionData<ndim>::setLocations(SEXP Rlocations)
{
	n_ = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	if(n_>0){
		locations_.reserve(n_);
		for(int i=0; i<n_; ++i)
			locations_.emplace_back(Point<ndim>(i, &REAL(Rlocations)[0], n_));
	}
}

template<UInt ndim>
void RegressionData<ndim>::setNrealizations(SEXP Rnrealizations) {
	nrealizations_ = INTEGER(Rnrealizations)[0];
}

template<UInt ndim>
void RegressionData<ndim>::setIncidenceMatrix(SEXP RincidenceMatrix)
{
	nRegions_ = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[0];
	UInt p = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[1];

	incidenceMatrix_.resize(nRegions_, p);

	for(auto i=0; i<nRegions_; ++i)
	{
		for(auto j=0; j<p; ++j)
		{
			incidenceMatrix_(i,j) = INTEGER(RincidenceMatrix)[i+nRegions_*j];
		}
	}
}

#endif

template<UInt ndim>
void RegressionData<ndim>::printObservations(std::ostream & out) const
{

	for(auto i=0;i<observations_.size(); i++)
	{
		out<<i<<"\t"<<observations_(i)<<std::endl;
	}
}

template<UInt ndim>
void RegressionData<ndim>::printCovariates(std::ostream & out) const
{

	for(auto i=0;i<covariates_.rows(); i++)
	{
		for(auto j=0; j<covariates_.cols(); j++)
		{
			out<<covariates_(i,j)<<"\t";
		}
		out<<std::endl;
	}
}

template<UInt ndim>
void RegressionData<ndim>::printLocations(std::ostream & out) const
{

	for(int i=0; i<locations_.size(); i++)
	{
		out<<locations_[i];
		//std::cout<<std::endl;
	}
}

template<UInt ndim>
void RegressionData<ndim>::printIncidenceMatrix(std::ostream & out) const
{
	for (auto i=0; i<incidenceMatrix_.rows(); i++)
	{
		for (auto j=0; j<incidenceMatrix_.cols(); j++)
		{
			out << incidenceMatrix_(i,j) << "\t";
		}
		out << std::endl;
	}
}

#endif
