#ifndef __EVALUATOR_IMP_HPP__
#define __EVALUATOR_IMP_HPP__

template <UInt ORDER>
void Evaluator<ORDER,2,2>::eval(Real* X, Real *Y, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside)
{
	static constexpr UInt NNodes = 3*ORDER;
	Element<NNodes,2,2> current_element;
	Point<2> current_point;
	Eigen::Matrix<Real,NNodes,1> coefficients;

	Element<NNodes,2,2> starting_element{mesh_.getElement(0)};
	for (int i = 0; i<length; ++i)
	{
		current_point = Point<2>({X[i],Y[i]});
		current_element = mesh_.findLocationWalking(current_point, starting_element);

		if(current_element.unassignedId() && redundancy == true)
			current_element = mesh_.findLocationNaive(current_point);

		if(current_element.unassignedId())
			isinside[i]=false;
		else
		{
			isinside[i]=true;

			for (int j=0; j<NNodes; ++j)
				coefficients[j] = coef[current_element[j].getId()];

			result[i] = current_element.evaluate_point(current_point, coefficients);
			starting_element = current_element;
		}
	}
}


template <UInt ORDER>
void Evaluator<ORDER,2,3>::eval(Real* X, Real *Y,  Real *Z, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside)
{

	static constexpr UInt NNodes = 3*ORDER;

	Element<NNodes,2,3> current_element;
	Point<3> current_point;

	Eigen::Matrix<Real,NNodes,1> coefficients;
	for (int i = 0; i<length; ++i)
	{
		current_point = Point<3>({X[i],Y[i],Z[i]});
		current_element = mesh_.findLocationNaive(current_point);

		if(current_element.unassignedId())
			isinside[i]=false;
		else
		{
			isinside[i]=true;

			for (int j=0; j<NNodes; ++j)
				coefficients[j] = coef[current_element[j].getId()];

			result[i] = current_element.evaluate_point(current_point, coefficients);

		}
	}
}


template <UInt ORDER>
void Evaluator<ORDER,3,3>::eval(Real* X, Real *Y,  Real *Z, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside)
{

	static constexpr UInt NNodes = 6*ORDER-2;
	Element<NNodes,3,3> current_element;

	Point<3> current_point;

	Eigen::Matrix<Real,NNodes,1> coefficients;
	Element<NNodes,3,3>starting_element = mesh_.getElement(0);

	for (int i = 0; i<length; ++i)
	{
		current_point = Point<3>({X[i],Y[i],Z[i]});
		current_element = mesh_.findLocationWalking(current_point, starting_element);

		if(current_element.unassignedId() && redundancy == true)
			current_element = mesh_.findLocationNaive(current_point);

		if(current_element.unassignedId())
			isinside[i]=false;
		else
		{
			isinside[i]=true;

			for (int j=0; j<NNodes; ++j)
				coefficients[j] = coef[current_element[j].getId()];

			result[i] = current_element.evaluate_point(current_point, coefficients);
			starting_element = current_element;
		}
	}
}


template <UInt ORDER>
void Evaluator<ORDER, 2, 2>::integrate(UInt** incidenceMatrix, UInt nRegions, UInt nElements, const Real *coef, Real* result)
{
	std::vector<Real> Delta(nRegions, 0);
	std::vector<Real> integral(nRegions, 0);
	static constexpr UInt NNodes = 3*ORDER;
	Element<NNodes, 2, 2> current_element;
	Eigen::Matrix<Real,NNodes,1> coefficients;

	for (int region=0; region<nRegions; ++region)
	{
		for (int elem=0; elem<nElements; ++elem)
		{
			if (incidenceMatrix[region][elem]==1) //elem is in region
			{
				current_element = mesh_.getElement(elem);
				for (int i=0; i<NNodes; ++i)
					coefficients[i]=coef[current_element[i].getId()];
				Delta[region] += current_element.getMeasure();
				integral[region] += current_element.integrate(coefficients);
			}
		}
		result[region]=integral[region]/Delta[region];
	}
}


template <UInt ORDER>
void Evaluator<ORDER, 2, 3>::integrate(UInt** incidenceMatrix, UInt nRegions, UInt nElements, const Real *coef, Real* result)
{
	std::vector<Real> Delta(nRegions, 0);
	std::vector<Real> integral(nRegions, 0);
	static constexpr UInt NNodes = 3*ORDER;
	Element<NNodes, 2, 3> current_element;
	Eigen::Matrix<Real,NNodes,1> coefficients;


	for (int region=0; region<nRegions; ++region)
	{
		for (int elem=0; elem<nElements; ++elem)
		{
			if (incidenceMatrix[region][elem]==1) //elem is in region
			{
				current_element = mesh_.getElement(elem);
				for (int i=0; i<NNodes; ++i)
					coefficients[i]=coef[current_element[i].getId()];
				Delta[region] += current_element.getMeasure();
				integral[region] += current_element.integrate(coefficients);
			}
		}
		result[region]=integral[region]/Delta[region];
	}
}


template <UInt ORDER>
void Evaluator<ORDER, 3, 3>::integrate(UInt** incidenceMatrix, UInt nRegions, UInt nElements, const Real *coef, Real* result)
{
	std::vector<Real> Delta(nRegions, 0);
	std::vector<Real> integral(nRegions, 0);
	static constexpr UInt NNodes = 6*ORDER-2;
	Element<NNodes, 3, 3> current_element;
	Eigen::Matrix<Real,NNodes,1> coefficients;

	for (int region=0; region<nRegions; ++region)
	{
		for (int elem=0; elem<nElements; ++elem)
		{
			if (incidenceMatrix[region][elem]==1) //elem is in region
			{
				current_element = mesh_.getElement(elem);
				for (int i=0; i<NNodes; ++i)
					coefficients[i]=coef[current_element[i].getId()];
				Delta[region] += current_element.getMeasure();
				integral[region] += current_element.integrate(coefficients);
			}
		}
		result[region]=integral[region]/Delta[region];
	}
}


#endif
