/*
 * mesh_object.cpp
 *
 *  Created on: Aug 16, 2015
 *      Author: eardi
 */

#include "mesh_objects.h"


void Point::print(std::ostream & out) const
{
	out<<"Point -"<< id_ <<"- "<<"("<<coord_[0]<<","<<coord_[1]<<","<<coord_[2]<<")"<<std::endl<<"------"<<std::endl;
}

void Edge::print(std::ostream & out) const
{
	out<<"Edge -"<< id_ <<"- "<<"("<<points_[0].getId()<<","<<points_[1].getId()<<")"<<std::endl;
}

std::vector<Real> point_diff(const Point &lhs, const Point &rhs){
		std::vector<Real> diff;
		for (int i=0; i<lhs.ndim; ++i)
				diff.push_back(lhs[i]-rhs[i]);
		return diff;
};
