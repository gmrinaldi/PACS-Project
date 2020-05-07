#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include "fdaPDE.h"
#include "mesh_objects.h"

struct IntegratorTriangleP2{
	static constexpr UInt ORDER = 1;
	//Number of nodes
	static constexpr UInt NNODES = 3;
	static constexpr std::array<Real,NNODES> WEIGHTS{{1./3, 1./3, 1./3}};
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<2>,NNODES> NODES{
		Point<2>({1./6,1./6}),
		Point<2>({2./3,1./6}),
		Point<2>({1./6,2./3})
	};
};

struct IntegratorTriangleP4{
	static constexpr UInt ORDER = 2;
	//Number of nodes
	static constexpr UInt NNODES = 6;
	static constexpr std::array<Real,NNODES> WEIGHTS{{0.223381589678011,0.223381589678011,0.223381589678011,0.109951743655322,0.109951743655322,0.109951743655322}};
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<2>,NNODES> NODES{
		Point<2>({0.445948490915965,0.445948490915965}),
		Point<2>({0.445948490915965,0.108103018168070}),
		Point<2>({0.108103018168070,0.445948490915965}),
		Point<2>({0.091576213509771,0.091576213509771}),
		Point<2>({0.091576213509771,0.816847572980459}),
		Point<2>({0.816847572980459,0.091576213509771})
	};
};

struct IntegratorTetrahedronP1{
	static constexpr UInt ORDER = 1;
	//Number of nodes
	static constexpr UInt NNODES = 1;
	static constexpr std::array<Real,NNODES> WEIGHTS{{1}};
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<3>,NNODES> NODES{
		Point<3>({1./4 ,1./4 ,1./4})
	};
};

struct IntegratorTetrahedronP2{
	static constexpr UInt ORDER = 1;
	//Number of nodes
	static constexpr UInt NNODES = 4;
	static constexpr std::array<Real,NNODES> WEIGHTS{{1./4,1./4,1./4,1./4}};
	//Point locations (in barycentric coordinates)
	static constexpr std::array<Point<3>,NNODES> NODES{
		Point<3>({0.585410196624969,0.138196601125011,0.138196601125011}),
		Point<3>({0.138196601125011,0.138196601125011,0.138196601125011}),
		Point<3>({0.138196601125011,0.138196601125011,0.585410196624969}),
		Point<3>({0.138196601125011,0.585410196624969,0.138196601125011})
	};
};


#endif
