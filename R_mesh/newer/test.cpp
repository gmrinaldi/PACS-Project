#include "matrix_assembler.h"
#include <iostream>
#include <chrono>

template <class T>
void doNotOptimizeAway(T&& datum) {
  asm volatile("" : "+r" (datum));
}

int main () {
  // Point<3> a(1,{0.109610,1,0.597550}),b(2,{0.070902,0.0770110,0.523130}),c(3,{0.090256, 0.0385055, 0.560340}),d(4,{0.085451, 0.0385055, 0.511565});
  // Point<3> a(1,{0,0,1}),b(2,{2,0,0}),c(3,{0,5,0});
  // Point<3> a(1,{0.109610,0.0000000,0.597550}),b(2,{0.070902,0.0770110,0.523130}),c(3,{0.100000,0.0000000,0.500000}),
  //           d(4,{0.085451,0.0385055,0.511565}),e(5,{0.104805, 0.0000000, 0.548775}),f(6,{0.090256, 0.0385055, 0.560340});
  //

  // Element<4,3,3> myel1(1,{a,b,c,d});
  // FiniteElement<1, 3, 3> fe;
  //
  // fe.updateElement(myel1);
  //
  // std::cout<<FiniteElement< 1, 3, 3>::Integrator::NODES[0]<<std::endl;
  // Element<3,2,3> myel1(1,{a,b,c});
  // FiniteElement<1, 2, 3> fe;
  // fe.updateElement(myel1);

  // Eigen::Matrix<Real, 6,6> result=Eigen::Matrix<Real, 6,6>::Zero();
  Point<2> a(1,{0,0}),b(2,{2,0}),c(3,{0,5}), d(4,{1,7}), e(5,{0,1}),f(6,{3,0});
  Element<3,2,2> myel1(1,{a,b,c});
  Element<3,2,2> myel2(1,{d,e,f});

  FiniteElement<1,2,2> fe;
  fe.updateElement(myel1);
  //
  // Eigen::Matrix<Real,2,2> K1, K2, K3;
  // K1 << 3, 2,
  //       2, 5;
  //
  // K2 << 3, 0,
  //       0, 1;
  //
  // K3 << 1, 4,
  //       4, 1;
  //
  // Diffusion<2> K_(K1);

  // Diffusion<2,true>::diff_matr_container K_({K1,K2,K3});
  //
  // Diffusion<2,true> K(K_);
  //
  // Eigen::Matrix<Real,2,1> b_;
  // b_<< 1, -3;
  //
  // Advection<2> b__(b_);

  Stiff OStiff;
  EOExpr<Stiff> stiff(OStiff);
  Mass OMass;
  EOExpr<Mass> mass(OMass);
  Grad OGrad;
  EOExpr<Grad> grad(OGrad);

  // Eigen::Matrix<Real,3,3> K;
  // K << 1,0,1,
  //      0,1,0,
  //      1,0,1;
  // Diffusion<3> K_(K);


  Eigen::Matrix<Real,3,3> matr=Eigen::Matrix<Real,3,3>::Zero();

  auto start = std::chrono::high_resolution_clock::now();

  // for (int u=0; u<100000; ++u)
    for (int i=0; i<1; ++i)
      matr+=(mass+stiff+grad)(fe, i)*FiniteElement<1, 2, 2>::Integrator::WEIGHTS[i];

  auto stop = std::chrono::high_resolution_clock::now();

  std::cout.precision(17);

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  std::cout << duration.count() << std::endl;

  std::cout << matr << std::endl;


  // std::cout<<result<<std::endl;

}
