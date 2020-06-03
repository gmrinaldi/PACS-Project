#ifndef __MESH_INPUT_HELPER_HPP__
#define __MESH_INPUT_HELPER_HPP__

#include <array>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <utility>
#include <type_traits>


template<UInt mydim>
class simplex{
public:

  using nodeIndices = std::array<UInt, mydim>;
  using const_iterator = typename nodeIndices::const_iterator;
  using const_reverse_iterator = typename nodeIndices::const_reverse_iterator;

  simplex()=delete;
  simplex(UInt elementID_, UInt subelementID_, std::array<UInt, mydim> nodes_) :
    elementID(elementID_), subelementID(subelementID_), nodes(nodes_) {}

  UInt i() const {return elementID;}
  UInt j() const {return subelementID;}
  const UInt& operator[](UInt i) const {return nodes[i];}

  friend bool operator==(const simplex& lhs, const simplex& rhs) {return std::equal(lhs.rbegin(),lhs.rend(), rhs.rbegin());}
  friend bool operator!=(const simplex& lhs, const simplex& rhs) {return !(lhs==rhs);}

  const_iterator begin() const {return nodes.begin();}
  const_iterator end() const {return nodes.end();}
  const_reverse_iterator rbegin() const {return nodes.rbegin();}
  const_reverse_iterator rend() const {return nodes.rend();}

private:
  UInt elementID;
  UInt subelementID;
  nodeIndices nodes;
};

template<UInt mydim>
class simplex_container{
  static_assert(mydim==2 || mydim==3,
    "ERROR! TRYING TO INSTANTIATE SIMPLEX_CONTAINER IN DIMENSION OTHER THAN 2 OR 3! See mesh_input_helper.h");

public:

  using OutputType=std::tuple<std::vector<UInt>, std::vector<bool>, std::vector<bool> >;
  using simplex_t = simplex<mydim>;
  using simplex_container_t = std::vector<simplex_t>;
  using const_iterator = typename simplex_container_t::const_iterator;

  simplex_container()=delete;

  template<std::size_t SIZE>
  simplex_container(const UInt* const elements_, UInt  num_elements_, UInt num_points_, const std::array<UInt, SIZE>& ORDERING) :
      elements(elements_), num_elements(num_elements_), num_points(num_points_) {this->fill_container(elements_, ORDERING);}

  OutputType assemble_output() const;
  std::vector<UInt> get_simplexes() const {return this->assemble_subs();};

  const simplex_t& operator[](UInt i) const {return simplexes[i];}
  const_iterator begin() const {return simplexes.begin();}
  const_iterator end() const {return simplexes.end();}

  bool is_repeated(UInt i) const {return duplicates[i];}

  UInt size() const {return simplexes.size();}
  UInt get_num_points() const {return num_points;}
  UInt get_num_elements() const {return num_elements;}

  std::vector<int> compute_neighbors() const;

private:
  simplex_container_t simplexes;
  std::vector<bool> duplicates;
  std::vector<UInt> distinct_indexes;
  const UInt num_elements;
  const UInt num_points;
  const UInt* const elements;

  template<std::size_t SIZE>
  void fill_container(const UInt* const, const std::array<UInt, SIZE>&);
  std::vector<UInt> compute_offsets(const UInt, std::vector<UInt>&);
  void bin_sort_(const UInt, std::vector<UInt>&);
  void bin_sort();
  void check_duplicates();
  void store_indexes();
  std::vector<bool> mark_boundary() const;
  std::vector<UInt> assemble_subs() const;

};

std::vector<UInt> order2extend(const simplex_container<2> &edge_container){
  std::vector<UInt> edges_extended(edge_container.size());
  UInt offset{edge_container.get_num_points()};
  {
    UInt i=0;
    for(auto const &curr : edge_container){
      offset += !edge_container.is_repeated(i);
      edges_extended[curr.i()+edge_container.get_num_elements()*curr.j()]=offset;
      ++i;
    }
  }
  return edges_extended;
}

std::vector<double> compute_midpoints(const double* const points, const std::vector<UInt>& edges, const UInt num_points){
  const UInt num_edges=edges.size()/2;
  std::vector<double> midpoints(3*num_edges);
  for (int i=0; i<num_edges; ++i)
    for (int j=0; j<3; ++j)
      midpoints[i+j*num_edges]=(points[edges[i]+j*num_points]+points[edges[i+num_edges]+j*num_points])/2;
  return midpoints;
}

// std::vector<UInt> split(simplex_container<2> &edge_container, std::vector<UInt> triangles){
//   std::vector<UInt> extended{order2extend(edge_container)};
//   triangles.insert(triangles.end(), extended.begin(), extended.end());
//
//   const UInt num_edges=triangles.size()/6;
//   std::vector<UInt> splitted_elements;
//   splitted_elements.reserve(12*num_edges);
//
//   for (auto const j : {0,1,2,3})
//     for (int i=0; i<num_edges; ++i)
//       splitted_elements.push_back(triangles[i+j*num_edges]);
//
//   for (auto const j : {5,3,4,4})
//     for (int i=0; i<num_edges; ++i)
//       splitted_elements.push_back(triangles[i+j*num_edges]);
//
//   for (auto const j : {4,5,3,5})
//     for (int i=0; i<num_edges; ++i)
//       splitted_elements.push_back(triangles[i+j*num_edges]);
//
//   return splitted_elements;
// }
//
// std::vector<UInt> split3D(simplex_container<2> &edge_container, std::vector<UInt> tetrahedrons){
//   std::vector<UInt> extended{order2extend(edge_container)};
//   tetrahedrons.insert(tetrahedrons.end(), extended.begin(), extended.end());
//
//   const UInt num_edges=tetrahedrons.size()/10;
//   std::vector<UInt> splitted_elements;
//   splitted_elements.reserve(32*num_edges);
//
//   for (auto const j : {0,4,5,6,4,4,5,5})
//     for (int i=0; i<num_edges; ++i)
//       splitted_elements.push_back(tetrahedrons[i+j*num_edges]);
//
//   for (auto const j : {4,1,7,9,5,5,6,7})
//     for (int i=0; i<num_edges; ++i)
//       splitted_elements.push_back(tetrahedrons[i+j*num_edges]);
//
//   for (auto const j : {5,7,2,8,6,7,9,9})
//     for (int i=0; i<num_edges; ++i)
//       splitted_elements.push_back(tetrahedrons[i+j*num_edges]);
//
//   for (auto const j : {6,9,8,3,9,9,8,8})
//     for (int i=0; i<num_edges; ++i)
//       splitted_elements.push_back(tetrahedrons[i+j*num_edges]);
//
//   return splitted_elements;
// }



#include "mesh_input_helper_imp.h"

#endif
