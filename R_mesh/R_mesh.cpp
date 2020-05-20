#include <Rcpp.h>
using namespace Rcpp;

#include <array>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <utility>
#include <type_traits>

using UInt=int;
using OutputType=std::tuple<std::vector<UInt>, std::vector<bool>, std::vector<bool>, std::vector<int> >;

// // Helper code to check is a type is an iterator type
// template <typename... >
// using void_t = void;
//
// template <class T, class = void>
// struct is_iterator : std::false_type {};
//
// template <class T>
// struct is_iterator<T, void_t<
//     typename std::iterator_traits<T>::iterator_category
// > > : std::true_type {};
//
//
// template <typename T, UInt ncols>
// class byRow{
// public:
//   static_assert(is_integral<T>::value, "Error! Value must be of integral type");
//
//   rowView(UInt nrows_, std::vector<value_t> &m_) :
//             nrows(nrows_), m(m_) {}
//
// private:
//   const UInt nrows;
//   std::vector<T> &m;
//   std::array<T,ncols> row;
//
// }

template<UInt mydim>
class simplex{
public:
  using simplex_const_it_t=decltype(std::declval<std::array<UInt, mydim> >().cbegin());
  using simplex_const_rev_it_t=decltype(std::declval<std::array<UInt, mydim> >().crbegin());

  simplex()=delete;
  simplex(UInt element_, UInt subelement_, std::array<UInt, mydim> nodes_) :
    element(element_), subelement(subelement_), nodes(nodes_) {}

  UInt i() const {return element;}
  UInt j() const {return subelement;}
  UInt operator[](UInt i) const {return nodes[i];}

  friend bool operator==(const simplex& lhs, const simplex& rhs) {return std::equal(lhs.rbegin(),lhs.rend(), rhs.rbegin());}
  friend bool operator!=(const simplex& lhs, const simplex& rhs) {return !(lhs==rhs);}

  simplex_const_it_t begin() const {return nodes.begin();}
  simplex_const_it_t end() const {return nodes.end();}
  simplex_const_rev_it_t rbegin() const {return nodes.rbegin();}
  simplex_const_rev_it_t rend() const {return nodes.rend();}

private:
  UInt element;
  UInt subelement;
  std::array<UInt, mydim> nodes;
};

template<UInt mydim>
class simplex_container{
public:
  static_assert(mydim==2 || mydim==3, "Error: this is intended for triangles or tetrahedrons only! See");
  using simplex_t=simplex<mydim>;
  using simplex_const_it_t=decltype(std::declval<std::vector<simplex_t> >().cbegin());

  simplex_container()=delete;

  simplex_container(const UInt* const elements_, UInt  num_elements_, UInt num_points_) :
      elements(elements_), num_elements(num_elements_), num_points(num_points_) {this->fill_container(elements_);}

  simplex_container(const UInt* const elements_, UInt  num_elements_, UInt num_points_, const std::vector<UInt> ORDERING) :
      elements(elements_), num_elements(num_elements_), num_points(num_points_) {this->fill_container(elements_, ORDERING);}

  OutputType assemble_output() const;
  std::vector<UInt> get_simplexes() const {return this->assemble_subs();};

  simplex_t operator[](UInt i) const {return simplexes[i];}
  simplex_const_it_t begin() const {return simplexes.begin();}
  simplex_const_it_t end() const {return simplexes.end();}

  bool is_repeated(UInt i) const {return duplicates[i];}

  UInt size() const {return simplexes.size();}
  UInt get_num_points() const {return num_points;}
  UInt get_num_elements() const {return num_elements;}

private:
  std::vector<simplex_t> simplexes;
  std::vector<bool> duplicates;
  std::vector<UInt> distinct_indexes;
  const UInt num_elements;
  const UInt num_points;
  const UInt* const elements;

  void fill_container(const UInt* const);
  void fill_container(const UInt* const, const std::vector<UInt>);
  std::vector<UInt> compute_offsets(const UInt, std::vector<UInt>&);
  void bin_sort_(const UInt, std::vector<UInt>&);
  void bin_sort();
  void check_duplicates();
  void store_indexes();
  std::vector<bool> mark_boundary() const;
  std::vector<UInt> assemble_subs() const;
  std::vector<int> compute_neighbors() const;
  UInt count_distinct() const;

};

// This function takes care of the initialization of the main container (simplexes)
template<UInt mydim>
void simplex_container<mydim>::fill_container(const UInt* const elements){
 simplexes.reserve((mydim+1)*num_elements);

 {
   std::array<UInt,mydim> curr;
   for(UInt i=0; i<num_elements; ++i){
     for(UInt j=0; j<mydim+1; ++j){
       // Rather ugly but necessary for keeping the right ordering of the edges/faces
       for(UInt k=0; k<mydim; ++k)
          curr[k]=elements[i+num_elements*((j+k+1)%(mydim+1))]-1;
       std::sort(curr.begin(), curr.end());
       simplexes.emplace_back(simplex_t(i,j,curr));
     }
   }
 }

 this->bin_sort();
 this->check_duplicates();
 this->store_indexes();

}

template<UInt mydim>
void simplex_container<mydim>::fill_container(const UInt* const elements, const std::vector<UInt> ORDERING){
 simplexes.reserve(num_elements*ORDERING.size()/mydim);

 {
   std::array<UInt,mydim> curr;
   for(UInt i=0; i<num_elements; ++i){
     for(UInt j=0; j<ORDERING.size()/mydim; ++j){
       for(UInt k=0; k<mydim; ++k)
        curr[k]=elements[i+num_elements*ORDERING[mydim*j+k]]-1;
       std::sort(curr.begin(), curr.end());
       simplexes.emplace_back(simplex_t(i,j,curr));
     }
   }
 }

  this->bin_sort();
  this->check_duplicates();
  this->store_indexes();

}

// This function counts the number of distinct edges/faces in the container
template<UInt mydim>
inline UInt simplex_container<mydim>::count_distinct() const {
  return (distinct_indexes.empty()) ? std::count(duplicates.begin(), duplicates.end(), false) : distinct_indexes.size();
}

template<UInt mydim>
OutputType simplex_container<mydim>::assemble_output() const {

  std::vector<UInt> subsimplexes{this->assemble_subs()};
  std::vector<bool> submarkers{this->mark_boundary()};
  std::vector<int> neighbors{this->compute_neighbors()};

  std::vector<bool> nodesmarkers(num_points);

  for(UInt k=0; k<mydim; ++k)
    for(UInt i=0; i<submarkers.size(); ){
      nodesmarkers[subsimplexes[i+k*submarkers.size()]]=submarkers[i];
      while(i<submarkers.size() && (!submarkers[i] || nodesmarkers[subsimplexes[i+k*submarkers.size()]]))
        ++i;
    }

  return std::make_tuple(std::move(subsimplexes),std::move(submarkers),std::move(nodesmarkers),std::move(neighbors));
}

template<UInt mydim>
std::vector<UInt> simplex_container<mydim>::assemble_subs() const {
  std::vector<UInt> subsimplexes;
  subsimplexes.reserve(mydim*this->count_distinct());

  for(UInt j=0; j<mydim; ++j)
    for(auto const &pos : distinct_indexes)
      subsimplexes.push_back(simplexes[pos][j]);

  return subsimplexes;
}

template<UInt mydim>
std::vector<bool> simplex_container<mydim>::mark_boundary() const {
  std::vector<bool> boundarymarkers;
  boundarymarkers.reserve(this->count_distinct());

  std::for_each(distinct_indexes.cbegin(), std::prev(distinct_indexes.cend()), [&] (UInt i) {
    boundarymarkers.push_back(!duplicates[i+1]);
  });

  boundarymarkers.push_back(distinct_indexes.back()+1==duplicates.size() || !duplicates[distinct_indexes.back()+1]);
  return boundarymarkers;
}

template<UInt mydim>
std::vector<int> simplex_container<mydim>::compute_neighbors() const {
  std::vector<int> neighbors(simplexes.size(), -1);

  auto rep_it=duplicates.cbegin();
  simplex_t prev{simplexes.front()};
  for (auto const &curr : simplexes){
    // Note: the first simplex cannot be a duplicate!
    if (*(rep_it++)){
      neighbors[curr.i()+curr.j()*num_elements]=prev.i();
      neighbors[prev.i()+prev.j()*num_elements]=curr.i();
    }
    prev=curr;
  }
  return neighbors;
}


template<UInt mydim>
void simplex_container<mydim>::bin_sort(){
  std::vector<UInt> positions;
  positions.reserve(simplexes.size());
  for(UInt i=0; i<simplexes.size(); ++i)
    positions.push_back(i);
  bin_sort_(mydim-1, positions);

  for(UInt i=0; i<positions.size(); ++i){
    UInt curr=i;
    while(i!=positions[curr]){
      UInt next=positions[curr];
      std::swap(simplexes[curr],simplexes[next]);
      positions[curr]=curr;
      curr=next;
    }
    positions[curr]=curr;
  }
}

// Recursive unction to sort container by ascending #(index+1) element of the arrays
template<UInt mydim>
void simplex_container<mydim>::bin_sort_(const UInt index, std::vector<UInt> &positions){
  // Note the scoping to avoid unnecessary storage!
  {
    std::vector<UInt> offsets{compute_offsets(index, positions)};
    for(UInt i=0; i<positions.size(); ++i){
      while(i!=offsets[i]){
        UInt next=offsets[i];
        std::swap(positions[i],positions[next]);
        std::swap(offsets[i],offsets[next]);
      }
    }
  }

  if(index)
    bin_sort_(index-1, positions);
}

template<UInt mydim>
std::vector<UInt> simplex_container<mydim>::compute_offsets(const UInt index, std::vector<UInt> &positions){
  std::vector<UInt> counts(num_points);
  for(auto const &pos : positions)
    ++counts[simplexes[pos][index]];

  UInt offset{0};
  for (auto &curr : counts){
    UInt count{curr};
    curr=offset;
    offset+=count;
  }

  std::vector<UInt> offsets;
  offsets.reserve(positions.size());
  for (auto const &pos : positions)
    offsets.push_back(counts[simplexes[pos][index]]++);
  return offsets;

}

template<UInt mydim>
void simplex_container<mydim>::check_duplicates(){
  duplicates.reserve(simplexes.size());
  // First face/edge cannot be a duplicate!
  duplicates.push_back(false);
  for(auto it=std::next(simplexes.cbegin()); it!=simplexes.cend(); ++it)
    duplicates.push_back(*std::prev(it) == *it);
}

template<UInt mydim>
void simplex_container<mydim>::store_indexes(){
  distinct_indexes.reserve(this->count_distinct());
  for(UInt i=0; i<duplicates.size(); ++i)
    if(!duplicates[i])
      distinct_indexes.push_back(i);
}

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

std::vector<UInt> split(simplex_container<2> &edge_container, std::vector<UInt> triangles){
  std::vector<UInt> extended{order2extend(edge_container)};
  triangles.insert(triangles.end(), extended.begin(), extended.end());

  const UInt num_edges=triangles.size()/6;
  std::vector<UInt> splitted_elements;
  splitted_elements.reserve(12*num_edges);

  for (auto const j : {0,1,2,3})
    for (int i=0; i<num_edges; ++i)
      splitted_elements.push_back(triangles[i+j*num_edges]);

  for (auto const j : {5,3,4,4})
    for (int i=0; i<num_edges; ++i)
      splitted_elements.push_back(triangles[i+j*num_edges]);

  for (auto const j : {4,5,3,5})
    for (int i=0; i<num_edges; ++i)
      splitted_elements.push_back(triangles[i+j*num_edges]);

  return splitted_elements;
}

std::vector<UInt> split3D(simplex_container<2> &edge_container, std::vector<UInt> tetrahedrons){
  std::vector<UInt> extended{order2extend(edge_container)};
  tetrahedrons.insert(tetrahedrons.end(), extended.begin(), extended.end());

  const UInt num_edges=tetrahedrons.size()/10;
  std::vector<UInt> splitted_elements;
  splitted_elements.reserve(32*num_edges);

  for (auto const j : {0,4,5,6,4,4,5,5})
    for (int i=0; i<num_edges; ++i)
      splitted_elements.push_back(tetrahedrons[i+j*num_edges]);

  for (auto const j : {4,1,7,9,5,5,6,7})
    for (int i=0; i<num_edges; ++i)
      splitted_elements.push_back(tetrahedrons[i+j*num_edges]);

  for (auto const j : {5,7,2,8,6,7,9,9})
    for (int i=0; i<num_edges; ++i)
      splitted_elements.push_back(tetrahedrons[i+j*num_edges]);

  for (auto const j : {6,9,8,3,9,9,8,8})
    for (int i=0; i<num_edges; ++i)
      splitted_elements.push_back(tetrahedrons[i+j*num_edges]);

  return splitted_elements;
}


// [[Rcpp::export]]
IntegerVector R_split_triangles(IntegerVector triangles, UInt num_points){

  const UInt n=triangles.length(), num_triangles=n/3;

  std::vector<UInt> CPP_triangles;
  CPP_triangles.reserve(n);
  for(UInt i=0; i<n; ++i)
    CPP_triangles.push_back(triangles[i]);

  // edges will contain arrays of this form: (node1, node2, triangleID, edgeID)
  // where triangleID tells which triangle the edge belongs to and edgeID tells
  // which node is in front of the edge (eg 1 if edge (2,3) )
  simplex_container<2> edges_list(&triangles[0],num_triangles,num_points);

  std::vector<UInt> splitted=split(edges_list,CPP_triangles);

  IntegerVector R_splitted(splitted.size());
  for(UInt i=0; i<splitted.size(); ++i)
    R_splitted[i]=splitted[i];

  return R_splitted;
}

// [[Rcpp::export]]
IntegerVector R_split_tetrahedrons(IntegerVector tetrahedrons, UInt num_points){

  const UInt n=tetrahedrons.length(), num_tetrahedrons=n/4;

  std::vector<UInt> CPP_tetrahedrons;
  CPP_tetrahedrons.reserve(n);
  for(UInt i=0; i<n; ++i)
    CPP_tetrahedrons.push_back(tetrahedrons[i]);


  // edges will contain arrays of this form: (node1, node2, tetrahedronID, edgeID)
  // where tetrahedronID tells which tetrahedron the edge belongs to and edgeID tells
  // which node is in front of the edge (eg 1 if edge (2,3) )
  simplex_container<2> edges_list(&tetrahedrons[0],num_tetrahedrons,num_points,{0,1,0,2,0,3,1,2,2,3,1,3});

  std::vector<UInt> splitted=split3D(edges_list,CPP_tetrahedrons);

  IntegerVector R_splitted(splitted.size());
  for(UInt i=0; i<splitted.size(); ++i)
    R_splitted[i]=splitted[i];

  return R_splitted;

}





// This function will be called by R create.mesh functions (see mesh.R script)
// It will compute edges, mark the ones on the boundary
// and possibly compute midpoints if asked to.
// triangles is a vector obtained unrolling a #triangles-by-3 matrix


// [[Rcpp::export]]
List R_mesh_helper_2_5D(IntegerVector triangles, NumericVector nodes, UInt num_points, bool extend=false){
  // triangle.size() is a multiple of 3 by construction!
  const UInt n=triangles.length(), num_triangles=n/3;

  // edges will contain arrays of this form: (node1, node2, triangleID, edgeID)
  // where triangleID tells which triangle the edge belongs to and edgeID tells
  // which node is in front of the edge (eg 1 if edge (2,3) )
  simplex_container<2> edges_list(&triangles[0],num_triangles,num_points);

  OutputType out{edges_list.assemble_output()};

  std::vector<UInt> &edges=std::get<0>(out);
  std::vector<bool> &edgesmarkers=std::get<1>(out);
  std::vector<bool> &nodesmarkers=std::get<2>(out);
  std::vector<int> &neighbors=std::get<3>(out);


  NumericVector R_edges(edges.size());
  LogicalVector R_edgesmarkers(edgesmarkers.size());
  NumericVector R_neighbors(neighbors.size());


  for (UInt i=0; i<edges.size(); ++i)
    R_edges[i]=edges[i]+1;

  for (UInt i=0; i<edgesmarkers.size(); ++i)
    R_edgesmarkers[i]=edgesmarkers[i];

  for (UInt i=0; i<neighbors.size(); ++i)
    R_neighbors[i]=neighbors[i]+(neighbors[i]>=0);

  if(extend){
    std::vector<UInt> extended{order2extend(edges_list)};
    IntegerVector R_extended(n+extended.size());
    for(UInt i=0; i<n; ++i)
      R_extended[i]=triangles[i];
    for(UInt i=0; i<extended.size(); ++i)
      R_extended[i+n]=extended[i];

    std::vector<double> midpoints{compute_midpoints(&nodes[0], edges, num_points)};
    NumericVector R_nodes(3*num_points+midpoints.size());
    for(UInt i=0; i<3*num_points; ++i)
      R_nodes[i]=nodes[i];
    for(UInt i=0; i<midpoints.size(); ++i)
      R_nodes[i+3*num_points]=midpoints[i];

    LogicalVector R_nodesmarkers(nodesmarkers.size()+midpoints.size()/3);
    for (UInt i=0; i<nodesmarkers.size(); ++i)
      R_nodesmarkers[i]=nodesmarkers[i];

    List L = List::create(Named("edges")=R_edges, _["edgesmarkers"]=R_edgesmarkers, _["nodesmarkers"]=R_nodesmarkers, _["neighbors"]=R_neighbors, _["triangles"]=R_extended, _["nodes"]=R_nodes);
    return L;
  }

  LogicalVector R_nodesmarkers(nodesmarkers.size());
  for (UInt i=0; i<nodesmarkers.size(); ++i)
    R_nodesmarkers[i]=nodesmarkers[i];

  List L = List::create(Named("edges")=R_edges, _["edgesmarkers"]=R_edgesmarkers, _["nodesmarkers"]=R_nodesmarkers, _["neighbors"]=R_neighbors);
  return L;

}

// This function will be called by R create.mesh functions (see mesh.R script)
// It will compute edges, mark the ones on the boundary
// and possibly compute midpoints if asked to.
// tetrahedrons is a vector obtained unrolling a #tetrahedrons-by-3 matrix


// [[Rcpp::export]]
List R_mesh_helper_3D(IntegerVector tetrahedrons, UInt num_points, bool extend=false){
  // tetrahedron.size() is a multiple of 3 by construction!
  const UInt n=tetrahedrons.length(), num_tetrahedrons=n/4;

  // edges will contain arrays of this form: (node1, node2, tetrahedronID, edgeID)
  // where tetrahedronID tells which tetrahedron the edge belongs to and edgeID tells
  // which node is in front of the edge (eg 1 if edge (2,3) )
  simplex_container<3> faces_list(&tetrahedrons[0],num_tetrahedrons,num_points);

  OutputType out{faces_list.assemble_output()};

  std::vector<UInt> &faces=std::get<0>(out);
  std::vector<bool> &facesmarkers=std::get<1>(out);
  std::vector<bool> &nodesmarkers=std::get<2>(out);
  std::vector<int> &neighbors=std::get<3>(out);

  LogicalVector R_nodesmarkers(nodesmarkers.size());
  NumericVector R_faces(faces.size());
  LogicalVector R_facesmarkers(facesmarkers.size());
  NumericVector R_neighbors(neighbors.size());

  for (UInt i=0; i<nodesmarkers.size(); ++i)
    R_nodesmarkers[i]=nodesmarkers[i];

  for (UInt i=0; i<faces.size(); ++i)
    R_faces[i]=faces[i]+1;

  for (UInt i=0; i<facesmarkers.size(); ++i)
    R_facesmarkers[i]=facesmarkers[i];

  for (UInt i=0; i<neighbors.size(); ++i)
    R_neighbors[i]=neighbors[i]+(neighbors[i]>=0);

  if(extend){
    simplex_container<2> edges_list(&tetrahedrons[0],num_tetrahedrons,num_points,{0,1,0,2,0,3,1,2,2,3,1,3});
    std::vector<UInt> edges{edges_list.get_simplexes()};
    NumericVector R_edges(edges.size());
    for (UInt i=0; i<edges.size(); ++i)
      R_edges[i]=edges[i]+1;

    std::vector<UInt> extended{order2extend(edges_list)};
    NumericVector R_extended(extended.size());
    for(UInt i=0; i<extended.size(); ++i)
      R_extended[i]=extended[i];

    List L = List::create(_["faces"]=R_faces, _["facesmarkers"]=R_facesmarkers, _["nodesmarkers"]=R_nodesmarkers, _["neighbors"]=R_neighbors,_["edges"]=R_edges, _["extra"]=R_extended);
    return L;
  }

  List L = List::create(_["nodesmarkers"]=R_nodesmarkers, _["faces"]=R_faces, _["facesmarkers"]=R_facesmarkers, _["neighbors"]=R_neighbors);
  return L;

}
