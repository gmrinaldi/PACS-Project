#include <Rcpp.h>
using namespace Rcpp;

#include <array>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <utility>

using UInt=int;
using OutputType=std::tuple<std::vector<UInt>, std::vector<bool>, std::vector<bool>, std::vector<int> >;

template <typename T>
class bin_list{
public:
  bin_list(const UInt size) {bins.resize(size);}
  std::list<T>& operator[](UInt i) {return bins[i];}

  template <typename Iterator>
  void unlist_into(Iterator);
private:
  std::vector<std::list<T> > bins;
};

template <typename T>
template <typename Iterator>
void bin_list<T>::unlist_into(Iterator it){
  // Note the use of it++! (returns unincremented it)
  for (auto const &bin : bins)
    for(auto const &i : bin)
      *(it++)=i;
}


template<UInt mydim>
class simplex_container{
public:
  static_assert(mydim==2 || mydim==3, "Error: this is intended for triangles or tetrahedrons only! See");
  using IntArray=std::array<UInt,mydim+2>;
  using SimplexConstIt=decltype(std::declval<std::vector<IntArray> >().cbegin());

  simplex_container()=delete;

  simplex_container(const UInt* const elements_, UInt  num_elements_, UInt num_points_) :
      elements(elements_), num_elements(num_elements_), num_points(num_points_) {this->fill_container(elements_);}

  simplex_container(const UInt* const elements_, UInt  num_elements_, UInt num_points_, const std::vector<UInt> ORDERING) :
      elements(elements_), num_elements(num_elements_), num_points(num_points_) {this->fill_container(elements_, ORDERING);}

  OutputType assemble_output() const;
  std::vector<int> get_simplexes() const {return this->assemble_subs();};

  IntArray operator[](UInt i) const {return simplexes[i];}
  SimplexConstIt begin() const {return simplexes.begin();}
  SimplexConstIt end() const {return simplexes.end();}

  bool is_repeated(UInt i) const {return duplicates[i];}

  UInt size() const {return simplexes.size();}
  UInt get_num_points() const {return num_points;}
  UInt get_num_elements() const {return num_elements;}

protected:
  std::vector<IntArray> simplexes;
  std::vector<bool> duplicates;
  std::vector<UInt> distinct_indexes;
  const UInt num_elements;
  const UInt num_points;
  const UInt* const elements;

  void fill_container(const UInt* const);
  void fill_container(const UInt* const, const std::vector<UInt>);
  void sort_by_index(const UInt);
  void bin_sort();
  bool simplex_cmp(const IntArray &, const IntArray &);
  void check_duplicates();
  void store_indexes();
  std::vector<bool> mark_boundary() const;
  std::vector<UInt> assemble_subs() const;
  std::vector<int> compute_neighbors() const;
  UInt count_distinct() const;

};


template<UInt mydim>
bool simplex_container<mydim>::simplex_cmp(const IntArray &lhs, const IntArray &rhs){
  for(UInt i=mydim+1; i>1; --i)
    if(lhs[i]!=rhs[i])
      return false;
  return true;
}

template<UInt mydim>
void simplex_container<mydim>::fill_container(const UInt* const elements){
 simplexes.resize((mydim+1)*num_elements);
 bin_list<IntArray> bins(num_points);
 IntArray simplex;
  for(UInt i=0; i<num_elements; ++i){
    for(UInt j=0; j<mydim+1; ++j){
      simplex={i,j};
      for(UInt k=0; k<mydim; ++k)
        simplex[k+2]=elements[i+num_elements*((j+k+1)%(mydim+1))]-1;
      std::sort(std::next(simplex.begin(),2),simplex.end());
      bins[simplex.back()].push_back(simplex);
    }
  }
  bins.unlist_into(simplexes.begin());

  this->bin_sort();
  this->check_duplicates();
  this->store_indexes();
}

template<UInt mydim>
void simplex_container<mydim>::fill_container(const UInt* const elements, const std::vector<UInt> ORDERING){
 simplexes.resize(num_elements*ORDERING.size()/mydim);
 bin_list<IntArray> bins(num_points);
 IntArray simplex;
  for(UInt i=0; i<num_elements; ++i){
    for(UInt j=0; j<ORDERING.size()/mydim; ++j){
      simplex={i,j};
      for(UInt k=0; k<mydim; ++k)
        simplex[k+2]=elements[i+num_elements*ORDERING[mydim*j+k]]-1;
      std::sort(std::next(simplex.begin(),2),simplex.end());
      bins[simplex.back()].push_back(simplex);
    }
  }
  bins.unlist_into(simplexes.begin());

  this->bin_sort();
  this->check_duplicates();
  this->store_indexes();
}


template<UInt mydim>
inline UInt simplex_container<mydim>::count_distinct() const {
  if (!distinct_indexes.empty())
    return distinct_indexes.size();
  UInt count=0;
  for (auto const &rep : duplicates)
    count+= !rep;
  return count;
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
  std::vector<UInt> subsimplexes(mydim*this->count_distinct());

  auto curr=subsimplexes.begin();
  for(UInt j=2; j<mydim+2; ++j){
    auto rep_it=duplicates.cbegin();
    for(auto const &pos : distinct_indexes)
      *(curr++)=simplexes[pos][j];
  }
  return subsimplexes;
}

template<UInt mydim>
std::vector<bool> simplex_container<mydim>::mark_boundary() const {
  std::vector<bool> boundarymarkers(this->count_distinct());

  auto boundary=boundarymarkers.begin();
  for (auto ind_it=distinct_indexes.cbegin(); std::next(ind_it)!=distinct_indexes.cend(); ++ind_it)
    *(boundary++) = !(duplicates[*ind_it + 1]);

  boundarymarkers.back() = (distinct_indexes.back()+1<duplicates.size())? !duplicates[distinct_indexes.back()+1] : true;
  return boundarymarkers;
}

template<UInt mydim>
std::vector<int> simplex_container<mydim>::compute_neighbors() const {
  std::vector<int> neighbors(simplexes.size(), -1);

  auto rep_it=duplicates.cbegin();
  IntArray prev{simplexes.front()};
  for (auto const &curr : simplexes){
    // Note: the first simplex cannot be a duplicate!
    if (*(rep_it++)){
      neighbors[curr[0]+curr[1]*num_elements]=prev[0];
      neighbors[prev[0]+prev[1]*num_elements]=curr[0];
    }
    prev=curr;
  }
  return neighbors;
}


// Function to sort container by ascending #(index+1) element of the arrays
template<UInt mydim>
void simplex_container<mydim>::sort_by_index(const UInt index){

  bin_list<IntArray> bins(num_points);
  for(auto const &curr : simplexes)
    bins[curr[index]].push_back(curr);

  bins.unlist_into(simplexes.begin());
}

template<UInt mydim>
void simplex_container<mydim>::bin_sort(){
  for(UInt i=mydim; i>1; --i)
    sort_by_index(i);
}

template<UInt mydim>
void simplex_container<mydim>::check_duplicates(){
  duplicates.resize(simplexes.size(),false);
  auto it=simplexes.cbegin();
  for(auto rep_it=std::next(duplicates.begin()); rep_it!=duplicates.end(); ++rep_it, ++it)
      *rep_it = simplex_cmp(*std::next(it), *it);
}

template<UInt mydim>
void simplex_container<mydim>::store_indexes(){
  distinct_indexes.resize(this->count_distinct());
  auto curr = distinct_indexes.begin();
  for(UInt i=0; i<duplicates.size(); ++i)
    if(!duplicates[i])
      *(curr++)=i;
}

std::vector<UInt> midpoints(simplex_container<2> &edge_container){
  std::vector<UInt> extended(edge_container.size());
  UInt offset{edge_container.get_num_points()};
  {
    UInt i=0;
    for(auto const &curr : edge_container){
      offset+= !edge_container.is_repeated(i);
      extended[curr[0]+edge_container.get_num_elements()*curr[1]]=offset;
      ++i;
    }
  }
  return extended;
}


// This function will be called by R create.mesh functions (see mesh.R script)
// It will compute edges, mark the ones on the boundary
// and possibly compute midpoints if asked to.
// triangles is a vector obtained unrolling a #triangles-by-3 matrix


// [[Rcpp::export]]
List R_mesh_helper_2_5D(IntegerVector triangles, UInt num_points, bool extend=false){
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
  LogicalVector R_nodesmarkers(nodesmarkers.size());
  NumericVector R_neighbors(neighbors.size());


  for (UInt i=0; i<edges.size(); ++i)
    R_edges[i]=edges[i]+1;

  for (UInt i=0; i<edgesmarkers.size(); ++i)
    R_edgesmarkers[i]=edgesmarkers[i];

  for (UInt i=0; i<nodesmarkers.size(); ++i)
    R_nodesmarkers[i]=nodesmarkers[i];

  for (UInt i=0; i<neighbors.size(); ++i)
    R_neighbors[i]=neighbors[i]+(neighbors[i]>=0);

  if(extend){
    std::vector<UInt> extended{midpoints(edges_list)};
    NumericVector R_extended(n);
    for(UInt i=0; i<n; ++i)
      R_extended[i]=extended[i];

    List L = List::create(Named("edges")=R_edges, _["edgesmarkers"]=R_edgesmarkers, _["nodesmarkers"]=R_nodesmarkers, _["neighbors"]=R_neighbors, _["extra"]=R_extended);
    return L;
  }

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

    std::vector<UInt> extended{midpoints(edges_list)};
    NumericVector R_extended(extended.size());
    for(UInt i=0; i<extended.size(); ++i)
      R_extended[i]=extended[i];

    List L = List::create(_["faces"]=R_faces, _["facesmarkers"]=R_facesmarkers, _["nodesmarkers"]=R_nodesmarkers, _["neighbors"]=R_neighbors,_["edges"]=R_edges, _["extra"]=R_extended);
    return L;
  }

  List L = List::create(_["nodesmarkers"]=R_nodesmarkers, _["faces"]=R_faces, _["facesmarkers"]=R_facesmarkers, _["neighbors"]=R_neighbors);
  return L;

}
