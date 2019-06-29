#include <Rcpp.h>
using namespace Rcpp;

typedef unsigned int uint;
typedef std::vector< uint > v_uint;

class aphylo_node;
class aphylo_edge;
class aphylo_edgelist;

// Basic structure: NODE -------------------------------------------------------
class aphylo_node {
  
private:
  uint id;
  uint type;
  std::string label;
  friend class aphylo_edge;
  friend class aphylo_edgelist;

public:
  ~aphylo_node() {};
  aphylo_node() {};
  // Generic init
  aphylo_node(uint id_, uint type_, std::string label_):
    id(id_), type(type_), label(label_) {};
  
  aphylo_node(uint id_):
    id(id_), type(0u), label("") {};
  
  // Copy constructor
  aphylo_node(const aphylo_node &n):
    id(n.id), type(n.type), label(n.label) {};
  
  // Assignmen operator
  aphylo_node& operator= (const aphylo_node& n) {
    
    id    = n.id;
    type  = n.type;
    label = n.label;
    
    return *this;
  };
  
  void print();
  
};

inline void aphylo_node::print() {
  Rprintf("aphylo_node: id(%d) type(%d) label(%s)\n",
          id, type, label);
}

// Basic structure: EDGE -------------------------------------------------------
class aphylo_edge {
  
private:
  aphylo_node source;
  aphylo_node target;
  uint type;
  friend class aphylo_edgelist;
  
public:
  ~aphylo_edge() {};
  
  aphylo_edge(aphylo_node source_, aphylo_node target_, uint type_):
    source(source_), target(target_), type(type_) {};
  
  aphylo_edge(uint source_, uint target_);
  
  // void print();
};

inline aphylo_edge::aphylo_edge(uint source_, uint target_) {

  aphylo_node s(source_);
  aphylo_node t(target_);
  source = s;
  target = t;
  type   = 0u;

}

// inline aphylo_edge::print() {
//   
//   Rprintf("")
//   
// }

// More complex structure: EDGELIST --------------------------------------------
class aphylo_edgelist {

private:
  std::vector< aphylo_edge > E;
  
public:
  ~aphylo_edgelist() {};
  aphylo_edgelist() {};
  aphylo_edgelist(const v_uint& source, const v_uint& target);
  
  void add_edge(aphylo_edge e) {
    E.push_back(e);
    return;
  };
  
  void add_edge(uint source, uint target) {
    // Should I use smart pointers here?
    aphylo_edge x(source, target);
    E.push_back(x);
    return;
  };
  
  void print();
  
};

inline aphylo_edgelist::aphylo_edgelist(
    const v_uint& source,
    const v_uint& target
  ) {
  
  int nedges = source.size();
  
  for (int i = 0; i < nedges; ++i)
    this->add_edge(source[i], target[i]);

}

inline void aphylo_edgelist::print() {
  
  Rprintf("This edgelist has the following contents:\n");
  for (int i = 0; i < E.size(); ++i)
    Rprintf("[%4d, %4d]\n", E[i].source.id, E[i].target.id);
  
  return;
}

// inline aphylo_node::aphylo_node(v_uint source, v_uint target) {
//   
// }


// [[Rcpp::export]]
int create_and_delete(uint id, uint type, std::string label) {
  
  aphylo_node x(id, type, label);
  x.print();
  
  aphylo_node y(id);
  y.print();
  
  return 1;
  
}

// [[Rcpp::export]]
int create2(std::vector< uint > source, std::vector< uint > target) {
  
  aphylo_edgelist E(source, target);
  E.print();
  
  return 1;
  
}

/*** R
# timesTwo(42)
create_and_delete(1, 2, "George")
x <- ape::rtree(10)
create2(x$edge[,1], x$edge[,2])
*/
