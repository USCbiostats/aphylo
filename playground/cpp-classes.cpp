#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

// Rcpp::plugins(cpp11)

typedef unsigned int uint;
typedef std::vector< uint > v_uint;

template <class T>
class NodeAnnotated;

class aphylo_edge;
class aphylo_tree;

// Basic structure: NODE -------------------------------------------------------
class Node {
private:
  uint id;
  template <class T>
  friend class NodeAnnotated;
  friend class aphylo_edge;
  friend class aphylo_tree;
public:
  ~Node() {};
  Node() {};
  Node(uint id_):
    id(id_) {};
  
  void print();
  
};

inline void Node::print() {
  Rprintf("a node with id %~04d.\n", id);
};

template < class T >
class NodeAnnotated : public Node {
  
private:
  std::vector< T > annotation;

public:
  ~NodeAnnotated() {};
  NodeAnnotated() {};
  // Generic init
  NodeAnnotated(uint id_, std::vector< T > annotation_):
    Node(id_), annotation(annotation_) {};
  
  // If a single annotation
  NodeAnnotated(uint id_, T annotation_): Node(id_) {
    annotation.push_back(annotation_);
  };
  
  void print(){
    Rprintf("id: %~04d, N annotations: %d) \n", id, annotation.size());
  };
  
};



// Basic structure: EDGE -------------------------------------------------------
class aphylo_edge {
  
private:
  int source;
  int target;
  uint type;
  friend class aphylo_tree;
  
public:
  ~aphylo_edge() {};
  
  aphylo_edge(int source_, int target_, uint type_):
    source(source_), target(target_), type(type_) {};
  
  aphylo_edge(uint source_, uint target_);
  
  // void print();
};

inline aphylo_edge::aphylo_edge(uint source_, uint target_) {

  int s(source_);
  int t(target_);
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
typedef std::unordered_map< unsigned int, Node > Nodes;
typedef std::vector< aphylo_edge > aphylo_edges;
class aphylo_tree {

private:
  Nodes nodes;
  aphylo_edges edges;
  
public:
  ~aphylo_tree() {};
  aphylo_tree() {};
  aphylo_tree(const v_uint& source, const v_uint& target);
  
  void add_edge(uint source, uint target) {
    
    // Checking if it exists
    if (!this->nodes.count(source)) 
      nodes[source] = *(new Node(source));
    
    if (!this->nodes.count(target)) {
      nodes[target] = *(new Node(target));
    }
    
    // Should I use smart pointers here?
    aphylo_edge x(source, target);
    edges.push_back(x);
    return;
  };
  
  void print();
  void print_nodes();

};



inline aphylo_tree::aphylo_tree(
    const v_uint& source,
    const v_uint& target
  ) {
  
  int nedges = source.size();
  
  for (int i = 0; i < nedges; ++i)
    this->add_edge(source[i], target[i]);

}

inline void aphylo_tree::print_nodes() {
  
  Nodes::iterator n;
  for (n = nodes.begin(); n != nodes.end(); ++n)
    n->second.print();
  
  return;
}

inline void aphylo_tree::print() {
  
  Rprintf("This edgelist has the following contents:\n");
  aphylo_edges::const_iterator e;
  for (e = edges.begin(); e != edges.end(); ++e)
    Rprintf("[%4d, %4d]\n", e->source, e->target);
  
  return;
}

// inline aphylo_node::aphylo_node(v_uint source, v_uint target) {
//   
// }


// [[Rcpp::export]]
int create_and_delete(uint id, uint type, std::string label) {
  
  Node x(id);
  x.print();
  
  NodeAnnotated< uint > y(id, type);
  y.print();
  
  return 1;
  
}

// [[Rcpp::export]]
int create2(std::vector< uint > source, std::vector< uint > target) {
  
  aphylo_tree E(source, target);
  E.print();
  E.print_nodes();
  
  return 1;
  
}

/*** R
# timesTwo(42)
create_and_delete(1, 2, "George")
set.seed(1);x <- ape::rtree(5)
create2(x$edge[,1], x$edge[,2])
*/
