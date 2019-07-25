#ifndef H_PRUNER
#include <memory>
#include "typedefs.hpp"
#include "treeiterator_bones.hpp"
#endif

#ifndef H_PRUNER_TREE_BONES
#define H_PRUNER_TREE_BONES

#ifndef MAX_TREE_SIZE
#define MAX_TREE_SIZE 20000
#endif

// Arbtrary set of arguments, this is the class that the creator function should
// inherit. Ideally it should have:
// - Model parameters can should be modifiable by a function
// - Additional data such as matrices and what not
// We start by declaring it, the user later on defines it

/** Arguments to be passed to Tree::fun.
 * 
 */
class TreeData;
class Tree;
class TreeIterator;


//! Tree class 
/** The Tree class is the core of pruner. The most relevant members are
 * - `parents`
 * - `offspring`
 * - `iterator`
 * - `args`
 * - `fun`
 */
class Tree {
  
private:
  bool is_dag_(int i = -1, int caller = -1, bool up_search = false);
  void postorder_(uint i);
  void postorder();
  TreeIterator iter;
  
protected:
  //! Each nodes' parents.
  vv_uint parents;
  //! Each nodes' offspring.
  vv_uint offspring;
  
  // Auxiliar variables
  //! List of already visited nodes (auxiliar)
  v_bool visited;
  //! Number of visits per node (auxilias)
  v_uint visit_counts;
  
  // Constant
  uint N_NODES;
  uint N_EDGES;
  uint N_TIPS;
  
  //! Postorder sequence
  /**
   * The POSTORDER vector shows the order in which TreeIterator runs through
   * the tree. This is set at the moment during which the tree is created.
   */
  v_uint POSTORDER;
  
  friend class TreeData;
  friend class TreeIterator;
  
public:
  
  //! Arbitrary set of arguments
  sptr_treedata args;
  
  //! Callable function during the the tree traversal
  /**
   * The idea of this function is that the users can specify what to do 
   * during a call to Tree::prune_postorder and Tree::prune_preorder conditional
   * on the current node. The argument of class TreeIterator allows them to
   * get that information by accessing the member function TreeIterator::id.
   */
  std::function<void(sptr_treedata, TreeIterator&)> fun;
  
  //! Evaluates the function by passing the arguments and the iterator
  void eval_fun() {
    
    if (fun)
      fun(this->args, this->iter);
    
    return;
  };
   
  // Creation ------------------------------------------------------------------
  ~Tree() {};
  Tree() {};
  //! Creating method using parents and offpring
  /**
   * The tree is initialized using an edgelist by listing each edges' parent and
   * offspring. The ids of the nodes must range 0 through `n - 1`.
   * 
   * @param parents_ Ids of the parents
   * @param offspring_ Ids of the offspring. Must be of the same length as `parents_`
   * @param out Return codes. 0 means success.
   */
  Tree(const v_uint & parents_, const v_uint & offspring_, uint & out);
  
  // Getter --------------------------------------------------------------------
  
  // As pointers
  const v_uint * get_parents_of(uint i);
  const v_uint * get_offspring_of(uint i);
  
  // As data
  vv_uint get_parents()   const {return this->parents;};
  vv_uint get_offspring() const {return this->offspring;};
  v_uint get_postorder()  const {return this->POSTORDER;};
  v_uint get_preorder()   const;
  
  //! Returns the numner of nodes.
  uint n_nodes()          const {return this->N_NODES;};
  //! Returns the numner of edges.
  uint n_edges()          const {return this->N_EDGES;};
  //! Return the number of tips defined as nodes with no offspring.
  uint n_tips() const;
  
  vv_uint get_edgelist()  const;
  
  void print(bool details = true) const;
  
  // Checker functions ---------------------------------------------------------
  bool is_dag();
  bool is_connected() const {return this->POSTORDER.size() == this->N_NODES;};
  
  // Setters -------------------------------------------------------------------
  void reset_visited() {
    std::fill(this->visited.begin(), this->visited.end(), false);
    std::fill(this->visit_counts.begin(), this->visit_counts.end(), 0u);
    return ;
  };
  
  // Pre-Post/order ------------------------------------------------------------
  
  //! Do the tree-traversal using the postorder
  /** The user defined function `fun` will be called at each step. This allows
   * users implementing their own model.
   */
  void prune_postorder();
  
  //! Do the tree-traversal using the postorder with a user specific sequence
  /** Same as `prune_postorder` but with a different sequence
   * @param seq A new sequence to apply.
   */
  void prune_postorder(v_uint & seq);
  
  //! Do the tree-traversal using the preorder
  /**
   * See Tree::prune_postorder.
   */
  void prune_preorder();
  //! Do the tree-traversal using the preorder with a user specific sequence
  /** Same as `prune_preorder` but with a different sequence
   * @param seq A new sequence to apply.
   */
  void prune_preorder(v_uint & seq);

  
};

/** \example 00-hello-world.cpp This shows how to create and print a Tree class
 * object.
 * 
 * \example 01-adding-function.cpp In this example we add a function to the tree
 * (including a pointer to the arguments). This includes the definition of the
 * class TreeData, which is stored as args in the Tree object.
 * 
 * \example 02-rcpp.cpp Shows how can this library be used within R.
 * 
 */

#endif
