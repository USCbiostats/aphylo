#ifndef H_PRUNER
#include "typedefs.hpp"
#endif

#ifndef H_PRUNER_TREEITERATOR_BONES
#define H_PRUNER_TREEITERATOR_BONES

template <typename Data_Type>
class Tree;

template <typename Data_Type = bool>
class TreeIterator {
  
private:
  
  uint current_node;
  uint pos_in_pruning_sequence;
  
  friend class Tree<Data_Type>;
public:
  
  Tree<Data_Type>* tree;
  
  ~TreeIterator() {};
  TreeIterator() {};
  TreeIterator(Tree<Data_Type> * tree);
  
  //! Begin of offpring `const_iterator` on the current node.
  v_uint::const_iterator begin_off() const;
  //! End of offpring `const_iterator` on the current node.
  v_uint::const_iterator end_off() const;
  
  //! Begin of parents `const_iterator` on the current node.
  v_uint::const_iterator begin_par() const;
  
  //! End of parents `const_iterator` on the current node.
  v_uint::const_iterator end_par() const;
  
  //! Sets the `current_node` to the next value as specified in Tree::POSTORDER.
  int up();
  
  //! Sets the `current_node` to the previous value as specified in Tree::POSTORDER.
  int down();
  
  int operator++();
  int operator--();
  
  void top();
  void bottom();
  
  //! Returns the id (index) of the current node
  uint id() const {return current_node;};
  
  //! Return the number of offsprings the current node has
  int n_offspring() const;
  
  //! Return the number of parents the current node has
  int n_parents() const;
  
  uint operator*() const {return current_node;};
  
  //! Check whether the current node is root
  bool is_root() const;
  
  //! Checks whether the current node is a tip (leaf) or not.
  bool is_tip() const;
  
  uint front() const;
  uint back() const;
  
};


#endif
