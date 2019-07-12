#ifndef H_PRUNER
#include "typedefs.hpp"
#endif

#ifndef H_PRUNER_TREEITERATOR_BONES
#define H_PRUNER_TREEITERATOR_BONES

class Tree;

class TreeIterator {
  
private:
  uint current_node;
  uint pos_in_pruning_sequence;
  Tree* tree;
  friend class Tree;
public:
  
  ~TreeIterator() {};
  TreeIterator() {};
  TreeIterator(Tree * tree);
  
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
  
  uint operator*() const {return current_node;};
  
  //! Check whether the current node is root
  bool is_root() const;
  
  //! Checks whether the current node is a leaf or not
  bool is_leaf() const;
  
  uint front() const;
  uint back() const;
  
};


#endif
