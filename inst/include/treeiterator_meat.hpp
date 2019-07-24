#ifndef H_PRUNER
#include "tree_bones.hpp"
#include "treeiterator_bones.hpp"
#endif

inline TreeIterator::TreeIterator(Tree * tree) {
  
  this->current_node = tree->POSTORDER.at(0);
  this->pos_in_pruning_sequence = 0u;
  this->tree = tree;
  
  return;
  
}

inline v_uint::const_iterator TreeIterator::begin_off() const {
  
  return (this->tree)->offspring.at(this->current_node).begin();
  
}

inline v_uint::const_iterator TreeIterator::end_off() const {
 
  return (this->tree)->offspring.at(this->current_node).end();
  
}

inline v_uint::const_iterator TreeIterator::begin_par() const {
  
  return (this->tree)->parents.at(this->current_node).begin();
  
}

inline v_uint::const_iterator TreeIterator::end_par() const {
  
  return (this->tree)->parents.at(this->current_node).end();
  
}

// Return codes:
// 0: At the requested point
// 1: End of road.
inline int TreeIterator::up() {
  
  if (++this->pos_in_pruning_sequence == this->tree->N_NODES) {
    --this->pos_in_pruning_sequence;
    return 1;
  }
  
  this->current_node = this->tree->POSTORDER[this->pos_in_pruning_sequence];
  return 0;
  
}

inline int TreeIterator::down() {
  
  if (this->pos_in_pruning_sequence == 0) {
    return 1;
  }
  
  this->current_node = this->tree->POSTORDER[--this->pos_in_pruning_sequence];
  return 0;
  
}

inline int TreeIterator::operator++() {return this->up();}
inline int TreeIterator::operator--() {return this->up();}

inline void TreeIterator::top() {
  this->current_node = this->tree->POSTORDER[this->tree->POSTORDER.size() - 1u];
  this->pos_in_pruning_sequence = this->tree->POSTORDER.size() - 1u;
  return;
}

inline void TreeIterator::bottom() {
  this->current_node = this->tree->POSTORDER[0u];
  this->pos_in_pruning_sequence = 0;
  return;
}


inline bool TreeIterator::is_root() const {
  return this->tree->parents[current_node].size() == 0u;
}

inline bool TreeIterator::is_leaf() const { 
  return this->tree->offspring[current_node].size() == 0u;
}

inline uint TreeIterator::front() const {
  
  return this->tree->POSTORDER.front();
  
}

inline uint TreeIterator::back() const {
  
  return this->tree->POSTORDER.back();
  
}
