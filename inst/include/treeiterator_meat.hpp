#ifndef H_PRUNER
#include "tree_bones.hpp"
#include "treeiterator_bones.hpp"
#endif

template <typename Data_Type>
inline TreeIterator<Data_Type>::TreeIterator(Tree<Data_Type> * tree) {
  
  this->current_node = tree->POSTORDER.at(0);
  this->pos_in_pruning_sequence = 0u;
  this->tree = tree;
  
  return;
  
}

template <typename Data_Type>
inline v_uint::const_iterator TreeIterator<Data_Type>::begin_off() const {
  
  return (this->tree)->offspring.at(this->current_node).begin();
  
}

template <typename Data_Type>
inline v_uint::const_iterator TreeIterator<Data_Type>::end_off() const {
 
  return (this->tree)->offspring.at(this->current_node).end();
  
}

template <typename Data_Type>
inline v_uint::const_iterator TreeIterator<Data_Type>::begin_par() const {
  
  return (this->tree)->parents.at(this->current_node).begin();
  
}

template <typename Data_Type>
inline v_uint::const_iterator TreeIterator<Data_Type>::end_par() const {
  
  return (this->tree)->parents.at(this->current_node).end();
  
}

// Return codes:
// 0: At the requested point
// 1: End of road.
template <typename Data_Type>
inline int TreeIterator<Data_Type>::up() {
  
  if (++this->pos_in_pruning_sequence == this->tree->POSTORDER.size()) {
    --this->pos_in_pruning_sequence;
    return 1;
  }
  
  this->current_node = this->tree->POSTORDER[this->pos_in_pruning_sequence];
  return 0;
  
}

template <typename Data_Type>
inline int TreeIterator<Data_Type>::down() {
  
  if (this->pos_in_pruning_sequence == 0) {
    return 1;
  }
  
  this->current_node = this->tree->POSTORDER[--this->pos_in_pruning_sequence];
  return 0;
  
}

template <typename Data_Type>
inline int TreeIterator<Data_Type>::operator++() {return this->up();}

template <typename Data_Type>
inline int TreeIterator<Data_Type>::operator--() {return this->up();}

template <typename Data_Type>
inline void TreeIterator<Data_Type>::top() {
  this->current_node = this->tree->POSTORDER[this->tree->POSTORDER.size() - 1u];
  this->pos_in_pruning_sequence = this->tree->POSTORDER.size() - 1u;
  return;
}

template <typename Data_Type>
inline void TreeIterator<Data_Type>::bottom() {
  this->current_node = this->tree->POSTORDER[0u];
  this->pos_in_pruning_sequence = 0;
  return;
}

template <typename Data_Type>
inline int TreeIterator<Data_Type>::n_offspring() const {
  return this->tree->n_offspring(this->current_node);
}

template <typename Data_Type>
inline int TreeIterator<Data_Type>::n_parents() const {
  return this->tree->n_parents(this->current_node);
}

template <typename Data_Type>
inline bool TreeIterator<Data_Type>::is_root() const {
  return this->tree->parents[current_node].size() == 0u;
}

template <typename Data_Type>
inline bool TreeIterator<Data_Type>::is_tip() const { 
  return this->tree->offspring[current_node].size() == 0u;
}

template <typename Data_Type>
inline uint TreeIterator<Data_Type>::front() const {
  
  return this->tree->POSTORDER.front();
  
}

template <typename Data_Type>
inline uint TreeIterator<Data_Type>::back() const {
  
  return this->tree->POSTORDER.back();
  
}
