#ifndef H_PRUNER
#include <vector>
#include "typedefs.hpp"
#include "tree_bones.hpp"
#include "treeiterator_bones.hpp"
#include "treeiterator_meat.hpp"
// #include <memory> // For the std::shared_ptr
#endif

#ifndef H_PRUNER_TREE_MEAT
#define H_PRUNER_TREE_MEAT

// Pruning ---------------------------------------------------------------------
 
// Function to get the pre and post order
template <typename Data_Type>
inline void Tree<Data_Type>::postorder() {
  
  // If it is not initialized
  if (this->POSTORDER.size() == 0u)
    POSTORDER.reserve(this->N_NODES);
  
  // We can start the algorithm from any place
  postorder_(0u);
  POSTORDER.shrink_to_fit();
  
  this->reset_visited();
  
  return;
}

template <typename Data_Type>
inline void Tree<Data_Type>::postorder_(uint i) {
  
  // Setting the state to visited
  this->visited[i] = true;
  
  // First, check the offspring
  for (uint j = 0u; j < this->offspring[i].size(); ++j) {
    
    // Nothing to do here
    if (this->visited[this->offspring[i][j]])
      continue;
    
    postorder_(this->offspring[i][j]);
    
  }
  
  // After visiting all of its childs, we need to add this node to the pruning
  // sequence and continue with its parent(s).
  POSTORDER.push_back(i);
  for (uint j = 0u; j < this->parents[i].size(); ++j) {
    
    // Nothing to do here
    if (this->visited[this->parents[i][j]])
      continue;
    
    postorder_(this->parents[i][j]);
    
  }
  
  return;
}

// Returns an edgelist in the form a vector of two uint vectors.
template <typename Data_Type>
inline vv_uint Tree<Data_Type>::get_edgelist() const {
  
  vv_uint res(2u);
  
  // We know the number of edges before hand, so we better save the space
  // up-front.
  res[0u].reserve(this->N_EDGES);
  res[1u].reserve(this->N_EDGES);
  
  for (uint i = 0u; i < this->N_NODES; ++i) {
    for (uint j = 0u; j < this->offspring[i].size(); ++j) {
      res[0u].push_back(i);
      res[1u].push_back(this->offspring[i][j]);
    }
  }
  
  return res;
  
}

template <typename Data_Type>
inline void Tree<Data_Type>::print(bool details) const {
  
  // Basic information
  printf(
    "tree structure with %i nodes and %i edges.\n",
    this->n_nodes(), this->n_edges()
  );
  
  // When details are true, information regarding the offspring and parents
  // is shown.
  if (!details)
    return;
  
  // What is sais below (just divide the sections :P)
  printf("List of offspring:\n");
  for (uint i = 0u; i < this->offspring.size(); ++i) {
    
    printf("%4i: [", i);
    
    if (this->offspring[i].size() == 0u)
      printf(" -");
    else {
      v_uint::const_iterator iter;
      for (iter = this->offspring[i].begin(); iter != this->offspring[i].end(); ++iter)
        printf(" %i", *iter);
    }
    
    printf(" ]\n");
    
  }
  
  // What is sais below (just divide the sections :P)
  printf("List of parents:\n");
  for (uint i = 0u; i < this->parents.size(); ++i) {
    
    printf("%4i: [", i);
    
    if (this->parents[i].size() == 0u)
      printf(" -");
    else {
      v_uint::const_iterator iter;
      for (iter = this->parents[i].begin(); iter != this->parents[i].end(); ++iter)
        printf(" %i", *iter);
    }
    
    printf(" ]\n");
    
  }
  
  return;
  
}

// Return codes:
// 0: OK
// 1: Sizes of parent and offspring differ
// 2: MAX_TREE_SIZE reached.
// 3: Disconnected tree
// 4: Not a Dag.
template <typename Data_Type>
inline Tree<Data_Type>::Tree(const v_uint & parents_, const v_uint & offspring_, uint & out) {
  
  // If different sizes, then it shouldn't work
  if (parents_.size() != offspring_.size()) {
    out = 1u;
    return;
  }
  
  // Checking ranges
  uint maxid = 0u, m = parents_.size();
  for (uint i = 0u; i < m; ++i) {
    if ((parents_[i] > MAX_TREE_SIZE) || (offspring_[i] > MAX_TREE_SIZE)) {
      out = 2u;
      return;
    }
    
    if (maxid < parents_[i])
      maxid = parents_[i];
    if (maxid < offspring_[i])
      maxid = offspring_[i];
  }
  
  // Resizing the vectors
  this->parents.resize(maxid + 1u);
  this->offspring.resize(maxid + 1u);
  
  this->visited.resize(maxid + 1u, false);
  this->visit_counts.resize(maxid + 1u, 0u);
  
  // Adding the data
  for (uint i = 0u; i < m; ++i) {
    this->offspring[parents_[i]].push_back(offspring_[i]);
    this->parents[offspring_[i]].push_back(parents_[i]);
  }
  
  
  // Constants
  this->N_NODES = (uint) maxid + 1u;
  this->N_EDGES = m;
  
  // Generating the postorder sequence
  this->postorder();
  
  // Initializing iterator 
  this->iter = TreeIterator<Data_Type>(this);
  
  // Some checks
  if (!this->is_connected()) {
    
    out = 3u;
    return;
    
  } else if (!this->is_dag()) {
    
    out = 4u;
    return;
    
  }
  
  // Marking tip nodes
  int status = 0;
  this->iter.bottom();
  this->TIPS.reserve(this->n_tips());
  while ( status == 0 ) {
    if (this->iter.is_tip())
      this->TIPS.push_back(*this->iter);
    status = this->iter.up();
  }
  this->iter.bottom();
  
  
  out = 0u;
  return;
  
}

// A recursive function to check whether the tree is a DAG or not. -------------
typedef v_uint::const_iterator v_uint_iter;

template <typename Data_Type>
inline bool Tree<Data_Type>::is_dag() {
  
  bool res = this->is_dag_();
  this->reset_visited();
  
  return res;
  
}

template <typename Data_Type>
inline bool Tree<Data_Type>::is_dag_(int i, int caller, bool up_search) {
  
  // For the first iteration
  if (i < 0) 
    i = 0u, caller = -1;
  
  // Yes, this is not a dag (came here multiple times)
  if (this->visited[i])
    return false;
  this->visited[i] = true;
  
  // Iterating through parents
  for (v_uint_iter n = this->parents[i].begin(); n != this->parents[i].end(); ++n) {
    
#ifdef DEBUG_TREE
    std::printf(
      "Tree<Data_Type>::is_dag() @ parents   (i, caller, *n, up_search): (%i, %i, %i, %i)\n",
      i, caller, *n, up_search
    );
#endif
    
    // Checking 1:1 cycles
    if ((int) *n == caller) {
#ifdef DEBUG_TREE
      std::printf("\tChecking 1:1 cycles.\n");
#endif
      if (up_search) return false;
      else continue; 
    }
    
    if (!(this->is_dag_((int) *n, i, true)))
      return false;
  }
  
  // Iterating through offspring
  for (v_uint_iter n = this->offspring[i].begin(); n != this->offspring[i].end(); ++n) {
    
#ifdef DEBUG_TREE
    std::printf(
      "Tree<Data_Type>::is_dag() @ offspring (i, caller, *n, up_search): (%i, %i, %i, %i)\n",
      i, caller, *n, up_search
    );
#endif
    
    // Checking 1:1 cycles
    if ((int) *n == caller) {
#ifdef DEBUG_TREE
      std::printf("\tChecking 1:1 cycles.\n");
#endif
      if (!up_search) return false;
      else continue; 
    }
    
    if (!(this->is_dag_((int) *n, i, false)))
      return false;
  }
  
  return true;
  
}

// The preorder is just the post order reversed
template <typename Data_Type>
inline v_uint Tree<Data_Type>::get_preorder() const {
  
  v_uint res = this->get_postorder();
  std::reverse(res.begin(), res.end());
  return res;
  
}

template <typename Data_Type>
inline uint Tree<Data_Type>::get_dist_tip2root_(uint i, uint count) {
  
  // If we are already in a root node, then return with the reached count
  if (parents[i].size() == 0u)
    return count;
  
  // We start by adding one step (since we are not at the root!)
  ++count;
  
  // Otherwise, let's see at its parents (who is closests to the root)
  v_uint counts(parents[i].size());
  uint j = 0;
  for (auto iter = parents[i].begin(); iter != parents[i].end(); ++iter) {
    
    // Looking how far do we reach
    counts[j] = get_dist_tip2root_(*iter, count);
    
    // If we reached to root in the next step. Otherwise we need to keep looking
    // until we get to that is closests.
    if (counts[j++] == count)
      return count;
    
  }
  
  // We keep the smallest one
  count = counts[0];
  for (auto iter = counts.begin() + 1; iter != counts.end(); ++iter)
    if (*iter < count)
      count = *iter;
  
  return count;
}

template <typename Data_Type>
inline v_uint Tree<Data_Type>::get_dist_tip2root() {
  
  if (this->DIST_TIPS2ROOT.size() == 0u) {
    
    // Making space available
    this->DIST_TIPS2ROOT.resize(this->n_tips());
    
    for (uint i = 0u; i < TIPS.size(); ++i)
      DIST_TIPS2ROOT[i] = get_dist_tip2root_(TIPS[i], 0u);
    
  }
  
  return this->DIST_TIPS2ROOT;
  
}

#define TOTAL(a) (a)->offspring.size()

template <typename Data_Type>
inline uint Tree<Data_Type>::set_postorder(const v_uint & POSTORDER_, bool check) { 
  
  // Checking the range of the data
  if (check) {
    uint min_idx = 999999u, max_idx = 0u;
    for (auto i = POSTORDER_.begin(); i != POSTORDER_.end(); ++i) {
      
      if (*i > max_idx) max_idx = *i;
      if (*i < min_idx) min_idx = *i;
      
    }
    
    if ((min_idx > max_idx) | (min_idx > TOTAL(this)) | (max_idx > TOTAL(this)))
      return 1u;
  }
  
  this->POSTORDER = POSTORDER_;
  
  return 0u;
}

template <typename Data_Type>
inline void Tree<Data_Type>::prune_postorder() {
  
  // Set the head in the first node of the sequence
  this->iter.bottom();
  int status = 0;
  while (status == 0) {
    
    this->eval_fun();
    status = this->iter.up();
    
  }
  
  return;
  
}

template <typename Data_Type>
inline void Tree<Data_Type>::prune_postorder(v_uint & seq) {

  // Let's reset the sequence
  v_uint OLDPOSTORDER = POSTORDER;
  POSTORDER.swap(seq);
    
  // Set the head in the first node of the sequence
  this->iter.bottom();
  int status = 0;
  while (status == 0) {
    
    this->eval_fun();
    status = this->iter.up();
    
  }
  
  // Going back to the previous order
  POSTORDER.swap(OLDPOSTORDER);
  
  return;
  
}

template <typename Data_Type>
inline void Tree<Data_Type>::prune_preorder() {
  
  // Set the head in the first node of the sequence
  this->iter.top();
  int status = 0;
  while (status == 0) {
    
    this->eval_fun();
    status = this->iter.down();
    
  }
  
  return;
  
}

template <typename Data_Type>
inline void Tree<Data_Type>::prune_preorder(v_uint & seq) {
  
  // Let's reset the sequence
  v_uint OLDPOSTORDER = POSTORDER;
  POSTORDER.swap(seq);
  
  // Set the head in the first node of the sequence
  this->iter.top();
  int status = 0;
  while (status == 0) {
    
    this->eval_fun();
    status = this->iter.down();
    
  }
  
  // Going back to the previous order
  POSTORDER.swap(OLDPOSTORDER);
  
  return;
  
}

template <typename Data_Type>
inline uint Tree<Data_Type>::n_tips() const {
  
  uint count = 0u;
  for (auto i = this->offspring.begin(); i != offspring.end(); ++i)
    if (i->size() == 0u)
      ++count;
  
  return count;
  
}

template <typename Data_Type>
inline int Tree<Data_Type>::n_offspring(uint i) const {
  
  if (i < this->offspring.size()) {
    return this->offspring.at(i).size();
  }
  return -1;
}

template <typename Data_Type>
inline int Tree<Data_Type>::n_parents(uint i) const {
  
  if (i < this->parents.size()) {
    return this->parents.at(i).size();
  }
  return -1;
}

#endif

