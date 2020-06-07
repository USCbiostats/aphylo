#include <RcppArmadillo.h>
#include <unordered_map>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

#ifndef APHYLO_EXPONENTIAL_H
#define APHYLO_EXPONENTIAL_H 1
#include "../include/barray/barray.hpp"
#include "../include/pruner/pruner.hpp"

template <typename T>
using Vec = std::vector< T >;

template <typename T, typename D, typename H>
using Map = std::unordered_map< T , D, H >;

template <typename T>
using BHash = barray::vecHasher< T >;

template <typename T>
bool equal(const Vec<T> & a, const Vec<T> & b, double eps = 1e-10) {
  
  if (a.size() != b.size())
    throw std::logic_error("a and b should be of the same size.");
  
  double diff = 0.0;
  for (unsigned int i = 0u; i < a.size(); ++i)
    diff += pow((double) (a[i] - b[i]), 2.0);
  
  return (sqrt(diff) < eps)?true:false;
  
}

// Brief declaration since banks should be friends of ingredients
class bank;

/**
 * This class holds all the info relevant for a given:
 * - Number of offspring.
 * - state of the parent.
 * - lengths of the branches.
 */
class ingredients {
private:
  
  static std::vector< double > empty_dbl;
  
protected:
  Vec< double > weights;
  Vec< double > params;
  Vec< Vec< double > > statmat;
  double numerator   = 0.0;
  double denominator = 0.0;
  double temp = 0.0;
  friend class bank;
public:
  ingredients() : weights(0u), params(0u), statmat(0u, empty_dbl) {};
  
  ingredients(
    const Vec< double > & W,
    const Vec< Vec< double > > & S,
    unsigned int npar
  ): weights(W), params(npar, 0.0), statmat(S) {};
  
  ~ingredients() {};
  
  double probability(
    const Vec< double > & par,
    const Vec< double > & target
    );
};

Vec<double> ingredients::empty_dbl = {};

inline double ingredients::probability(
    const Vec< double > & par, const Vec< double > & target
) {
  
  // Checking if it has changed
  numerator = 0.0;
  for (unsigned int i = 0u; i < par.size(); ++i)
    numerator += par[i] * target[i];
  numerator = exp(numerator);
   
  // Need to update the denominator
  if (!equal(par, params)) {
    
    denominator = 0.0;
    
    for (unsigned int i = 0u; i < weights.size(); ++i) {
      temp = 0.0;
      for (unsigned int j = 0u; j < par.size(); ++j) {
        temp += par[j] * statmat[i][j];
      }
      denominator += exp(temp) * weights[i];
    }
    
    // Saving the results
    params = par;
    
  }
  
  // Returning the probability
  return numerator/denominator;
  
}

/**
 * The bank is a collection of ingredients!
 */

typedef std::unordered_map< Vec<double>, ingredients, BHash< double > > ingredients_map;

class bank {
protected:
  
  Vec< std::unordered_map< Vec<bool>, ingredients_map, BHash< bool > >> data;
  // Probably will need to add a list of statistics that needs to be
  // passed to the actual counter function.
  // Vec< phylocounters::PhyloCounter > x;
  
public:
  
  bank() : data(0u) {};
  ~bank() {};
  
  void add(
      const unsigned int & noff,
      const Vec<bool> & state,
      const Vec<double> & blength
  );
  
  void print() const;
  
};

/**
 * Adds a new entry to the bank.
 */
inline void bank::add(
    const unsigned int & noff,
    const Vec<bool> & state,
    const Vec<double> & blength
) {
  
  if (noff == 0u)
    throw std::logic_error("Not supported for zero offspring!.");
  
  
  // Is it in the boundary and is non zero?
  unsigned int noff_less1 = noff - 1u;
  if ((noff_less1 < data.size()) && (!data[noff_less1].empty())) {
    
    // Is it there? 
    auto state_loc = data[noff_less1].find(state); 
    if (state_loc != data[noff_less1].end()) {
      
      // Is it in the needed branch, if already in there, meaning that it is not
      // at the end fo it, we can safely return.
      if (state_loc->second.find(blength) != data[noff_less1][state_loc->first].end()) {
        return;
      }
      
    } else {
      ingredients_map empty_ingredient(0u);
      data[noff_less1][state] = empty_ingredient;
    }
      
  } else {
    // Not even containing it
    if (data.size() < noff)
      data.resize(noff);
    
    // Adding the first map
    Map< Vec<bool>, ingredients_map, BHash<bool>> empty_map(0u);
    data[noff_less1] = empty_map;
    
    // Adding the second map
    ingredients_map empty_ingredient(0u);
    data[noff_less1][state] = empty_ingredient;
  }
  
  // Here is the bulk of the work, we need to calculate the powerset of the
  // thing
  phylocounters::PhyloArray Array(state.size(), blength.size());
  Array.data = new phylocounters::NodeData(blength, state);
  
  // Generating the support counter
  barray::Support< phylocounters::PhyloArray, Vec< unsigned int > > S(&Array);
  
  // Adding a few statistics
  S.add_counter(phylocounters::overall_gains);
  S.add_counter(phylocounters::overall_loss);
  
  // Computing and retrieving
  S.calc(0u, true);
  
  delete Array.data;
  Array.data = nullptr;
  
  // Now, iterating through the data to generate the ingredients
  Vec< double > W(0u);
  Vec< Vec< double > > StatMat(0u);
  for (auto iter = S.support.stats.begin(); iter != S.support.stats.end(); ++iter) {
    StatMat.push_back(iter->first);
    W.push_back(iter->second);
  }
  
  data[noff_less1][state][blength] = ingredients(W, StatMat, StatMat[0u].size());
  
  std::cout << "Adding a new set of stats " << std::endl;
  
  return;
  
}

inline void bank::print() const {
  
  unsigned int n = 0u;
  
  // Iterating on number of offspring
  for (auto iter1 = data.begin(); iter1 != data.end(); ++iter1) {
    
    std::cout << "Noffspring: " << ++n << " total records: " << iter1->size() << std::endl;
    if (iter1->size() == 0u)
      continue;
    else
      
    for (auto iter2 = iter1->begin(); iter2 != iter1->end(); ++iter2) {
      
      std::cout << "  Listing records for state: [";
      for (auto states_i = iter2->first.begin(); states_i != iter2->first.end(); ++states_i) 
        std::cout << ((*states_i)? 1 : 0) << ", ";

      std::cout << "],  total records " << iter2->second.size() << std::endl;
      for (auto states_i = iter2->second.begin(); states_i != iter2->second.end(); ++states_i) {
        
        std::cout << "    Listing records for state: [";
        for (auto len_i = states_i->first.begin(); len_i != states_i->first.end(); ++len_i) 
          std::cout << *len_i << ", ";
        
        std::cout << "] the stats have: " << states_i->second.weights.size() <<
          " many unique stats" << std::endl;
        
      }
        
    }
  }
  
  return;
  
}

#endif

// [[Rcpp::export]]
int testing_a_bank(
    const std::vector< int > & noff,
    const std::vector< std::vector< bool > > & states,
    const std::vector< std::vector< double > > & lenghts
    ) {
  
  Rcpp::XPtr< bank > ptr(new bank(), true);
  for (unsigned int i = 0u; i < noff.size(); ++i)
    ptr->add( (unsigned int) noff[i], states[i], lenghts[i]);
  
  ptr->print();
  
  return 0;
}

/***R

# Generating some data
nstates <- 20
nfuns   <- 4
set.seed(5544)
noff    <- sample.int(2, nstates, TRUE) + 1
states  <- replicate(nstates, as.logical(sample.int(2, nfuns, TRUE) - 1), simplify = FALSE)
lenghts <- lapply(noff, runif)

testing_a_bank(noff, states, lenghts)

*/
