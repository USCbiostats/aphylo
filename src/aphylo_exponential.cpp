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

/**
 * Compares two vectors of type T and returns `true` if the norm between
 * the two is smaller than `eps`.
 */
template <typename T>
bool equal(const Vec<T> & a, const Vec<T> & b, double eps = 1e-10) {
  
  if (a.size() != b.size())
    throw std::logic_error("a and b should be of the same size.");
  
  double diff = 0.0;
  for (unsigned int i = 0u; i < a.size(); ++i)
    diff += pow((double) (a[i] - b[i]), 2.0);
  
  return (sqrt(diff) < eps)?true:false;
  
}

inline void as_phylo_key(
    Vec<double> & res,
    const unsigned int & noff,
    const Vec<bool> & states,
    const Vec<double> & blengths
) {
  
  // Vec<double> res(1u + blengths.size(), 0.0);
  res[0u] = noff * pow(10, states.size());
  for (unsigned int i = 0u; i < states.size(); ++i)
    if (states[i])
      res[0u] += pow(10, i);
    
  std::copy(blengths.begin(), blengths.end(), res.begin() + 1);
      
  return;
  
}

// [[Rcpp::export]]
NumericVector aspk(
    const unsigned int & noff,
    const Vec<bool> & states,
    const Vec<double> & blengths
    ) {
  
  Vec<double> res(1 + blengths.size(), 0.0);
  as_phylo_key(res, noff, states, blengths);
  
  return Rcpp::wrap(res);
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
  
  unsigned int nqueries = 0u;
  
  friend class bank;
  
public:
  ingredients() :
  weights(0u), params(0u), statmat(0u, empty_dbl) {};
  
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
  
  // Getters
  const Vec<double> * get_weights() const {return &weights;};
  const Vec<double> * get_params() const {return &params;};
  const Vec<Vec<double>> * get_statmat() const {return &statmat;};
  
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
  
  ingredients_map data;
  Vec<double> tmpkey;
  // Probably will need to add a list of statistics that needs to be
  // passed to the actual counter function.
  // Vec< phylocounters::PhyloCounter > counters;
  
public:
  
  bank() : data(0u), tmpkey(2u, 0.0) {};
  ~bank() {
    // for (auto iter = counters.begin(); iter != counters.end(); ++iter) {
    //   delete iter->data;
    //   iter->data = nullptr;
    // }
  };

  // Query functions
  unsigned int size() const {return data.size();};

  const ingredients_map::const_iterator begin() const {return data.begin();};
  const ingredients_map::const_iterator end() const {return data.end();};

  // Manipulation functions
  void add(
      const unsigned int & noff,
      const Vec<bool> & state,
      const Vec<double> & blength
  );
  
  // void add_counter(const phylocounters::PhyloCounter & counter);

  // Vec< unsigned int > * counter_data_ptr(unsigned int i);
  
  // Misc
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
  
  // Generating double vector key for the hash map
  if (tmpkey.size() != (1u + blength.size()))
    tmpkey.resize(1u+blength.size());
  as_phylo_key(tmpkey, noff, state, blength);
  
   
  // Is it in the boundary and is non zero?
  auto key_iter = data.find(tmpkey);
  if (key_iter != data.end()) {
    key_iter->second.nqueries++;
    return;
  }
  
  // Here is the bulk of the work, we need to calculate the powerset of the
  // thing
  phylocounters::PhyloArray Array(state.size(), blength.size());
  Array.data = new phylocounters::NodeData(blength, state);
  
  // Generating the support counter
  barray::Support< phylocounters::PhyloArray, Vec< unsigned int > > S(&Array);
  
  // Adding a few statistics
  // Adding some model terms
  S.add_counter(phylocounters::overall_gains);
  S.add_counter(phylocounters::overall_loss);
  
  // Longest branch gains a function
  phylocounters::PhyloCounter count2 = phylocounters::longest;
  count2.data = new Vec< unsigned int >({});
  
  // // Co-evol of functions
  // Vec< phylocounters::PhyloCounter > coevolve
  // for (unsigned int i = 0u; i < state.size(); ++i) {
  //   
  // }

  // Computing and retrieving
  S.calc(0u, true);
  
  delete Array.data;
  delete count2.data;
  Array.data = nullptr;
  count2.data = nullptr;
  
  // Now, iterating through the data to generate the ingredients
  Vec< double > W(0u);
  Vec< Vec< double > > StatMat(0u);
  for (auto iter = S.support.stats.begin(); iter != S.support.stats.end(); ++iter) {
    StatMat.push_back(iter->first);
    W.push_back(iter->second);
  }
  
  data[tmpkey] = ingredients(W, StatMat, StatMat[0u].size());
  
  return;
  
}

// inline Vec< unsigned int > * bank::counter_data_ptr(unsigned int i) {
//   return counters.at(i).data;
// }
// 
// inline void bank::add_counter(const phylocounters::PhyloCounter & counter) {
//   counters.push_back(counter);
//   counters[counters.size() - 1u].data = new Vec< unsigned int >({0u, 0u});
//   return;
// };

inline void bank::print() const {
  
  unsigned int n = 0u;
  unsigned int n_unique = 0u;
  
  // Iterating on number of offspring
  for (auto iter1 = data.begin(); iter1 != data.end(); ++iter1) {
    std::cout << "Key: [";
    ++n;
    for (auto iter2 = iter1->first.begin(); iter2 != iter1->first.end(); ++iter2)
      std::cout << *iter2 <<", ";
    std::cout << "], # elements: " << iter1->second.weights.size() <<
      " # queried: " << iter1->second.nqueries << std::endl;
    
    // Overall counters
    n        += iter1->second.nqueries;
    n_unique += iter1->second.weights.size();
    
  }
  std::cout << "Total entries parsed: " << n << " with " << 
    n_unique << " unique entries compared to " << data.size() <<
      " ingredients." << std::endl;
  
  return;
  
}

#endif

// [[Rcpp::export]]
List testing_a_bank(
    const std::vector< int > & noff,
    const std::vector< std::vector< bool > > & states,
    const std::vector< std::vector< double > > & lenghts
    ) {
  
  // Creating the bank
  Rcpp::XPtr< bank > ptr(new bank(), true);
  
  for (unsigned int i = 0u; i < noff.size(); ++i)
    ptr->add( (unsigned int) noff[i], states[i], lenghts[i]);
  
  ptr->print();
  
  // Collecting the results
  Rcout << "Getting the result" << std::endl;
  List res(ptr->size());
  unsigned int i = 0u;
  for (auto iter = ptr->begin(); iter != ptr->end(); ++iter) {
    
    const ingredients * tmp = &iter->second;
    NumericMatrix statmat(
        tmp->get_statmat()->size(),
        tmp->get_params()->size()
        );
    
    unsigned int nrow = 0u;
    for (auto entry = tmp->get_statmat()->begin(); entry != tmp->get_statmat()->end(); ++entry) {
      for (unsigned int ncol = 0u; ncol < statmat.ncol(); ++ncol)
        statmat(nrow, ncol) = entry->operator[](ncol);
      nrow++;
    }
    
    res[i++] = List::create(
      _["key"]     = wrap(iter->first),
      _["weights"] = wrap((*tmp->get_weights())),
      _["statmat"] = clone(statmat)
    );
    
  }
  
  return res;
}

/***R

# Generating some data
nstates <- 100
nfuns   <- 2
set.seed(5544)
noff    <- sample.int(2, nstates, TRUE) + 1
states  <- replicate(nstates, as.logical(sample.int(2, nfuns, TRUE) - 1), simplify = FALSE)

# Random branch length
lenghts <- lapply(noff, runif)
testing_a_bank(noff, states, lenghts)

# Ranked branch length
lenghts <- lapply(lenghts, order)
testing_a_bank(noff, states, lenghts)

# No differences on branch length
lenghts <- lapply(noff, function(i) rep(1, i))
testing_a_bank(noff, states, lenghts)

*/
