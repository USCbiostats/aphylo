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

template <typename T, typename D, typename H = std::hash<T>>
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

/**
 * Creates a map to 
 */
inline void as_phylo_key(
    Vec<double> & res,
    const Vec<bool> & states,
    const Vec<double> & blengths
) {
  
  // Checking the size
  if (res.size() != (blengths.size() + 1u))
    res.resize(1u + blengths.size());
  
  // Vec<double> res(1u + blengths.size(), 0.0);
  res[0u] = 1.0 * pow(10, states.size());
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
  as_phylo_key(res, states, blengths);
  
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
  
  // To compute the normalizing constant
  Vec< double >        weights;
  Vec< Vec< double > > statmat;
  double               numerator = 0.0;
  
  // To recompute the probabilities
  double        denominator = 0.0;
  Vec< double > params;
  
  
  // Other variables
  double temp           = 0.0;
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
  
  Vec< double > prob_all(const Vec<double> & par);
  
  // Getters
  const Vec<double> * get_weights() const {return &weights;};
  const Vec<double> * get_params() const {return &params;};
  const Vec<Vec<double>> * get_statmat() const {return &statmat;};
  
};

Vec<double> ingredients::empty_dbl = {};

inline double ingredients::probability(
    const Vec< double > & par,
    const Vec< double > & target
) {
  
  // Numerator is always re calculated
  numerator = 0.0;
  for (unsigned int i = 0u; i < par.size(); ++i)
    numerator += par[i] * target[i];
  numerator = exp(numerator);
   
   // std::cout << "Computing probs for target [";
   // for (auto v = target.begin(); v!=target.end(); ++v)
   //   std::cout << *v << ", ";
   // std::cout << std::endl;
   
  // Need to update the denominator
  if (!equal(par, params)) {
    
    // std::cout << "Computing probs for vector [";
    // for (auto v = par.begin(); v!=par.end(); ++v)
    //   std::cout << *v << ", ";
    // std::cout << std::endl;
    
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
  
  // printf("numerator: %.4f, denominator %.4f\n",numerator, denominator);
  
  // Returning the probability
  return numerator/denominator;
  
}

inline Vec<double> ingredients::prob_all(const Vec<double> & par) {
  
  // Reserving space, and iterating through all possible cases
  Vec<double> res(this->weights.size(), 0.0);
  for (unsigned int i = 0u; i < res.size(); ++i){
    res[i] = probability(par, statmat[i]) * weights[i];
  }
  
  return res;
  
}

/**
 * The bank is a collection of ingredients!
 */

typedef std::unordered_map< Vec<double>, ingredients, BHash< double > > ingredients_map;

class bank {
protected:
  
  std::unordered_map< unsigned int, ingredients_map> data;
  Vec<double> tmpkey;
  Vec<std::string> counters;
  Vec<Vec<unsigned int>> counters_parameters;
  // Probably will need to add a list of statistics that needs to be
  // passed to the actual counter function.
  // Vec< phylocounters::PhyloCounter > counters;
  
public:
  
  bank() : data(0u), tmpkey(2u, 0.0), counters(0u), counters_parameters(0u) {};
  ~bank() {
    // for (auto iter = counters.begin(); iter != counters.end(); ++iter) {
    //   delete iter->data;
    //   iter->data = nullptr;
    // }
  };

  /**@brief Computes the probability as a wrapper of the ingredient.
   * 
   * This function checks whether the requested type of event exists and
   * returns the probability associated to observing that event.
   */
  double probability(
      const Vec<double> & target_stats,
      const Vec<double> & par,
      const unsigned int & noff,
      const Vec<bool> & state,
      const Vec<double> & blength,
      bool check = true
  );
  
  /**@brief Computes the probability of observing a given array.
   * 
   * This is more computationally cost than the other version of this
   * function.
   * 
   */
  double probability(
      const phylocounters::PhyloArray & Array,
      const Vec<double> & par,
      bool check = true
  );
  
  // Query functions
  unsigned int size() const;

  ingredients_map * get_ingredients(unsigned int i) {
    if (data.find(i) == data.end())
      throw std::range_error("The offspring is out of range!");
    return &(data[i]);
  };
  
  const Map< unsigned int, ingredients_map>::const_iterator begin() const {return data.begin();};
  const Map< unsigned int, ingredients_map>::const_iterator end() const {return data.end();};

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

inline double bank::probability(
    const Vec<double> & target_stats,
    const Vec<double> & par,
    const unsigned int & noff,
    const Vec<bool> & state,
    const Vec<double> & blength,
    bool check
) {
  
  // Concatenating to search
  as_phylo_key(tmpkey, state, blength);
  
  // Finding the corresponding record
  if (check) {
    auto record = data.find(noff);
    if (record == data.end()) 
      throw std::range_error("Number of offspring out of range.");
    
    if (record->second.find(tmpkey) == record->second.end())
      throw std::range_error("Combination of (state, branch length) not found.");
    
  }
  
  return data[noff][tmpkey].probability(par, target_stats);

}

inline double bank::probability(
    const phylocounters::PhyloArray & Array,
    const Vec<double> & par,
    bool check
) {
  
  // Generating the counter function
  barray::StatsCounter< phylocounters::PhyloArray, Vec<unsigned int> >
    S(&Array);
  
  // Adding stats
  S.add_counter(phylocounters::overall_gains);
  S.add_counter(phylocounters::overall_loss);
  
  // Longest branch gains a function
  phylocounters::PhyloCounter count2 = phylocounters::longest;
  count2.data = new Vec< unsigned int >({});
  S.add_counter(count2);
  
  S.count_all();
  delete count2.data;
  
  return 
    probability(
      S.current_stats,
      par, 
      Array.M,
      Array.data->states,
      Array.data->blengths,
      true
    );
  
}

inline unsigned int bank::size() const {
  
  unsigned int res = 0u;
  for (auto iter = data.begin(); iter != data.end(); ++iter)
    res += iter->second.size();
    
  return res;
  
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
  as_phylo_key(tmpkey, state, blength);
  
  // Any record there?
  if (data.find(noff) != data.end()) {
    auto key_iter = data[noff].find(tmpkey);
    if (key_iter != data[noff].end()) {
      key_iter->second.nqueries++;
      return;
    }
    
  } else // Need to initialize it
    data[noff] = ingredients_map(0u);

  
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
  
  S.add_counter(count2);
  
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
  
  data[noff][tmpkey] = ingredients(W, StatMat, StatMat[0u].size());
  
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
    
    std::cout << "#offspring: " << iter1->first << std::endl;
    
    // Iterating over the map
    for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2) {
      std::cout << "Key: [";
      for (auto key = iter2->first.begin(); key != iter2->first.end(); ++key)
        std::cout << *key <<", ";
      
      std::cout << "], # elements: " << iter2->second.weights.size() <<
        " # queried: " << iter2->second.nqueries << std::endl;
      
      // Overall counters
      n        += iter2->second.nqueries + 1u;
      n_unique += iter2->second.weights.size();
      
    }

  }
  std::cout << "Total entries parsed: " << n << " with " << 
    n_unique << " unique entries compared to " << this->size() <<
      " ingredients." << std::endl;
  
  return;
  
}

#endif

// [[Rcpp::export]]
List testing_a_bank(
    const std::vector< int > & noff,
    const std::vector< std::vector< bool > > & states,
    const std::vector< std::vector< double > > & blengths,
    const std::vector< double > & par
    ) {
  
  // Creating the bank
  Rcpp::XPtr< bank > ptr(new bank(), true);
  
  for (unsigned int i = 0u; i < noff.size(); ++i)
    ptr->add( (unsigned int) noff[i], states[i], blengths[i]);
  
  // ptr->print();
/*  
  // Example computing transition probabilities
  Vec< double > probs_empty(noff.size());
  Vec< double > probs_nicer(noff.size());
  for (unsigned int i = 0u; i < noff.size(); ++i) {
    
    // Creating the situation (for now, just counting on empty trees)
    phylocounters::PhyloArray tree(states[0u].size(), noff[i]);
    tree.data = new phylocounters::NodeData(blengths[i], states[i]);
    
    probs_empty[i] = ptr->probability(tree, par);
    
    // Picking just the first case
    Vec<double> key(1u + blengths[i].size());
    as_phylo_key(key, states[i], blengths[i]);
    probs_nicer[i] = ptr->probability(
      ptr->get_ingredients(noff[i])->at(key).get_statmat()->at(0u),
      par,
      noff[i],
      states[i],
      blengths[i]
    );
    
    delete tree.data;
    tree.data = nullptr;
    
  }
  */
  // Computing each ingredients entire set of probabilities. These 
  // should add up to one.
  ingredients_map * ingr_ptr2 = ptr->get_ingredients(2);
  ingredients_map * ingr_ptr3 = ptr->get_ingredients(3);
  Vec< Vec<double> > probs(0u);
  for (auto iter = ingr_ptr2->begin(); iter != ingr_ptr2->end(); ++iter) 
    probs.push_back(iter->second.prob_all(par));
  
  for (auto iter = ingr_ptr3->begin(); iter != ingr_ptr3->end(); ++iter) 
    probs.push_back(iter->second.prob_all(par));
  
  // Collecting the results
  List res(ptr->size());
  unsigned int i = 0u;
  for (auto iter0 = ptr->begin(); iter0 != ptr->end(); ++iter0) {
    
    for (auto iter = iter0->second.begin(); iter != iter0->second.end(); ++iter) {
    
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
      
      res[i] = List::create(
        _["key"]     = wrap(iter->first),
        _["weights"] = wrap((*tmp->get_weights())),
        _["statmat"] = clone(statmat),
        _["probs"]   = wrap(probs[i])
      );
      
      i++;
    }
    
  }
  
  return res;
}

/***R

# Sanity check: We expect that all probabilities add up to one once we calculate
# them over the entire sample space (including weights)
check_probabilities <- function(x) {
  
  sums <- sapply(x, function(i) sum(i$probs))
  
  sum(abs(sums - 1) < 1e-10)/length(sums)

  
}

# Generating some data
nstates <- 400
nfuns   <- 4
set.seed(5544)

par     <- runif(3)

noff    <- sample.int(2, nstates, TRUE) + 1
states  <- replicate(nstates, as.logical(sample.int(2, nfuns, TRUE) - 1), simplify = FALSE)


# Random branch length
lenghts <- lapply(noff, runif)
ans0 <- testing_a_bank(noff, states, lenghts, par)

# Ranked branch length
lenghts <- lapply(lenghts, order)
ans1 <- testing_a_bank(noff, states, lenghts, par)

# No differences on branch length
lenghts <- lapply(noff, function(i) rep(1, i))
ans2 <- testing_a_bank(noff, states, lenghts, par)

# These tell us the proportion of each case that has total probability near to 1
check_probabilities(ans0)
check_probabilities(ans1)
check_probabilities(ans2)

# rm(list = ls())
*/
