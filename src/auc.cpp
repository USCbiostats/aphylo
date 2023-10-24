#include <Rcpp.h>
using namespace Rcpp;

//' Area Under the Curve and Receiving Operating Curve
//' 
//' The AUC values are computed by approximation using the area of the polygons formed
//' under the ROC curve.
//' @param pred A numeric vector with the predictions of the model. Values must
//' range between 0 and 1.
//' @param labels An integer vector with the labels (truth). Values should be either
//' 0 or 1.
//' @param nc Integer. Number of cutoffs to use for computing the rates and AUC.
//' @param nine_na Logical. When `TRUE`, 9 is treated as `NA`.
//' @return A list:
//' - `tpr` A vector of length `nc` with the True Positive Rates.
//' - `tnr` A vector of length `nc` with the True Negative Rates.
//' - `fpr` A vector of length `nc` with the False Positive Rates.
//' - `fnr` A vector of length `nc` with the False Negative Rates.
//' - `auc` A numeric value. Area Under the Curve.
//' - `cutoffs` A vector of length `nc` with the cutoffs used.
//' @export
//' @examples
//' set.seed(8381)
//' x   <- rdrop_annotations(raphylo(50), .3)
//' ans <- aphylo_mcmc(x ~ mu_d + mu_s + Pi)
//' ans_auc <- auc(predict(ans, loo = TRUE), x[,1,drop=TRUE])
//' print(ans_auc)
//' plot(ans_auc)
// [[Rcpp::export]]
List auc(
    const NumericVector & pred,
    const IntegerVector & labels,
    int nc = 200,
    bool nine_na = true
) {
  
  // Checking the dimension
  unsigned int nobs = pred.size();
  if (nobs != labels.size())
    stop(
      "Number of observed predicted points %i does not match the number of expected points %i.",
      nobs, labels.size()
    );
  
  // Tagging values
  std::vector< unsigned int > locations;
  locations.reserve(nobs);
  int n0 = 0, n1 = 0;
  
  for (unsigned int j = 0; j < nobs; ++j)
  {
    
    if (NumericVector::is_na(pred[j]))
      continue;
    
    if (!IntegerVector::is_na(labels[j]))
    {
      
      if (labels[j] == 0) ++n0;
      else if (labels[j] == 1) ++n1;
      else if (nine_na && labels[j] == 9) continue;
      else 
        stop("Only values 0/1 are supported.");
      
    } else
      continue;
    
    locations.push_back(static_cast<unsigned int>(j));
    
  }
  locations.shrink_to_fit();

  // Rates
  NumericVector TPR(nc);
  NumericVector TNR(nc);
  NumericVector FPR(nc);
  NumericVector FNR(nc);
  
  // Generating the sequence
  NumericVector s(nc);
  for (int i = 0; i < nc; ++i) 
    s[i] = ((double) i)/((double) nc);
  
  // Evaluating the sequence
  double auc = 0.0;
  for (int i = 0; i < nc; ++i) {
    
    TPR[i] = 0.0;
    TNR[i] = 0.0;
    FPR[i] = 0.0;
    FNR[i] = 0.0;
    
    for (
        std::vector< unsigned int >::const_iterator j = locations.begin();
        j != locations.end(); ++j) {
      
      if (pred[*j] > s[i]) {
        if (labels[*j] == 1) ++TPR[i];
        else ++FPR[i];
      } else {
        if (labels[*j] == 1) ++FNR[i];
        else ++TNR[i];
      }
      
    }
    
    TPR[i] /= n1;
    TNR[i] /= n0;
    FPR[i] /= n0;
    FNR[i] /= n1;
    
    if (i == 0) auc += (1+TPR[0])*(1 - FPR[0])/2.0;
    else auc += (TPR[i] + TPR[i - 1])*(FPR[i - 1] - FPR[i])/2.0;
    
  }
  
  List ans = List::create(
    _["tpr"] = TPR,
    _["tnr"] = TNR,
    _["fpr"] = FPR,
    _["fnr"] = FNR,
    _["auc"] = auc,
    _["n_used"] = locations.size(),
    _["cutoffs"] = s
  );
  
  ans.attr("class") = "aphylo_auc";
  
  return ans;
  
}

