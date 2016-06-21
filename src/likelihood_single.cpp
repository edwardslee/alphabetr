#include <Rcpp.h>
using namespace Rcpp;

//' Calculate likelihood of a single TCR clone
//'
//' @param est a
//' @param err a
//' @param numb_cells a
//' @param numb_wells a
//' @param numb_sample a
//'
//' @export
// [[Rcpp::export]]
double likelihood_single(double est, double err, NumericVector numb_wells, NumericVector numb_cells, NumericVector numb_sample)
{
  double obj = 0;            // value of the likelihood objective function
  int n = numb_wells.size(); // number of wells with the clones in each column

  for (int i = 0; i < n; i++) {
    double appear = numb_wells[i];      // number of wells in column w/ clone
    double sample_size = numb_cells[i]; // sample size per well in column
    double p_none = 0;                  // prob of not seeing clone in well

    for (int j = 0; j < sample_size; j++) {
      int k = j + 1;
      p_none += (2 * pow(err, k) - pow(err, 2*k)) * R::choose(sample_size, k) *
        pow(est, k) * pow((1 - est), (sample_size - k));
    }
    p_none += pow((1 - est), sample_size);
    p_none = -appear * log(1 - p_none) - (numb_sample[i] - appear) * log(p_none);

    obj = obj + p_none;
  } // end for - i
  return obj;
} // end function likelihood_single
