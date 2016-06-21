#include <Rcpp.h>
using namespace Rcpp;
//' Calculate likelihood of dual clone
//'
//' @param est a
//' @param err a
//' @param numb_cells a
//' @param numb_wells a
//' @param binomials a
//'
//' @export
// [[Rcpp::export]]
NumericVector dual_discrim_dual_likelihood(double est, double err, NumericVector numb_cells, NumericMatrix numb_wells, List binomials) {
  NumericVector log_like = NumericVector::create(0.0);                // log likelihood value
  int numb_ss = numb_wells.nrow();     // number of distinct sample size values
  /* For each sample size, we look at the wells of that sample size and
  * calculate the likelihood of the well apperances given a frequency estimate
  * "est"
  */
  for (int samp = 0; samp < numb_ss; samp++) {
    NumericVector well_1(1, numb_wells(samp, 0));    // numb of wells with chains b, a_1
    NumericVector well_2(1, numb_wells(samp, 1));    // numb of wells with chains b, a_2
    NumericVector well_3(1, numb_wells(samp, 2));    // numb of wells with chains a_1, a_2
    NumericVector well_d(1, numb_wells(samp, 3));    // numb of wells with chains b, a_1, a_2
    NumericVector well_o(1, numb_wells(samp, 4));    // numb of wells not in the above
    NumericVector well_T(1, numb_wells(samp, 0) + numb_wells(samp, 1) +
      numb_wells(samp, 2) + numb_wells(samp, 3) + numb_wells(samp, 4));

    double size = numb_cells[samp]; // number of cells per well

    double prob_2  = 0; // prob of a well to have just two of three chains
    double prob_3  = 0; // prob of a well to have all three chains
    double prob_0  = 0; // prob of a well to have just one or none of the chains

    NumericVector binom = as<NumericVector>(binomials[samp]);

    for (int cell = 1; cell < size + 1; cell++) {
      prob_2 += binom[cell - 1] * pow(est, cell) * pow(1 - est, size - cell) *
        pow(err, cell) * pow(1 - pow(err, cell), 2);
      prob_3 += binom[cell - 1] * pow(est, cell) * pow(1 - est, size - cell) *
        pow(1 - pow(err, cell), 3);
    }
    prob_0 = 1 - 3 * prob_2 - prob_3;

    log_like = log_like + lfactorial(well_T) - lfactorial(well_1) - lfactorial(well_2) -
      lfactorial(well_3) - lfactorial(well_d) - lfactorial(well_o) +
      (well_1 + well_2 + well_3) * log(prob_2) + well_d * log(prob_3) +
      well_o * log(prob_0);
  } // end for - int
  return -log_like;
} // end dual_discrim_dual_likelihood
