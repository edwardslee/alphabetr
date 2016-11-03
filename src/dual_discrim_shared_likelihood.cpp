#include <Rcpp.h>
using namespace Rcpp;
//' Calculate likelihood of two beta-sharing candidate alpha-beta pairs deriving from a dual clone
//'
//' \code{dual_discrim_shared_likelihood()} is used within \code{\link{dual_top}}
//' to calculate the likelihood that two alpha-beta pairs identified by
//' \code{\link{bagpipe}} sharing the same beta chain derive from a two
//' distinct clones sharing the same beta chain dual-alpha clone (instead of
//' a single dual-alpha clone)
//'
//' @param est1 Frequency estimate of the first alpha-beta clone
//' @param est2 Frequency estimate of the second alpha-beta clone
//' @param err Mean drop error rate
//' @param numb_cells Vector containing the number of cells per well
//' @param numb_wells Vector containing the number of wells with the sample
//'    sizes given by \code{numb_cells}
//' @param binomials Calculations of the needed binomial coefficients; this is
//'    faster in R than in Rcpp (from my own tests)
//' @param multinomials Calculations of the needed multinomial coefficients;
//'    this is way faster in R due to vectorization
//'
//' @return A numeric containing the negative log likelihood
//' @export
// [[Rcpp::export]]
NumericVector dual_discrim_shared_likelihood(double est1, double est2, double err, NumericVector numb_cells, NumericMatrix numb_wells, List binomials, List multinomials) {
  NumericVector log_like = NumericVector::create(0.0);                // log likelihood value
  int numb_ss = numb_wells.nrow();     // number of distinct sample size values

  for (int samp = 0; samp < numb_ss; samp++) {    //
    NumericVector well_1(1, numb_wells(samp, 0));    // numb of wells with chains b, a_1
    NumericVector well_2(1, numb_wells(samp, 1));    // numb of wells with chains b, a_2
    NumericVector well_3(1, numb_wells(samp, 2));    // numb of wells with chains a_1, a_2
    NumericVector well_d(1, numb_wells(samp, 3));    // numb of wells with chains b, a_1, a_2
    NumericVector well_o(1, numb_wells(samp, 4));    // numb of wells not in the above
    NumericVector well_T(1, numb_wells(samp, 0) + numb_wells(samp, 1) +
      numb_wells(samp, 2) + numb_wells(samp, 3) + numb_wells(samp, 4));

    int size = numb_cells[samp]; // number of cells per well

    // recording the probabilities of each individual case
    double prob_ba1   = 0;    // case when only beta and alpha1 found in wells
    double prob_ba2   = 0;    // case when only beta and alpha2 found in wells
    double prob_a1a2  = 0;    // case when only alpha1 and alpha2 found in wells
    double prob_ba1a2 = 0;    // case when beta alpha1 and alpha2 found in wells
    double prob_0     = 0;    // every other case

    NumericVector binom = as<NumericVector>(binomials[samp]);
    NumericMatrix multinom = as<NumericMatrix>(multinomials[samp]);
    for (int cell = 1; cell < size + 1; cell++) {
      prob_ba1 += binom[cell - 1] * pow(est1, cell) * pow(1 - est1 - est2, size - cell) * (1 - pow(err, cell)) * (1 - pow(err, cell));
      prob_ba2 += binom[cell - 1] * pow(est2, cell) * pow(1 - est1 - est2, size - cell) * (1 - pow(err, cell)) * (1 - pow(err, cell));
    }
    for (int n1 = 1; n1 < size; n1++) {
      for (int n2 = 1; n2 <= size - n1;n2 ++) {
        prob_ba1 += multinom(n1 - 1, n2 - 1) * pow(est1, n1) * pow(est2, n2) * pow(1 - est1 - est2, size - n1 - n2) *
          (pow(err, n2) * (pow(err, n1) - 1) * (pow(err, n2) - 1));

      }
    }
    for (int n2 = 1; n2 < size; n2++) {
      for (int n1 = 1; n1 <= size - n2; n1++) {
        prob_ba2 += multinom(n1 - 1, n2 - 1) * pow(est2, n2) * pow(est1, n1) * pow(1 - est1 - est2, size - n1 - n2) *
          (pow(err, n1) * (pow(err, n1) - 1) * (pow(err, n2) - 1));
      }
    }

    for (int n1 = 1; n1 < size; n1++) {
      for (int n2 = 1; n2 <= size - n1; n2++) {
        prob_a1a2 += multinom(n1 - 1, n2 - 1) * pow(est1, n1) * pow(est2, n2) * pow(1 - est1 - est2, size - n1 - n2) *
          pow(err, n1) * (1 - pow(err, n1)) * pow(err, n2) * (1 - pow(err, n2));
      }
    }

    for (int n1 = 1; n1 < size; n1++) {
      for (int n2 = 1; n2 <= size - n1; n2++) {
        prob_ba1a2 += multinom(n1 - 1, n2 - 1) * pow(est1, n1) * pow(est2, n2) * pow(1 - est1 - est2, size - n1 - n2) *
          (pow(err, n1) * (1 - pow(err, n2)) * (1 - pow(err, n2)) + (1 - pow(err, n1)) * (1 - pow(err, n1)) *
          (1 - pow(err, n2)) * (1 - pow(err, n2)) +  pow(err, n2) * (1 - pow(err, n1)) *(1 - pow(err, n1)));
      }
    }

    prob_0 = 1 - prob_ba1 - prob_ba2 - prob_a1a2 - prob_ba1a2;

    log_like += lfactorial(well_T) - lfactorial(well_1) - lfactorial(well_2) -
      lfactorial(well_3) - lfactorial(well_d) - lfactorial(well_o) +
      well_1 * log(prob_ba1) + well_2 * log(prob_ba2) + well_3 * log(prob_a1a2) +
      well_d * log(prob_ba1a2) + well_o * log(prob_0);
  } // end for - int
  return -log_like;
}
