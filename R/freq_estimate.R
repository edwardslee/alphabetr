#' Estimation frequencies of clones identified by alphabetr
#'
#' \code{freq_estimate()} estimates the frequencies of clones with confidence
#'    intervals by using a maximum likelihood approach. The function looks at
#'    the wells that a chains of a clone appear in and determines the most
#'    likely frequency that explains the data.
#'
#' @param alpha Matrix recording which alpha chains appear in each well of the
#'    data. See  .
#' @param beta Matrix recording which beta chains appear in the each well of the
#'    data. See .
#' @param pair A matrix where each row is a beta/alpha pair, column 1 is the
#'    beta index, and column 2 is the alpha index.
#' @param error The error or "dropped" chain rate due to PCR or sequencing
#'    errors in the experimental pipeline.
#' @param cells The number of cells per well in each column of the plates.
#'    Should be a vector of 12 elements.
#'
#' @return A data frame with frequency estimates and confidence intervals
#' @examples
#' TCR <- create_clones()
#' data <- create_data()
#' results <- bagpipe(alpha = data$alpha, beta = data$beta, )
#' freq <- freq_estimate(data$alpha, data$beta, results, 0.15, cellstor)
#' @export


freq_estimate <- function(alpha, beta, pair, error = .15, cells) {
  #-------- CI function - single --------#
  ci_single <- function(x, MLE) {
    likelihood_single(est = MLE$minimum, err = error, numb_wells = well_clone,
                      numb_cells = sample_size_well, numb_sample = numb_sample) -
      likelihood_single(est = x, err = error, numb_wells = well_clone,
                        numb_cells = sample_size_well, numb_sample = numb_sample) + 1.96
  } # end function - ci.function
  #-------- CI function - single --------#

  #-------- CI function - single --------#
  ci_dual <- function(x, MLE) {
    likelihood_dual(est = MLE$minimum, err = error, numb_wells = well_clone,
                    numb_cells = cells) -
      likelihood_dual(est = x, err = error, numb_wells = well_clone,
                      numb_cells = cells) + 1.96
  } # end function - ci.function
  #-------- CI function - single --------#

  # Find the number of wells and then determine the number of plates
  numb_wells <- nrow(alpha)
  if(nrow(beta) != numb_wells)
    stop("Different number of wells in alpha and beta data sets")
  numb_plates <- numb_wells / 96

  # Find the number of unique alpha and beta chains
  numb_alph <- ncol(alpha)
  numb_beta <- ncol(beta)

  abundance.results <- matrix(0, nrow = nrow(pair), ncol = 8)         # preallocating matrix to record results

  for (clone in 1:nrow(pair)) {
    # indices of the chains of the clone
    ind_beta  <- pair[clone, 1]
    ind_alph1 <- pair[clone, 2]
    ind_alph2 <- pair[clone, 3]

    # Determing which wells that contain the clone a clone is counted to be in a
    # well if their component chains are found in the same well
    sample_size_well <- cells[, 1]     # number of cells per well
    numb_sample <- cells[, 2]          # number of wells w/ sample size

    numb_distinct <- length(sample_size_well)
    well_clone <- rep(0, numb_distinct)
    if (ind_alph1 == ind_alph2 || ncol(pair) == 3) {
      for (size in 1:numb_distinct) {
        ind1 <- cumsum(numb_sample[1:size])[size] - numb_sample[size] + 1
        ind2 <- cumsum(numb_sample[1:size])[size]
        well_clone[size] <- sum(beta[ind1:ind2, ind_beta] == 1 &
                                alpha[ind1:ind2, ind_alph1] == 1)
      }
    } else {
      for (size in 1:numb_distinct) {
        ind1 <- cumsum(numb_sample[1:size])[size] - numb_sample[size] + 1
        ind2 <- cumsum(numb_sample[1:size])[size]
        well_clone[size] <- sum(beta[ind1:ind2, ind_beta] == 1 &
                                  alpha[ind1:ind2, ind_alph1] == 1 &
                                  alpha[ind1:ind2, ind_alph2] == 1)
      }
    }


    all_wells <- sum(well_clone) == nrow(alpha)

    if ((ind_alph1 == ind_alph2 | ncol(pair) == 3) & !all_wells) {
      mle <- optimize(likelihood_single, interval = c(0, .5), maximum = FALSE,
                      err = error, numb_wells = well_clone,
                      numb_cells = sample_size_well, numb_sample = numb_sample)
    } else if(!all_wells) {
      mle <- optimize(likelihood_dual, interval = c(0, .5), maximum = FALSE,
                      err = error, numb_wells = well_clone,
                      numb_cells = sample_size_well, numb_sample = numb_sample)
    }

    if (!all_wells) {
      # Calculate upper bound for the 95% Confidence interval
      CI_upper <- 0
      if (ci_single(x = mle$minimum, MLE = mle)*ci_single(x = .9, MLE = mle) < 0){
        CI_upper  <- uniroot(ci_single, MLE = mle, lower = mle$minimum,
                             upper = .9, tol = .Machine$double.eps)
        CI_upper <- CI_upper$root
      }
      # Calculate lower bound for the 95% confidence interval
      CI_lower <- 0
      unbound_ci <- likelihood_single(0.01*mle$minimum, err = error,
                                      numb_wells = well_clone,
                                      numb_cells = sample_size_well, numb_sample = numb_sample) == -Inf
      if (unbound_ci) {
        CI.lower <- 0
      } else {
        # browser()
        if (ci_single(1e-16, mle) * mle$minimum < 0) {
          CI_lower <- uniroot(ci_single, MLE = mle, lower = 1e-16,
                              upper = mle$minimum, tol = .Machine$double.eps)
          CI_lower <- CI_lower$root

        } else {
          CI_lower <- 0
        }
      }
      abundance.results[clone, 1] <- ind_beta
      abundance.results[clone, 2] <- ind_alph1
      abundance.results[clone, 3] <- ind_alph2
      abundance.results[clone, 4] <- mle$minimum
      abundance.results[clone, 5] <- CI_upper
      abundance.results[clone, 6] <- CI_lower
      abundance.results[clone, 7] <- CI_upper - CI_lower
      if (ind_alph1 == ind_alph2) {
        abundance.results[clone, 8] <- pair[clone, 4]
      } else {
        abundance.results[clone, 8] <- -1
      }
    } else {
      abundance.results[clone, 1] <- ind_beta
      abundance.results[clone, 2] <- ind_alph1
      abundance.results[clone, 3] <- ind_alph2
      abundance.results[clone, 4] <- NA
      abundance.results[clone, 5] <- NA
      abundance.results[clone, 6] <- NA
      abundance.results[clone, 7] <- NA
      if (ind_alph1 == ind_alph2) {
        abundance.results[clone, 8] <- pair[clone, 4]
      } else {
        abundance.results[clone, 8] <- -1
      }
    }
    colnames(abundance.results) <- c("beta", "alpha1", "alpha2", "MLE", "CI_up",
                                     "CI_low", "CI_length", "pct_replicates")

  } # end for clone
  as.data.frame(abundance.results)
}
