#' Estimation of frequencies of clones identified by \code{alphabetr}
#'
#' \code{freq_estimate()} estimates the frequencies of clones with confidence
#'    intervals by using a maximum likelihood approach. The function looks at
#'    the wells that a chains of a clone appear in and determines the most
#'    likely frequency that explains the data.
#'
#' @param alpha Matrix recording which alpha chains appear in each well of the
#'    data. See \code{\link{create_data}}.
#' @param beta Matrix recording which beta chains appear in the each well of the
#'    data. See \code{\link{create_data}}.
#' @param pair A matrix where each row is a beta/alpha pair, column 1 and 2 are
#'    the beta indices, and column 3 and 4 are the alpha indices, and column 5
#'    is the proportion of replicates the clone was found in (or equal to -1 if
#'    the clone is dual)
#' @param error The mean error "dropped" chain rate due to PCR or sequencing
#'    errors.
#' @param numb_cells The number of cells per well in each column of the plates.
#'    Should be a vector of 12 elements.
#'
#' @return A data frame with frequency estimates and confidence intervals
#' @examples
#'  \dontrun{
#'  # obtained from the output of bagpipe()
#'  pairs <- pairs[pairs[, 5] > 0.3, ]
#'  freq  <- freq_estimate(alpha = dat$alpha,
#'                         beta = dat$beta,
#'                         pair = pairs,
#'                         numb_cells = matrix(c(50, 480), ncol = 2))
#'  }
#' @export


freq_estimate <- function(alpha, beta, pair, error = .15, numb_cells) {
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
                    numb_cells = numb_cells) -
      likelihood_dual(est = x, err = error, numb_wells = well_clone,
                      numb_cells = numb_cells) + 1.96
  } # end function - ci.function
  #-------- CI function - single --------#

  # Find the number of wells and then determine the number of plates
  numb_wells <- nrow(alpha)
  if(nrow(beta) != numb_wells)
    stop("Different number of wells in alpha and beta data sets")
  if (nrow(pair) == 0)
    stop("There no alpha-beta pairs to check (the pairs argument is empty).")
  numb_plates <- numb_wells / 96
  # if (all(names(pair)[1:4] != c("beta1", "beta2", "alpha1", "alpha2")) | ncol(pairs) == 5)
  #   stop("The pair argument should be a 5 col matrix, first 2 col for betas,
  # next 2 col for alphas, and 5th col for proportion of replicates.
  # Update alphabetr if bagpipe() isn't outputting this format.")

  # Find the number of unique alpha and beta chains
  numb_alph <- ncol(alpha)
  numb_beta <- ncol(beta)

  freq_results <- matrix(0, nrow = nrow(pair), ncol = 9)         # preallocating matrix to record results

  for (clone in 1:nrow(pair)) {
    # indices of the chains of the clone
    ind_beta1 <- pair[clone, 1]
    ind_beta2 <- pair[clone, 2]
    ind_alph1 <- pair[clone, 3]
    ind_alph2 <- pair[clone, 4]

    # Determing which wells that contain the clone a clone is counted to be in a
    # well if their component chains are found in the same well
    sample_size_well <- numb_cells[, 1]     # number of cells per well
    numb_sample <- numb_cells[, 2]          # number of wells w/ sample size

    numb_distinct <- length(sample_size_well)
    well_clone <- rep(0, numb_distinct)

    # logicals whether the clone has dual alphas and/or dual betas
    dual_alph <- ind_alph1 != ind_alph2
    dual_beta <- ind_beta1 != ind_beta2

    # determining the number of wells of each sample size that all chains appear in
    for (size in 1:numb_distinct) {
      ind1 <- cumsum(numb_sample[1:size])[size] - numb_sample[size] + 1
      ind2 <- cumsum(numb_sample[1:size])[size]
      well_clone[size] <- sum(beta[ind1:ind2, ind_beta1]  == 1 &
                              beta[ind1:ind2, ind_beta2]  == 1 &
                              alpha[ind1:ind2, ind_alph1] == 1 &
                              alpha[ind1:ind2, ind_alph2] == 1)
    }

    # checking to see if all the clones' chains appear in all of the wells;
    # method falls apart when this occurs and # of cells in wells will need to
    # be reduced in a future experiment in estimate these clones
    all_wells <- sum(well_clone) == nrow(alpha)

    # calculating the MLE for the frequency point estimate
    # a different C++ function is called depending on what chain(s) are duals
    if (dual_alph & dual_beta & !all_wells) {
      mle <- optimize(likelihood_dualdual, interval = c(0, .5), maximum = FALSE,
                      err = error, numb_wells = well_clone,
                      numb_cells = sample_size_well, numb_sample = numb_sample)
    } else if (dual_alph & !dual_beta & !all_wells) {
      mle <- optimize(likelihood_dual, interval = c(0, .5), maximum = FALSE,
                      err = error, numb_wells = well_clone,
                      numb_cells = sample_size_well, numb_sample = numb_sample)
    } else if (dual_beta & !dual_alph & !all_wells) {
      mle <- optimize(likelihood_dual, interval = c(0, .5), maximum = FALSE,
                      err = error, numb_wells = well_clone,
                      numb_cells = sample_size_well, numb_sample = numb_sample)
    } else if (!dual_alph & !dual_beta & !all_wells) {
      mle <- optimize(likelihood_single, interval = c(0, 0.5), maximum = FALSE,
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
        if (ci_single(1e-16, mle) * mle$minimum < 0) {
          CI_lower <- uniroot(ci_single, MLE = mle, lower = 1e-16,
                              upper = mle$minimum, tol = .Machine$double.eps)
          CI_lower <- CI_lower$root

        } else {
          CI_lower <- 0
        }
      }
      freq_results[clone, 1] <- ind_beta1
      freq_results[clone, 2] <- ind_beta2
      freq_results[clone, 3] <- ind_alph1
      freq_results[clone, 4] <- ind_alph2
      freq_results[clone, 5] <- mle$minimum
      freq_results[clone, 6] <- CI_upper
      freq_results[clone, 7] <- CI_lower
      freq_results[clone, 8] <- CI_upper - CI_lower
      if (!dual_alph & !dual_beta) {
        freq_results[clone, 9] <- pair[clone, 5]
      } else {
        freq_results[clone, 9] <- -1
      }
    } else {
      freq_results[clone, 1] <- ind_beta1
      freq_results[clone, 2] <- ind_beta2
      freq_results[clone, 3] <- ind_alph1
      freq_results[clone, 4] <- ind_alph2
      freq_results[clone, 5] <- NA
      freq_results[clone, 6] <- NA
      freq_results[clone, 7] <- NA
      freq_results[clone, 8] <- NA
      if (!dual_alph & !dual_beta) {
        freq_results[clone, 9] <- pair[clone, 5]
      } else {
        freq_results[clone, 9] <- -1
      }
    }
    colnames(freq_results) <- c("beta1", "beta2", "alpha1", "alpha2",
                                "MLE", "CI_up", "CI_low", "CI_length",
                                "pct_replicates")

  } # end for clone
  as.data.frame(freq_results)
}
