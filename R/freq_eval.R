#' Calculate the precision, CV, and accuracy of frequency estimates
#'
#' \code{freq_eval()} will evaluated how well \code{\link{freq_estimate}} performed
#' by calculating the precision and CV of the frequency estimates for the top
#' clones and by determining the proportion of the top clones whose true clonal
#' frequency lies in the 95-percent CI determined by \code{\link{freq_estimate}}
#'
#' @param freq The output of \code{\link{freq_estimate}}
#' @param number_skewed The number of clones represent the top proportion of
#'    the T cell population by frequency (this is the same \code{number_skewed}
#'    argument used when \code{\link{create_clones}} is called)
#' @param TCR The clonal structure of the simulated T cell population. This is
#'    obtained by subsetting the \code{TCR} element of the output of
#'    \code{\link{create_clones}}
#' @param numb_clones Total number of distinct clones in the parent population
#' @param prop_top The proportion of the population in frequency represented by
#'    the number of clones specified by \code{skewed}.
#'
#' @return A list with the precision, cv, and accuracy of the frequency estimation.
#' @export
freq_eval <- function(freq, number_skewed, TCR, numb_clones, prop_top) {
  numb_clones <- nrow(TCR)
  last_term <- (1 - prop_top)/(numb_clones - number_skewed) * 1.1
  lhs11 <- 1
  lhs12 <- number_skewed - 1
  lhs21 <- (1 + number_skewed - 1)
  lhs22 <- (number_skewed - 1)/2 + (number_skewed - 1)^2/2
  A <- matrix(c(lhs11, lhs21, lhs12, lhs22), nrow = 2)
  b <- matrix(c(last_term, prop_top), nrow = 2)
  solveAb <- solve(A,b)
  prob1 <- solveAb[1] + 0:(number_skewed-1) * solveAb[2]
  prob2 <- rep((1 - prop_top)/(numb_clones - number_skewed), numb_clones - number_skewed)
  dist_vector <- c(prob1, prob2)

  TCR <- cbind(TCR, dist_vector)

  top_clones <- TCR[1:(number_skewed), ]
  correct_top <- vector()
  for (i in 1:nrow(freq)) {
    ind_beta  <- freq[i, "beta1"]
    ind_alph1 <- freq[i, "alpha1"]
    ind_alph2 <- freq[i, "alpha2"]

    if (ind_alph1 == ind_alph2) {
      if (any(top_clones[, 1] == ind_beta &
             (top_clones[, 3] == ind_alph1 | top_clones[, 4] == ind_alph1)))
        correct_top <- c(correct_top, i)
    } else {
      if (any(top_clones[, 1] == ind_beta &
              ((top_clones[, 3] == ind_alph1 & top_clones[, 4] == ind_alph2) |
               (top_clones[, 3] == ind_alph2 & top_clones[, 4] == ind_alph1))
             )
      )
        correct_top <- c(correct_top, i)
    }
  }

  freq <- freq[correct_top,, drop = FALSE ]
  if (nrow(freq) != 0) {
    freq$answer  <- 0
    freq$correct <- 0
  } else {
    freq <- as.data.frame(matrix(nrow = 0, ncol = 11))
    names(freq) <- c("beta1", "beta2", "alpha1", "alpha2", "MLE",
                     "CI_up", "CI_low", "CI_length", "pct_replicates", "answer",
                     "correct")
  }

  for (i in seq_len(nrow(freq))) {
    ind_beta  <- freq[i, "beta1"]
    ind_alph1 <- freq[i, "alpha1"]
    ind_alph2 <- freq[i, "alpha2"]
    if (!is.na(freq[i, "MLE"])) { #check to see if freq was estimatable at all
      if (ind_alph1 == ind_alph2) {
        ind <- which(top_clones[, 1] == ind_beta &
                       (top_clones[, 3] == ind_alph1 | top_clones[, 4] == ind_alph1))
        freq[i, "answer"] <- TCR[ind[1], 4]  #in situation where two dual clones have the same beta and alph1 but diff alph2s
        if (TCR[ind, 5] <= freq[i, "CI_up"] & TCR[ind, 5] >= freq[i, "CI_low"]) {
          freq[i, "correct"] <- 1
        }
      } else {
        ind <- which(top_clones[, 1] == ind_beta &
                       ((top_clones[, 3] == ind_alph1 & top_clones[, 4] == ind_alph2) |
                        (top_clones[, 3] == ind_alph2 & top_clones[, 4] == ind_alph1)))
        ind <- ind[1]
        if (TCR[ind, 5] <= freq[i, "CI_up"] & TCR[ind, 5] >= freq[i, "CI_low"]) {
          freq[i, "correct"] <- 1
          freq[i, "answer"] <- TCR[ind, 5]
        }
      } # end if else - dual or not
    } # end if - check if estimatable
  } # end for - i
  numb_top <- nrow(freq)
  freq <- dplyr::mutate(
            dplyr::filter(freq, correct == 1),
              precision = CI_length/MLE, CV = CI_length/MLE/3.92
            )
  if (nrow(freq) == 0) {
    list(precision = NA, cv = NA, correct = 0)
  } else {
    list(precision = mean(freq$precision), cv = mean(freq$CV), correct = sum(freq$correct)/numb_top)
  }
}
