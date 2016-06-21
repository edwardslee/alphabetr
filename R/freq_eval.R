#' @export
freq_eval <- function(freq, number_skewed, TCR_sizes, numb_clones, prop_top) {
  numb_clones <- nrow(TCR_sizes)
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

  TCR_sizes <- cbind(TCR_sizes, dist_vector)

  top_clones <- TCR_sizes[1:(number_skewed), ]
  correct_top <- vector()
  for (i in 1:nrow(freq)) {
    ind_beta  <- freq[i, "beta"]
    ind_alph1 <- freq[i, "alpha1"]
    ind_alph2 <- freq[i, "alpha2"]

    if (ind_alph1 == ind_alph2) {
      if (any(top_clones[, 1] == ind_beta &
              (top_clones[, 2] == ind_alph1 | top_clones[, 3] == ind_alph1)))
        correct_top <- c(correct_top, i)
    } else {
      if (any(top_clones[, 1] == ind_beta &
              ((top_clones[, 2] == ind_alph1 & top_clones[, 3] == ind_alph2) |
              (top_clones[, 2] == ind_alph2 & top_clones[, 3] == ind_alph1))
             )
      )
        correct_top <- c(correct_top, i)
    }
  }

  freq <- freq[correct_top,, drop = FALSE ]
  freq$answer  <- 0
  freq$correct <- 0

  for(i in 1:nrow(freq)) {
    ind_beta  <- freq[i, "beta"]
    ind_alph1 <- freq[i, "alpha1"]
    ind_alph2 <- freq[i, "alpha2"]
    if(!is.na(freq[i, "MLE"])) { #check to see if freq was estimatable at all
      if (ind_alph1 == ind_alph2) {
        ind <- which(top_clones[, 1] == ind_beta &
                       (top_clones[, 2] == ind_alph1 | top_clones[, 3] == ind_alph1))
        freq[i, "answer"] <- TCR_sizes[ind, 4]
        if (TCR_sizes[ind, 4] <= freq[i, "CI_up"] & TCR_sizes[ind, 4] >= freq[i, "CI_low"]) {
          freq[i, "correct"] <- 1
        }
      } else {
        ind <- which(top_clones[, 1] == ind_beta &
                       ((top_clones[, 2] == ind_alph1 & top_clones[, 3] == ind_alph2) |
                          (top_clones[, 2] == ind_alph2 & top_clones[, 3] == ind_alph1)))
        if (TCR_sizes[ind, 4] <= freq[i, "CI_up"] & TCR_sizes[ind, 4] >= freq[i, "CI_low"]) {
          freq[i, "correct"] <- 1
          freq[i, "answer"] <- TCR_sizes[ind, 4]
        }
      } # end if else - dual or not
    } # end if - check if estimatable
  } # end for - i
  numb_top <- nrow(freq)
  freq <- dplyr::mutate(
            dplyr::filter(freq, correct == 1),
              precision = CI_length/MLE, CV = CI_length/MLE/3.92
            )
  list(precision = mean(freq$precision), cv = mean(freq$CV), correct = sum(freq$correct)/numb_top)
}
