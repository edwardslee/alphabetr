#' Calculate dual depths and false dual rates for simulated alphabetr experiments
#'
#' \code{dual_eval()} is used in simulation situations to compare the duals
#' determined by \code{\link{dual_top}} and \code{\link{dual_tail}} (which can
#' be combined with \code{rbind()}) to the duals in the simulated T cell
#' population.
#'
#' @param duals A 4 column matrix (col 1 + 2 = beta indices, col 3 + 4 = alpha
#'    indices) containing the indices of dual-alpha clones. The output of
#'    \code{\link{dual_top}} and \code{\link{dual_tail}} are in this form (and
#'    the outputs of these two functions can combined by using \code{rbind()})
#' @param pair The output of \code{\link{bagpipe}}
#' @param TCR The clonal structure of the simulated T cell population. This is
#'    obtained by subsetting the \code{TCR} element of the output of
#'    \code{\link{create_clones}}
#' @param number_skewed The number of clones represent the top proportion of
#'    the T cell population by frequency (this is the same \code{number_skewed}
#'    argument used when \code{\link{create_clones}} is called)
#' @param TCR_dual The dual clones of the simulated T cell population. This is
#'    obtained by subsetting the \code{dual_alph} element of the output of
#'    \code{\link{create_clones}}
#'
#' @return A data.frame with the following columns:
#'    \itemize{
#'      \item \code{fdr}, the false dual rate
#'      \item \code{numb_cand_duals}, the number of duals identified
#'      \item \code{adj_depth_top}, the adjusted dual depth of top clones
#'      \item \code{abs_depth_top}, the absolute dual depth of top clones
#'      \item \code{numb_correct_top}, the number of correctly identified dual
#'         clones in the top
#'      \item \code{numb_duals_ans_top}, the number of top dual clones in the
#'         simulated T cell population
#'      \item \code{numb_poss_top}, the number of top dual clones whose beta and
#'         both alpha chains were identified by \code{bagpipe()}
#'      \item \code{numb_unestimated_top}, number of top dual clones whose
#'         frequencies could not be calculated (usually because the clones
#'         appeared in every well of the data)
#'      \item \code{adj_depth_tail}, the adjusted dual depth of tail clones
#'      \item \code{abs_depth_tail}, the absolute dual depth of tail clones
#'      \item \code{numb_correct_tail}, the number of correctly identified tail
#'         clones
#'      \item \code{numb_duals_ans_tail}, the number of dual tail clones in the
#'         simulated T cell population
#'      \item \code{numb_poss_tail}, the number of tail dual cloens whose beta
#'         and both alpha chains were identified by \code{bagpipe()}
#'      \item \code{numb_unestimated_tail}, the number of tail clones whose
#'         frequencies could not be calculated
#'    }
#'
#' @export
dual_eval <- function(duals, pair, TCR, number_skewed, TCR_dual) {
  # determine top (common) duals and tail (rare) duals
  dual_top <- TCR[which(TCR[1:number_skewed, 3] != TCR[1:number_skewed, 4]), , drop = FALSE]

  # record the top dual indices, remove them to get the dual tail indices
  indices_dual <- vector(length = nrow(dual_top))
  for (i in 1:nrow(dual_top)) {
    ind_beta  <- dual_top[i, "beta1"]
    ind_alph1 <- dual_top[i, "alpha1"]
    ind_alph2 <- dual_top[i, "alpha2"]
    indices_dual[i] <- which(TCR_dual[, 1] == ind_beta &
                               ((TCR_dual[, 3] == ind_alph1 & TCR_dual[, 4] == ind_alph2) |
                                (TCR_dual[, 3] == ind_alph2 & TCR_dual[, 4] == ind_alph1))
    )
  } # end for - i; remove dual clones
  dual_tail <- TCR_dual[-indices_dual, ]    # remove top duals, leaving the tail

  # if the freq estimation fails, then can't do the dual thing
  unest_pairs <- pair[is.na(pair[, "MLE"]), 1:4]

  tail_dual_record <- matrix(nrow = 0, ncol = 3)
  if (nrow(duals) > 0) { # if any duals are determined
    #----- number correct ------#
    # determine the number of correct top and tail dual clones in the list
    numb_correct_top  <- 0
    numb_correct_tail <- 0

    for (i in 1:nrow(duals)) {
      ind_beta  <- duals[i, 1]
      ind_alph1 <- duals[i, 3]
      ind_alph2 <- duals[i, 4]

      # check if clone is a top dual; need to match beta and both alphas
      top_cond  <- any(dual_top[, 1] == ind_beta &
                         ((dual_top[, 3] == ind_alph1 & dual_top[, 4] == ind_alph2) |
                          (dual_top[, 3] == ind_alph2 & dual_top[, 4] == ind_alph1)))
      # check if clone is a tail dual; need to match beta and both alphas
      tail_cond <- any(dual_tail[, 1] == ind_beta &
                         ((dual_tail[, 3] == ind_alph1 & dual_tail[, 4] == ind_alph2) |
                          (dual_tail[, 3] == ind_alph2 & dual_tail[, 4] == ind_alph1)))

      if (top_cond)  numb_correct_top  <- numb_correct_top + 1
      if (tail_cond) {
        numb_correct_tail <- numb_correct_tail + 1
        tail_dual_record  <- rbind(tail_dual_record,
                                   matrix(c(ind_beta, ind_alph1, ind_alph2),
                                          nrow = 1, ncol = 3)
                                   )
      }
    }
    numb_correct <- numb_correct_top + numb_correct_tail
    #----- number correct ------#

    #------ number possible --------#
    numb_poss_top   <- 0
    numb_poss_tail  <- 0
    numb_unest_top  <- 0
    numb_unest_tail <- 0

    for (i in 1:nrow(dual_top)) {
      ind_beta  <- dual_top[i, 1]
      ind_alph1 <- dual_top[i, 3]
      ind_alph2 <- dual_top[i, 4]

      if (any(pair[, 1] == ind_beta & pair[, 3] == ind_alph1) &
          any(pair[, 1] == ind_beta & pair[, 4] == ind_alph2)) {
        numb_poss_top <- numb_poss_top + 1
        # then check to see if the top dual has any unestimateable freq
        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 3] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 3] == ind_alph2)) {
          numb_unest_top <- numb_unest_top + 1
        }
      }
    }

    for (i in 1:nrow(dual_tail)) {
      ind_beta  <- dual_tail[i, 1]
      ind_alph1 <- dual_tail[i, 3]
      ind_alph2 <- dual_tail[i, 4]

      if (any(pair[, 1] == ind_beta & pair[, 3] == ind_alph1) &
          any(pair[, 1] == ind_beta & pair[, 4] == ind_alph2)) {
        numb_poss_tail <- numb_poss_tail + 1
        # then check to see if the top dual has any unestimateable freq
        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 3] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 3] == ind_alph2)) {
          numb_unest_tail <- numb_unest_tail + 1
        }
      }
    }

    # calculating top depths
    if (numb_poss_top == 0) {
      adj_depth_top <- NA
    } else {
      adj_depth_top <- numb_correct_top / numb_poss_top
    }
    abs_depth_top <- numb_correct_top / nrow(dual_top)

    # calculating tail depths
    if (numb_poss_tail == 0) {
      adj_depth_tail <- NA
    } else {
      adj_depth_tail <- numb_correct_tail / numb_poss_tail
    }
    abs_depth_tail <- numb_correct_tail / nrow(dual_tail)

    # false dual rate
    fdr <- (nrow(duals) - numb_correct) / nrow(duals)

    return(tibble::tibble(
      fdr = fdr,
      numb_cand_duals = nrow(duals),
      adj_depth_top = adj_depth_top,
      abs_depth_top = abs_depth_top,
      numb_correct_top = numb_correct_top,
      numb_duals_ans_top = nrow(dual_top),
      numb_poss_top = numb_poss_top,
      numb_unestimated_top = numb_unest_top,
      adj_depth_tail = adj_depth_tail,
      abs_depth_tail = abs_depth_tail,
      numb_correct_tail = numb_correct_tail,
      numb_duals_ans_tail = nrow(dual_tail),
      numb_poss_tail = numb_poss_tail,
      numb_unestimated_tail = numb_unest_tail))
  } else { # no duals determined
    # since no duals were determined, there are no correct top or tail duals
    numb_correct_top  <- 0
    numb_correct_tail <- 0

    #------ number possible --------#
    numb_poss_top   <- 0
    numb_poss_tail  <- 0
    numb_unest_top  <- 0
    numb_unest_tail <- 0

    for (i in 1:nrow(dual_top)) {
      ind_beta  <- TCR_dual[i, 1]
      ind_alph1 <- TCR_dual[i, 3]
      ind_alph2 <- TCR_dual[i, 4]

      # check if these indices represent a top dual
      if (any(dual_top[, 1] == ind_beta & dual_top[, 2] == ind_alph1) &
          any(dual_top[, 1] == ind_beta & dual_top[, 3] == ind_alph2)) {
        numb_poss_top <- numb_poss_top + 1
        # then check to see if the top dual has any unestimateable freq
        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 3] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 3] == ind_alph2)) {
          numb_unest_top <- numb_unest_top + 1
        }
      }

      # check if these indices represent a tail dual
      if (any(dual_tail[, 1] == ind_beta & dual_tail[, 2] == ind_alph1) &
          any(dual_tail[, 1] == ind_beta & dual_tail[, 3] == ind_alph2)) {
        numb_poss_tail <- numb_poss_tail + 1
        # then check to see if the top dual has any unestimateable freq
        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 3] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 3] == ind_alph2)) {
          numb_unest_tail <- numb_unest_tail + 1
        }
      }
    } # endfor - i

    # calculating top depths
    if (numb_poss_top == 0) {
      adj_depth_top <- NA
    } else {
      adj_depth_top <- 0
    }
    abs_depth_top <- 0

    # calculating tail depths
    if (numb_poss_top == 0) {
      adj_depth_tail <- NA
    } else {
      adj_depth_tail <- 0
    }
    abs_depth_tail <- 0

    # false dual rate
    fdr <- NA

    data.frame(
      fdr = fdr,
      numb_cand_duals = nrow(duals),
      adj_depth_top = adj_depth_top,
      abs_depth_top = abs_depth_top,
      numb_correct_top = numb_correct_top,
      numb_duals_ans_top = nrow(dual_top),
      numb_poss_top = numb_poss_top,
      numb_unestimated_top = numb_unest_top,
      adj_depth_tail = adj_depth_tail,
      abs_depth_tail = abs_depth_tail,
      numb_correct_tail = numb_correct_tail,
      numb_duals_ans_tail = nrow(dual_tail),
      numb_poss_tail = numb_poss_tail,
      numb_unestimated_tail = numb_unest_tail)
  } # end if else
}
