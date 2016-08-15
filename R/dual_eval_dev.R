#' @export
dual_eval_dev <- function(duals, pair, TCR_sizes, number_skewed, TCR_dual) {
  # determine top (common) duals and tail (rare) duals
  dual_top <- TCR_sizes[which(TCR_sizes[1:number_skewed, 2] != TCR_sizes[1:number_skewed, 3]), , drop = FALSE]

  # record the top dual indices, remove them to get the dual tail indices
  indices_dual <- vector(length = nrow(dual_top))
  for (i in 1:nrow(dual_top)) {
    ind_beta  <- dual_top[i, 1]
    ind_alph1 <- dual_top[i, 2]
    ind_alph2 <- dual_top[i, 3]
    indices_dual[i] <- which(TCR_dual[, 1] == ind_beta &
                               ((TCR_dual[, 2] == ind_alph1 & TCR_dual[, 3] == ind_alph2) |
                                (TCR_dual[, 2] == ind_alph2 & TCR_dual[, 3] == ind_alph1))
    )
  } # end for - i; remove dual clones
  dual_tail <- TCR_dual[-indices_dual, ]    # remove top duals, leaving the tail

  # if the freq estimation fails, then can't do the dual thing
  unest_pairs <- pair[is.na(pair[, "MLE"]), 1:3]

  tail_dual_record <<- matrix(nrow = 0, ncol = 3)
  if (nrow(duals) > 0) { # if any duals are determined
    #----- number correct ------#
    # determine the number of correct top and tail dual clones in the list
    numb_correct_top  <- 0
    numb_correct_tail <- 0

    for (i in 1:nrow(duals)) {
      ind_beta  <- duals[i, 1]
      ind_alph1 <- duals[i, 2]
      ind_alph2 <- duals[i, 3]

      # check if clone is a top dual; need to match beta and both alphas
      top_cond  <- any(dual_top[, 1] == ind_beta &
                      ((dual_top[, 2] == ind_alph1 & dual_top[, 3] == ind_alph2) |
                       (dual_top[, 2] == ind_alph2 & dual_top[, 3] == ind_alph1)))
      # check if clone is a tail dual; need to match beta and both alphas
      tail_cond <- any(dual_tail[, 1] == ind_beta &
                      ((dual_tail[, 2] == ind_alph1 & dual_tail[, 3] == ind_alph2) |
                       (dual_tail[, 2] == ind_alph2 & dual_tail[, 3] == ind_alph1)))

      if (top_cond)  numb_correct_top  <- numb_correct_top + 1
      if (tail_cond) {
        numb_correct_tail <- numb_correct_tail + 1
        tail_dual_record <<- rbind(tail_dual_record,
                                   matrix(c(ind_beta, ind_alph1, ind_alph2), nrow = 1, ncol = 3))
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
      ind_alph1 <- dual_top[i, 2]
      ind_alph2 <- dual_top[i, 3]

      if (any(pair[, 1] == ind_beta & pair[, 2] == ind_alph1) &
          any(pair[, 1] == ind_beta & pair[, 3] == ind_alph2)) {
        numb_poss_top <- numb_poss_top + 1
        # then check to see if the top dual has any unestimateable freq
        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph2)) {
          numb_unest_top <- numb_unest_top + 1
        }
      }
    }

    for (i in 1:nrow(dual_tail)) {
      ind_beta  <- dual_tail[i, 1]
      ind_alph1 <- dual_tail[i, 2]
      ind_alph2 <- dual_tail[i, 3]

      if (any(pair[, 1] == ind_beta & pair[, 2] == ind_alph1) &
          any(pair[, 1] == ind_beta & pair[, 3] == ind_alph2)) {
        numb_poss_tail <- numb_poss_tail + 1
        # then check to see if the top dual has any unestimateable freq
        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph2)) {
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
    #------ number possible --------#
    numb_poss_top   <- 0
    numb_poss_tail  <- 0
    numb_unest_top  <- 0
    numb_unest_tail <- 0

    for (i in 1:nrow(dual_top)) {
      ind_beta  <- TCR_dual[i, 1]
      ind_alph1 <- TCR_dual[i, 2]
      ind_alph2 <- TCR_dual[i, 3]

      # check if these indices represent a top dual
      if (any(dual_top[, 1] == ind_beta & dual_top[, 2] == ind_alph1) &
          any(dual_top[, 1] == ind_beta & dual_top[, 3] == ind_alph2)) {
        numb_poss_top <- numb_poss_top + 1
        # then check to see if the top dual has any unestimateable freq
        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph2)) {
          numb_unest_top <- numb_unest_top + 1
        }
      }

      # check if these indices represent a tail dual
      if (any(dual_tail[, 1] == ind_beta & dual_tail[, 2] == ind_alph1) &
          any(dual_tail[, 1] == ind_beta & dual_tail[, 3] == ind_alph2)) {
        numb_poss_tail <- numb_poss_tail + 1
        # then check to see if the top dual has any unestimateable freq
        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph2)) {
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
  } # end if else
}
