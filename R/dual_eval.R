#' @export
dual_eval <- function(duals, prop, pair, TCR_sizes, number_skewed, TCR_dual) {
  if (nrow(duals) > 0) {
    unest_pairs <- pair[is.na(pair[, "MLE"]), 1:3, drop = FALSE]
    #-------------------------#
    #        Top clones       #
    #-------------------------#
    if (prop == "top") {
      dual_top <- TCR_sizes[which(TCR_sizes[1:number_skewed, 2] != TCR_sizes[1:number_skewed, 3]), , drop = FALSE]

      #------------ number correct -------------#
      # determine the number of correct dual top clones
      numb_correct <- 0
      for (i in 1:nrow(duals)) {
        ind_beta  <- duals[i, 1]
        ind_alph1 <- duals[i, 2]
        ind_alph2 <- duals[i, 3]

        if (any(TCR_dual[, 1] == ind_beta &
                ((TCR_dual[, 2] == ind_alph1 & TCR_dual[, 3] == ind_alph2) |
                 (TCR_dual[, 2] == ind_alph2 & TCR_dual[, 3] == ind_alph1))
        )) {
          numb_correct <- numb_correct + 1
        }
      } # end for - i
      #------------ number correct -------------#

      #------------ number possible ------------#
      numb_possible <- 0
      numb_dual_unest <- 0
      for (i in 1:nrow(dual_top)) {
        ind_beta  <- dual_top[i, 1]
        ind_alph1 <- dual_top[i, 2]
        ind_alph2 <- dual_top[i, 3]

        if (any(pair[, 1] == ind_beta & pair[, 2] == ind_alph1) &
            any(pair[, 1] == ind_beta & pair[, 3] == ind_alph2)) {
          numb_possible <- numb_possible + 1
        }

        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph2)
        ) {
          numb_dual_unest <- numb_dual_unest + 1
        }
      }
      adj_depth <- numb_correct / numb_possible
      abs_depth <- numb_correct / nrow(dual_top)
      fdr <- (nrow(duals) - numb_correct) / nrow(duals)
      #------------ number possible ------------#

      return(list(adj_depth = adj_depth, abs_depth = abs_depth, fdr = fdr,
                  numb_correct = numb_correct, numb_duals_ans = nrow(dual_top),
                  numb_cand_duals = nrow(duals), numb_duals_poss = numb_possible,
                  numb_duals_unest = numb_dual_unest))


    #---------------------------#
    #        Tail clones        #
    #---------------------------#
    } else if (prop == "tail") {
      # remove dual clones from the matrix of dual clones
      dual_top <- TCR_sizes[which(TCR_sizes[1:number_skewed, 2] !=
                            TCR_sizes[1:number_skewed, 3]),, drop = FALSE]
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
      dual_tail <- TCR_dual[-indices_dual, ]

      #------------ number possible ------------#
      # determine the number of possible
      numb_possible <- 0
      numb_dual_unest <- 0
      # poss_dual <- matrix(nrow = 0, ncol = 3)
      for (i in 1:nrow(dual_tail)) {
        ind_beta  <- dual_tail[i, 1]
        ind_alph1 <- dual_tail[i, 2]
        ind_alph2 <- dual_tail[i, 3]
        if (any(pair[, 1] == ind_beta & (pair[, 2] == ind_alph1 | pair[, 3] == ind_alph1)) &
            any(pair[, 1] == ind_beta & (pair[, 2] == ind_alph2 | pair[, 3] == ind_alph2))
        ){
          # poss_dual <- rbind(poss_dual, c(ind_beta, ind_alph1, ind_alph2))
          numb_possible <- numb_possible + 1
        }
        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph2)
        ) {
          numb_dual_unest <- numb_dual_unest + 1
        }
      } # end for - i; numb_possible
      #------------ number possible ------------#

      #------------ number correct -------------#
      # determine the number of correct dual top clones
      numb_correct <- 0
      for (i in 1:nrow(duals)) {
        ind_beta  <- duals[i, 1]
        ind_alph1 <- duals[i, 2]
        ind_alph2 <- duals[i, 3]

        if (any(TCR_dual[, 1] == ind_beta &
                ((TCR_dual[, 2] == ind_alph1 & TCR_dual[, 3] == ind_alph2) |
                 (TCR_dual[, 2] == ind_alph2 & TCR_dual[, 3] == ind_alph1))
        )) {
          numb_correct <- numb_correct + 1
        }
      } # end for - i
      #------------ number correct -------------#

      adj_depth <- numb_correct / numb_possible
      abs_depth <- numb_correct / nrow(dual_tail)
      fdr <- (nrow(duals) - numb_correct) / nrow(duals)

      return(list(adj_depth = adj_depth, abs_depth = abs_depth, fdr = fdr,
                  numb_correct = numb_correct, numb_duals_ans = nrow(dual_tail),
                  numb_cand_duals = nrow(duals), numb_duals_poss = numb_possible,
                  numb_duals_unest = numb_dual_unest))
    } # end if - top/tail
  } else {
    unest_pairs <- pair[is.na(pair[, "MLE"]), 1:3]
    if (prop == "top") {
      dual_top <- TCR_sizes[which(TCR_sizes[1:number_skewed, 2] !=
                                  TCR_sizes[1:number_skewed, 3]),, drop = FALSE]
      numb_possible <- 0
      numb_dual_unest <- 0
      for (i in 1:nrow(dual_top)) {
        ind_beta  <- dual_top[i, 1]
        ind_alph1 <- dual_top[i, 2]
        ind_alph2 <- dual_top[i, 3]

        if (any(pair[, 1] == ind_beta & pair[, 2] == ind_alph1) &
            any(pair[, 1] == ind_beta & pair[, 3] == ind_alph2)) {
          numb_possible <- numb_possible + 1
        }
        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph2)
        ) {
          numb_dual_unest <- numb_dual_unest + 1
        }
      }
      return(list(adj_depth = 0, abs_depth = 0, fdr = 0,
                  numb_correct = 0, numb_duals_ans = nrow(dual_top),
                  numb_cand_duals = nrow(duals), numb_duals_poss = numb_possible,
                  numb_duals_unest = numb_dual_unest))
    } else if (prop == "tail") {
      # remove dual clones from the matrix of dual clones
      dual_top <- TCR_sizes[which(TCR_sizes[1:number_skewed, 2] != TCR_sizes[1:number_skewed, 3]), , drop = FALSE]
      indices_dual <- vector(length = nrow(dual_top))
      for (i in 1:nrow(dual_top)) {
        ind_beta  <- dual_top[i, 1]
        ind_alph1 <- dual_top[i, 2]
        ind_alph2 <- dual_top[i, 3]
        indices_dual[i] <- which(TCR_dual[, 1] == ind_beta &
                    ((TCR_dual[, 2] == ind_alph1 & TCR_dual[, 3] == ind_alph2) |
                     (TCR_dual[, 2] == ind_alph2 & TCR_dual[, 3] == ind_alph1)))
      } # end for - i; remove dual clones
      dual_tail <- TCR_dual[-indices_dual, ]

      numb_possible <- 0
      numb_dual_unest <- 0
      # poss_dual <- matrix(nrow = 0, ncol = 3)
      for (i in 1:nrow(dual_tail)) {
        ind_beta  <- dual_tail[i, 1]
        ind_alph1 <- dual_tail[i, 2]
        ind_alph2 <- dual_tail[i, 3]
        if (any(pair[, 1] == ind_beta & (pair[, 2] == ind_alph1 | pair[, 3] == ind_alph1)) &
            any(pair[, 1] == ind_beta & (pair[, 2] == ind_alph2 | pair[, 3] == ind_alph2))
        ){
          # poss_dual <- rbind(poss_dual, c(ind_beta, ind_alph1, ind_alph2))
          numb_possible <- numb_possible + 1
        }
        if (any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph1) |
            any(unest_pairs[, 1] == ind_beta & unest_pairs[, 2] == ind_alph2)
        ) {
          numb_dual_unest <- numb_dual_unest + 1
        }
      }

      return(list(adj_depth = 0, abs_depth = 0, fdr = 0,
                  numb_correct = 0, numb_duals_ans = nrow(dual_tail),
                  numb_cand_duals = nrow(duals), numb_duals_poss = numb_possible,
                  numb_duals_unest = numb_dual_unest))
    }
  }
}