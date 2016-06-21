#' @export
dual_top <- function(alpha, beta, pair, error, cells) {
  number_plates <- nrow(alpha)/96  # number of plates
  max_beta <- ncol(beta)           # determine maximum beta index

  sample_size_well <- cells[, 1]     # number of cells per well
  numb_sample <- cells[, 2]          # number of wells w/ sample size

  # Determine the wells with the small sample sizes, which is necessary to
  # determine top clone duals properly
  small <- which(sample_size_well < 50)       # find the sample sizes < 50 cells
  if (length(small) == 0 ) return(matrix(nrow = 0, ncol = 3))

  if (length(small) > 1) {
    if (any(small == 1)) {
      small_wells <- c(1:numb_sample[small[1]])
      for (i in 2:length(small)) {
        j <- small[i]
        small_wells <- c(small_wells,
                         (sum(numb_sample[1:(j-1)]) + 1):sum(numb_sample[1:j])
        )
      }
    } else {
      small_wells <- vector()
      for (i in 1:length(small)) {
        j <- small[i]
        small_wells <- c(small_wells,
                         (sum(numb_sample[1:(j-1)]) + 1):sum(numb_sample[1:j])
        )
      }
    }
  } else {
    small_wells <- 1:numb_sample[small]
  }

  # look at only the wells with the small sample size
  alpha <- alpha[small_wells, ]
  beta <- beta[small_wells, ]

  # pre-allocate matrix to record the indices of the candidate dual TCR clones,
  # the shared and dual likelihoods, and the estimated frequency of the
  # candidate dual
  rec <- matrix(nrow = 0, ncol = 6)
  colnames(rec) <- c("beta", "alpha1", "alpha2", "shared_LL",
                     "dual_LL", "dual_freq")

  #
  freq <- pair
  freq <- freq[freq[, 8] > .8,]
  freq <- freq[order(freq[, "MLE"], decreasing = TRUE), ]
  freq <- freq[!is.na(freq[, "MLE"]), ]
  for (clon in 1:max_beta) {
    # browser()
    x <- freq[freq[, "beta"] == clon, , drop = FALSE]  # find clones with the beta index
    numb_cand <- nrow(x)                            # find number of alphas associated with beta
    if (numb_cand > 1) {                            # if more than 1 alpha
      combos <- combn(numb_cand, 2)                   # find all combos of pairs
      for (ind in 1:ncol(combos)) {                   # check each combo
        ind_beta <- clon
        ind_alph1 <- x[combos[1, ind], "alpha1"]
        ind_alph2 <- x[combos[2, ind], "alpha1"]
        f1 <- x[combos[1, ind], "MLE"]                # freq est of beta alpha1
        f2 <- x[combos[2, ind], "MLE"]                # freq est of beta alpha2
        exp_plates <- 0                               # expected number of plates


        # Determing which wells that contain the clone a clone is counted to be in a
        # well if their component chains are found in the same well
        sample_size_well <- cells[, 1]     # number of cells per well
        numb_sample <- cells[, 2]          # number of wells w/ sample size

        sample_size_well <- sample_size_well[small]
        numb_sample <- numb_sample[small]

        numb_distinct <- length(sample_size_well)
        well_clone <- rep(0, numb_distinct)

        well_clone <- matrix(ncol = 5, nrow = numb_distinct)
        colnames(well_clone) <- c("ba1", "ba2", "a1a2", "ba1a2", "none")


        for (size in 1:numb_distinct) {
          ind1 <- cumsum(numb_sample[1:size])[size] - numb_sample[size] + 1
          ind2 <- cumsum(numb_sample[1:size])[size]
          well_clone[size, "ba1"] <- sum(beta[ind1:ind2, ind_beta] == 1 &
                                           alpha[ind1:ind2, ind_alph1] == 1 &
                                           alpha[ind1:ind2, ind_alph2] == 0)
          well_clone[size, "ba2"] <- sum(beta[ind1:ind2, ind_beta] == 1 &
                                           alpha[ind1:ind2, ind_alph1] == 0 &
                                           alpha[ind1:ind2, ind_alph2] == 1)
          well_clone[size, "a1a2"] <- sum(beta[ind1:ind2, ind_beta] == 0 &
                                            alpha[ind1:ind2, ind_alph1] == 1 &
                                            alpha[ind1:ind2, ind_alph2] == 1)
          well_clone[size, "ba1a2"] <- sum(beta[ind1:ind2, ind_beta] == 1 &
                                             alpha[ind1:ind2, ind_alph1] == 1 &
                                             alpha[ind1:ind2, ind_alph2] == 1)
          well_clone[size, "none"] <- numb_sample[size] - sum(well_clone[size, 1:4])
        }


        binomial_coeff <- list()
        for (i in 1:numb_distinct) {
          binomial_coeff[[i]] <- choose(sample_size_well[i], 1:sample_size_well[i])
        }

        multinomial_coeff <- list()
        for(i in 1:numb_distinct) {
          cells_per_well <- sample_size_well[i]
          multi_mat <- matrix(0, nrow = cells_per_well, ncol = cells_per_well)
          for (j in 1:(cells_per_well-1)) {
            for (k in 1:(cells_per_well - j)) {
              # print(c(i, j, k))
              multi_mat[j, k] <- multicool::multinom(c(j, k, sample_size_well[i] - j - k),
                                                     counts = TRUE)
            }
          }
          multinomial_coeff[[i]] <- multi_mat
        }

        shared_LL <- dual_discrim_shared_likelihood(f1, f2, .15, numb_wells = well_clone, numb_cells = sample_size_well, binomials = binomial_coeff, multinomials = multinomial_coeff)
        dual_LL   <- optimize(dual_discrim_dual_likelihood, interval = c(0, .4), err = .15, numb_wells = well_clone,
                              numb_cells = sample_size_well, binomials = binomial_coeff)
        dual_LL_LL   <- dual_LL$objective
        dual_LL_freq <- dual_LL$minimum
        rec <- rbind(rec, c(clon, ind_alph1, ind_alph2, shared_LL, dual_LL_LL, dual_LL_freq))
      }
    }
  }

  rec <- as.data.frame(rec)
  rec <- dplyr::mutate(rec, diff = shared_LL - dual_LL)
  filt_rec <- dplyr::filter(rec, diff > 10)
  filt_rec <- dplyr::filter(filt_rec, shared_LL > 40 & shared_LL < 100)

  filt_rec[, 1:3, drop = FALSE]

}
