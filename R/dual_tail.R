#' Discriminate between beta sharing clones and dual TCR clones.
#'
#' \code{dual_tail()} distinguishes between clones that share a common beta
#'    chain and dual TCR clones with two productive alpha chains. The procedure
#'    tests the null hypothesis that two candidate alpha, beta pairs with the
#'    same beta represent two separate clones by using the frequency estimates
#'    to calculate the number of wells that both clones are expected to be in.
#'    This is compared to the actual number of wells that both clones appear in,
#'    and if the actual number is greater than the expected number, than the
#'    pairs are chosen to represent a dual TCR clone.
#' @param data_alpha Matrix recording which alpha chains appear in each well of the
#'    data. See  .
#' @param data_beta Matrix recording which beta chains appear in the each well of the
#'    data. See .
#' @param freq_results Output of \code{freq_estimate()}
#' @param numb_wells Vector containing the number of cells sampled in the wells
#'    of each column of the plates.
#'
#' @return A n x 3 matrix where n is the number of candidate clones, column 1
#'    is the beta index of the clone, and column 2-3 are the alpha indices of
#'    the clone
#' @examples asdfas
#' @export
dual_tail <- function(alpha, beta, freq_results, cells, population) {
  freq_results <- freq_results[order(freq_results[, "MLE"], decreasing = TRUE), ]
  freq_results <- freq_results[!is.na(freq_results[, "MLE"]), ]

  dual_procedure <- function(data_alph, data_beta, freq, cells) {
    # preallocate matrix to record the expected
    # col 1: beta index,      col 2: alpha1,             col 3: alpha2
    # col 4: expt # of wells, col 5: actual # of wells
    rec <- matrix(nrow = 0, ncol = 6)
    colnames(rec) <- c("beta1", "beta2",
                       "alpha1", "alpha2",
                       "expt_wells", "act_wells")

    # checking the clones with each beta chain
    #   -calculate the expected number of plates for every pair of alpha associated
    #       associated with beta under assumption that they are independent clones
    #   -if beta, alpha_k1, alpha_k2 show up more often than expected, then they
    #       must be dual
    for (clon in 1:max_beta) {
      x <- freq[freq[, "beta1"] == clon, , drop = FALSE]  # find clones with the beta index
      numb_cand <- nrow(x)                            # find number of alphas associated with beta
      if (numb_cand > 1) {                            # if more than 1 alpha
        combos <- combn(numb_cand, 2)                   # find all combos of pairs
        for (ind in 1:ncol(combos)) {                   # check each combo
          alpha1 <- x[combos[1, ind], "alpha1"]
          alpha2 <- x[combos[2, ind], "alpha1"]
          f1 <- x[combos[1, ind], "MLE"]                # freq est of beta alpha1
          f2 <- x[combos[2, ind], "MLE"]                # freq est of beta alpha2
          exp_plates <- 0                               # expected number of plates
          for (samp in seq_along(sample_size_well)) {                  # find the number of expected well apperances in each column
            number_cells <- sample_size_well[samp]
            exp_plates <- numb_sample[samp] *
              (1 - (1 - f1)^number_cells - (1 - f2)^number_cells +
                 (1 - (f1 + f2))^number_cells) + exp_plates
          }
          dat_plates <- sum(data_beta[, clon] == 1 &
                              data_alph[, alpha1] == 1 &
                              data_alph[, alpha2] == 1)
          rec <- rbind(rec, c(clon, clon, alpha1, alpha2, exp_plates, dat_plates))
        }
      }
    } # end for going through all beta chains

    rec_tail <- rec
    if (nrow(rec) > 2) {
      # make decisions on if two clones are independent or actually a dual TCR by
      # looking at the ratio of actual # of wells vs expected # of wells
      # -do this by clustering the ratios into "high" and "low" ratios using k-means
      # -the "high" cluster are more likely to be duals
      well_ratio <- rec[, "act_wells"]/rec[, "expt_wells"] # ratio of actual to expected
      if (length(unique(well_ratio)) > 1) {
        km_out <- kmeans(well_ratio, 2, nstart = 200)      # kmeans of ratios in 2 clusters
        # figure out which cluster represents "high"
        clus_ind1 <- which(km_out$cluster == 1)            # one cluster
        clus_ind2 <- which(km_out$cluster == 2)            # the other cluster
        # determine which cluster has the higher ratios, making it "high" cluster
        if (mean(well_ratio[clus_ind1]) > mean(well_ratio[clus_ind2])) {
          clus_ind <- clus_ind1
        } else {
          clus_ind <- clus_ind2
        } # end if - else
      } # end if
    } # end if - nrow

    list(results = rec, cluster = clus_ind)
  } # end function - dual_procedure

  # parameters
  max_beta <- ncol(beta)              # determine maximum beta index
  sample_size_well <- cells[, 1]      # number of cells per well
  numb_sample <- cells[, 2]           # number of wells w/ sample size

  # perform the dual discrimination for the tail
  tail_dual <- dual_procedure(alpha, beta, freq_results, cells)
  tail_rec <- tail_dual$results
  tail_rec <- tail_rec[tail_dual$cluster, 1:4, drop = FALSE]

  colnames(tail_rec) <- c("beta1", "beta2", "alpha1", "alpha2")
  as.data.frame(tail_rec)
}
