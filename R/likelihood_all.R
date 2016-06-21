
# freq_estimate_new(dat_alph_all, dat_beta_all, aa, error = error_rate, cells = numb_cells)
likelihood_all <- function(f, err, alpha, beta, cells, ind_beta, ind_alph) {
  obj <- 0
  prob_a  <- 0
  prob_b  <- 0
  prob_ab <- 0

  # Determing which wells that contain the clone a clone is counted to be in a
  # well if their component chains are found in the same well
  sample_size_well <- cells[, 1]     # number of cells per well
  numb_sample <- cells[, 2]          # number of wells w/ sample size


  numb_distinct <- length(sample_size_well)
  well_clone <- rep(0, numb_distinct)

  well_clone <- matrix(ncol = 4, nrow = numb_distinct)
  colnames(well_clone) <- c("a", "b", "ab", "none")


  # browser()
  for (size in 1:numb_distinct) {
    ind1 <- cumsum(numb_sample[1:size])[size] - numb_sample[size] + 1
    ind2 <- cumsum(numb_sample[1:size])[size]
    well_clone[size, "a"] <- sum(beta[ind1:ind2, ind_beta] == 0 &
                                   alpha[ind1:ind2, ind_alph] == 1)
    well_clone[size, "b"] <- sum(beta[ind1:ind2, ind_beta] == 1 &
                                   alpha[ind1:ind2, ind_alph] == 0)
    well_clone[size, "ab"] <- sum(beta[ind1:ind2, ind_beta] == 1 &
                                    alpha[ind1:ind2, ind_alph] == 1)
    well_clone[size, "none"] <- numb_sample[size] - sum(well_clone[size, 1:3])
  }


  browser()
  for(samp in seq_len(nrow(well_clone))) {
    numb_a  <- well_clone[samp, "a"]
    numb_b  <- well_clone[samp, "b"]
    numb_ab <- well_clone[samp, "ab"]
    numb_0  <- well_clone[samp, "none"]
    number_cells_well <- numb_cells[samp]
    for (i in 1:number_cells_well) {
      prob_a  <- choose(number_cells_well, i) * f^i * (1-f)^(number_cells_well - i) * (1 - err^i) * err^i
      prob_ab <- choose(number_cells_well, i) * f^i * (1-f)^(number_cells_well - i) * (1 - err^i)^2
    }
    prob_b <- prob_a
    prob_0 <- 1 - prob_a - prob_b - prob_ab

    p_a <- 0
    p_b <- 0
    p_ab <- 0
    p_0 <- 0

    if (numb_a != 0) {
      p_a <- numb_a * log(prob_a)
    }
    if (numb_b != 0) {
      p_b <- numb_b * log(prob_b)
    }
    if (numb_ab != 0) {
      p_ab <- numb_ab * log(prob_ab)
    }
    if (numb_0 != 0) {
      p_0 <- numb_0 * log(prob_0)
    }

    obj <- obj - p_a - p_b - p_ab - p_0
  }
  obj
}
