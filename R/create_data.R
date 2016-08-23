#' @export
create_data <- function(TCR, plates = 5, error_drop = c(.15, .01),
                        error_seq = c(.05, .01), error_mode = c("constant", "constant"),
                        skewed = 15, prop_top = .5, dist = "linear",
                        numb_cells = matrix(rep(200, 12), rep(5*plates, 12), ncol = 2)
                        ){
  wells <- 96 * plates    # 96 well plates

  # collecting some parameters
  numb_clones <- nrow(TCR)      # number of clones
  numb_alph <- max(TCR[, 3:4])  # number of unique alpha chains in the population
  numb_beta <- max(TCR[, 1:2])  # number of unique beta chains in the population

  # creating the skewed distribution from which the clones will be sampled; the first number.skewed clones will comprise 50% of the population, the rest of the clones the other 50%
  # The sampling is done by sampling indices from a skewed distribution, and then choosing the corresponding clone/row from the input data; the input data should already have the clones in random order
  # linear distribution: first n clones repesent 50%, p_k = p0 - k * step
  if (dist == "linear") {
    last_term <- (1 - prop_top)/(numb_clones - skewed) * 1.1
    lhs11 <- 1
    lhs12 <- skewed - 1
    lhs21 <- (1 + skewed - 1)
    lhs22 <- (skewed - 1)/2 + (skewed - 1)^2/2
    A <- matrix(c(lhs11, lhs21, lhs12, lhs22), nrow = 2)
    b <- matrix(c(last_term, prop_top), nrow = 2)
    solveAb <- solve(A,b)
    prob1 <- solveAb[1] + 0:(skewed-1) * solveAb[2]
    prob2 <- rep((1 - prop_top)/(numb_clones - skewed), numb_clones - skewed)
    dist_vector <- c(prob1, prob2)
  }


  #------------- error model ----------------#
  # Create columns to attach to TCR (the clones matrix) for the errors in the exp
  # col 5 is be the drop rate (i.e. the rate at which the chain won't show up at all)
  # col 6 is be the substitution sequencing rates (i.e. rate of false in frame seq)
  # col 7 is the number of different possible substitute sequences

  # drop rates: can be a constant drop rate or drop rates distributed on LN dist
  if (error_mode[1] == "constant") {
    err_drop <- error_drop[1]
    TCR_drop <- cbind(TCR, err_drop)
  } else if (error_mode[1] == "lognormal") {
    if (length(error_drop) != 2) stop("The error_drop argument must be length 2.")

    # input mean and sd of the error drop rates
    err_mean <- error_drop[1]
    err_sd   <- error_drop[2]

    # need to calculate the input parameters of the association normal dist; see vign
    mu <- log(err_mean) - log((err_sd / err_mean) ^ 2 + 1) / 2
    sd <- sqrt(log((err_sd / err_mean) ^ 2 + 1))

    drop_vec <- rlnorm(numb_clones, meanlog = mu, sdlog = sd)
    while (any(drop_vec > .9)) {
      drop_vec <- rlnorm(numb_clones, meanlog = mu, sdlog = sd)
    }
    TCR_drop <- cbind(TCR, drop_vec)
  } else {
    stop("Invalid error model specification. Choose either 'constant' or 'lognormal'")
  }


  err_seq <- error[2]

  err_seq_alph <- vector(length = numb_alph)
  err_seq_beta <- vector(length = numb_beta)
  err_num_alph <- vector(length = numb_alph)
  err_num_beta <- vector(length = numb_beta)

  if (error_mode[2] == "constant") {
    err_seq_alph <- rep(error_seq[1],  numb_alph)
    err_seq_beta <- rep(error_seq[1],  numb_beta)
    err_num_alph <- rep(2,  numb_alph)
    err_num_beta <- rep(2,  numb_beta)
  } else if (error_mode[2] == "lognormal") {
    # input mean and sd of the error drop rates
    err_mean <- error_seq[1]
    err_sd   <- error_seq[2]

    # need to calculate the input parameters of the association normal dist; see vign
    mu <- log(err_mean) - log((err_sd / err_mean) ^ 2 + 1) / 2
    sd <- sqrt(log((err_sd / err_mean) ^ 2 + 1))

    err_seq_alph <- rlnorm(numb_alph, meanlog = mu, sdlog = sd)
    err_seq_beta <- rlnorm(numb_beta, meanlog = mu, sdlog = sd)

  } else {
    stop("Invalid error model specification. Choose either 'constant' or 'lognormal'")
  }

  #-----------------------------------#
  #   Creation of the fake data set   #
  #-----------------------------------#
  # Create a matrix to represent which chains are present in which well
  # Entry (i, j) = 0 if chain j is not present in well i, (i,j) = 1 if chain j
  # is present in well i
  data_alph <- matrix(0, nrow = wells, ncol = numb_alph)
  data_beta <- matrix(0, nrow = wells, ncol = numb_beta)

  distinct_ss <- numb_cells[, 1]  # vector of the distinct sample sizes in well
  ind_well <- 1                   # index of well we will be sampling cells into
  for (sampl in seq_along(distinct_ss)) {
    sample_size <- distinct_ss[sampl]   # number of cells in the wells
    numb_wells  <- numb_cells[sampl, 2] # number of wells with sample_size cells
    for (n_well in 1:numb_wells) {
      rand <- sample(numb_clones, size = sample_size, prob = dist_vector,
                     replace = TRUE)
      samp_clones <- TCR_drop[rand, ]

      samp_beta <- matrix(nrow = 2 * nrow(samp_clones), ncol = 2)
      dual_beta <- which(samp_clones[, 1] != samp_clones[, 2])
      ind <- 1
      for (i in seq_len(nrow(samp_clones))) {
          if (i %in% dual_beta) {
            samp_beta[ind:(ind + 1), 1] <- samp_clones[i, 1:2]
            samp_beta[ind:(ind + 1), 2] <- samp_clones[i, 5]
            ind <- ind + 2
          } else {
            samp_beta[ind, 1] <- samp_clones[i, 1]
            samp_beta[ind, 2] <- samp_clones[i, 5]
            ind <- ind + 1
          }
      }
      # col 2 of samp beta has the drop rate; beta_i isn't dropped when runif > err
      samp_beta <- samp_beta[runif(nrow(samp_beta)) > samp_beta[, 2], 1]

      samp_alph <- matrix(nrow = 2 * nrow(samp_clones), ncol = 2)
      dual_alph <- which(samp_clones[, 3] != samp_clones[, 4])
      ind <- 1
      for (i in seq_len(nrow(samp_clones))) {
        if (i %in% dual_alph) {
          samp_alph[ind:(ind + 1), 1] <- samp_clones[i, 3:4]
          samp_alph[ind:(ind + 1), 2] <- samp_clones[i, 5]
          ind <- ind + 2
        } else {
          samp_alph[ind, 1] <- samp_clones[i, 3]
          samp_alph[ind, 2] <- samp_clones[i, 4]
          ind <- ind + 1
        }
      }
      # col 2 of samp beta has the drop rate; beta_i isn't dropped when runif > err
      samp_alph <- samp_alph[runif(nrow(samp_alph)) > samp_alph[, 2], 1]

      alph_drop <- rep(samp_clones[, 4], each = 2)

      for (x in choice_alph)  if(runif(1) > error) data_alph[ind_well, x] <- 1
      for (x in choice_beta)  if(runif(1) > error) data_beta[ind_well, x] <- 1

      # Internal QA controls
      # sample_sizes[rand] <- sample_sizes[rand] + 1
      # sample_in_wells[[ind_well]] <- rand

      ind_well <- ind_well + 1
    } # end for - wells
  } # end for - sampl

  #-----------------------------------#
  # END Creation of the fake data set #
  #-----------------------------------#

  TCR <- TCR[order(TCR[, 1]), ]

  return(list(alpha = data_alph, beta = data_beta, ans = TCR))#, tracker = sample_sizes, wells_tracker = sample_in_wells))#, ans1 = TCR.single, ans2 = TCR.dual))

}
# end main function
