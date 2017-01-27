#' Simulate sequencing data obtained single-cell sequencing
#'
#' \code{create_data_singlecells()} simulates a single-cell sequencing
#' experiment by sampling clones from a clonal structure specified by the user
#' and using the same error models and frequency distributions used in
#' \code{\link{create_data}}. These functions are almost identical except this
#' one simulates the sampling and sequencing of single T cells.
#'
#' @param TCR The specified clonal structure, which can be created from
#'    \code{\link{create_clones}}.
#' @param plates The number of plates of data. The number of single-cells is 96
#'    times \code{plates}.
#' @param error_drop A vector of length 2 with the mean of the drop error rate
#'    and the sd of the drop error rate.
#' @param error_seq A vector of length 2 with the mean of the in-frame error
#'    rate and the sd of the in-frame error rate.
#' @param error_mode A vector of two strings determining the "mode" of the error
#'    models. The first element sets the mode of the drop errors, and the second
#'    element sets the mode of the in-frame errors. The two modes available are
#'    "constant" for a constant error rate and "lognormal" for error rates
#'    drawn from a lognormal distribution. If the mode is set to "constant" the
#'    sd specified in \code{error_drop} and/or \code{error_seq} will be ignored.
#' @param skewed Number of clones represent the top proportion of the population
#'    by frequency (which is specified by \code{prop_top}).
#' @param prop_top The proportion of the population in frequency represented by
#'    the number of clones specified by \code{skewed}.
#' @param dist The distribution of frequency of the top clones. Currently only
#'    "linear" is available.
#'
#' @return A list of length 3. The first element is a matrix representing the
#'    data of the alpha chains ($alpha), and the second element is a matrix representing
#'    the data of beta chains ($beta). The matrix represents the sequencing data by
#'    representing the wells of the data by rows and the chain indices by
#'    column. Entry [i, j] of the matrix represents if chain j is found in
#'    well i (yes == 1, no == 0). e.g. if alpha chain 25 is found in well 10,
#'    then [10, 25] of the alpha matrix will be 1.
#'
#'    The third element of the list ($drop) is a matrix that records the index
#'    of the \strong{clone} sampled in the well (col 1), records if a drop error
#'    occurred (col 2), and record if an in-frame error occurred (col 3).
#'
#' @examples
#'  # see the help for create_clones() for details of this function call
#'  clones <- create_clones(numb_beta = 1000,
#'                       dual_alpha = .3,
#'                       dual_beta  = .06,
#'                       alpha_sharing = c(0.80, 0.15, 0.05),
#'                       beta_sharing  = c(0.75, 0.20, 0.05))
#'
#'  # creating a data set with 480 single cells, lognormal error rates, 10 clones
#'  # making up the top 60% of the population in frequency, and a constant
#'  # sampling strategy of 50 cells per well for 480 wells (five 96-well plates)
#'  dat <- create_data_singlecells(clones$TCR, plate = 5,
#'                                 error_drop = c(.15, .01),
#'                                 error_seq  = c(.05, .001),
#'                                 error_mode = c("lognormal", "lognormal"),
#'                                 skewed = 10,
#'                                 prop_top = 0.6,
#'                                 dist = "linear")
#'
#'
#' @export
create_data_singlecells <- function(TCR, plates = 5, error_drop = c(.15, .01),
                        error_seq = c(.05, .01), error_mode = c("constant", "constant"),
                        skewed = 15, prop_top = .5, dist = "linear")
{
  wells <- 96 * plates    # 96 well plates

  # collecting some parameters
  numb_clones <- nrow(TCR)      # number of clones
  numb_alph <- max(TCR[, 3:4])  # number of unique alpha chains in the population
  numb_beta <- max(TCR[, 1:2])  # number of unique beta chains in the population

  # creating the skewed distribution from which the clones will be sampled; the first number.skewed clones will comprise 50% of the population, the rest of the clones the other 50%
  # The sampling is done by sampling indices from a skewed distribution, and then choosing the corresponding clone/row from the input data; the input data should already have the clones in random order
  # linear distribution: first n clones repesent 50%, p_k = p0 - k * step
  if (dist == "linear") {
    last_term <- (1 - prop_top) / (numb_clones - skewed) * 1.1
    lhs11 <- 1
    lhs12 <- skewed - 1
    lhs21 <- (1 + skewed - 1)
    lhs22 <- (skewed - 1) / 2 + (skewed - 1) ^ 2 / 2
    A <- matrix(c(lhs11, lhs21, lhs12, lhs22), nrow = 2)
    b <- matrix(c(last_term, prop_top), nrow = 2)
    solveAb <- solve(A,b)
    prob1 <- solveAb[1] + 0:(skewed - 1) * solveAb[2]
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
    # log of 0 is -Inf, so need to manually specify when err_mean is 0
    if (err_mean == 0) {
      mu <- 0
      sd <- 0
    }

    # sample drop rates from lognormal distribution; ensure none are greater than 90%
    drop_vec <- stats::rlnorm(numb_clones, meanlog = mu, sdlog = sd)
    while (any(drop_vec > .9)) {
      drop_vec <- stats::rlnorm(numb_clones, meanlog = mu, sdlog = sd)
    }
    TCR_drop <- cbind(TCR, drop_vec)
  } else {
    stop("Invalid error model specification. Choose either 'constant' or 'lognormal'")
  }

  # in-frame error rates: can be constant or occur on a distribution
  if (error_mode[2] == "constant") {
    # when constant, all chains have same rate of being sequenced as one of two distinct erroneous sequences
    err_seq_alph <- rep(error_seq[1],  numb_alph)
    err_seq_beta <- rep(error_seq[1],  numb_beta)
    err_num_alph <- rep(2,  numb_alph)
    err_num_beta <- rep(2,  numb_beta)

    # if the err rate is 0, then there are no false in frame chains
    if (error_seq[1] != 0) {
      # summing up the number of false in-frame sequences
      numb_false_alph <- sum(err_num_alph)
      numb_false_beta <- sum(err_num_beta)

      # lists record the erroneous chains associated with each chain
      false_alph <- vector(mode = "list", length = numb_false_alph)
      false_beta <- vector(mode = "list", length = numb_false_beta)
      for (i in 1:numb_false_alph) {
        false_alph[[i]] <- numb_alph + (2 * i - 1):(2 * i)
      }
      for (i in 1:numb_false_beta) {
        false_beta[[i]] <- numb_beta + (2 * i - 1):(2 * i)
      }
    } else {
      # when error rate is 0, set everything to 0
      numb_false_alph <- 0
      numb_false_beta <- 0

      # lists are empty as well
      false_alph <- vector(mode = "list", length = numb_false_alph)
      false_beta <- vector(mode = "list", length = numb_false_beta)
    }
  } else if (error_mode[2] == "lognormal") {
    # input mean and sd of the error drop rates
    err_mean <- error_seq[1]
    err_sd   <- error_seq[2]

    # if err_mean is 0, then no need to bother with this
    # if err_mean > 0, then we give error rates and # of erroneous in frame seqs
    if (err_mean != 0) {
      # need to calculate the input parameters of the association normal dist; see vign
      mu <- log(err_mean) - log((err_sd / err_mean) ^ 2 + 1) / 2
      sd <- sqrt(log((err_sd / err_mean) ^ 2 + 1))

      # sample error rates from LN dist; each chain can be erroneously
      # sequenced into 1-4 distinct different wrong chains
      err_seq_alph <- stats::rlnorm(numb_alph, meanlog = mu, sdlog = sd)
      err_seq_beta <- stats::rlnorm(numb_beta, meanlog = mu, sdlog = sd)
      # err_num_alph <- sample(1:4, size = numb_alph, replace = TRUE)
      # err_num_beta <- sample(1:4, size = numb_beta, replace = TRUE)
      err_num_alph <- rep(3, numb_alph)
      err_num_beta <- rep(3, numb_beta)

      # number of false chains
      numb_false_alph <- sum(err_num_alph)
      numb_false_beta <- sum(err_num_beta)

      # recording which errorenous chains are associated with which true chains
      false_alph <- vector(mode = "list", length = numb_alph)
      false_beta <- vector(mode = "list", length = numb_beta)

      ind <- 1
      for (i in 1:numb_alph) {
        false_alph[[i]] <- numb_alph + ind:(ind + err_num_alph[i] - 1)
        ind <- ind + err_num_alph[i]
      }
      ind <- 1
      for (i in 1:numb_beta) {
        false_beta[[i]] <- numb_beta + ind:(ind + err_num_beta[i] - 1)
        ind <- ind + err_num_beta[i]
      }
    } else { # when err_mean == 0
      # no false sequences, no false sequences associated with any chain
      numb_false_alph <- 0
      numb_false_beta <- 0
      err_seq_alph <- rep(0, numb_alph) # stats::rlnorm(numb_alph, meanlog = mu, sdlog = sd)
      err_seq_beta <- rep(0, numb_beta) # stats::rlnorm(numb_beta, meanlog = mu, sdlog = sd)
      false_alph <- vector(mode = "list", length = numb_alph)
      false_beta <- vector(mode = "list", length = numb_beta)
    }
  } else {
    stop("Invalid error model specification. Choose either 'constant' or 'lognormal'")
  }

  #-----------------------------------#
  #   Creation of the fake data set   #
  #-----------------------------------#
  # Create a matrix to represent which chains are present in which well
  # Entry (i, j) = 0 if chain j is not present in well i, (i,j) = 1 if chain j
  # is present in well i
  data_alph <- vector(mode = "list", length = wells)
  data_beta <- vector(mode = "list", length = wells)
  data_clone <- matrix(nrow = wells, ncol = 3)

  for (n_well in 1:wells) {
    rand <- sample(numb_clones, size = 1, prob = dist_vector, replace = TRUE)
    data_clone[n_well, 1] <- rand
    samp_clones <- TCR_drop[rand, , drop = FALSE]

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
    # col 2 of samp beta has the drop rate; beta_i isn't dropped when stats::runif > err
    samp_beta <- samp_beta[!is.na(samp_beta[, 1]), , drop = FALSE]
    nodrop_beta <- samp_beta
    samp_beta <- samp_beta[stats::runif(nrow(samp_beta)) > samp_beta[, 2], 1]

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
        samp_alph[ind, 2] <- samp_clones[i, 5]
        ind <- ind + 1
      }
    }
    samp_alph <- samp_alph[!is.na(samp_alph[, 1]), , drop = FALSE ]
    nodrop_alph <- samp_alph
    samp_alph <- samp_alph[stats::runif(nrow(samp_alph)) > samp_alph[, 2], 1]

    #checking if dropped
    if (length(samp_beta) < nrow(nodrop_beta) | length(samp_alph) < nrow(nodrop_alph)) {
      data_clone[n_well, 2] <- 1
    } else {
      data_clone[n_well, 2] <- 0
    }


    switch_alph <- which(stats::runif(length(samp_alph)) < err_seq_alph[samp_alph])
    for (a in switch_alph) {
      ind_alph <- samp_alph[a]
      samp_alph[a] <- sample(false_alph[[ind_alph]], size = 1)
    }

    switch_beta <- which(stats::runif(length(samp_beta)) < err_seq_beta[samp_beta])
    for (b in switch_beta) {
      ind_beta <- samp_beta[b]
      samp_beta[b] <- sample(false_beta[[ind_beta]], size = 1)
    }

    if (length(switch_beta) > 0 | length(switch_alph) > 0) {
      data_clone[n_well, 3] <- 1
    }

    data_alph[[n_well]] <- samp_alph
    data_beta[[n_well]] <- samp_beta
  } # end for - wells

  #-----------------------------------#
  # END Creation of the fake data set #
  #-----------------------------------#

  TCR <- TCR[order(TCR[, 1]), ]

  list(alpha = data_alph, beta = data_beta, drop = data_clone)
}
# end main function
