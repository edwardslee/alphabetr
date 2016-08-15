#' @export
create_data <- function(TCR, plates = 5, error = .15, skewed = 15,
                        prop_top = .5, dist = "linear",
                        numb_cells = matrix(rep(200, 12), rep(5*plates, 12), ncol = 2)
                        ){
  wells <- 96 * plates    # 96 well plates

  # collecting some parameters
  numb_clones <- nrow(TCR)     # number of clones
  numb_alph <- max(TCR[, 2:3]) # number of unique alpha chains in the population
  numb_beta <- max(TCR[, 1])   # number of unique beta chains in the population

  # creating the skewed distribution from which the clones will be sampled; the first number.skewed clones will comprise 50% of the population, the rest of the clones the other 50%
  # The sampling is done by sampling indices from a skewed distribution, and then choosing the corresponding clone/row from the input data; the input data should already have the clones in random order
  # linear distribution: first n clones repesent 50%, p_k = p0 - k * step
  # geometric distribution: first n clones represent 50%, p_k = p0 * .7^k

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

  # sample_sizes <- rep(0, length = nrow(TCR))
  # sample_in_wells <- list(length = wells)

  #-----------------------------------#
  #   Creation of the fake data set   #
  #-----------------------------------#
  # Create a matrix to represent which chains are present in which well
  # Entry (i, j) = 0 if chain j is not present in well i, (i,j) = 1 if chain j
  # is present in well i
  data_alph <- matrix(0, nrow = wells, ncol = numb_alph)
  data_beta <- matrix(0, nrow = wells, ncol = numb_beta)

  # dropdrop1 <- 1
  # dropdrop2 <- 1

  #  matrix(rep(200, 12),
  #         rep(5*plates, 12), ncol = 2)
  distinct_ss <- numb_cells[, 1]  # vector of the distinct sample sizes in well
  ind_well <- 1                   # index of well we will be sampling cells into
  for (sampl in seq_along(distinct_ss)) {
    sample_size <- distinct_ss[sampl]   # number of cells in the wells
    numb_wells  <- numb_cells[sampl, 2] # number of wells with sample_size cells
    for (n_well in 1:numb_wells) {
      rand <- sample(numb_clones, size = sample_size, prob = dist_vector,
                     replace = TRUE)

      # Collecting the alpha chains that are chosen
      choice_alph    <- as.vector(TCR[rand,2:3])
      half_len <- length(choice_alph)/2
      remove_ind <- vector()
      for (x in half_len:1) {
        i1 <- x
        i2 <- x + half_len
        if (choice_alph[i1] == choice_alph[i2]) remove_ind <- c(remove_ind, i2)
      }
      if (length(remove_ind > 0)) choice_alph <- choice_alph[-remove_ind]
      choice_beta <- as.vector(TCR[rand, 1])

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
