create_data_singlecells <- function(TCR, plates = 5, error = .15, skewed = 15,
                        prop_top = .5, dist = "linear") {
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

  data_alph <- vector(mode = "list", length = wells)
  data_beta <- vector(mode = "list", length = wells)

  for (well in 1:wells) {
    rand <- sample(numb_clones, size = 1, prob = dist_vector, replace = TRUE)
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
    choice_alph <- choice_alph[order(choice_alph, decreasing = FALSE)]
    choice_beta <- as.vector(TCR[rand, 1])

    data_alph[[well]] <- choice_alph[runif(length(choice_alph)) > error]
    data_beta[[well]] <- choice_beta[runif(length(choice_beta)) > error]
  }

  TCR <- TCR[order(TCR[, 1]), ]

  list(alpha = data_alph, beta = data_beta, ans = TCR)
}
