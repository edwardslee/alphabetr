#' Create a synthetic set of clones with a specific underlying clonal structure
#'
#' \code{create_clones()} creates a set of (beta1, beta2, alpha1, alpha2)
#'    quadruples that represent the indices of the chains of clones. The function
#'    will take a fixed number of unique beta chains that are in the T cell
#'    population, and then use the degree of beta and alpha sharing to determine
#'    the number of unique alpha chains in the populations. These chains will
#'    then be randomly assigned to each other, with a proportion of them being
#'    dual TCR clones (i.e. alpha1 != alpha2 and/or beta1 != beta2), forming
#'    our random list of clones with their chain indices.
#'
#' @param numb_beta The number of unique betas in the clonal population
#' @param dual_beta The proportion of clone that are dual TCRbeta clones, i.e.
#'    has two distinct beta chains
#' @param dual_alpha The proportion of clones that are dual TCRalpha clones,
#'    i.e. has two distinct alpha chains
#' @param beta_sharing A vector where the ith position represents the
#'    proportion of beta chains that are shared by i clones; beta chains can be
#'    shared by up to 5 clones
#' @param alpha_sharing A vector where the ith position represents the
#'    proportion of alpha chains that are shared by i clones; alpha chains can
#'    be shared by up to 7 clones
#'
#' @return A list of four different matrices. Each matrix is has dimensions
#'    n x 4, where n is the total number of clones and each row represents the
#'    chains of a clone. Column 1 and column 2 are the beta index/indices of the
#'    beta chain(s) used by the clone. Column 3 and 4 are the alpha index/indices
#'    of the alpha chain(s) used by the clone. If a clone has a single beta
#'    chain, then col 1 and col 2 will be equal. If a clone has a single alpha
#'    chain, then col 3 and col 4 will be equal.
#'
#' @examples
#'  # Creating a population containing 1000 beta chains; 10% of clones with
#'  # dual-beta TCRs and 30% of clones with dual TCRs; 75% beta shared by one
#'  # clone, 20% by two clones, 5% by three clones; 80% alpha chains shared by
#'  # one clone, 15% by two clones, and 5% by three clones
#'
#'  clones <- create_clones(numb_beta = 1000,
#'                          dual = .3,
#'                          alpha_sharing = c(0.80, 0.15, 0.05),
#'                          beta_sharing  = c(0.75, 0.20, 0.05))
#'
#' @export
create_clones <- function(numb_beta, dual_beta, dual_alpha, alpha_sharing, beta_sharing) {
  # sharing for paper:
  # alpha_sharing = c(.816, .085, .021, .007, .033, .005, .033)
  # beta_sharing = c(.859, .076, .037, .019, .009)

  # Checking sharing inputs
  if (length(alpha_sharing) == 0 | any(is.na(alpha_sharing)) | any(alpha_sharing < 0)){
    stop("Please include a valid alpha chain sharing vector")
  }
  if (length(beta_sharing) == 0 | any(is.na(beta_sharing)) | any(beta_sharing < 0)){
    stop("Please include a valid beta chain sharing vector")
  }
  if (sum(alpha_sharing) != 1) {
    stop("Alpha sharing proportion do not add up to 100%")
  }
  if (sum(beta_sharing) != 1) {
    stop("Beta sharing proportion do not add up to 100%")
  }
  if (length(numb_beta) != 1) {
    stop("Length of the numb_beta argument must be 1")
  }

  if (length(alpha_sharing) < 7) {
    # add enough zeros to alpha_sharing so that it is 6 elements long
    len_zero <- 7 - length(alpha_sharing)
    alpha_sharing <- c(alpha_sharing, rep(0, len_zero))
  }
  if (length(beta_sharing) < 5) {
    # add enough zeros to beta_sharing so that it is 6 elements long
    len_zero <- 5 - length(beta_sharing)
    beta_sharing <- c(beta_sharing, rep(0, len_zero))
  }



  # n_beta_i is the number of beta chains that will be shared by i clones
  n_beta2 <- ceiling(beta_sharing[2] * numb_beta)
  n_beta3 <- ceiling(beta_sharing[3] * numb_beta)
  n_beta4 <- ceiling(beta_sharing[4] * numb_beta)
  n_beta5 <- ceiling(beta_sharing[5] * numb_beta)
  n_beta1 <- numb_beta - n_beta2 - n_beta3 - n_beta4 - n_beta5

  # this vector will randomly allocate the # of clones sharing a beta chain
  clones_per_beta <- sample(c(rep(1, n_beta1),
                              rep(2, n_beta2),
                              rep(3, n_beta3),
                              rep(4, n_beta4),
                              rep(5, n_beta5)
  ))

  ind_beta <- sample(1:numb_beta)

  beta_mat <- matrix(c(ind_beta, clones_per_beta), ncol = 2)
  beta_vec <- vector(length = sum(clones_per_beta))
  ind <- 1
  for (i in seq_len(nrow(beta_mat))) {
    numb_clones <- clones_per_beta[i]
    beta_vec[ind:(ind + numb_clones - 1)] <- rep(ind_beta[i], numb_clones)
    ind <- ind + numb_clones
  }
  beta_vec <- sample(beta_vec)
  numb_dual <- ceiling(dual_beta * sum(clones_per_beta) / (1 + dual_beta))

  assign_dual <- sample(seq(from = 1, to = sum(clones_per_beta) - 1, by = 2), size = numb_dual)


  # check to see if the duals we assigned have different indices
  check_dual <- 1
  while(check_dual == 1) {
    check_dual <- 0
    for (i in assign_dual) {
      if (beta_vec[i] == beta_vec[i + 1]) {
        check_dual <- 1
        assign_dual <- sample(seq(from = 1, to = sum(clones_per_beta) - 1, by = 2), size = numb_dual)
      }
    }
  }


  clones <- matrix(nrow = sum(clones_per_beta) - numb_dual, ncol = 4)
  colnames(clones) <- c("beta1", "beta2", "alpha1", "alpha2")

  ind_beta_vec <- 1
  ind_clones <- 1
  while(ind_clones <= nrow(clones)) {
    if (any(assign_dual == ind_beta_vec)) {
      clones[ind_clones, 1] <- beta_vec[ind_beta_vec]
      clones[ind_clones, 2] <- beta_vec[ind_beta_vec + 1]
      ind_beta_vec <- ind_beta_vec + 2
      ind_clones <- ind_clones + 1
    } else {
      clones[ind_clones, 1] <- beta_vec[ind_beta_vec]
      clones[ind_clones, 2] <- beta_vec[ind_beta_vec]
      ind_beta_vec <- ind_beta_vec + 1
      ind_clones <- ind_clones + 1
    }
  }


  numb_alph_used <- nrow(clones) + ceiling(dual_alpha * nrow(clones))

  numb_alph <- floor(numb_alph_used / sum(1:length(alpha_sharing) * alpha_sharing))
  n_alph2 <- floor(alpha_sharing[2] * numb_alph)
  n_alph3 <- floor(alpha_sharing[3] * numb_alph)
  n_alph4 <- floor(alpha_sharing[4] * numb_alph)
  n_alph5 <- floor(alpha_sharing[5] * numb_alph)
  n_alph6 <- floor(alpha_sharing[6] * numb_alph)
  n_alph7 <- floor(alpha_sharing[7] * numb_alph)

  n_alph1   <- numb_alph_used - (2 * n_alph2 + 3 * n_alph3 +
                                   4 * n_alph4 + 5 * n_alph5 +
                                   5 * n_alph6 + 7 * n_alph7)

  numb_alph <- n_alph1 + n_alph2 + n_alph3 + n_alph4 +
    n_alph5 + n_alph6 + n_alph7

  # this vector will randomly allocate the # of clones sharing a beta chain
  clones_per_alph <- sample(c(rep(1, n_alph1),
                              rep(2, n_alph2),
                              rep(3, n_alph3),
                              rep(4, n_alph4),
                              rep(5, n_alph5),
                              rep(6, n_alph6),
                              rep(7, n_alph7)
  ))

  ind_alph <- sample(1:numb_alph)

  alph_mat <- matrix(c(ind_alph, clones_per_alph), ncol = 2)
  alph_vec <- vector(length = sum(clones_per_alph))
  ind <- 1
  for (i in seq_len(nrow(alph_mat))) {
    numb_clones <- clones_per_alph[i]
    alph_vec[ind:(ind + numb_clones - 1)] <- rep(ind_alph[i], numb_clones)
    ind <- ind + numb_clones
  }
  alph_vec <- sample(alph_vec)

  numb_dual <- ceiling(dual_alpha * nrow(clones))
  assign_dual <- sample(seq(from = 1, to = sum(clones_per_alph) - 1, by = 2), size = numb_dual)

  # check to see if the duals we assigned have different indices
  check_dual <- 1
  while(check_dual == 1) {
    check_dual <- 0
    for (i in assign_dual) {
      if (alph_vec[i] == alph_vec[i + 1]) {
        check_dual <- 1
        assign_dual <- sample(seq(from = 1, to = sum(clones_per_alph) - 1, by = 2), size = numb_dual)
      }
    }
  }

  ind_alph_vec <- 1
  ind_clones <- 1
  while(ind_clones <= nrow(clones)) {
    if (any(assign_dual == ind_alph_vec)) {
      clones[ind_clones, 3] <- alph_vec[ind_alph_vec]
      clones[ind_clones, 4] <- alph_vec[ind_alph_vec + 1]
      ind_alph_vec <- ind_alph_vec + 2
      ind_clones <- ind_clones + 1
    } else {
      clones[ind_clones, 3] <- alph_vec[ind_alph_vec]
      clones[ind_clones, 4] <- alph_vec[ind_alph_vec]
      ind_alph_vec <- ind_alph_vec + 1
      ind_clones <- ind_clones + 1
    }
  }

  ordered   <- clones[order(clones[, 1], decreasing = FALSE), ]
  dual_alph <- clones[clones[, 3] != clones[, 4], ]
  dual_beta <- clones[clones[, 1] != clones[, 2], ]
  colnames(ordered)   <- c("beta1", "beta2", "alpha1", "alpha2")
  colnames(dual_alph) <- c("beta1", "beta2", "alpha1", "alpha2")
  colnames(dual_beta) <- c("beta1", "beta2", "alpha1", "alpha2")

  list(TCR = clones,
       ordered = ordered,
       dual_alpha = dual_alph,
       dual_beta = dual_beta)
} # end function
