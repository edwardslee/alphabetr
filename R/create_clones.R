#' Create a synthetic set of clones with a specific underlying clonal structure
#'
#' \code{create_clones()} creates a set of (beta, alpha1, alpha2) triples that
#'    represents the indices of the chains of clones. The function will take a
#'    fixed number of unique beta chains that are in the T cell population, and
#'    then use the degree of beta and alpha sharing to determine the number of
#'    unique alpha chains in the populations. These chains will then be randomly
#'    assigned to each other, with a proportion of them being dual TCR clones
#'    (i.e. alpha1 != alpha2), forming our random list of clones with their
#'    chain indices.
#'
#' @param numb_beta The number of unique betas in the clonal population
#' @param dual The proportion of clones that are dual TCR clones, i.e. has two
#'    distinct alpha chains
#' @param beta_sharing A vector where the ith position represents the
#'    proportion of beta chains that are shared by i clones; beta chains can be
#'    shared by up to 5 clones
#' @param alpha_sharing A vector where the ith position represents the
#'    proportion of alpha chains that are shared by i clones; alpha chains can
#'    be shared by up to 6 clones
#'
#' @return A n by 3 matrix, where n is the total number of clones, column 1 is
#'    the beta index of the clone, column 2 is the alpha1 index of the clone,
#'    and column 3 is the alpha2 index of the clone. If the clone is a single
#'    TCR clone, then col 1 and col2 will have the same value. If the clone is
#'    a dual TCR clone, then the indices in col 1 and col 2 will be different.
#'
#' @examples
#'    # Creating a population containing 1000 beta chains; 30% of clones with
#'    # dual TCRs; 75% beta shared by one clone, 20% by two clones, 5% by three
#'    # clones; 80% alpha chains shared by one clone, 15% by two clones, and 5%
#'    # by three clones
#'    TCR_pairings <- create_clones(numb_beta = 1000, dual = .3,
#'                               alpha_sharing = c(0.80, 0.15, 0.05),
#'                               beta_sharing = c(0.75, 0.20, 0.05))
#'
#' @export
create_clones <- function(numb_beta = 1000, dual, alpha_sharing, beta_sharing) {
  # sharing for paper:
  #   alpha_sharing = c(.816, .085, .021, .007, .033, .005, .033)
  #   beta_sharing = c(.859, .076, .037, .019, .009)

  # Checking sharing inputs
  if (length(alpha_sharing) == 0 | any(is.na(alpha_sharing)) | any(alpha_sharing < 0)){
    stop("Please include a valid alpha chain sharing vector")
  }
  if (length(beta_sharing) == 0 | any(is.na(beta_sharing)) | any(beta_sharing < 0)){
    stop("Please include a valid beta chain sharing vector")
  }
  if(sum(alpha_sharing) != 1) {
    stop("Alpha sharing proportion do not add up to 100%")
  }
  if(sum(beta_sharing) != 1) {
    stop("Beta sharing proportion do not add up to 100%")
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

  # The try_again while loop ensures indices are set up correctly because
  # current implementation can mess up at times; will need to fix in future
  # versions
  try_again <- 1
  while (try_again == 1) {
    # shared_mat summarizes what the T cell population looks like;
    # dim(shared_mat): rows = # of beta chains, cols = 4
    #   col 1 = beta indices
    #   col 2 = # of clones that share beta_i
    #   col 3 = # of dual TCRs among clones that share beta_i
    #   col 4 = # of alpha chains associated wtih beta_i
    shared_mat <- matrix(0, ncol = 4, nrow = numb_beta)
    shared_mat[,1] <- 1:numb_beta

    # choosing the number of clones that share beta_i;
    # numb_betaX is the number of beta chains shared by X clones
    numb_beta1 <- floor(beta_sharing[1]*numb_beta)
    numb_beta2 <- floor(beta_sharing[2]*numb_beta)
    numb_beta3 <- floor(beta_sharing[3]*numb_beta)
    numb_beta4 <- floor(beta_sharing[4]*numb_beta)
    numb_beta5 <- floor(beta_sharing[5]*numb_beta)

    numb_beta1 <- numb_beta - (numb_beta2 + numb_beta3 + numb_beta4 + numb_beta5)
    if((numb_beta1 + numb_beta2 + numb_beta3 + numb_beta4 + numb_beta5) != numb_beta) {
      stop("Every beta index not being represented. Tweak beta sharing numbers")
    }
    clones_per_beta <- sample(c(rep(1, numb_beta1),
                                rep(2, numb_beta2),
                                rep(3, numb_beta3),
                                rep(4, numb_beta4),
                                rep(5, numb_beta5))
                              )
    shared_mat[, 2] <- clones_per_beta   # row i, col 2 of shared_mat represents how many clones share beta_i
    rm(clones_per_beta)                     # clear memory

    # choosing the number of dual alphas for each beta
    numb_clones <- sum(shared_mat[, 2])     # total number of clones
    numb_dual <- floor(dual * numb_clones)  # number of dual clones
    # col 3 = # of duals per beta; 1 means dual, 0 means single
    dual_assign <- sample(c(rep(1, numb_dual),
                            rep(0, numb_clones - numb_dual))
                          )
    ind_dual <- 1
    for (bet in 1:numb_beta) {                       # for each beta_i, determine how many of the clones sharing beta_i are dual clones
      shared_number <- shared_mat[bet, 2]
      shared_mat[bet, 3] <- sum(dual_assign[ind_dual:(ind_dual + shared_number - 1)])
      ind_dual <- ind_dual + shared_number
    } # end for - bet
    rm(shared_number, ind_dual, dual_assign)   # clearing memory

    # determining various important numbers; number of alpha chains per beta and other tidbits
    shared_mat[, 4] <- shared_mat[, 2] + shared_mat[, 3]  # the number of clones sharing beta_i + the number of dual TCRs = number of times that alphas associate with beta_i (not necessarily distinct)
    numb_alph_used  <- sum(shared_mat[, 4])                   # total number of times that alpha chains are used in a TCR pairing
    numb_clones     <- sum(shared_mat[, 2])                   # number of distinct clones in the population (def: two clones are distinct if at least one alpha chain is different OR beta chains are different)
    numb_dual       <- sum(shared_mat[, 3])                   # number of dual clones in the population

    # determining number of shared alphas
      shared_alph2 <- floor(alpha_sharing[2] * numb_alph_used)
      shared_alph3 <- floor(alpha_sharing[3] * numb_alph_used)
      shared_alph4 <- floor(alpha_sharing[4] * numb_alph_used)
      shared_alph5 <- floor(alpha_sharing[5] * numb_alph_used)
      shared_alph6 <- floor(alpha_sharing[6] * numb_alph_used)
      shared_alph7 <- floor(alpha_sharing[7] * numb_alph_used)

    # number of (distinct) alpha chains that are shared by more than one clone
    numb_shared_alphs <- shared_alph2 + shared_alph3 + shared_alph4 +
                         shared_alph5 + shared_alph6 + shared_alph7
    # number of distinct alphas; take the number of times alphas are used then
    # subtract the number of times shared alphas are used
    numb_unq_alph   <- numb_alph_used - (shared_alph2 + 2*shared_alph3 +
                        3*shared_alph4 + 4*shared_alph5 + 5*shared_alph6 +
                        6*shared_alph7)
    # number of alphas involved in single TCR clone
    numb_sing_alph <- numb_unq_alph - numb_shared_alphs
    if (numb_sing_alph < 0) stop("Impossible alpha sharing scenario. Reduce the degree of sharing")
    # matrix recording how many clones share alpha_i;
    # col 1 = i, col 2 = # of clones sharing alpha_i
    clones_per_alph <- matrix(c(1:numb_unq_alph, rep(1, numb_unq_alph)), ncol = 2)
    clones_per_alph[, 2] <- sample(c(rep(1,numb_sing_alph), rep(2, shared_alph2),
                                      rep(3, shared_alph3), rep(4, shared_alph4),
                                      rep(5, shared_alph5), rep(6, shared_alph6),
                                      rep(7, shared_alph7)
                                     ), replace = FALSE)

    # creating a vector of random alpha indices (each index repeated the number of times they're shared) for  random assignment of alphas to betas to form alpha/beta TCR pairs
    alpha <- vector(length = numb_alph_used)
    alph_added <- 1
    for(j in 1:nrow(clones_per_alph)) {            # go through each row of clones_per_alph;
      alph_ind  <- clones_per_alph[j, 1]            # record the alpha index
      numb.alph <- clones_per_alph[j, 2]               # record the number of times alph_ind is shared
      alpha[alph_added:(alph_added + numb.alph - 1)] <- alph_ind   # record alph_ind in numb.alph entries
      alph_added <- alph_added + numb.alph          # update next index
    } # end for - j
    alpha <- sample(alpha)        # randomly rerranged these alphas to form a random assignment (used below)

    # assiging random alphas to betas
    TCR <- matrix(ncol=3, nrow = numb_clones)               # recording the pairs, three columns, beta = col 1, alpha = col 2 and 3; if single TCR, then col 2 == col 3
    ind_beta <- 1                                           # keeps track of which number of clones we are in the TCR matrix as we loop through all the betas
    ind_alph <- 1                                           # keeps track of where we are in assigning the alphas (ie which index/indices to use next in the vector alpha)
    ind_TCR  <- 1                                           # ? i dunno yet this code is confusing
    for (betclone in 1:numb_beta) {                             # go through all of the beta chains
      numb.shared <- shared_mat[betclone, 2]                     # determine the number of clones that share this beta chain
      numb.dual   <- shared_mat[betclone, 3]                     # determine the number of dual TCR clones with this beta chain
      TCR[ind_TCR:(ind_TCR+numb.shared-1), 1] <- betclone        # record this beta chain's index in the TCR matrix
      # assigning alpha to the beta chain
      if (numb.dual > 0) {                                          # if there are dual clones associated with beta i
        last_alph_dual <- ind_alph + 2*numb.dual - 1                                      # the alphas to be used are ind_alph to last_alph_dual
        if (length(unique(alpha[ind_alph:last_alph_dual])) != length(alpha[ind_alph:last_alph_dual])) try_again <- 1
        TCR[ind_TCR:(ind_TCR+numb.dual - 1),2:3] <- alpha[ind_alph:last_alph_dual]    # assign the alphas to the dual TCR beta clones
        ind_alph <- ind_alph + 2 * numb.dual                                            # the next alpha after last_alph_dual is ind_alph + 2*numb_dual
        ind_TCR <- ind_TCR + numb.dual                                                    # the next TCR index is ind_TCR + numb_dual
        leftover_single <- numb.shared - numb.dual                                            # if numb_dual != number of clones sharing our beta clone, then we have single TCR clones
        if (leftover_single > 0) {                                                            # if there are single TCR clones
          last_ind_TCR <- ind_TCR + leftover_single - 1
          last_ind_alph <- ind_alph + leftover_single - 1
          TCR[ind_TCR:last_ind_TCR, 2:3] <- alpha[ind_alph:last_ind_alph]             # assign alpha to single TCR clone
          ind_TCR <- ind_TCR + leftover_single                                              # increase index by the number of single clones assigned
          ind_alph <- ind_alph + leftover_single                                          # increase index by the number of single clones assigned
        } # end if
      } else {
        leftover_single <- numb.shared
        last_ind_TCR <- ind_TCR + leftover_single - 1
        last_ind_alph <- ind_alph + leftover_single - 1
        TCR[ind_TCR:last_ind_TCR, 2:3] <- alpha[ind_alph:last_ind_alph]             # assign alpha to single TCR clone
        ind_TCR <- ind_TCR + leftover_single                                              # increase index by the number of single clones assigned
        ind_alph <- ind_alph + leftover_single                                          # increase index by the number of single clones assigned
      } # end if else
    } # end for - betclone
    rm(alpha, ind_beta, ind_alph, ind_TCR)     # memory clearing

    # collecting the single TCR and dual TCR clones
    TCR_single <- matrix(ncol = 2, nrow = (numb_clones - numb_dual))
    TCR_dual   <- matrix(ncol = 3, nrow = numb_dual)
    ind_single <- 1
    ind_dual <- 1
    for (clon in 1:numb_clones) {
      if (TCR[clon, 2] == TCR[clon, 3]) {
        if (ind_single > (numb_clones - numb_dual)) {
          try_again <- 1
        } else {
          TCR_single[ind_single, ] <- TCR[clon, 1:2]
          ind_single <- ind_single + 1
          try_again <- 0
        }
      } else {
        if (ind_dual > numb_dual) {
          try_again <- 1
        } else {
          TCR_dual[ind_dual, ] <- TCR[clon, ]
          ind_dual <- ind_dual + 1
          try_again <- 0
        }
      } # end if-else
    } # end for
    # randomizing the clones (i.e. randomize the rows of the TCR matrix)
  }
  TCR_dual  <- TCR_dual[complete.cases(TCR_dual), ]
  TCR_sizes <- TCR[sample(1:numb_clones, replace = FALSE),]
  colnames(TCR)       <- c("beta", "alpha1", "alpha2")
  colnames(TCR_sizes) <- c("beta", "alpha1", "alpha2")
  colnames(TCR_dual)  <- c("beta", "alpha1", "alpha2")
  list(ordered = TCR, random = TCR_sizes, dual = TCR_dual)
} # end function
