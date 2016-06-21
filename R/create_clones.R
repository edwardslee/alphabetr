#' Create a synthetic list of clones
#'
#' \code{create_clones()} creates a list of (beta, alpha1, alpha2) sets that
#'    represents the indices of the chains of clones. The function will take a
#'    fixed number of unique beta chains that are in the T cell population, and
#'    then use the degree of beta and alpha sharing to determine the number of
#'    unique alpha chains in the populations. These chains will then be randomly
#'    assigned to each other, with a proportion of them being dual TCR clones
#'    (i.e. alpha1 != alpha2), forming our random list of clones with their
#'    chain indices.
#'
#' @param numb_beta The number of unique betas in the population
#' @param dual The proportion of clones that are dual TCR clones, i.e. has two
#'    distinct alpha chains
#'
#' @export
create_clones <- function(numb_beta = 1000, dual = .1) {
  # createClones is used to create the alpha/beta TCR pairs and incorporates dual TCRs (ie clones with one beta and two alpha chains),
  # clones sharing the same beta chain, and clones sharing the same alpha chains.


  try_again <- 1
  while (try_again == 1) {
    # shared_mat summarizes what the T cell population looks like; dim(shared_mat): rows = # of beta chains, cols = 4
    #   col 1 = beta indices
    #   col 2 = # of clones that share beta_i
    #   col 3 = # of dual TCRs among clones that share beta_i
    #   col 4 = # of alpha chains associated wtih beta_i
    shared_mat <- matrix(0, ncol = 4, nrow = numb_beta) # preallocating matrix
    shared_mat[,1] <- 1:numb_beta                       # assigning the indices of beta chains to the first column

    # choosing the number of clones that share beta_i; numb_betaX is the number of beta chains shared by X clones
    numb_beta1 <- floor(0.859*numb_beta)       # 40% of the distinct beta chains will be found in 1 clone only
    numb_beta2 <- floor(0.076*numb_beta)       # 20% of the distinct beta chains will be shared by 2 clones
    numb_beta3 <- floor(0.037*numb_beta)       # 20% of the distinct beta chains will be shared by 3 clones
    numb_beta4 <- floor(0.019*numb_beta)       # 10% of the distinct beta chains will be shared by 4 clones
    numb_beta5 <- floor(0.009*numb_beta)       # 10% of the distinct beta chains will be shared by 5 clones

    # numb_beta1 <- floor(.40*numb_beta)       # 40% of the distinct beta chains will be found in 1 clone only
    # numb_beta2 <- floor(.20*numb_beta)       # 20% of the distinct beta chains will be shared by 2 clones
    # numb_beta3 <- floor(.20*numb_beta)       # 20% of the distinct beta chains will be shared by 3 clones
    # numb_beta4 <- floor(.10*numb_beta)       # 10% of the distinct beta chains will be shared by 4 clones
    # numb_beta5 <- floor(.10*numb_beta)       # 10% of the distinct beta chains will be shared by 5 clones

    # numb_beta1 <- floor(.70*numb_beta)       # 40% of the distinct beta chains will be found in 1 clone only
    # numb_beta2 <- floor(.10*numb_beta)       # 20% of the distinct beta chains will be shared by 2 clones
    # numb_beta3 <- floor(.10*numb_beta)       # 20% of the distinct beta chains will be shared by 3 clones
    # numb_beta4 <- floor(.05*numb_beta)       # 10% of the distinct beta chains will be shared by 4 clones
    # numb_beta5 <- floor(.05*numb_beta)       # 10% of the distinct beta chains will be shared by 5 clones
    numb_beta1 <- numb_beta - (numb_beta2 + numb_beta3 + numb_beta4 + numb_beta5)
    if((numb_beta1 + numb_beta2 + numb_beta3 + numb_beta4 + numb_beta5) != numb_beta) return(print("ERROR in beta chains")) # check that every beta will get an assignment of the number of clones that share it
    cutoff1 <- numb_beta1
    cutoff2 <- numb_beta1 + numb_beta2
    cutoff3 <- numb_beta1 + numb_beta2 + numb_beta3
    cutoff4 <- numb_beta1 + numb_beta2 + numb_beta3 + numb_beta4
    cutoff5 <- numb_beta
    clones_per_beta <- sample(c(rep(1, numb_beta1),   rep(2, numb_beta2), rep(3, numb_beta3), rep(4, numb_beta4), rep(5, numb_beta5)))
    # clones_per_beta <- sample(c(rep(1, numb_beta1),   rep(1, numb_beta2), rep(1, numb_beta3), rep(1, numb_beta4), rep(1, numb_beta5)))
    shared_mat[, 2] <- clones_per_beta   # row i, col 2 of shared_mat represents how many clones share beta_i
    rm(clones_per_beta)                     # clear memory

    ## stochastic method of assigning
    #number.shared <- runif(numb_beta)                      # choose a random number to determine how many clones will share this beta
    #shared_mat[,2] <- ifelse(number.shared >= 0   & number.shared <  .4, 1, shared_mat[,2])   # 40% chance of just 1 clone having this beta chain
    #shared_mat[,2] <- ifelse(number.shared >= .4  & number.shared <  .6, 2, shared_mat[,2])   # 20% chance of 2 clones having this beta chain
    #shared_mat[,2] <- ifelse(number.shared >= .6  & number.shared <  .8, 3, shared_mat[,2])   # 20% chance of 3 clones having this beta chain
    #shared_mat[,2] <- ifelse(number.shared >= .8  & number.shared < .95, 4, shared_mat[,2])   # 15% chance of 4 clones having this beta chain
    #shared_mat[,2] <- ifelse(number.shared >= .95  & number.shared <=  1, 5, shared_mat[,2])   # 5% chance of 5 clones having this beta chain
    #rm(number.shared)                                                                               # clearning number.shared; memory management

    # choosing the number of dual alphas for each beta
    numb_clones <- sum(shared_mat[, 2])                                                # calculate the number of clones by adding up the number of clones that share each distinct beta
    numb_dual <- floor(dual * numb_clones)                                              # calculate the number of dual clones in the population
    dual_assign <- sample(c(rep(1, numb_dual), rep(0, numb_clones - numb_dual)))     # col 3 = # of duals per beta; 1 means dual, 0 means single
    ind_dual <- 1
    for (bet in 1:numb_beta) {                       # for each beta_i, determine how many of the clones sharing beta_i are dual clones
      shared.number <- shared_mat[bet, 2]
      shared_mat[bet, 3] <- sum(dual_assign[ind_dual:(ind_dual + shared.number - 1)])
      ind_dual <- ind_dual + shared.number
    } # end for - bet
    rm(shared.number, ind_dual, dual_assign)   # clearing memory

    # determining various important numbers; number of alpha chains per beta and other tidbits
    shared_mat[, 4] <- shared_mat[, 2] + shared_mat[, 3]  # the number of clones sharing beta_i + the number of dual TCRs = number of times that alphas associate with beta_i (not necessarily distinct)
    numb_alph_used <- sum(shared_mat[, 4])                   # total number of times that alpha chains are used in a TCR pairing
    numb_clones     <- sum(shared_mat[, 2])                   # number of distinct clones in the population (def: two clones are distinct if at least one alpha chain is different OR beta chains are different)
    numb_dual       <- sum(shared_mat[, 3])                   # number of dual clones in the population

    # determining number of shared alphas
    # shared_alph.2 <- floor(.09*numb_alph_used)                 # 15% of alpha chains are shared by 2 different clones
    # shared_alph.3 <- floor(.03*numb_alph_used)                 # 10% of alpha chains are shared by 3 different clones
    # shared_alph.4 <- floor(.03*numb_alph_used)                 #  5% of alpha chains are shared by 4 different clones
    # shared_alph.5 <- floor(.03*numb_alph_used)                 #  5% of alpha chains are shared by 4 different clones
    # shared_alph.6 <- floor(.03*numb_alph_used)                 #  5% of alpha chains are shared by 4 different clones
    # shared_alph.7 <- floor(.03*numb_alph_used)                 #  5% of alpha chains are shared by 4 different clones
      shared_alph.2 <- floor(0.085*numb_alph_used)                 # 15% of alpha chains are shared by 2 different clones
      shared_alph.3 <- floor(0.021*numb_alph_used)                 # 10% of alpha chains are shared by 3 different clones
      shared_alph.4 <- floor(0.007*numb_alph_used)                 #  5% of alpha chains are shared by 4 different clones
      shared_alph.5 <- floor(0.033*numb_alph_used)                 #  5% of alpha chains are shared by 4 different clones
      shared_alph.6 <- floor(0.005*numb_alph_used)                 #  5% of alpha chains are shared by 4 different clones
      shared_alph.7 <- floor(0.033*numb_alph_used)                 #  5% of alpha chains are shared by 4 different clones

    number.of.shared_alphs <- shared_alph.2 + shared_alph.3 + shared_alph.4 + shared_alph.5 + shared_alph.6 + shared_alph.7                             # number of (distinct) alpha chains that are shared by more than one clone
    number.distinct.alpha   <- numb_alph_used - (shared_alph.2 + 2*shared_alph.3 + 3*shared_alph.4
                                                    + 4*shared_alph.5 + 5*shared_alph.6 + 6*shared_alph.7)   # number of distinct alphas; take the number of times alphas are used then subtract the number of times shared alphas are used
    number.of.single.alphas <- number.distinct.alpha - number.of.shared_alphs                              # number of alphas involved in single TCR clone
    clones.per.alpha <- matrix(c(1:number.distinct.alpha, rep(1, number.distinct.alpha)), ncol = 2)         # preallocating matrix recording how many clones share alpha_i; col 1 = i, col 2 = # of clones sharing alpha_i
    clones.per.alpha[, 2]  <- sample(c(rep(1,number.of.single.alphas), rep(2, shared_alph.2),
                                       rep(3, shared_alph.3), rep(4, shared_alph.4),
                                       rep(5, shared_alph.5), rep(6, shared_alph.6),
                                       rep(7, shared_alph.7)
    ), replace = FALSE)                                                   # randomly assigning 1-4 clones sharing alpha_i

    # creating a vector of random alpha indices (each index repeated the number of times they're shared) for  random assignment of alphas to betas to form alpha/beta TCR pairs
    alpha <- vector(length = numb_alph_used)
    alph_added <- 1
    for(j in 1:nrow(clones.per.alpha)) {            # go through each row of clones.per.alpha;
      alph_ind  <- clones.per.alpha[j, 1]            # record the alpha index
      numb.alph <- clones.per.alpha[j, 2]               # record the number of times alph_ind is shared
      alpha[alph_added:(alph_added + numb.alph - 1)] <- alph_ind   # record alph_ind in numb.alph entries
      alph_added <- alph_added + numb.alph          # update next index
    } # end for - j
    alpha <- sample(alpha)        # randomly rerranged these alphas to form a random assignment (used below)

    # assiging random alphas to betas
    TCR <- matrix(ncol=3, nrow = numb_clones)                   # recording the pairs, three columns, beta = col 1, alpha = col 2 and 3; if single TCR, then col 2 == col 3
    ind_beta <- 1                                               # keeps track of which number of clones we are in the TCR matrix as we loop through all the betas
    ind_alph <- 1                                              # keeps track of where we are in assigning the alphas (ie which index/indices to use next in the vector alpha)
    ind_TCR <- 1                                                # ? i dunno yet this code is confusing
    for (betclone in 1:numb_beta) {                             # go through all of the beta chains
      numb.shared <- shared_mat[betclone, 2]                     # determine the number of clones that share this beta chain
      numb.dual   <- shared_mat[betclone, 3]                     # determine the number of dual TCR clones with this beta chain
      TCR[ind_TCR:(ind_TCR+numb.shared-1),1] <- betclone        # record this beta chain's index in the TCR matrix
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
      if (TCR[clon,2] == TCR[clon,3]) {
        if (ind_single > (numb_clones - numb_dual)) {
          try_again <- 1
        } else {
          TCR_single[ind_single,] <- TCR[clon,1:2]
          ind_single <- ind_single + 1
          try_again <- 0
        }
      } else {
        if (ind_dual > numb_dual) {
          try_again <- 1
        } else {
          TCR_dual[ind_dual,] <- TCR[clon,]
          ind_dual <- ind_dual + 1
          try_again <- 0
        }
      } # end if-else
    } # end for
    # randomizing the clones (i.e. randomize the rows of the TCR matrix)
  }
  TCR_dual <- TCR_dual[complete.cases(TCR_dual), ]
  TCR_sizes <- TCR[sample(1:numb_clones, replace = FALSE),]
  list(ordered = TCR, random = TCR_sizes, dual = TCR_dual)
} # end function
