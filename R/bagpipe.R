#' @useDynLib alphabetr
#' @importFrom Rcpp sourceCpp
NULL

#' Identify candidate alpha/beta pairs.
#'
#' \code{bagpipe()} uses the alphabetr resampling procedure on sequencing data
#'    to identify candidate alpha/beta pairs. The procedure takes a subsample
#'    of the data without replacement, calculates association scores using ASOIDHAOISJ,
#'    and then for each well uses the Hungarian algorithm to determine the most
#'    likely pairings for the chains found in the well. Each time this is done
#'    is a replicate, and the number of replicates is specified as an option.
#'    A threshold is then used to filter the candidate pairs that appear in
#'    proportion of the replicates larger than the threshold, resulting in the
#'    final list of candidate pairs. Bagpipe is an acronym for
#'    \strong{b}ootstrapping \strong{a}lphabetr \strong{g}enerated
#'    \strong{p}a\strong{i}rs \strong{p}rocedur\strong{e} (based on older
#'    versions that utilized bootstrapping)
#'
#' @param alpha Matrix recording which alpha chains appear in each well of the
#'    data. See  .
#' @param beta Matrix recording which beta chains appear in the each well of the
#'    data. See .
#' @param replicates The number of times the resampling procedure is repeated,
#'    i.e. the number of replicates. At least 100 replicates is recommended.
#' @param frac The fraction of the wells resampled in each replicate. Default
#'    is 75\% of the wells
#' @param thres A vector of thresholds (from 0 to 1) to use to filter the
#'    results from all of the candidate alpha/beta pairs from the resampling
#'    procedure.
#' @param bootstrap Legacy option. Calls a bootstrapping strategy (which
#'    resamples with replacement) instead of sampling a subset without
#'    replacement.
#'
#' @return A n x 2 matrix where
#' @examples
#' set.seed(110)
#' createClones()
#' createData()
#' bagpipe(alpha = data_alph, beta = data_beta)
#' @export

bagpipe <- function(alpha, beta, replicates = 100, frac = 0.75,
                    thres = c(.3, .6, .9), bootstrap = FALSE) {
  # If user asks for bootstrapping, then the resampled size must be equal
  # to the original sample size
  if (bootstrap == TRUE) frac <- 1.0
  # Ensure that there are the same number of wells in the alpha and beta data
  if (nrow(alpha) != nrow(beta))
    stop("Different number of wells in the alpha and beta data inputs.")
  if (replicates < 1) stop("The replicates arguments need to be at least 1.")

  numb_plates <- nrow(alpha) / 96
  size_rep    <- floor(frac * nrow(alpha))  # number of wells to sample

  for (plat in 1:replicates) {
    # choose a random subset of the wells of the data without replacement
    # (unless the user wants to bootstrap and use replacement)
    ind <- sample(1:(96 * numb_plates), size_rep, replace = bootstrap)
    data_alph <- alpha[ind, ]
    data_beta <- beta[ind, ]

    # calculate association scores between
    pref <- chain_scores(data_a = data_alph, data_b = data_beta)

    # number of alpha/beta chains
    numb_alph  <- ncol(data_alph) # columns of data_alph is the number of alpha chains
    numb_beta  <- ncol(data_beta) # columns of data_beta is the number of beta chains

    # matrices of scores
    score_alph  <- pref$ascores   # matrices of scores that alpha (row i) has for beta (col j)
    score_beta  <- pref$bscores   # matrices of scores that beta (row i) has for beta (col j)
    rm(pref)

    # This section loops through each well, looks at the alpha and beta chains
    # in the well and use the Hungarian algorithm to choose the set of one-to-
    # one beta-alpha pairings that maximizes the sum of the association scores

    # Entry (i, j) of track_pairs keeps track the number of wells beta_i and
    # alpha_j are chosen in; a threshold is then calculated to filter out the
    # beta_i/alpha_j pairs that are most likely not true pairs
    track_pairs <- matrix(0, nrow = numb_beta, ncol = numb_alph)
    #track_assign <- list(length = nrow(data_beta))
    score_mat <- score_beta + t(score_alph)
    for (well in 1:nrow(data_beta)) {
      # find the alpha and beta chain indices in the well
      ind_alph <- which(data_alph[well, ] == 1)
      ind_beta <- which(data_beta[well, ] == 1)
      # Hungarian algorithm is posed differently depending on whether there are
      # more alpha or beta chains in the well
      if (length(ind_beta) <= length(ind_alph)) {
        # use the hungarian algorithm to obtain the beta_i/alpha_j pairings that
        # maximize the sum of the scores for the well
        assign <- clue::solve_LSAP(score_mat[ind_beta, ind_alph, drop = FALSE],
                                   maximum = TRUE)
        assign <- matrix(c(ind_beta, ind_alph[assign]), ncol = 2)
        #track_assign[[well]] <- assign
        track_pairs[assign] <- track_pairs[assign] + 1
      } else {
        assign <- clue::solve_LSAP(t(score_mat[ind_beta, ind_alph, drop=FALSE]),
                                   maximum = TRUE)
        assign <- matrix(c(ind_beta[assign], ind_alph), ncol = 2)
        #track_assign[[well]] <- assign
        track_pairs[assign] <- track_pairs[assign] + 1
      }
    }

    # Determine a threshold for the # number of wells that a candidate pair must
    # be picked in in order to make the cut
    threshold <- mean(track_pairs[track_pairs > 0])

    # arr.ind = TRUE will return a n by 2 matrix, which col 1 is the beta index
    # and col 2 is the alpha index; finding the beta, alpha pairs that are above
    # the threshold
    app <- which(track_pairs > threshold, arr.ind = TRUE)
    app <- app[order(app[, 1]), ]

    # Sorting the resulting candidate pairs into a list pair where the ith
    # element of the list contains the alpha indices pairing with beta_i
    list_pair <- list(length = numb_beta)
    for (i in 1:numb_beta) {
      # find the candidate pairs with the beta index and save those alpha
      # indices to the ith element of the list
      ind_beta <- which(app[, 1] == i)
      list_pair[[i]] <- unique(as.vector(app[ind_beta, 2]))
    } # end for i

    # recording the beta_i\alpha_j pairs from the output of the replicate to the
    # variable "list_pairX" where X is the X_th replicate
    result_name <- paste("list_pair", plat, sep = "")
    assign(result_name, list_pair)
  } # end for - plat

  # combine all the results of the replicates into one master list
  # The vector of the i_th element of list_master represents the indices of the
  # alpha chains paired with beta_i
  list_master <- mapply(c, list_pair1, list_pair2, SIMPLIFY = FALSE)
  rm(list_pair1, list_pair2)
  for (i in 3:replicates) {
    list_name   <- paste("list_pair", i, sep = "")
    list_master <- mapply(c, list_master, eval(parse(text = list_name)),
                          SIMPLIFY = FALSE)
    rm(list = list_name)
  }

  # determine the number of candidate pairs by adding up all the number of
  # unique alpha chain indices paired with each beta chain
  cand_pairs <- sum(sapply(lapply(list_master, unique), length))

  jack_results <- matrix(0, nrow = cand_pairs, ncol = 3)
  ind <- 1
  for (ind_beta in 1:length(list_master)) {
    if (length(list_master[[ind_beta]]) > 0) {
      alphas <- table(list_master[[ind_beta]])
      indices <- as.numeric(names(alphas))
      prop_replicates  <- as.vector(alphas) / replicates
      i1 <- ind
      i2 <- ind + length(unique(indices)) - 1
      jack_results[i1:i2, 1] <- ind_beta
      jack_results[i1:i2, 2] <- indices
      jack_results[i1:i2, 3] <- prop_replicates
      ind <- ind + length(unique(indices))
    }
  }

  jack_results <- jack_results[, c(1, 1, 2, 2, 3), drop = FALSE]
  colnames(jack_results) <- c("beta1", "beta2", "alpha1", "alpha2", "prop_replicates")
  jack_results
} # end function
