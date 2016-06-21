#' Calculate association scores between alpha and beta chain pairs.
#'
#' \code{chain_scores()} calculates association scores between every pair of
#'    alpha and beta chains based on the number of concurrent well appearances
#'    each alpha and beta pair makes, scaled inversely by the number of unique
#'    chains in that well. See Lee et. al. 2016 for more information about this
#'    procedure.
#'
#' @param data_a Matrix recording which alpha chains appear in each well of the
#'    data. See .
#' @param beta Matrix recording which beta chains appear in the each well of the
#'    data. See .
#' @param numb_wells Vector containing the number of cells sampled in the wells
#'    of each column of the plates.
#'
#' @return A list containing the alpha and beta association scores. Accessed
#'    with list$alpha and list$beta respectively.
#' @examples
#'    set.seed(110)
#'    createClones()
#'    createData()
#'    scores <- chain_scores(data_a = data_alph, data_b = data_beta,
#'                           numb_cells = cells_vec)
#'    scores <- scores$alpha + scores$beta
#' @export
chain_scores <- function(data_a, data_b) {
  # Determining the number of unique alpha and beta chains and number of wells
  # Each column of data_a and data_b represents an alpha/beta chain respectively
  # Each row is a well; entry (i,j) determines if chain j is present in well i
  numb_alph <- ncol(data_a)
  numb_beta <- ncol(data_b)
  numb_well <- nrow(data_a)
  if (numb_well != nrow(data_b)) stop()

  # Creating the matrices that will record the association scores between
  # each pair of alpha and beta chain. In scores_alph matrix, entry (i,j)
  # represents how associated beta_j is associated with alpha_i (analogous
  # situation for scores_beta matrix). These two matrices will be added to get
  # a composite score (i.e. S_ij = scores_alph[i, j] + scores_beta[j, i])
  scores_alph <- matrix(0, nrow = numb_alph, ncol = numb_beta)
  scores_beta <- matrix(0, nrow = numb_beta, ncol = numb_alph)

  # loop through each well, find which alphas and betas are present, and
  # calculate association scores based on concurrent well appearances
  for (well in 1:numb_well) {
    # find which alphas and betas are in the well
    well_alph <- which(data_a[well, ] == 1)
    well_beta <- which(data_b[well, ] == 1)

    # Each concurrent apperance is scaled by the number of unique partner chains
    # found in the well (e.g. for a given beta, each concurrent alpha appearance
    # is scaled inversely by the number of unique alphas found in that well)
    scale_alph <- length(well_alph)
    scale_beta <- length(well_beta)
    scores_alph[well_alph, well_beta] <- scores_alph[well_alph, well_beta] +
                                          1/scale_alph
    scores_beta[well_beta, well_alph] <- scores_beta[well_beta, well_alph] +
                                          1/scale_beta
  } # end for - well

  list(ascores = scores_alph, bscores = scores_beta)
} # end function
