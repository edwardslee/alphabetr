#' Read in alphabetr sequencing data into the binary matrix form needed by bagpipe()
#'
#' \code{read_alphabetr()} will read in two different forms of a csv file to
#' convert sequencing data using the alphabetr approach into the binary
#' matrices required by \code{\link{bagpipe}}. The csv file(s) can have one of
#' two forms. (1) A single csv file with three columns: column 1 containing
#' whether the sequence is "TCRA" or "TCRB"; column 2 containing the well
#' number; and column 3 containing the CDR3 sequence
#' (2) Two CSV files, one for TCRA and one for TCRB, with two columns: column
#' 1 containing the well number and column 2 containing the CDR3 sequence
#'
#' @param data To read in a 3-column csv file containing both TCRA and TCRB
#'    sequencing information
#' @param data_alpha To read in a 2-column csv file containing TCRA sequencing
#'    information. Must be used in conjunction with the \code{data_beta}
#'    argument and cannot be used with the \code{data} argument.
#' @param data_beta To read in a 2-column csv file containing TCRB sequencing
#'    information. Must be used in conjunction with the \code{data_alpha}
#'    argument and cannot eb used with the \code{data} argument.
#'
#' @return A list of two binary matrices that represent the sequencing data and
#'    two character vectors that give the CDR3 sequences associated with each
#'    chain index.
#'
#' @examples
#' \dontrun{
#'   dat <- read_alphabetr(data = "alphabetr_data.csv")
#'
#'   # saving the alpha and beta binary matrices
#'   data_alpha <- dat$alpha
#'   data_beta  <- dat$beta
#'
#'   # finding the cdr3 sequences of alpha_2 and beta_4 respectively
#'   cdr3_alpha2 <- dat$alpha_lib[2]
#'   cdr3_beta4  <- dat$beta_lib[4]
#' }
#' @export
read_alphabetr <- function(data = NULL, data_alpha = NULL, data_beta = NULL) {
  if (!is.null(data) & !is.null(data_alpha) & !is.null(data_beta))
    stop("Must choose to use just the data arg or the data_alph/data_beta arg. See documentation.")

  if (!is.null(data)) {
    # reading in the data
    df <- read.csv(data, stringsAsFactors = FALSE)
    if (ncol(df) != 3)
      stop("Must supply a csv with 3 colums if using the data arguments.
            col 1 specifies 'alpha' or 'beta', col 2 specifies well #,
            col 3 specifies the CDR3 sequence.")

    # subsetting for only alpha or only beta chains
    df_alph <- df[df[, 1] == "alpha" | df[, 1] == "TCRA", ]
    df_beta <- df[df[, 1] == "beta" | df[, 1] == "TCRB", ]

    # determining all of the distinct alpha and beta seq found in the data
    alphs <- unique(df_alph[, 3])
    betas <- unique(df_beta[, 3])

    # giving each sequence an index and making lookup tables for them
    ind_alph <- 1:length(alphs)
    ind_beta <- 1:length(betas)
    names(ind_alph) <- alphs
    names(ind_beta) <- betas

    # initializing the matrices used by bagpipe()
    mat_alph <- matrix(0, nrow = max(df[, 2]), ncol = length(alphs))
    mat_beta <- matrix(0, nrow = max(df[, 2]), ncol = length(betas))

    # determining the well index and the chain index for each sequence in df_alph
    well_vec <- df_alph[, 2]
    alph_vec <- ind_alph[df_alph[, 3]]
    for (i in seq_along(well_vec)) {
      mat_alph[well_vec[i], alph_vec[i]] <- 1
    }

    # determining the well index and the chain index for each sequence in df_beta
    well_vec <- df_beta[, 2]
    beta_vec <- ind_beta[df_beta[, 3]]
    for (i in seq_along(well_vec)) {
      mat_beta[well_vec[i], beta_vec[i]] <- 1
    }
  } else {
    # reading in the data
    df_alph <- read.csv(data_alpha, stringsAsFactors = FALSE)
    df_beta <- read.csv(data_beta, stringsAsFactors = FALSE)
    if (ncol(df_alph) != 2)
      stop("Must supply a csv with 2 cols if using the data_alph argument.
            col 1 specifies well #, col 2 specifies the CDR3 sequence.")
    if (ncol(df_beta) != 2)
      stop("Must supply a csv with 2 cols if using the data_beta argument.
            col 1 specifies well #, col 2 specifies the CDR3 sequence.")

    # determining all of the distinct alpha and beta seq found in the data
    alphs <- unique(df_alph[, 3])
    betas <- unique(df_beta[, 3])

    # giving each sequence an index and making lookup tables for them
    ind_alph <- 1:length(alphs)
    ind_beta <- 1:length(betas)
    names(ind_alph) <- alphs
    names(ind_beta) <- betas

    # initializing the matrices used by bagpipe()
    mat_alph <- matrix(0, nrow = max(df[, 2]), ncol = length(alphs))
    mat_beta <- matrix(0, nrow = max(df[, 2]), ncol = length(betas))

    # determining the well index and the chain index for each sequence in df_alph
    well_vec <- df_alph[, 2]
    alph_vec <- ind_alph[df_alph[, 3]]
    for (i in seq_along(well_vec)) {
      mat_alph[well_vec[i], alph_vec[i]] <- 1
    }

    # determining the well index and the chain index for each sequence in df_beta
    well_vec <- df_beta[, 2]
    beta_vec <- ind_beta[df[, 3]]
    for (i in seq_along(well_vec)) {
      mat_beta[well_vec[i], beta_vec[i]] <- 1
    }
  }
  list(alpha = mat_alph, beta = mat_beta, alpha_lib = alphs, beta_lib = betas)
}
