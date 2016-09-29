read_alphabetr <- function(data = NULL, data_alpha = NULL, data_beta = NULL) {
  if (!is.null(data) & !is.null(data_alpha) & !is.null(data_beta))
    stop("Must choose to use just the data arg or the data_alph/data_beta arg. See documentation.")

  if (!is.null(data)) {
    # reading in the data
    df <- read.csv(data)
    if (ncol(df) != 3)
      stop("Must supply a csv with 3 colums if using the data arguments.
            col 1 specifies 'alpha' or 'beta', col 2 specifies well #,
            col 3 specifies the CDR3 sequence.")

    # subsetting for only alpha or only beta chains
    df_alph <- df[df[, 1] == "alpha", ]
    df_beta <- df[df[, 1] == "beta", ]

    # determining all of the distinct alpha and beta seq found in the data
    alphs <- unique(df_alph[, 3])
    betas <- unique(df_beta[, 3])

    # giving each sequence an index and making lookup tables for them
    ind_alph <- 1:length(alphs)
    ind_beta <- 1:length(betas)
    names(ind_alph) <- alphs
    names(ind_beta) <- betas

    # initializing the matrices used by bagpipe()
    mat_alph <- matrix(nrow = max(df[, 1]), ncol = length(alphs))
    mat_beta <- matrix(nrow = max(df[, 1]), ncol = length(betas))

    # determining the well index and the chain index for each sequence in df_alph
    well_vec <- df_alph[, 2]
    alph_vec <- ind_alph[df[, 3]]
    mat_alph[well_vec, alph_vec] <- 1

    # determining the well index and the chain index for each sequence in df_beta
    well_vec <- df_beta[, 2]
    beta_vec <- ind_beta[df[, 3]]
    mat_alph[well_vec, beta_vec] <- 1
  } else {
    # reading in the data
    df_alph <- read.csv(data_alpha)
    df_beta <- read.csv(data_beta)
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
    mat_alph <- matrix(nrow = max(df[, 1]), ncol = length(alphs))
    mat_beta <- matrix(nrow = max(df[, 1]), ncol = length(betas))

    # determining the well index and the chain index for each sequence in df_alph
    well_vec <- df_alph[, 2]
    alph_vec <- ind_alph[df[, 3]]
    mat_alph[well_vec, alph_vec] <- 1

    # determining the well index and the chain index for each sequence in df_beta
    well_vec <- df_beta[, 2]
    beta_vec <- ind_beta[df[, 3]]
    mat_alph[well_vec, beta_vec] <- 1
  }
}
