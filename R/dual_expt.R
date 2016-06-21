# rec <- matrix(nrow = 0, ncol = 9)
# colnames(rec) <- c("beta", "alpha1", "alpha2", "numb_2", "numb_3", "numb_0", "act_2", "act_3", "act_0")
# for (i in 1:100) {
#   print(i)
# dual_prop <- .1
# number_pl <- 5
# error_rate <- .15
# number_skewed <- 25
# percent_top <- .5
# distribution <- "linear"
# # numb_cells_vec <- rep(300, 12) #c(20, 50, 50, rep(100,3), rep(200, 3), rep(300, 3))
# # numb_cells_vec <- c(1, 20, 50, 50, rep(100,3), rep(200, 3), rep(300, 3))
# # numb_cells <- matrix(c(numb_cells_vec, 96, rep(8*(number_pl-1), 12)), ncol = 2)
# numb_cells_vec <- c(20, 20, 50, 50, rep(100,3), rep(200, 3), rep(300, 3))
# numb_cells <- matrix(c(numb_cells_vec, c(96, rep(8*(number_pl -1), 12))), ncol = 2)
#
# # numb_cells <- matrix(c(20, 20, 50, 50, rep(100,2), rep(200, 2), rep(300, 2), 24, 24, 8, 8, 8,8,8,8,8,8), ncol = 2)
# retry <- 1
# input_beta <- 1685
# while (retry == 1) {
#   TCRpairings <- create_clones(numb_beta = input_beta, dual = dual_prop)
#   TCR_clones  <- TCRpairings$random
#   TCR_answer  <- TCRpairings$ordered
#   TCR_dual    <- TCRpairings$dual
#   clone.size <- nrow(TCR_clones)
#
#
#   data <- create_data(TCR_clones, plates = number_pl, error = error_rate,
#                       skewed = number_skewed, prop_top = percent_top,
#                       dist = distribution,
#                       numb_cells = numb_cells)
#
#   # data <- createData(TCR_clones, plates = number_pl, error = error_rate,
#   #                   skewed = number_skewed, percent.top = percent_top,
#   #                   dist = distribution, number.cells.vector = numb_cells_vec)
#   dat_alph <- data$alpha
#   dat_beta <- data$beta
#   TCR <- data$ans
#
#   dual_top <- TCR.sizes[which(TCR.sizes[1:number_skewed, 2] != TCR.sizes[1:number_skewed, 3]), , drop = FALSE]; dual_top[order(dual_top[,1]), ]
#   if (nrow(dual_top) > 1) retry <- 0
# }
# cells.vec <- numb_cells_vec
# number.plates <- number_pl
# skewed.number <- number_skewed
# err <- error_rate
# pct.top <- percent_top
# dual.prop <- .1
#
# alpha = dat_alph
# beta = dat_beta
# # pair = bb
# error = 0.15
# cells = numb_cells
#
# number_plates <- nrow(alpha)/96  # number of plates
# max_beta <- ncol(beta)           # determine maximum beta index
#
# sample_size_well <- cells[, 1]     # number of cells per well
# numb_sample <- cells[, 2]          # number of wells w/ sample size
#
# # Determine the wells with the small sample sizes, which is necessary to
# # determine top clone duals properly
# small <- which(sample_size_well < 50)       # find the sample sizes < 50 cells
# if (length(small) > 1) {
#   if (any(small == 1)) {
#     small_wells <- c(1:numb_sample[small[1]])
#     for (i in 2:length(small)) {
#       j <- small[i]
#       small_wells <- c(small_wells,
#                        (sum(numb_sample[1:(j-1)]) + 1):sum(numb_sample[1:j])
#       )
#     }
#   } else {
#     small_wells <- vector()
#     for (i in 1:length(small)) {
#       j <- small[i]
#       small_wells <- c(small_wells,
#                        (sum(numb_sample[1:(j-1)]) + 1):sum(numb_sample[1:j])
#       )
#     }
#   }
# } else {
#   small_wells <- 1:numb_sample[small]
# }
#
# # look at only the wells with the small sample size
# alpha <- alpha[small_wells, ]
# beta <- beta[small_wells, ]
#
# # pre-allocate matrix to record the indices of the candidate dual TCR clones,
# # the shared and dual likelihoods, and the estimated frequency of the
# # candidate dual
#
#
# #
# # freq <- pair
# # freq <- freq[freq[, 4] > .9,]
# # freq <- freq[order(freq[, "MLE"], decreasing = TRUE), ]
#
# dual_top <- TCR.sizes[which(TCR.sizes[1:number_skewed, 2] != TCR.sizes[1:number_skewed, 3]), , drop = FALSE]; dual_top <- dual_top[order(dual_top[,1]), ]
#
# for (clon in seq_along(dual_top[, 1])) {
#   ind_beta  <- dual_top[clon, 1]
#   ind_alph1 <- dual_top[clon, 2]
#   ind_alph2 <- dual_top[clon, 3]
#   rank <- which(TCR.sizes[, 1] == ind_beta & ((TCR.sizes[, 2] == ind_alph1 & TCR.sizes[, 3] == ind_alph2) | (TCR.sizes[, 2] == ind_alph2 & TCR.sizes[, 3] == ind_alph1)))
#   ff <- prob1[rank]
#
#   # Determing which wells that contain the clone a clone is counted to be in a
#   # well if their component chains are found in the same well
#   sample_size_well <- cells[, 1]     # number of cells per well
#   numb_sample <- cells[, 2]          # number of wells w/ sample size
#
#   sample_size_well <- sample_size_well[small]
#   numb_sample <- numb_sample[small]
#
#   numb_distinct <- length(sample_size_well)
#   well_clone <- rep(0, numb_distinct)
#
#   well_clone <- matrix(ncol = 5, nrow = numb_distinct)
#   colnames(well_clone) <- c("ba1", "ba2", "a1a2", "ba1a2", "none")
#
#
#   for (size in 1:numb_distinct) {
#     ind1 <- cumsum(numb_sample[1:size])[size] - numb_sample[size] + 1
#     ind2 <- cumsum(numb_sample[1:size])[size]
#     well_clone[size, "ba1"] <- sum(beta[ind1:ind2, ind_beta] == 1 &
#                                      alpha[ind1:ind2, ind_alph1] == 1 &
#                                      alpha[ind1:ind2, ind_alph2] == 0)
#     well_clone[size, "ba2"] <- sum(beta[ind1:ind2, ind_beta] == 1 &
#                                      alpha[ind1:ind2, ind_alph1] == 0 &
#                                      alpha[ind1:ind2, ind_alph2] == 1)
#     well_clone[size, "a1a2"] <- sum(beta[ind1:ind2, ind_beta] == 0 &
#                                       alpha[ind1:ind2, ind_alph1] == 1 &
#                                       alpha[ind1:ind2, ind_alph2] == 1)
#     well_clone[size, "ba1a2"] <- sum(beta[ind1:ind2, ind_beta] == 1 &
#                                        alpha[ind1:ind2, ind_alph1] == 1 &
#                                        alpha[ind1:ind2, ind_alph2] == 1)
#     well_clone[size, "none"] <- numb_sample[size] - sum(well_clone[size, 1:4])
#   }
#
#
#   binomial_coeff <- list()
#   for (i in 1:numb_distinct) {
#     binomial_coeff[[i]] <- choose(sample_size_well[i], 1:sample_size_well[i])
#   }
#
#   multinomial_coeff <- list()
#   for(i in 1:numb_distinct) {
#     cells_per_well <- sample_size_well[i]
#     multi_mat <- matrix(0, nrow = cells_per_well, ncol = cells_per_well)
#     for (j in 1:(cells_per_well-1)) {
#       for (k in 1:(cells_per_well - j)) {
#         # print(c(i, j, k))
#         multi_mat[j, k] <- multicool::multinom(c(j, k, sample_size_well[i] - j - k),
#                                                counts = TRUE)
#       }
#     }
#     multinomial_coeff[[i]] <- multi_mat
#   }
#
#   # browser()
#   dual_LL <-  dual_prob_binom(ff, err = .15, numb_wells = well_clone,
#                               numb_cells = sample_size_well, binomials = binomial_coeff)
#   # browser()
#
#   act <- colSums(well_clone)
#
#   numb_2 <- dual_LL[1, 1] * 128
#   numb_3 <- dual_LL[1, 2] * 128
#   numb_0 <- dual_LL[1, 3] * 128
#   act_2 <- sum(act[1:3])
#   act_3 <- act[4]
#   act_0 <- act[5]
#
#   rec <- rbind(rec, c(ind_beta, ind_alph1, ind_alph2, numb_2, numb_3, numb_0, act_2, act_3, act_0))
# }
# }
# rec <- as.data.frame(rec)
# rec <- rec %>% mutate(diff2 = numb_2 - act_2, diff3 = numb_3 - act_3, diff0 = numb_0 - act_0)
