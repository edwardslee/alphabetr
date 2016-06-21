# # test_freq
# library(dplyr)
#
#
# dual_prop <- .1
# number_pl <- 1
# error_rate <- .15
# number_skewed <- 10
# percent_top <- .5
# distribution <- "linear"
# # numb_cells_vec <- rep(300, 12) #c(20, 50, 50, rep(100,3), rep(200, 3), rep(300, 3))
# numb_cells_vec <- rep(300, 12)#c(20, 50, 50, rep(100,3), rep(200, 3), rep(300, 3))
# numb_cells <- matrix(c(numb_cells_vec, rep(8*number_pl, 12)), ncol = 2)
#
# total <- 0
# correct <- 0
# for (i in 1:1000) {
#   TCRpairings <- create_clones(numb_beta = 1685, dual = dual_prop)
#   TCR_clones  <- TCRpairings$random
#   TCR_answer  <- TCRpairings$ordered
#   TCR_dual    <- TCRpairings$dual
#   clone.size <- nrow(TCR_clones)
#
#   data <- create_data(TCR_clones, plates = number_pl, error = error_rate,
#                       skewed = number_skewed, prop_top = percent_top,
#                       dist = distribution,
#                       numb_cells = numb_cells)
#   dat_alph <- data$alpha
#   dat_beta <- data$beta
#   TCR <- data$ans
#
#   number.clones <- nrow(TCR.sizes)
#   number.skewed <- number_skewed
#   rhs1 <- .5                                              # the top "number.skewed" clones represent 50% of the population
#   rhs2 <- rhs1/(number.clones - number.skewed)            # the last clone in the top 50% has the same freq as the clones in the second 50%
#   lhs11 <- (1+number.skewed-1)                            #
#   lhs21 <- -((number.skewed-1)/2 + (number.skewed-1)^2/2) #
#   lhs12 <- 1
#   lhs22 <- -(number.skewed - 1)
#   A <- matrix(c(lhs11, lhs12, lhs21, lhs22), nrow = 2)
#   b <- matrix(c(rhs1, rhs2), nrow = 2)
#   solveAb <- solve(A, b)
#   prob1 <- solveAb[1] - 0:(number.skewed-1) * solveAb[2]
#   prob2 <- rep(.5 / (number.clones - number.skewed), number.clones - number.skewed)
#
#
#   aa <- TCR.sizes[1:number_skewed, ]
#   bb <- freq_estimate(dat_alph, dat_beta, aa, error = error_rate, cells = numb_cells)
#   bb$answer <- prob1
#   numb_na <- bb %>% filter(is.na(MLE)) %>% summarize(n()) %>% as.vector()
#   corr <- bb %>% filter(!is.na(MLE) & answer < CI_up & answer > CI_low) %>% summarize(n())
#
#   total <- total + number_skewed - numb_na
#   correct <- correct + corr
# }
#
