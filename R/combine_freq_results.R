#' @export
combine_freq_results <- function(single, dual) {
  remove_dual <- vector(length = 2*nrow(dual))
  for (dd in 1:nrow(dual)) {
    ind_beta  <- dual[dd, "beta1"]
    ind_alph1 <- dual[dd, "alpha1"]
    ind_alph2 <- dual[dd, "alpha2"]

    clone1 <- which(single[, "beta1"] == ind_beta & single[, "alpha1"] == ind_alph1)
    clone2 <- which(single[, "beta1"] == ind_beta & single[, "alpha1"] == ind_alph2)
    remove_dual[2*dd - 1] <- clone1
    remove_dual[2*dd]     <- clone2
  }
  single <- single[-remove_dual, ]
  freq_results <- rbind(single, dual)
  freq_results <- freq_results[order(freq_results[, "MLE"], decreasing = TRUE),, drop = FALSE]
  rownames(freq_results) <- 1:nrow(freq_results)
  freq_results
}
