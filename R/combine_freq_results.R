#' Combines the frequency estimation results from single TCR clones and dual TCR clones
#'
#' \code{combine_freq_results()} combines the results of the frequency estimation
#' performed on single TCR clones (from the output of \code{\link{bagpipe}}) and
#' the frequency estimation performed on dual clones. The code will find the
#' rows of the single TCR frequency results that are represented by the dual
#' clones and replace them with the appropriate dual clone entry.
#'
#' @param single Frequency estimation results of single TCR clones (usually from
#'   the first time \code{\link{freq_estimate}} is called)
#' @param dual Frequency estimation results of dual TCR-alpha clones
#'
#' @return A data.frame with the same structure as the output of
#'    \code{\link{freq_estimate}}. If two single "clones" in the \code{single}
#'    data.frame is represented by a dual clone in \code{dual}, then it is
#'    removed and replaced with one row represented by the dual clone.
#'
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
  freq_results <- freq_results[order(freq_results[, "MLE"], decreasing = TRUE), , drop = FALSE]
  rownames(freq_results) <- 1:nrow(freq_results)
  freq_results
}
