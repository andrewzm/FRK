## Diagnostic functions
.RMSPE <- function(pred, true) sqrt(mean((pred - true)^2))
.Coverage <- function(lower, upper, true) mean((true > lower) & (true < upper))
.IS90 <- function(lower, upper, true, alpha = 0.1) {
  
  ## lower and upper denote the lower and upper bounds of the prediction interval
  ## Interval score for all locations:
  ISs <- (upper - lower) + 2/alpha * (lower - true) * (true < lower) +
    2/alpha * (true - upper) * (true > upper)
  
  ## Average interval score:
  mean(ISs)
}

.CRPS_wrapper <- function(true, samples) {
  if(is.null(samples)) {
    return(NA)
  } else {
    return(scoringRules::crps_sample(y = true, dat = samples))
  }
}

.diagnostic_stats <- function(lower, upper, pred = rowMeans(cbind(lower, upper)), samples, true, ignoreCRPS = FALSE) {
  
    out <- data.frame(RMSPE = .RMSPE(true, pred),
                      CRPS = mean(.CRPS_wrapper(true = true, samples = samples)),
                      IS90 = .IS90(lower, upper, true),
                      Coverage = .Coverage(lower, upper, true)) 
  
}



#' Summary of diagnostic functions for the FRK model
#' 
#' This function provides a summary of several diagnostic functions used to 
#' assess performance of the model. The diagnostics are root mean squared prediction error, 
#' the interval score (assuming a 90% interval), coverage, and, if a Monte Carlo sample is supplied, 
#' the continuous ranked probability score. The function assumes the same format that 
#' the \code{predict()} function returns when \code{method = 'TMB'}.
#' 
#' @param M SRE object.
#' @param newdata Dataframe returned by \code{predict}.
#' @param MC List of Monte Carlo samples returned by \code{predict}. This argument
#' is only used to compute the continuous ranked probability score; set \code{MC = NULL}
#' if this is not desired. 
#' @param true_value Named list, where the names indicate the quantities of interest.
#' Must contain at least one element named 'Y', 'mu', or 'prob'.
#' @param out_sample Are we interested in out of sample locations? If \code{out_sample = TRUE}, 
#' only the unobserved BAUs are used for computing the diagnostics. 
#' @return A data.frame containing diagnostics for each quantity of interest.
.diagnostic_stats_FRK <- function(M, newdata, MC = NULL, true_value, out_sample = TRUE) {
  
  ## True value must be a named list, where the names indicate the quantities of interest.
  if (class(true_value) != "list" || !any(names(true_value) %in% c("Y", "mu", "prob"))) {
    stop("true_value must be a named list, with at least one element named 'Y', 'mu', or 'prob'.")
  }
  
  ## If we are interested in out of sample diagnostics, focus on unobserved locations only.
  if (out_sample == TRUE) {
    obsidx   <- apply(M@Cmat, 1, function(x) which(x == 1)) # Observed BAUs indices
    newdata <- newdata[-obsidx, ]
    MC <- lapply(MC, function(x)x[-obsidx,])
    true_value <- lapply(true_value, function(x) x[-obsidx])
  }
  
  ## Compute the summary statistics (and output as a dataframe)
  out <- data.frame(RMPSE = numeric(), CRPS = numeric(), IS90 = numeric(), Cov90 = numeric())
  
  for (i in c("Y", "mu", "prob")) {
    if(i %in% names(true_value)) {
      temp <- .diagnostic_stats(lower = newdata[, paste0(i, "_percentile_05")], 
                               upper = newdata[, paste0(i, "_percentile_95")],
                               true  = true_value[[i]],
                               pred  = newdata[, paste0("p_", i)], 
                               samples = MC[[paste0(i, "_samples")]])
      out  <- rbind(out, temp)
      rownames(out)[nrow(out)] <- i
    }
  }

  return(out)
}
