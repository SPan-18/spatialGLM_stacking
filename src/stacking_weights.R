# find stacking weights wrt log score

CVXR_stacking_weights <- function(lpd_point, solver = "ECOS", verbose = TRUE){
  if(verbose){ message("Using ", solver, " solver to get stacking weights.") }
  if(sum(!is.finite(lpd_point)) > 0){
    minelpd <- min(lpd_point[is.finite(lpd_point)])
    lpd_point[!is.finite(lpd_point)] <- minelpd - 1
  }

  lpd_m = mean(lpd_point)
  lpd_point = lpd_point - lpd_m ## rescale the log-density for numerical stability
  exp_lpd_point <- exp(lpd_point)
  G <- ncol(lpd_point)

  w <- Variable(G, nonneg = TRUE)
  obj <- Maximize(sum(log(exp_lpd_point %*% w)))
  # constr <- list(sum(w) == 1, w >= 0)
  constr <- list(sum(w) == 1)
  prob <- Problem(objective = obj, constraints = constr)
  result <- solve(prob, solver = solver)
  if(verbose){ message("Solver status: ", result$status) }
  wts <- as.numeric(result$getValue(w))
  return(wts)
}

# CVXR_stacking_weights <- function(lpd_point, solver = "ECOS"){
#   message("Using ", solver, " solver to get stacking weights.")
#   lpd_m = mean(lpd_point)
#   lpd_point = lpd_point - lpd_m ## rescale the log-density for numerical stability
#   exp_lpd_point <- exp(lpd_point)
#   G <- ncol(lpd_point)
#   
#   w <- Variable(G)
#   obj <- Maximize(sum(log(exp_lpd_point %*% w)))
#   constr <- list(sum(w) == 1, w >= 0)
#   prob <- Problem(objective = obj, constraints = constr)
#   result <- solve(prob, solver = solver)
#   message("Solver status: ", result$status)
#   wts <- as.numeric(result$getValue(w))
#   return(wts)
# }
