# find stacking weights wrt log score

CVXR_stacking_weights <- function(lpd_point, solver = "ECOS"){
  message("Using ", solver, " solver to get stacking weights.")
  lpd_m = mean(lpd_point)
  lpd_point = lpd_point - lpd_m ## rescale the log-density for numerical stability
  exp_lpd_point <- exp(lpd_point)
  G <- ncol(lpd_point)
  
  w <- Variable(G)
  obj <- Maximize(sum(log(exp_lpd_point %*% w)))
  constr <- list(sum(w) == 1, w >= 0)
  prob <- Problem(objective = obj, constraints = constr)
  result <- solve(prob, solver = solver)
  
  wts <- structure(
    result$getValue(w)[,1],
    names = paste0("model", 1:G),
    class = c("stacking_weights")
  )
  flush.console()
  return(wts)
}