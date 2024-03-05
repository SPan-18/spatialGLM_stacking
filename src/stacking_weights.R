# find stacking weights wrt log score

CVXR_stacking_weights <- function(lpd_point, solver = "ECOS"){
  message("Using ", solver, " solver to get stacking weights.")
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
  
  # solver1 <- solver
  # solver2 <- setdiff(c("MOSEK", "ECOS"), solver1)
  # result <- tryCatch(solve(prob, solver = solver1),
  #                    error = function(ex){
  #                      warning("Had an error running solver ", solver1, "-", ex, 
  #                              "Using solver ", solver2)
  #                      solve(prob, solver = solver2)
  #                    })
  # 
  # wts <- structure(
  #   result$getValue(w)[,1],
  #   names = paste0("model", 1:G),
  #   class = c("stacking_weights")
  # )
  
  # solvefun <- function(){
  #   solve(prob, solver = solver)
  # }
  # 
  # expr <- substitute(solvefun())
  # result <- retry_expr(expr)
  result <- solve(prob, solver = solver)
  message("Solver status: ", result$status)
  wts <- as.numeric(result$getValue(w))
  return(wts)
}


# retry_expr <- function(expr, max.attempts = 3) {
#   x <- try(eval(expr))
#   
#   if(inherits(x, "try-error") && max.attempts > 0) {
#     return(retry_expr(expr, max.attempts - 1))
#   }
#   
#   x
# }

# elpd_mat <- read.table("../test/elpd_mat.txt")
# elpd_mat <- as.matrix(elpd_mat)
# w1 <- loo::stacking_weights(elpd_mat)
# 
# w2 <- CVXR_stacking_weights(elpd_mat, solver = "MOSEK")


