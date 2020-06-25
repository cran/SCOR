#' Optimizing Different Estimators Of Hyper Volume Under Manifold
#'
#' As we know `SCOptim` is efficient in estimating maximizing Hyper Volume Under Manifolds Estimators,
#' we made some pre-functions that optimizes specific Problems of EHUM,SHUM and ULBA.
#'
#' Optimization of EHUM, SHUM and ULBA using SCOptim.
#'
#' @param beta_start The initial guess for optimum \eqn{\beta} by user
#' @param labels Sample Sizes vector of that has number of elements in each category. It works like the labels of data matrix.
#' @param x_mat The Data Matrix
#' @param p This parameter exists for the case of optimized_SHUM only.p decides whether to use \eqn{s_n(x)} or \eqn{\phi_n(x)}. p = 1 stands for \eqn{\phi_n(x)} and p = 0 stands for \eqn{s_n(x)}
#' @param rho Step Decay Rate with default value 2
#' @param phi Lower Bound Of Global Step Size. Default value is \eqn{10^{-6}}
#' @param max_iter Max Number Of Iterations In each Run. Default Value is 50,000.
#' @param s_init Initial Global Step Size. Default Value is 2.
#' @param tol_fun Termination Tolerance on the function value. Default Value is \eqn{10^{-6}}
#' @param tol_fun_2 Termination Tolerance on the difference of solutions in two consecutive runs. Default Value is \eqn{10^{-6}}
#' @param minimize Binary Command to set SCOptim on minimization or maximization. FALSE is for minimization which is set default.
#' @param time Time Allotted for execution of SCOptim
#' @param print Binary Command to print optimized value of objective function after each iteration. FALSE is set fault
#' @param lambda Sparsity Threshold. Default value is \eqn{10^{-3}}
#' @param parallel Binary Command to ask SCOptim to perform parallel computing. Default is set at TRUE.
#'
#' @return Optimum Values Of HUM Estimates
#'
#' @examples
#'
#' \donttest{
#' R <- optimized_SHUM(rep(1, 12), colnames(AL), AL, parallel = FALSE)
#' estimate_SHUM(R, colnames(AL), AL)
#' # This run will take about 10 mins on average based on computational capacity of the system
#' # Optimum value of HUM estimate noticed for this case : 0.8440681
#' }
#'
#' R <- optimized_EHUM(rep(1, 12), colnames(AL), AL, parallel = FALSE)
#' estimate_EHUM(R, colnames(AL), AL)
#' # Optimum value of HUM estimate noticed for this case : 0.8403805
#'
#' R <- optimized_ULBA(rep(1, 12), colnames(AL), AL, parallel = FALSE)
#' estimate_ULBA(R, colnames(AL), AL)
#' # Optimum value of HUM estimate noticed for this case : 0.9201903
#'
#' @name optimized_HUM
NULL

#' @rdname optimized_HUM
#' @export
optimized_EHUM <-
  function(beta_start,
           labels,
           x_mat,
           rho = 2,
           phi = 1e-3,
           max_iter = 5e+04,
           s_init = 2,
           tol_fun = 1e-6,
           tol_fun_2 = 1e-6,
           minimize = FALSE,
           time = 3.6e+04,
           print = FALSE,
           lambda = 1e-3,
           parallel = TRUE)
  {
    stopifnot(nrow(x_mat) == length(beta_start))
    stopifnot(ncol(x_mat) == length(labels))

    # create a function that calculates EHUM, to perform SCOptim with other parameters fixed except beta
    func <- function(beta_start)
      return(estimate_EHUM(beta_start, labels, x_mat))

    R <-
      SCOptim(
        beta_start,
        func,
        rho,
        phi,
        max_iter,
        s_init,
        tol_fun,
        tol_fun_2,
        minimize,
        time,
        print,
        lambda,
        parallel
      )
    #Thus the performing function of the vaue obtained through SCOptim gives maximum value
    return(R)
  }

#' @rdname optimized_HUM
#' @export
optimized_SHUM <-
  function(beta_start,
           labels,
           x_mat,
           p = 0,
           rho = 2,
           phi = 1e-3,
           max_iter = 5e+04,
           s_init = 2,
           tol_fun = 1e-6,
           tol_fun_2 = 1e-6,
           minimize = FALSE,
           time = 3.6e+04,
           print = FALSE,
           lambda = 1e-3,
           parallel = TRUE)
  {
    stopifnot(nrow(x_mat) == length(beta_start))
    stopifnot(ncol(x_mat) == length(labels))

    # create a function that calculates SHUM, to perform SCOptim with other parameters fixed except beta
    func <- function(beta_start)
      return(estimate_SHUM(beta_start, labels, x_mat, p))

    R <-
      SCOptim(
        beta_start,
        func,
        rho,
        phi,
        max_iter,
        s_init,
        tol_fun,
        tol_fun_2,
        minimize,
        time,
        print,
        lambda,
        parallel
      )
    #Thus the performing function of the vaue obtained through SCOptim gives maximum value
    return(R)
  }

#' @rdname optimized_HUM
#' @export
optimized_ULBA <-
  function(beta_start,
           labels,
           x_mat,
           rho = 2,
           phi = 1e-3,
           max_iter = 5e+04,
           s_init = 2,
           tol_fun = 1e-6,
           tol_fun_2 = 1e-6,
           minimize = FALSE,
           time = 3.6e+04,
           print = FALSE,
           lambda = 1e-3,
           parallel = TRUE)
  {
    stopifnot(nrow(x_mat) == length(beta_start))
    stopifnot(ncol(x_mat) == length(labels))

    # create a function to calculate ULBA, to perform SCOptim with other parameters fixed except beta
    func <- function(beta_start)
      return(estimate_ULBA(beta_start, labels, x_mat))

    R <-

      SCOptim(
        beta_start,
        func,
        rho,
        phi,
        max_iter,
        s_init,
        tol_fun,
        tol_fun_2,
        minimize,
        time,
        print,
        lambda,
        parallel
      )
    #Thus the performing function of the vaue obtained through SCOptim gives maximum value
    return(R)
  }
