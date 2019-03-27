#' @export

run_CV_clinical <- function(
  M, K_fold_repetitions, K, X, Z, times, delta, boost_intercepts_continually=boost_intercepts_continually
) {
  return(run_CV_cyclic_helper_mu(M_fixed=1, M, K_fold_repetitions, K, X, Z, times, delta, boost_intercepts_continually=boost_intercepts_continually))
}
