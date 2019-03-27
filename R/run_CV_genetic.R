#' @export

run_CV_genetic <- function(
  M, K_fold_repetitions, K, X, Z, times, delta, boost_intercepts_continually=boost_intercepts_continually
) {
  return(run_CV_cyclic_helper_y0(M, M_fixed=1, K_fold_repetitions, K, X, Z, times, delta, boost_intercepts_continually=boost_intercepts_continually))
}
