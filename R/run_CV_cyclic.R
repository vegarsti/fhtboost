run_CV_cyclic <- function(
  M1, M2, K_fold_repetitions, K, X, Z, times, delta, boost_intercepts_continually=boost_intercepts_continually
) {
  for (M_fixed in 1:M2) {
    run_CV_cyclic_helper_y0(M, M_fixed, K_fold_repetitions, K, X, Z, times, delta, boost_intercepts_continually=boost_intercepts_continually)
  }
}
