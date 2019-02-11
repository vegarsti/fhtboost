#' @export

plot_wiener_processes <- function() {
  set.seed(1)
  N <- 5
  y0 <- rep(10, N)
  mu <- -1
  time_increment <- 0.01
  times <- seq(0, 10, time_increment)
  increments <- matrix(rnorm(N*length(times), 0, sqrt(time_increment)), ncol=N)
  noise_matrix <- apply(increments,2,cumsum)
  y0_matrix <- matrix(rep(y0, length(times)), ncol=N, byrow=T)
  mu_matrix <- matrix(rep(mu*times,N),ncol=N)
  w <- y0_matrix + mu_matrix + noise_matrix
  matplot(times, w, type='l', ylab="Y(t)", xlab="t")
  abline(h=0)

  hitting_times <- sapply(apply(w, 2, function(col) which(col < 0)), min)
}
