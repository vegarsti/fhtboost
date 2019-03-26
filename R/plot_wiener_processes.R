#' @export

plot_wiener_processes <- function(y0, mu) {
  set.seed(1)
  N <- 100
  time_increment <- 0.01
  times <- seq(0, 40, by=time_increment)
  increments <- matrix(
    rnorm(
      N*length(times), mean=0, sd=sqrt(time_increment)
    ), ncol=N)
  noise_matrix <- apply(increments, 2, cumsum)
  y0_matrix <- matrix(rep(y0, length(times)*N), ncol=N, byrow=T)
  mu_matrix <- matrix(rep(mu*times, N), ncol=N)
  w <- y0_matrix + mu_matrix + noise_matrix
  w_before_max_time <- w[1:1300, ]
  below_0_2 <- apply(w_before_max_time, 2, min) < 0
  below_0 <- apply(w, 2, min) < 0
  dead <- sum(below_0)
  red_colors <- which(below_0)
  colors <- rep(rgb(0, 0, 0, alpha = 0.2), N)
  colors[red_colors] <- rgb(1, 0, 0, alpha = 0.2)
  #pdf("lol.pdf", width=12, height=6)
  matplot(times, w, type='l', ylab="Y(t)", xlab="t", ylim=c(-1, 20), xlim=c(0, 40), col=colors)
  abline(h=0, lty=3)
  #dev.off()
}
