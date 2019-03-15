maximum_likelihood_intercepts_blah <- function(times, delta) {
  new_initial <- c(1.6, -0.2)
  FHT_only_intercepts(new_initial, times, delta)

  func_with_log_transform <- function(initial, times, delta) {
    new_initial <- c(initial[1], -exp(initial[2]))
    return(FHT_only_intercepts(new_initial, times, delta))
  }

  make_optim_function <- function(times, delta) {
    optimizable_function <- function(param_vector) {
      return(FHT_only_intercepts(new_initial, times, delta))
    }
  }

  make_y0_gradient <- function(times, delta) {
    optimizable_function <- function(param_vector) {
      y0 <- param_vector[1]
      mu <- param_vector[2]
      sigma2 <- 1
      return(
        -sum(loss_function_derivative_y0(y0, mu, sigma2, times, delta))
      )
    }
  }

  make_mu_gradient <- function(times, delta) {
    optimizable_function <- function(param_vector) {
      y0 <- param_vector[1]
      mu <- param_vector[2]
      sigma2 <- 1
      return(
        -sum(loss_function_derivative_mu(y0, mu, sigma2, times, delta))
      )
    }
  }

  y0_gradient <- make_y0_gradient(times, delta)
  mu_gradient <- make_mu_gradient(times, delta)

  initial <- c(0, 0)
  optim_result <- optim(par=initial, fn=func_with_log_transform,
                        gr=function(param_vector) {
                          return(y0_gradient(param_vector), mu_gradient(param_vector))
                        },
                        times, delta
  )

  #mean_time <- mean(times) # 47
  initial <- c(0, 0)

  ## BEST SO FAR -1.25775 16.55883

  initial <- c(log(47), log(1))

  gradient_func

  optim_result <- optim(par=initial, fn=func_with_log_transform, gr=NULL, times, delta)

  nlm_result <- nlm(FHT_only_intercepts, initial, times, delta_n)
  max_num_retries <- 1000
  retries <- 0

  ### CREATE GRID
  min_x <- -5
  max_x <- 5
  min_y <- -5
  max_y <- 5
  N <- 100
  x_values <- seq(from=min_x, to=max_x, length.out=N)
  y_values <- seq(from=min_y, to=max_y, length.out=N)
  grid <- sample(expand.grid(x_values, y_values))

  feasible_points <- c()
  minimum_values <- c()
  num_grid_points <- dim(grid)[1]
  for (i in 1:num_grid_points) {
    initial <- as.numeric(grid[i, ])
    nlm_result <- tryCatch(
      {
        nlm(FHT_only_intercepts, initial, times, delta)
      },
      error=function(cond) {
        return(nlm(FHT_only_intercepts, initial, times, delta))
      },
      warning=function(cond) {
        return(nlm(FHT_only_intercepts, initial, times, delta))
      },
      finally={}
    )

    if (nlm_result$iterations > 0) {
      minimum_values <- rbind(nlm_result$minimum, minimum_values)
      feasible_points <- rbind(nlm_result$estimate, feasible_points)
    }
  }


  for (i in 1:max_num_retries) {
    if (all(nlm_result != initial)) {
      break
    } else {
      retries <- retries + 1
      nlm_result <- tryCatch({
        nlm(FHT_only_intercepts, initial, times, delta)$estimate
      },
      error=function(cond) {
        return(initial)
      },
      warning=function(cond) {
        return(nlm(FHT_only_intercepts, initial, times, delta)$estimate)
      },
      finally={}
      )
    }
  }
  if (retries == max_num_retries) {
    stop("Couldn't find intercepts numerically.")
  }
  return(nlm_result)
}
