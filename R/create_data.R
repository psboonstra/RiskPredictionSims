
#' Simulate binary outcomes from a logistic model with center-level random
#' intercepts and year-within-center random intercepts
#'
#' @param num_years Number of years to simulate (integer, scalar)
#' @param max_num_centers Maximum possible number of centers to simulate. Since
#'   runs, i.e. observations, within a center within a year are generated from a
#'   Poisson, a given center may have no run in a given year (integer, scalar)
#' @param alpha Global intercept. (double, scalar)
#' @param gamma_var Variance of the center-level random intercepts (double,
#'   scalar, positive)
#' @param omega_var Variance of the year-within-center random intercepts
#'   (double, scalar, positive)
#' @param omega_cor AR(1) correlation of the year-within-center random
#'   intercepts (double, scalar, in (0, 1))
#' @param beta Regression coefficient for the run-specific risk factor, which is
#'   N(0,1) (double, scalar)
#' @param mean_runs_per_center Optional vector as long as `max_number_centers`
#'   corresponding to the mean of the Poisson distribution. If provided, the
#'   number of runs from each center from each year will be independently
#'   generated from a Poisson distribution with this mean. If not provided, this
#'   is created from a shifted half-t distribution to mimic right-skewness, i.e.
#'   many centers with small numbers of runs (double, vector, positive)
#' @param seed Random seed for reproducibility (integer, scalar, positive)
#'
#' @return A data.frame
#' @export
#'
#' @importFrom mnormt rmnorm
#' @importFrom purrr map2
#' @importFrom stats rt rpois rnorm plogis rbinom
#'
#'
#' @examples create_data()
create_data <- function(num_years = 5,
                        max_num_centers = 400,
                        alpha = -1.75,
                        gamma_var = 0.5^2,
                        omega_var = 0.15^2,
                        omega_cor = 0.9,
                        beta = 1,
                        mean_runs_per_center = NULL,
                        seed = sample(.Machine$integer.max, 1)) {

  set.seed(seed)
  if(!is.null(mean_runs_per_center)) {
    stopifnot(length(mean_runs_per_center) == max_num_centers)
  } else {
    mean_runs_per_center <- 1 + 3 * abs(rt(max_num_centers, 2))
  }

  omega_cov <- omega_var * omega_cor ^ abs(matrix(1:num_years - 1, nrow = num_years, ncol = num_years, byrow = TRUE) - (1:num_years - 1))

  runs_per_center <-
    replicate(n = num_years, rpois(max_num_centers, mean_runs_per_center))
  year_index <- unlist(map2(seq_len(num_years), colSums(runs_per_center), ~rep(.x, .y)))
  center_index <- unlist(map2(rep(seq_len(max_num_centers), times = num_years), as.numeric(runs_per_center), ~rep(.x, .y)))

  gamma <- rnorm(max_num_centers, 0, sqrt(gamma_var))
  omega <- rmnorm(max_num_centers, varcov = omega_cov)

  x <- rnorm(sum(runs_per_center), 0, sd = 1.0)

  true_fixeff <- alpha + x;
  true_raneff <- gamma[center_index] + omega[cbind(center_index, year_index)];
  true_prob <- plogis(true_fixeff + true_raneff);
  y <- rbinom(length(true_prob), 1, true_prob);

  all_dat <-
    data.frame(y = y,
               x = x,
               center_id = center_index,
               year_id = year_index,
               true_fixeff = true_fixeff,
               true_raneff = true_raneff,
               true_prob = true_prob)
}
