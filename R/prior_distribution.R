#' @title Multivariate Uniform Distribution on a ball
#' @description This function generates random samples from multivariate uniform distribution on a ball in \strong{R}^{d}, equipped with \eqn{L^{2}} norm (\emph{i.e.,} Euclidean distance), centered in \strong{0}, with radius R.
#' @param n number of desired samples.
#' @param d positive integer, representing the dimension of the observations.
#' @param R positive real value, the radius of the ball in \strong{R}^{d}.
#' @details This function generates samples from the multivariate uniform distribution whose density is \deqn{\pi(c, R) = \Gamma(d/2 + 1)/ \pi^{d/2} * 1/(R)^{d} 1_{B_{d}(R)}(c),} where \eqn{1_{B_{d}(R)}} is a centered \eqn{L^{2}} ball with radius R.
#' @importFrom stats runif
#' @return a matrix of \code{n} samples of length d.
#' @examples
#' ##generating 10000 samples from uniform distribution on a unit ball in \strong{R}^{2}
#' result <- runiform_ball(10000, 2, 1)
#' plot(result)
#' @export

runiform_ball = function(n, d, R){
  simulations = matrix(NA, nrow = n, ncol = d)
  for (i in 1:n){
      simul = runif(d, min = -R, max = R)
      while(sqrt(sum(simul^2)) > R){
      simul = runif(d, min = -R, max = R)
      }
      simulations[i,] = simul
    }
  return (simulations)
}



