% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PACBO.R
\name{PACBO}
\alias{PACBO}
\title{PACBO}
\usage{
PACBO(mydata, R, coeff = 2, K_max = 50, scaling = FALSE,
  var_ind = FALSE, N_iterations = 500, plot_ind = FALSE, axis_ind = c(1,
  2))
}
\arguments{
\item{mydata}{a matrix where each row corresponds to an observation of length d.}

\item{R}{a positive real value that should be larger than the maximum Euclidean distance of all the observations in \code{mydata}. We recommend to set R equaling to this maximum Euclidean distance.}

\item{coeff}{a positive real value, enforcing large number of cells. The default, 2, should be convenient for most users. A larger value brings more cells for the clustering.}

\item{K_max}{a positive integer indicating the maximum number of cells allowed for the clustering.}

\item{scaling}{logical indicating whether the matrix \code{mydata} should be centered and scaled. The centering is done by subtracting the column means of \code{mydata} from their corresponding columns; the scaling is done by dividing the (centered) columns of \code{mydata} by their standard deviations. We recommend to set it to \code{TRUE} only when the maximum Euclidean distance of all the observations in \code{mydata} is smaller than 1.}

\item{var_ind}{logical indicating whether predicted centers of cells will be calculated sequentially. If \code{TRUE}, at each round, predicted centers of cells will be calculated on the basis of the past observations and past predicted centers. Setting this to \code{FALSE} will largely save execution time.}

\item{N_iterations}{a positive integer indicating the number of iterations of algorithm.}

\item{plot_ind}{logical indicating whether clusters should be plotted.}

\item{axis_ind}{numeric indicating which axes are to be plotted if d >= 2. The default is the first two coordinates of observations.}
}
\value{
Returns a list including
\item{predicted_centers}{a matrix of predicted centers of cells, where each row corresponds to a center.}
\item{nb_of_clusters}{positive integer indicating the estimation of the number of cells for the dataset.}
\item{labels}{labels for observations in \code{mydata}.}
}
\description{
This function performs clustering on online datasets. The number of cells is data-driven and need not to be chosen in advance by the user.
}
\details{
The PACBO algorithm is introduced and fully described in Le Li, Benjamin Guedj, Sebastien Loustau (2016), "PAC-Bayesian Online Clustering" (\url{https://arxiv.org/abs/1602.00522}). It relies on PAC-Bayesian approach, allowing for a dynamic (\emph{i.e.,} time-dependent) estimation of the number of clusters, up to \code{K_max} clusters. Its implementation is done via an RJMCMC-flavored algorithm.
}
\examples{
## generating 4 clusters of 100 points in \\strong{R}^{5}.
set.seed(100)
Nb <- 4
d <- 5
T <- 100
proportion = rep(1/Nb, Nb)
Mean_vectors <- matrix(runif(d*Nb,min=-10, max=10),nrow=Nb,ncol=d, byrow=TRUE)
mydata <- matrix(replicate(T, rmnorm(1, mean= Mean_vectors[sample(1:Nb, 1, prob = proportion),],
varcov = diag(1,d))), nrow = T, byrow=T)
R <- max(sqrt(rowSums(mydata^2)))
##run the algorithm.
result <- PACBO(mydata, R, plot_ind = TRUE)
}
\author{
Le Li <le@iadvize.com>
}
\references{
Le Li, Benjamin Guedj and Sebastien Loustau (2016), PAC-Bayesian Online Clustering, arXiv preprint: \url{https://arxiv.org/abs/1602.00522}.
}

