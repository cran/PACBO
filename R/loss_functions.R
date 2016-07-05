#' @title Instantaneous Loss
#' @description This function computes the instantaneous loss (\emph{i.e.,} squared Euclidean distantce) of an observation to its nearest center within a set of centers.
#' @param centers a matrix containing m centers of length d, where each row corresponds to d coordinates of a center.
#' @param instant_observation a vector of length d.
#' @details  Given a set \emph{C} of m centers of length d (\emph{i.e.,} \emph{C} = \eqn{{c_{1}, c_{2}, \dots, c_{m}}}) and a vector observation x of length d, this function computes the squared euclidean distance of x to its nearest center within \emph{C}, \emph{i.e.,} \deqn{L(x,C) = min_{1<= i <= m}|x-c_{i}|_{2}^{2}.}
#' @return The squared Euclidean distance of \code{instant_observation} to its nearest center within \code{centers}.
#' @examples
#' ## generating 4 centers of length 3.
#' centers <- matrix(1:12, nrow = 4, ncol = 3)
#' instant_observation <- c(2,6,10)
#' instantaneous_loss(centers, instant_observation)
#' @keywords internal
#' @export


instantaneous_loss = function(centers, instant_observation){
  if (class(centers) != 'matrix'){
    return (sum((centers-instant_observation)^2))
  }else{
    replicate_instant_observation = t(replicate(length(centers[,1]), instant_observation))
    distances = apply((centers - replicate_instant_observation)^2, 1, sum)
    return (min(distances))
  }
}


#' @title Cumulative Loss
#' @description This function computes the sum of squared Euclidean distance of each observation to its nearest center.
#' @param centers a matrix containing m centers of length d, where each row corresponds to coordinates of a center.
#' @param observations a matrix containing T observations of length d, where each row of the matrix is an observation of length d.
#' @details Given a set \emph{C} of m centers of length d (\emph{i.e.,} \emph{C} = \eqn{{c_{1}, c_{2}, \dots, c_{m}}}) and a set \emph{X} of T observations of length d (\emph{i.e.,} \emph{X} = \eqn{{x_{1}, x_{2}, \dots, x_{T}}}), this function computes the sum of squared euclidean distance of each observation in X to its nearest center in \emph{C}, \emph{i.e.,} \deqn{S_{T}(C) =\sum_{t=1,2,\dots,T} min_{1<= i <= m}|x_{t}-c_{i}|_{2}^{2}.}
#' @return The sum of squared Euclidean distance of each of T observations in matrix \code{observations} to its nearest center within \code{centers}.
#' @examples
#' ## generating 4 centers of length 3.
#' centers <- matrix(1:12, nrow = 4, ncol = 3)
#' ## generating 10 observations of length 3.
#' observations <- matrix(rmnorm(10, mean = c(0,3,5), varcov = diag(1,3)), nrow = 10)
#' cumulative_loss(centers, observations)
#' @keywords internal
#' @export





cumulative_loss = function(centers, observations){
  if (class(observations) != 'matrix'){
    return (instantaneous_loss(centers, observations))
  }else{
    if (class(centers) != 'matrix'){
      number_of_obs = length(observations[,1])
      centers_transform = t(replicate(number_of_obs, centers))
      distances = apply((centers_transform - observations)^2, 1, sum)
      return (sum(distances))

    }else{

      d = length(observations[1,])
      lth_centers=length(centers[,1])
      size_obs=length(observations[,1])
      C_mtx=replicate(size_obs,centers)
      Data_mtx=array(apply(observations,1,function(x) matrix(rep(x,lth_centers),nrow                           =lth_centers, byrow=T)),dim=c(lth_centers,d,size_obs))
      Diff=apply((C_mtx-Data_mtx)^2, c(1,3), sum)
      return(sum(apply(Diff,2,min)))
      }
  }
}


#' @title Label Caculation Function
#' @description This function gives cluster label to each of the observations.
#' @param centers a matrix containing m centers of length d, where each row corresponds to coordinates of a center.
#' @param observations a matrix containing T observations of length d, where each row of the matrix is an observation of length d.
#' @details Given a set \emph{C} of m centers of length d (\emph{i.e.,} \emph{C} = \eqn{{c_{1}, c_{2}, \dots, c_{m}}}), and a set \emph{X} of T observations of length d (\emph{i.e.,} \emph{X} = \eqn{{x_{1}, x_{2}, \dots, x_{T}}}), this function associates each observation \eqn{x_{t}} to the nearest centers to form different cells, and gives a corresponding cluster label to it.
#' @return cluster label of each of the observations.
#' @examples
#' ## generating 4 centers of length 3.
#' centers <- matrix(1:12, nrow = 4, ncol = 3)
#' ## generating 10 observations of length 3.
#' observations <- matrix(rmnorm(10, mean = c(2,6,10), varcov = diag(1,3)), nrow = 10)
#' labels_function(centers, observations)
#' @keywords internal
#' @export



labels_function = function(centers, observations){
  if (class(observations) != 'matrix'){
    if (class(centers) != 'matrix'){
      return (c(1))
    }else{
      replicate_instant_observation = t(replicate(length(centers[,1]), observations))
      distances = apply((centers - replicate_instant_observation)^2, 1, sum)
      return (which(distances == min(distances)))
    }
  }else{
    if (class(centers) != 'matrix'){
      number_of_obs = length(observations[,1])
      centers_transform = t(replicate(number_of_obs, centers))
      distances = apply((centers_transform - observations)^2, 1, sum)
      return (rep(1, number_of_obs))

    }else{
      d = length(observations[1,])
      lth_centers=length(centers[,1])
      size_obs=length(observations[,1])
      C_mtx=replicate(size_obs,centers)
      Data_mtx=array(apply(observations,1,function(x) matrix(rep(x,lth_centers),nrow                           =lth_centers, byrow=T)),dim=c(lth_centers,d,size_obs))
      Diff=apply((C_mtx-Data_mtx)^2, c(1,3), sum)
      return(apply(Diff, 2, function(x) which(x == min(x))))
    }
  }
}



