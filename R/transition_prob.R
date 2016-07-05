#' @title State Transition Probability
#' @description This function computes the transition probability from current state x to the proposal state z belonging to the neighbor of x (\emph{i.e.,} \eqn{z \in {x-1, x, x+1}}).
#' @param x positive integer indicating the current state.
#' @param z positive integer indicating the proposal state.
#' @param K_max positive integer indicating the maximum value allowed for \code{x} and \code{z}.
#' @param prob probability of jumping from state \code{x} to state \code{x} if \code{x} doesn't equal to K_max nor 1.
#' @details This function computes the transition probability from current state x to the proposal state z. If x doesn't equal to \code{K_max} nor 1 and z belongs to the neighbor of x (\emph{i.e., {x-1,x,x+1}}), then the transition probability is \code{prob}, else 0. If x equals to 1 or K_max, then the transition probability from x to z belongs to {1,2} or {K_max-1, K_max} is 0.5, else 0.
#' @return A transition probability from state \code{x} to state \code{z}.
#' @examples
#' transition_probability(2, 3, 50, 1/3)
#' @keywords internal
#' @export


transition_probability = function(x, z, K_max, prob = 1/3){
  if((x > K_max)|(z > K_max)|(x < 1)|(z < 1)){
    q = 0
  }else{ if(x == 1){
    q = 0.5*((z == 1)|(z == 2))
  }else{ if(x == K_max){
    q = 0.5*((z == K_max)|(z == (K_max-1)))
  }else{if(x == z){
    q = prob
  }else{q = (1-prob)/2 * (( x < z & z < (x+2))|(z < x & (z > (x-2))))}

  }
  }
  }
  return(q)
}
