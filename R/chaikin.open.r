#' Chaikin corner cutting algorithm for open segments
#' Used recursively in a loop
#' @param x Original line segment
#' @return A matrix of points
#' @export 

chaikin.open = function(x)
# https://gist.github.com/mstrimas/ac50a38a7e656a2b3a173f3a6b31a760
# https://observablehq.com/@pamacha/chaikins-algorithm
# https://observablehq.com/@infowantstobeseen/chaikins-curves
{
  n_pts = nrow(x)
  qr = matrix(NA, nrow = 2 * (n_pts - 1), ncol = 2)
  qr[seq(1, nrow(qr), by = 2), ] <- 0.75 * x[-n_pts, ] + 0.25 * x[-1, ]
  qr[seq(2, nrow(qr), by = 2), ] <- 0.75 * x[-1, ] + 0.25 * x[-n_pts, ]
  return(qr)
}


