#' Observation Indexing and Counting in Defined Partitions
#'
#' These functions identify and count the number of observations within
#' specified intervals defined by the change points in `p`. The intervals are
#' assumed to be left-closed and right-open.
#'
#' @param x A numeric vector of observations.
#' @param p A numeric vector specifying the minimum value \eqn{\tau_{(0)}} and
#' change points \eqn{\tau_{(1)}, \ldots, \tau_{(k)}} that define the partitions,
#' i.e., \eqn{p = (\tau_{(0)}, \tau_{(1)}, \ldots, \tau_{(k)})}. It must be sorted
#' in increasing order.
#'
#' @return
#' - `index_each_interval()`: A list of integer vectors, where each element
#'   contains the indices of observations belonging to the corresponding interval.
#' - `n_each_interval()`: A numeric vector indicating the count of observations
#'   in each interval.
#'
#' @details The intervals follow the structure:
#' \deqn{\mathcal{R}_1 = [\tau_{(0)}, \tau_{(1)}), \mathcal{R}_2 = [\tau_{(1)}, \tau_{(2)}),
#' ..., \mathcal{R}_{k} = [\tau_{(k-1)}, \tau_{(k)})}
#' except for the last partition, which is closed: \eqn{\mathcal{R}_{k+1} = [\tau_{(k)}, \infty)}.
#'
#' @examples
#' x <- c(7, 7, 6, 4, 1, 5, 2, 6, 3, 5)
#'
#' # Get indices of observations within each interval
#' index_each_interval(x, c(1, 3, 5))
#'
#' # Count observations in each interval
#' n_each_interval(x, c(1, 3, 5))
#'
#' @export
index_each_interval <- function (x, p)
{
  k <- length(p)
  nj <- vector("list", length = k)
  for (j in 1:k)
  {
    if (j < k)
    {
      nj[[j]] <- which(x >= p[j] & x < p[j + 1])
    } else {
      nj[[j]] <- which(x >= p[j])
    }
  }

  nj
}

#' @rdname index_each_interval
#' @export
n_each_interval <- function (x, p)
{
  lengths(index_each_interval(x, p))
}
