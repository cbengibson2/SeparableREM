#' Militarized Interstate Dispute conflict edgelist
#'
#' A directed event sequence of militarized conflict acts between sovereign
#' states, drawn from the Militarized Interstate Dispute (MID) dataset.
#' Each row records one country using military force against another.
#'
#' @format A data frame with 1,708 rows and 3 variables:
#' \describe{
#'   \item{ev}{Integer event index (1–1708).}
#'   \item{id_snd}{Integer sender country ID.}
#'   \item{id_rec}{Integer receiver country ID.}
#' }
#'
#' @source Militarized Interstate Dispute (MID) dataset.
"edgelist"
