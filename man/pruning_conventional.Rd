% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pruning.R
\name{pruning_conventional}
\alias{pruning_conventional}
\title{Conventional pruning at specified height or stem length}
\usage{
pruning_conventional(
  qsm,
  threshold_m = 3,
  method = c("length", "height"),
  remove = FALSE
)
}
\arguments{
\item{qsm}{An object of class \code{QSM}.}

\item{threshold_m}{\code{numeric}, tree height or stem length in meters below
which all first order branches should be removed.}

\item{method}{\code{character}, whether the threshold refers to tree height
(\code{"height"}) or stem length (\code{"length"}).}

\item{remove}{\code{boolean}, whether the to be pruned cylinders should be
removed (\code{TRUE}) or labelled (\code{FALSE}).}
}
\value{
\code{QSM}, with removed or labelled cylinders (column \code{cylinder$prune}).
}
\description{
\code{pruning_conventional} prunes all first order and their subsequent
branches below the specified tree height or stem length.
}
\examples{
# load qsm
file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
qsm <- readQSM(file_path)

# prune specified branches, without removing
not_removed <- pruning_conventional(qsm, threshold_m = 3, method = "length", remove = FALSE)

# plot qsm
plot(not_removed, col_var = "prune")

# prune specified branches, with removing
removed <- pruning_conventional(qsm, threshold_m = 3, method = "length", remove = TRUE)

# plot qsm
plot(removed)
}
\seealso{
\code{\link{pruning_selective}}, \code{\link{pruning_whorlwise}},
\code{\link{pruning_tips}}
}
