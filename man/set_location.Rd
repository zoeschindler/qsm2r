% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{set_location}
\alias{set_location}
\title{Set location of the stem base}
\usage{
set_location(qsm, location = c(0, 0, 0))
}
\arguments{
\item{qsm}{An object of class \code{QSM}.}

\item{location}{\code{numeric}, \code{xyz}-coordinates of the new stem base.}
}
\value{
An object of class \code{QSM}.
}
\description{
\code{set_location} changes the location of the QSM so that \code{location}
is moved to the point \code{(0|0|0)}.
}
\examples{
# load qsm
file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
qsm <- readQSM(file_path)

# check old location
get_location(qsm)

# shift stem base to (0|0|0)
qsm <- set_location(qsm, location = c(0,0,0))

# check new location
get_location(qsm)
}
\seealso{
\code{\link{get_location}}
}
