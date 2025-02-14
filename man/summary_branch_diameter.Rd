% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize.R
\name{summary_branch_diameter}
\alias{summary_branch_diameter}
\title{Summarize branches in diameter classes}
\usage{
summary_branch_diameter(qsm)
}
\arguments{
\item{qsm}{An object of class \code{QSM}.}
}
\value{
A \code{data.table} containing stats like volume (l), length (m) and surface
area (m^2).
}
\description{
\code{summary_branch_diameter} summarizes the branches into 1 cm diameter
classes.
}
\examples{
# load qsm
file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
qsm <- readQSM(file_path)

# summarize
summary_branch_diameter(qsm)
}
