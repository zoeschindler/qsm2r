% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize.R
\name{summary_cylinder_azimuth}
\alias{summary_cylinder_azimuth}
\title{Summarize cylinders in azimuth classes}
\usage{
summary_cylinder_azimuth(qsm)
}
\arguments{
\item{qsm}{An object of class \code{QSM}.}
}
\value{
A \code{data.table} containing stats like volume (l), length (m) and surface
area (m^2).
}
\description{
\code{summary_cylinder_azimuth} summarizes the cylinders into 10 degree
azimuth classes.
}
\examples{
# load qsm
file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
qsm <- readQSM(file_path)

# summarize
summary_cylinder_azimuth(qsm)
}
