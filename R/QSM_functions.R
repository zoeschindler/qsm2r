################################################################################
# MAIN FUNCTIONS
################################################################################

writeQSM <- function(qsm, file_path) {
  # TODO: which datatype? rds?
  # TODO: change readQSM depending on output datatype?
  return(invisible(file_path))
}

################################################################################

#' Get location of the stem base
#'
#' @description
#' \code{get_location} obtains the location of the stem base.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' A list containing the \code{xyz}-coordinates of the stem base.
#'
#' @seealso \code{\link{set_location}}
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # get stem base location
#' stem_base <- get_location(qsm)
get_location <- function(qsm) {

  # get starting coordinates of stem base
  base_id <- qsm@cylinder$cyl_id[qsm@cylinder$BranchOrder == 0 & qsm@cylinder$PositionInBranch == 1]
  location <- list(
    "x" = qsm@cylinder$start_X[qsm@cylinder$cyl_id == base_id],
    "y" = qsm@cylinder$start_Y[qsm@cylinder$cyl_id == base_id],
    "z" = qsm@cylinder$start_Z[qsm@cylinder$cyl_id == base_id])

  # show results
  cat("x:", round(location[["x"]], 2), "\n")
  cat("y:", round(location[["y"]], 2), "\n")
  cat("z:", round(location[["z"]], 2), "\n")

  # return results
  return(invisible(location))
}

################################################################################

#' Set location of the stem base
#'
#' @description
#' \code{set_location} changes the location of the QSM so that \code{location}
#' is moved to the point \code{(0|0|0)}.
#'
#' @param qsm An object of class \code{QSM}.
#' @param location Coordinates of a point which is supposed to be moved to
#' \code{(0|0|0)}. A list or a data.frame with items or columns \code{xyz}.
#'
#' @return
#' An object of class \code{QSM}.
#'
#' @seealso \code{\link{get_location}}
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # get stem base location
#' stem_base <- get_location(qsm)
#'
#' # shift stem base to (0|0|0)
#' qsm <- set_location(qsm, stem_base)
set_location <- function(qsm, location) {

  # recenter cylinders from location to (0|0|0)
  qsm@cylinder$start_X <- qsm@cylinder$start_X - location$x
  qsm@cylinder$start_Y <- qsm@cylinder$start_Y - location$y
  qsm@cylinder$start_Z <- qsm@cylinder$start_Z - location$z

  # return results
  return(qsm)
}

################################################################################

get_stemtaper <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################

summary_branchorder <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################

summary_cylinder_diameter <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################

summary_cylinder_height <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################

summary_cylinder_zenith <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################

summary_cylinder_azimuth <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################

summary_branch_diameter <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################

summary_branch_height <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################

summary_branch_zenith <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################

summary_branch_azimuth <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################
