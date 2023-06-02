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
#' @export
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
#' @export
set_location <- function(qsm, location) {

  # recenter cylinders from location to (0|0|0)
  qsm@cylinder$start_X <- qsm@cylinder$start_X - location$x
  qsm@cylinder$start_Y <- qsm@cylinder$start_Y - location$y
  qsm@cylinder$start_Z <- qsm@cylinder$start_Z - location$z

  # return results
  return(qsm)
}

################################################################################

summary_stemtaper <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################

summary_spreads <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################

#' Summarize cylinders in diameter classes
#'
#' @description
#' \code{summary_cylinder_diameter} summarizes the cylinders into 1 cm diameter
#' classes.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' A \code{data.table} containing stats like volume (l), length (m) and surface
#' area (m^2).
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # summarize
#' summary_cylinder_diameter(qsm)
#' @export
summary_cylinder_diameter <- function(qsm) {

  # prepare classes
  var_class <- qsm@cylinder$radius
  n_class <- ceiling(max(200 * var_class))
  diff_class <- 0.005 # 1 cm diameter = 0.005 m radius classes

  # derive & save stats
  stats <- summarize_cylinders(qsm, n_class, diff_class, var_class)
  classes <- data.table::data.table(
    "diam_start_cm" = ((1:n_class) - 1)*diff_class*200,
    "diam_end_cm" = (1:n_class)*diff_class*200)
  out <- cbind(classes, stats)

  # return results
  return(out)
}

################################################################################

#' Summarize cylinders in height classes
#'
#' @description
#' \code{summary_cylinder_height} summarizes the cylinders into 1 m height
#' classes.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' A \code{data.table} containing stats like volume (l), length (m) and surface
#' area (m^2).
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # summarize
#' summary_cylinder_height(qsm)
#' @export
summary_cylinder_height <- function(qsm) {

  # make subsets
  radius <- qsm@cylinder$radius
  length <- qsm@cylinder$length

  # get number of classes
  max_height <- ceiling(qsm@overview$TreeHeight)

  # prepare data
  base <- qsm@cylinder$start_Z[qsm@cylinder$BranchOrder == 0 & qsm@cylinder$PositionInBranch == 1]
  start <- qsm@cylinder$start_Z - base # B
  end <- qsm@cylinder$start_Z + qsm@cylinder$length * qsm@cylinder$axis_Z - base # T

  # prepare storage
  out <- data.table::data.table()

  # loop through classes
  for (class in 1:max_height) {

    # prepare storage
    curr <- data.table::data.table(
      "height_start_m" = class - 1,
      "height_end_m" = class,
      "volume_l" = 0,
      "area_m2" = 0,
      "length_m" = 0)

    # get indices of different cases
    sb <- start >= (class - 2) & start < (class - 1)  # start below
    si <- start >= (class - 1) & start < class        # start inside
    sa <- start >= class & start < (class + 1)        # start above
    eb <- end >= (class - 2) & end < (class - 1)      # end below
    ei <- end >= (class - 1) & end < class            # end inside
    ea <- end >= class & end < (class + 1)            # end above

    # case 1: start and end inside
    idx <- si & ei
    curr$volume_l <- curr$volume_l + 1000*sum(pi*length[idx]*radius[idx]**2)
    curr$area_m2  <- curr$area_m2 + sum(2*pi*radius[idx]*length[idx])
    curr$length_m <- curr$length_m + sum(length[idx])

    # case 2: start inside, end above
    idx <- si & ea
    rel_length <- (class - start[idx])/(end[idx] - start[idx])
    curr$volume_l <- curr$volume_l + 1000*sum(pi*rel_length*length[idx]*radius[idx]**2)
    curr$area_m2  <- curr$area_m2 + sum(2*pi*radius[idx]*length[idx]*rel_length)
    curr$length_m <- curr$length_m + sum(length[idx] * rel_length)

    # case 3: start inside, end below
    idx <- si & eb
    rel_length <- (start[idx] - class + 1)/(start[idx] - end[idx])
    curr$volume_l <- curr$volume_l + 1000*sum(pi*rel_length*length[idx]*radius[idx]**2)
    curr$area_m2  <- curr$area_m2 + sum(2*pi*radius[idx]*length[idx]*rel_length)
    curr$length_m <- curr$length_m + sum(length[idx] * rel_length)

    # case 4: start below, end inside
    idx <- sb & ei
    rel_length <- (end[idx] - class + 1)/(end[idx] - start[idx])
    curr$volume_l <- curr$volume_l + 1000*sum(pi*rel_length*length[idx]*radius[idx]**2)
    curr$area_m2  <- curr$area_m2 + sum(2*pi*radius[idx]*length[idx]*rel_length)
    curr$length_m <- curr$length_m + sum(length[idx] * rel_length)

    # case 5: start above, end inside
    idx <- sa & ei
    rel_length <- (class - end[idx])/(start[idx] - end[idx])
    curr$volume_l <- curr$volume_l + 1000*sum(pi*rel_length*length[idx]*radius[idx]**2)
    curr$area_m2  <- curr$area_m2 + sum(2*pi*radius[idx]*length[idx]*rel_length)
    curr$length_m <- curr$length_m + sum(length[idx] * rel_length)

    # save stats
    out <- rbind(out, curr)
  }

  # return results
  return(out)
}

################################################################################

#' Summarize cylinders in zenith classes
#'
#' @description
#' \code{summary_cylinder_zenith} summarizes the cylinders into 10 degree zenith
#' classes.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' A \code{data.table} containing stats like volume (l), length (m) and surface
#' area (m^2).
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # summarize
#' summary_cylinder_zenith(qsm)
#' @export
summary_cylinder_zenith <- function(qsm) {

  # prepare classes
  var_class <- 180/pi*acos(qsm@cylinder$axis_Z)
  n_class <- 18
  diff_class <- 10 # 10 degree classes

  # derive & save stats
  stats <- summarize_cylinders(qsm, n_class, diff_class, var_class)
  classes <- data.table::data.table(
    "zenith_start_deg" = as.integer(((1:n_class) - 1) * diff_class),
    "zenith_end_deg" = as.integer((1:n_class) * diff_class))
  out <- cbind(classes, stats)

  # return results
  return(out)
}

################################################################################

#' Summarize cylinders in azimuth classes
#'
#' @description
#' \code{summary_cylinder_azimuth} summarizes the cylinders into 10 degree
#' azimuth classes.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' A \code{data.table} containing stats like volume (l), length (m) and surface
#' area (m^2).
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # summarize
#' summary_cylinder_azimuth(qsm)
#' @export
summary_cylinder_azimuth <- function(qsm) {

  # prepare classes
  var_class <- 180/pi*atan2(qsm@cylinder$axis_Y,qsm@cylinder$axis_X) + 180
  n_class <- 36
  diff_class <- 10 # 10 degree classes

  # derive & save stats
  stats <- summarize_cylinders(qsm, n_class, diff_class, var_class)
  classes <- data.table::data.table(
    "azimuth_start_deg" = as.integer(((1:n_class) - 1) * diff_class),
    "azimuth_end_deg" = as.integer((1:n_class) * diff_class))
  out <- cbind(classes, stats)

  # return results
  return(out)
}

################################################################################

#' Summarize branches in branch orders
#'
#' @description
#' \code{summary_branch_order} summarizes the branches into branch order
#' classes.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' A \code{data.table} containing stats like volume (l), length (m) and surface
#' area (m^2).
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # summarize
#' summary_branch_order(qsm)
#' @export
summary_branch_order <- function(qsm) {

  # prepare storage
  out <- c()

  # loop through branch orders
  for (n in 1:max(qsm@branch$order)) {

    # get indices of subset
    idx <- qsm@branch$order == n

    # derive & save stats
    curr <- data.table::data.table(
      "branch_order" = n,
      "volume_l" = sum(qsm@branch$volume[idx]),
      "area_m2" = sum(qsm@branch$area[idx]),
      "length_m" = sum(qsm@branch$length[idx]),
      "count_branch" = nrow(qsm@branch[idx,]))
    out <- rbind(out, curr)
  }

  # return results
  return(out)
}

################################################################################

#' Summarize branches in angle classes
#'
#' @description
#' \code{summary_branch_angle} summarizes the branches into 10 degree angle
#' classes.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' A \code{data.table} containing stats like volume (l), length (m) and surface
#' area (m^2).
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # summarize
#' summary_branch_angle(qsm)
#' @export
summary_branch_angle <- function(qsm) {

  # prepare classes
  idx_all <- qsm@branch$order > 0
  idx_1st <- qsm@branch$order == 1
  var_class <- qsm@branch$angle
  n_class <- 18
  diff_class <- 10 # 10 degree classes

  # derive & save stats
  stats <- summarize_branches(qsm, n_class, diff_class, var_class, idx_all, idx_1st)
  classes <- data.table::data.table(
    "angle_start_deg" = ((1:n_class) - 1)*diff_class,
    "angle_end_deg" = (1:n_class)*diff_class)
  out <- cbind(classes, stats)

  # return results
  return(out)
}

################################################################################


#' Summarize branches in diameter classes
#'
#' @description
#' \code{summary_branch_diameter} summarizes the branches into 1 cm diameter
#' classes.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' A \code{data.table} containing stats like volume (l), length (m) and surface
#' area (m^2).
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # summarize
#' summary_branch_diameter(qsm)
#' @export
summary_branch_diameter <- function(qsm) {

  # prepare classes
  idx_all <- qsm@branch$order > 0
  idx_1st <- qsm@branch$order == 1
  var_class <- qsm@branch$diameter
  n_class <- ceiling(max(100 * var_class[idx_all]))
  diff_class <- 0.01 # 1 cm diameter classes

  # derive & save stats
  stats <- summarize_branches(qsm, n_class, diff_class, var_class, idx_all, idx_1st)
  classes <- data.table::data.table(
    "diam_start_cm" = ((1:n_class) - 1)*diff_class*100,
    "diam_end_cm" = (1:n_class)*diff_class*100)
  out <- cbind(classes, stats)

  # return results
  return(out)
}

################################################################################

#' Summarize branches in height classes
#'
#' @description
#' \code{summary_branch_height} summarizes the branches into 1 m height
#' classes.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' A \code{data.table} containing stats like volume (l), length (m) and surface
#' area (m^2).
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # summarize
#' summary_branch_height(qsm)
#' @export
summary_branch_height <- function(qsm) {

  # prepare classes
  idx_all <- qsm@branch$order > 0
  idx_1st <- qsm@branch$order == 1
  var_class <- qsm@branch$height
  n_class <- ceiling(qsm@overview$TreeHeight)
  diff_class <- 1 # 1 m height classes

  # derive & save stats
  stats <- summarize_branches(qsm, n_class, diff_class, var_class, idx_all, idx_1st)
  classes <- data.table::data.table(
    "height_start_m" = ((1:n_class) - 1)*diff_class,
    "height_end_m" = (1:n_class)*diff_class)
  out <- cbind(classes, stats)

  # return results
  return(out)
}

################################################################################

#' Summarize branches in zenith classes
#'
#' @description
#' \code{summary_branch_zenith} summarizes the branches into 10 degree zenith
#' classes.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' A \code{data.table} containing stats like volume (l), length (m) and surface
#' area (m^2).
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # summarize
#' summary_branch_zenith(qsm)
#' @export
summary_branch_zenith <- function(qsm) {

  # prepare classes
  idx_all <- qsm@branch$order > 0
  idx_1st <- qsm@branch$order == 1
  var_class <- qsm@branch$zenith
  n_class <- 18
  diff_class <- 10 # 10 degree classes

  # derive & save stats
  stats <- summarize_branches(qsm, n_class, diff_class, var_class, idx_all, idx_1st)
  classes <- data.table::data.table(
    "zenith_start_deg" = ((1:n_class) - 1)*diff_class,
    "zenith_end_deg" = (1:n_class)*diff_class)
  out <- cbind(classes, stats)

  # return results
  return(out)
}

################################################################################

#' Summarize branches in azimuth classes
#'
#' @description
#' \code{summary_branch_azimuth} summarizes the branches into 10 degree azimuth
#' classes.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' A \code{data.table} containing stats like volume (l), length (m) and surface
#' area (m^2).
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # summarize
#' summary_branch_azimuth(qsm)
#' @export
summary_branch_azimuth <- function(qsm) {

  # prepare classes
  idx_all <- qsm@branch$order > 0
  idx_1st <- qsm@branch$order == 1
  var_class <- qsm@branch$azimuth + 180
  n_class <- 36
  diff_class <- 10 # 10 degree classes

  # derive & save stats
  stats <- summarize_branches(qsm, n_class, diff_class, var_class, idx_all, idx_1st)
  classes <- data.table::data.table(
    "azimuth_start_deg" = ((1:n_class) - 1)*diff_class,
    "azimuth_end_deg" = (1:n_class)*diff_class)
  out <- cbind(classes, stats)

  # return results
  return(out)
}

################################################################################
