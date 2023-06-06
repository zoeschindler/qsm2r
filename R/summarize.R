################################################################################
# MAIN FUNCTIONS
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

summarize_branches <- function(qsm, n_class, diff_class, var_class, idx_all, idx_1st = NULL) {

  # prepare storage
  out <- c()

  # loop through classes
  for (n in 1:n_class) {

    # get indices of subset
    class_start <- (n - 1)*diff_class
    class_end <- n*diff_class
    idx_all_class <- var_class >= class_start & var_class < class_end & idx_all
    idx_1st_class <- var_class >= class_start & var_class < class_end & idx_1st

    # derive & save stats
    curr <- data.table::data.table(
      "volume_all_l" = sum(qsm@branch$volume[idx_all_class]), # in l
      "volume_1st_l" = sum(qsm@branch$volume[idx_1st_class]), # in l
      "area_all_m2" = sum(qsm@branch$area[idx_all_class]), # in m2
      "area_1st_m2" = sum(qsm@branch$area[idx_1st_class]), # in m2
      "length_all_m" = sum(qsm@branch$length[idx_all_class]), # in m
      "length_1st_m" = sum(qsm@branch$length[idx_1st_class]), # in m
      "count_all_branch" = nrow(qsm@branch[idx_all_class,]),
      "count_1st_branch" = nrow(qsm@branch[idx_1st_class,]))
    out <- rbind(out, curr)
  }

  # return results
  return(out)
}

################################################################################

summarize_cylinders <- function(qsm, n_class, diff_class, var_class) {

  # prepare storage
  out <- c()

  # loop through classes
  for (n in 1:n_class) {

    # get indices of subset
    class_start <- (n - 1)*diff_class
    class_end <- n*diff_class
    idx <- var_class >= class_start & var_class < class_end

    # derive & save stats
    curr <- data.table::data.table(
      "volume_l" = 1000*sum(pi*qsm@cylinder$length[idx]*qsm@cylinder$radius[idx]**2),
      "area_m2" = sum(2*pi*qsm@cylinder$radius[idx]*qsm@cylinder$length[idx]),
      "length_m" = sum(qsm@cylinder$length[idx]),
      "count_cylinder" = nrow(qsm@cylinder[idx,]))
    out <- rbind(out, curr)
  }

  # return results
  return(out)
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
