################################################################################
# MAIN FUNCTIONS
################################################################################

updateQSM_basics <- function(qsm) {

  # prepare data
  overview <- qsm@overview

  # create subsets
  trunk_cyl  <- qsm@cylinder[qsm@cylinder$branch == 1,]
  branch_cyl <- qsm@cylinder[qsm@cylinder$branch > 1,]

  # derive basic stats
  overview$TotalVolume    <- sum(pi * qsm@cylinder$length * qsm@cylinder$radius ** 2) * 1000
  overview$TrunkVolume    <- sum(pi * trunk_cyl$length * trunk_cyl$radius ** 2) * 1000
  overview$BranchVolume   <- sum(pi * branch_cyl$length * branch_cyl$radius ** 2) * 1000
  #
  overview$TreeHeight     <- max(c(
    qsm@cylinder$start_Z, # base
    qsm@cylinder$start_Z + qsm@cylinder$axis_Z * qsm@cylinder$length)) - # tip
    qsm@cylinder$start_Z[qsm@cylinder$BranchOrder == 0 & qsm@cylinder$PositionInBranch == 1]
  #
  overview$TrunkLength    <- sum(trunk_cyl$length)
  overview$BranchLength   <- sum(branch_cyl$length)
  overview$TotalLength    <- sum(qsm@cylinder$length)
  #
  overview$NumberBranches <- length(unique(qsm@cylinder$branch)) - 1
  overview$MaxBranchOrder <- max(qsm@cylinder$BranchOrder)
  #
  overview$TrunkArea      <- sum(2 * pi * trunk_cyl$radius * trunk_cyl$length)
  overview$BranchArea     <- sum(2 * pi * branch_cyl$radius * branch_cyl$length)
  overview$TotalArea      <- sum(2 * pi * qsm@cylinder$radius * qsm@cylinder$length)
  #
  overview$DBHqsm         <- trunk_cyl$radius[min(which(cumsum(trunk_cyl$length) > 1.3))] * 2
  # overview$DBHcyl is calculated from point cloud and can't be reproduced

  # overwrite old data
  qsm@overview <- overview

  # return results
  return(qsm)
}

################################################################################

updateQSM_crown <- function(qsm, method = "TreeQSM") {

  # requires up-to-date branch and basic overview data

  # prepare output
  overview <- qsm@overview

  # check method argument
  if (!method %in% c("TreeQSM", "qsm2r")) {
    stop("method must be either 'TreeQSM' or 'qsm2r'")
  }

  ### PART 1: create point cloud from cylinders

  # prepare data
  cylinder <- qsm@cylinder

  # create empty columns
  cylinder$orth_1_X <- cylinder$orth_1_Y <- cylinder$orth_1_Z <- NA

  # get unit vectors orthogonal to axis
  cylinder[,c("orth_1_X", "orth_1_Y", "orth_1_Z")] <- random_orth_norm(
    cylinder$axis_X, cylinder$axis_Y, cylinder$axis_Z)

  # create subsets
  trunk_cyl  <- cylinder[cylinder$branch == 1,]
  branch_cyl <- cylinder[cylinder$branch > 1,]

  # duplicate rows accordingly
  trunk_cyl  <- trunk_cyl[rep(1:nrow(trunk_cyl), each = 4 * 12),]
  branch_cyl <- branch_cyl[rep(1:nrow(branch_cyl), each = 4),]

  # add column length_mod
  trunk_cyl$length_mod  <- rep(rep(seq(0.25, 1, 0.25), each = 12), length.out = nrow(trunk_cyl))
  branch_cyl$length_mod <- 0.5

  # add column angle_mod
  trunk_cyl$angle_mod  <- rep(seq(0, 2*pi - 2*pi/12, length.out = 12), length.out = nrow(trunk_cyl))
  branch_cyl$angle_mod <- rep(seq(0, 2*pi - 2*pi/4,  length.out = 4), length.out = nrow(branch_cyl))

  # combine data frames
  cylinder <- rbind(trunk_cyl, branch_cyl)
  rm(trunk_cyl, branch_cyl)

  # get circle centers along cylinder axis
  cylinder$center_X <- cylinder$start_X + cylinder$length_mod * cylinder$length * cylinder$axis_X
  cylinder$center_Y <- cylinder$start_Y + cylinder$length_mod * cylinder$length * cylinder$axis_Y
  cylinder$center_Z <- cylinder$start_Z + cylinder$length_mod * cylinder$length * cylinder$axis_Z

  # create empty columns
  cylinder$orth_2_X <- cylinder$orth_2_Y <- cylinder$orth_2_Z <- NA

  # get unit vector orthogonal to orth_1 and axis
  cylinder[,c("orth_2_X", "orth_2_Y", "orth_2_Z")] <- two_vector_orth_norm(
    cylinder$axis_X, cylinder$axis_Y, cylinder$axis_Z,
    cylinder$orth_1_X, cylinder$orth_1_Y, cylinder$orth_1_Z)

  # prepare storage
  points <- data.table::data.table()

  # derive surface points
  points$X <- cylinder$center_X +
    cylinder$radius * cylinder$orth_1_X * cos(cylinder$angle_mod) +
    cylinder$radius * cylinder$orth_2_X * sin(cylinder$angle_mod)
  points$Y <- cylinder$center_Y +
    cylinder$radius * cylinder$orth_1_Y * cos(cylinder$angle_mod) +
    cylinder$radius * cylinder$orth_2_Y * sin(cylinder$angle_mod)
  points$Z <- cylinder$center_Z +
    cylinder$radius * cylinder$orth_1_Z * cos(cylinder$angle_mod) +
    cylinder$radius * cylinder$orth_2_Z * sin(cylinder$angle_mod)

  # add tips and starts
  points <- rbind(
    points,
    data.table::data.table(
      "X" = cylinder$start_X,
      "Y" = cylinder$start_Y,
      "Z" = cylinder$start_Z),
    data.table::data.table(
      "X" = cylinder$start_X + cylinder$length * cylinder$axis_X,
      "Y" = cylinder$start_Y + cylinder$length * cylinder$axis_Y,
      "Z" = cylinder$start_Z + cylinder$length * cylinder$axis_Z))

  # delete duplicates
  points <- unique(round(points, 6))

  ### PART 2: vertical profiles

  # I don't care about spreads (yet)

  ### PART 3: crown area & diameter

  # prepare data
  cylinder <- qsm@cylinder

  # get coordinates of convex hull
  chull_idx <- grDevices::chull(points[,c("X", "Y")])
  chull_idx <- c(chull_idx, chull_idx[1])
  chull_n <- length(chull_idx)
  chull_x <- points$X[chull_idx]
  chull_y <- points$Y[chull_idx]

  # get centroid of the polygon
  chull_area <- polygon_area(chull_x, chull_y, chull_n)
  chull_centroid <- polygon_centroid(chull_x, chull_y, chull_n, chull_area)

  # choosing method
  if (method == "TreeQSM") {

    # calculate vector from center to cylinder tips
    vec_to_tips <- data.table::data.table(
      "vec_X" = (cylinder$start_X + cylinder$length * cylinder$axis_X) - chull_centroid$X,
      "vec_Y" = (cylinder$start_Y + cylinder$length * cylinder$axis_Y) - chull_centroid$Y)

    # sort according to angle
    vec_to_tips$angle <- atan2(vec_to_tips$vec_Y, vec_to_tips$vec_X) + pi
    vec_to_tips <- vec_to_tips[order(vec_to_tips$angle),]

    # get lengths of the vectors
    vec_to_tips$length <- sqrt(vec_to_tips$vec_X**2 +  vec_to_tips$vec_Y**2) # TODO: check this

    # prepare storage
    crown_widths <- c()

    # loop through 10 degree steps
    for (section in 1:18) {

      # get length in one direction
      idx_1 <- which(vec_to_tips$angle >= (section - 1)*pi/18 & vec_to_tips$angle < section*pi/18)
      len_1 <- ifelse(length(idx_1) > 0, max(vec_to_tips$length[idx_1]), 0)

      # get length in opposite direction
      idx_2 <- which(vec_to_tips$angle >= ((section - 1)*pi/18 + pi) & vec_to_tips$angle < (section*pi/18 + pi))
      len_2 <- ifelse(length(idx_2) > 0, max(vec_to_tips$length[idx_2]), 0)

      # add crown width
      crown_widths <- c(crown_widths, len_1 + len_2)
    }

    # derive average crown width from section crown widths
    avg_crown_width <- mean(crown_widths)

    # choosing method
  } else if (method == "qsm2r") {

    # Pretzsch? Forest Dynamics, Growth and Yield?
    # derive average crown width from crown projection area
    avg_crown_width <- 2 * sqrt(abs(chull_area) / pi)
  }

  # calculate maximum spread in the data
  max_crown_width <- max(dist(cbind(chull_x, chull_y)))

  # derive points of alpha shape
  alpha_value <- max(0.5, avg_crown_width/10)
  alpha_points <- unique(points[,c("X","Y")])
  alpha_shape <- alphahull::ashape(x = alpha_points$X, y = alpha_points$Y, alpha = alpha_value) # takes ages
  alpha_idx <- c(alpha_shape$alpha.extremes, alpha_shape$alpha.extremes[1]) # points are unordered

  # derive area of alpha shape
  # https://stackoverflow.com/questions/74054828/unordered-lines-to-closed-polygon
  points_unsorted <- alpha_points[alpha_idx,]
  points_unsorted$X <- points_unsorted$X - chull_centroid$X
  points_unsorted$Y <- points_unsorted$Y - chull_centroid$Y
  points_unsorted$angle <- atan2(points_unsorted$X, points_unsorted$Y)
  points_unsorted <- unique(points_unsorted[order(points_unsorted$angle),])
  alpha_x <- c(points_unsorted$X, points_unsorted$X[1])
  alpha_y <- c(points_unsorted$Y, points_unsorted$Y[1])
  alpha_n <- nrow(points_unsorted) + 1
  alpha_area <- abs(polygon_area(alpha_x, alpha_y, alpha_n))

  ### PART 4: crown base, length, ratio & volume

  # first major branch:
  # - lowest branch whose diameter > min(0.05 * dbh, 5cm)
  # - horizontal relative reach > median reach of 1st-ord.branches (or maximum 10)

  # prepare data
  branch <- qsm@branch
  cylinder <- qsm@cylinder

  # get cylinder tips
  cylinder$end_X <- cylinder$start_X + cylinder$length * cylinder$axis_X
  cylinder$end_Y <- cylinder$start_Y + cylinder$length * cylinder$axis_Y
  cylinder$end_Z <- cylinder$start_Z + cylinder$length * cylinder$axis_Z

  # get first order branches
  branch_1_id <- branch$bra_id[branch$order == 1] # IDs of first order branches
  branch_1 <- branch[branch$order == 1,]
  branch_1$rel_reach <- 0

  # loop through first branches
  for (bra_id in branch_1_id) {
    cyl_sub <- cylinder[cylinder$branch == bra_id,]
    if (nrow(cyl_sub) > 0) {

      # calculate relative reach
      first <- cyl_sub[1,]
      last <- cyl_sub[nrow(cyl_sub),]
      reach <- sqrt((first$start_X - last$end_X)**2 + (first$start_Y - last$end_Y)**2)
      rel_reach <- reach / qsm@overview$DBHqsm * 2
      branch_1$rel_reach[branch_1$bra_id == bra_id] <- rel_reach
    }
  }

  # derive threshold
  thresh_reach <- min(10, median(branch_1$rel_reach))
  thresh_diam  <- min(0.05, 0.05 * qsm@overview$DBHqsm)

  # find lowest branch meeting the demands (reach & diameter)
  branch_1$thresh <- branch_1$diameter > thresh_diam & branch_1$rel_reach > thresh_reach
  branch_1 <- branch_1[branch_1$thresh,]
  branch_selected <- branch_1[which.min(branch_1$height),]

  # there is a fitting branch
  if (nrow(branch_selected) > 0) {

    # get child cylinders of the selected branch
    child_bra_id <- find_childs_recursive_branch(cylinder, branch_selected$bra_id, TRUE)
    child_cyl <- cylinder[cylinder$branch %in% child_bra_id,] # B

    # get crown base height
    base_height <- min(child_cyl[,c("start_Z","end_Z")])
    crown_base_height <- base_height -
      cylinder$start_Z[cylinder$BranchOrder == 0 & cylinder$PositionInBranch == 1]

    # get crown length and ratio
    crown_length <- qsm@overview$TreeHeight - crown_base_height
    crown_ratio  <- crown_length / qsm@overview$TreeHeight

    # get points above the crown base
    points_crown_base <- as.matrix(unique(points[points$Z >= base_height,]))

    # crown volume from convex hull (3d)
    chull_3d <- cxhull::cxhull(points_crown_base)
    crown_volume_conv <- chull_3d$volume

    # crown volume from alpha shape (3d)
    alpha_value_3d <- max(0.5, avg_crown_width/5)
    alpha_shape_3d <- alphashape3d::ashape3d(x = points_crown_base, alpha = alpha_value_3d) # TODO: prevent holes in the middle?
    crown_volume_alpha <- alphashape3d::volume_ashape3d(alpha_shape_3d)

  # there is no fitting branch
  } else {
    message("tree has no crown")
    crown_base_height <- overview$TreeHeight
    crown_length <- 0
    crown_ratio <- 0
    crown_volume_conv <- 0
    crown_volume_alpha <- 0
  }

  ### PART 5: updating data

  # overwrite old data
  overview$CrownAreaConv    <- abs(chull_area)
  overview$CrownAreaAlpha   <- abs(alpha_area)
  overview$CrownDiamAve     <- avg_crown_width
  overview$CrownDiamMax     <- max_crown_width
  overview$CrownBaseHeight  <- crown_base_height
  overview$CrownLength      <- crown_length
  overview$CrownRatio       <- crown_ratio
  overview$CrownVolumeConv  <- crown_volume_conv
  overview$CrownVolumeAlpha <- crown_volume_alpha

  # overwrite old data
  qsm@overview <- overview

  # return results
  return(qsm)
}

################################################################################

updateQSM_branch <- function(qsm) {

  # prepare storage
  branch <- data.table::data.table()

  # create subsets
  radius <- data.table::data.table(
    "cyl_id" = qsm@cylinder$cyl_id,
    "radius" = qsm@cylinder$radius)
  length <- data.table::data.table(
    "cyl_id" = qsm@cylinder$cyl_id,
    "length" = qsm@cylinder$length)
  axis <- data.table(
    "cyl_id" = qsm@cylinder$cyl_id,
    "axis_X" = qsm@cylinder$axis_X,
    "axis_Y" = qsm@cylinder$axis_Y,
    "axis_Z" = qsm@cylinder$axis_Z)

  # loop through branches
  for (branch_id in unique(qsm@cylinder$branch)) {

    # get indices of current branch rows
    sub_idx <- qsm@cylinder$cyl_id[qsm@cylinder$branch == branch_id]
    first_cyl <- qsm@cylinder$cyl_id[qsm@cylinder$branch == branch_id & qsm@cylinder$PositionInBranch == 1]

    # derive basic stats
    branch_order <- qsm@cylinder$BranchOrder[qsm@cylinder$cyl_id == first_cyl]
    branch_diameter <-  2 * radius$radius[radius$cyl_id == first_cyl]
    branch_volume <- sum(1000 * pi * length$length[length$cyl_id %in% sub_idx] * radius$radius[radius$cyl_id %in% sub_idx] ** 2)
    branch_area <- sum(2 * pi * length$length[length$cyl_id %in% sub_idx] * radius$radius[radius$cyl_id %in% sub_idx])
    branch_length <- sum(length$length[length$cyl_id %in% sub_idx])
    branch_height <- qsm@cylinder$start_Z[first_cyl] - qsm@cylinder$start_Z[1]

    # if first cylinder was added, use second cylinder to compute angle
    if (qsm@cylinder$added[qsm@cylinder$cyl_id == first_cyl] == 1 & length(sub_idx) > 1) {
      first_considered <- sub_idx[2]
    } else {
      first_considered <- first_cyl
    }
    parent_cyl <- qsm@cylinder$parent[qsm@cylinder$cyl_id == first_cyl]
    branch_angle <- ifelse(parent_cyl > 0, 180 / pi *
                             acos(as.numeric(axis[axis$cyl_id == first_considered,2:4]) %*%
                                    as.numeric(axis[axis$cyl_id == parent_cyl,2:4])), 0)

    # derive branch azimuth and zenith
    branch_azimuth <- 180 / pi * atan2(axis$axis_Y[axis$cyl_id == first_cyl], axis$axis_X[axis$cyl_id == first_cyl])
    branch_zenith <- 180 / pi * acos(axis$axis_Z[axis$cyl_id == first_cyl])

    # get parent cylinder
    first_cyl_idx <- qsm@cylinder$cyl_id[qsm@cylinder$branch == branch_id & qsm@cylinder$PositionInBranch == 1]
    parent_cyl_idx <- qsm@cylinder$parent[qsm@cylinder$cyl_id == first_cyl_idx]
    if (parent_cyl_idx > 0) {
      branch_parent <- qsm@cylinder$branch[qsm@cylinder$cyl_id == parent_cyl_idx]
    } else {
      branch_parent <- 0
    }

    # append branch data
    branch_curr <- data.table(
      "bra_id" = branch_id,
      "order" = as.integer(branch_order),
      "parent" = as.integer(branch_parent),
      "diameter" = branch_diameter,
      "volume" = branch_volume,
      "area" = branch_area,
      "length" = branch_length,
      "angle" = branch_angle,
      "height" = branch_height,
      "azimuth" = branch_azimuth,
      "zenith" = branch_zenith
    )
    branch <- rbind(branch, branch_curr)
  }

  # overwrite old data
  qsm@branch <- branch

  # return results
  return(qsm)
}

################################################################################

#' Update QSM
#'
#' @description
#' \code{updateQSM} updates the \code{overview} and \code{branch} data of a
#' \code{QSM} object based on its \code{cylinder} data.
#'
#' @param qsm An object of class \code{QSM}.
#'
#' @return
#' An updated object of class \code{QSM}.
#'
#' @seealso \code{\link{checkQSM}}
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # delete some cylinders
#' qsm@cylinder <- qsm@cylinder[qsm@cylinder$BranchOrder <= 3,]
#'
#' # update qsm
#' updateQSM(qsm)
updateQSM <- function(qsm, method = "TreeQSM") {

  # update basic statistics based on cylinder
  message("updating basic statistics ...")
  qsm <- updateQSM_basics(qsm)

  # update branch based on cylinder
  message("updating branch statistics ...")
  qsm <- updateQSM_branch(qsm)

  # this function needs to be called last
  # update crown statistics based on cylinder
  message("updating crown statistics ...")
  qsm <- updateQSM_crown(qsm, method)

  # return results
  return(qsm)
}

################################################################################
