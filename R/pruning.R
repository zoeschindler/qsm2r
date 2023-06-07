################################################################################
# HELPER FUNCTIONS
################################################################################

pruning_remove <- function(qsm, cylinder, remove = TRUE) {

  # remove cylinders & prune column
  if (remove) {
    cylinder <- cylinder[!cylinder$prune,]
    cylinder$prune <- NULL
    message(" please use 'updateQSM()' to update the QSM")
  }

  # return modified qsm
  qsm@cylinder <- cylinder
  return(qsm)
}

################################################################################
# MAIN FUNCTIONS
################################################################################

#' Pruning according to cylinder ID
#'
#' @description
#' \code{pruning_cylinders} prunes the tree at the specified cylinders and also
#' removes subsequent cylinders.
#'
#' @param qsm An object of class \code{QSM}.
#' @param cyl_id \code{integer}, cylinder IDs at which the tree should be pruned.
#' @param remove \code{boolean}, whether the to be pruned cylinders should be
#' removed (\code{TRUE}) or labelled (\code{FALSE}).
#'
#' @return
#' \code{QSM}, with removed or labelled cylinders (column \code{cylinder$prune}).
#'
#' @seealso \code{\link{pruning_branches}}
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # prune at specified cylinders, without removing
#' not_removed <- pruning_cylinders(qsm, cyl_id = c(33,55,77), remove = FALSE)
#'
#' # plot qsm
#' plot(not_removed, col_var = "prune")
#'
#' # prune at specified cylinders, with removing
#' removed <- pruning_cylinders(qsm, cyl_id = c(33,55,77), remove = TRUE)
#'
#' # plot qsm
#' plot(removed)
#' @export
pruning_cylinders <- function(qsm, cyl_id, remove = FALSE) {

  # access cylinder data
  cylinder <- qsm@cylinder

  # get all children cylinders
  cyl_childs <- find_childs_recursive_cylinder(cylinder, cyl_id)

  # label all to be removed branches
  cylinder$prune <- FALSE
  cylinder$prune[cylinder$cyl_id %in% cyl_childs] <- TRUE

  # remove cylinders
  qsm <- pruning_remove(qsm, cylinder, remove)

  # return modified qsm
  return(qsm)
}

################################################################################

#' Pruning according to branch ID
#'
#' @description
#' \code{pruning_branches} prunes the tree at the specified branches and also
#' removes subsequent branches.
#'
#' @param qsm An object of class \code{QSM}.
#' @param branch_id \code{integer}, branch IDs of to be pruned branches.
#' @param remove \code{boolean}, whether the to be pruned cylinders should be
#' removed (\code{TRUE}) or labelled (\code{FALSE}).
#'
#' @return
#' \code{QSM}, with removed or labelled cylinders (column \code{cylinder$prune}).
#'
#' @seealso \code{\link{pruning_cylinders}}
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # prune specified branches, without removing
#' not_removed <- pruning_branches(qsm, branch_id = c(6,9), remove = FALSE)
#'
#' # plot qsm
#' plot(not_removed, col_var = "prune")
#'
#' # prune specified branches, with removing
#' removed <- pruning_branches(qsm, branch_id = c(6,9), remove = TRUE)
#'
#' # plot qsm
#' plot(removed)
#' @export
pruning_branches <- function(qsm, branch_id, remove = FALSE) {

  # access cylinder data
  cylinder <- qsm@cylinder

  # get all children cylinders
  branch_childs <- find_childs_recursive_branch(cylinder, branch_id)

  # label all to be removed branches
  cylinder$prune <- FALSE
  cylinder$prune[cylinder$branch %in% branch_childs] <- TRUE

  # remove cylinders
  qsm <- pruning_remove(qsm, cylinder, remove)

  # return modified qsm
  return(qsm)
}

################################################################################

#' Conventional pruning at specified height or stem length
#'
#' @description
#' \code{pruning_conventional} prunes all first order and their subsequent
#' branches below the specified tree height or stem length.
#'
#' @param qsm An object of class \code{QSM}.
#' @param threshold_m \code{numeric}, tree height or stem length in meters below
#' which all first order branches should be removed.
#' @param method \code{character}, whether the threshold refers to tree height
#' (\code{"height"}) or stem length (\code{"length"}).
#' @param remove \code{boolean}, whether the to be pruned cylinders should be
#' removed (\code{TRUE}) or labelled (\code{FALSE}).
#'
#' @return
#' \code{QSM}, with removed or labelled cylinders (column \code{cylinder$prune}).
#'
#' @seealso \code{\link{pruning_selective}}, \code{\link{pruning_whorlwise}}
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # prune specified branches, without removing
#' not_removed <- pruning_conventional(qsm, threshold_m = 3, method = "length", remove = FALSE)
#'
#' # plot qsm
#' plot(not_removed, col_var = "prune")
#'
#' # prune specified branches, with removing
#' removed <- pruning_conventional(qsm, threshold_m = 3, method = "length", remove = TRUE)
#'
#' # plot qsm
#' plot(removed)
#' @export
pruning_conventional <- function(qsm, threshold_m = 3, method = c("length", "height"), remove = FALSE) {

  # check input validity
  if (length(method) > 1 | !any(method %in% c("length", "height"))) {
    stop("method must be 'length' or 'height")
  }

  # first cylinders of main branches
  cylinder <- qsm@cylinder
  sub <- cylinder[cylinder$BranchOrder == 1 & cylinder$PositionInBranch == 1,]

  # method length
  if (method == "length") {

    # get stem length
    sub$method_m  <- apply(sub, 1, function(cyl) {
      parent_pos <- cylinder$PositionInBranch[cylinder$cyl_id == cyl["parent"]]
      stem_below <- cylinder[cylinder$BranchOrder == 0 & cylinder$PositionInBranch <= parent_pos,]
      return(sum(stem_below$length))
    })

  } else if (method == "height") {
    #  method height

    # get height above stem base
    stem_base_z <- cylinder$start_Z[cylinder$BranchOrder == 0 & cylinder$PositionInBranch == 1]
    sub$method_m <- sub$start_Z - stem_base_z
  }

  # get branch id's below threshold
  branch_id <- sub$branch[sub$method_m <= threshold_m]

  # get children of branches
  branch_id <- find_childs_recursive_branch(cylinder, branch_id, TRUE)

  # remove branches
  qsm <- pruning_branches(qsm, branch_id = branch_id, remove = remove)

  # return modified qsm
  return(qsm)
}

################################################################################

#' Selective pruning according to specified branch diameter and angle
#'
#' @description
#' \code{pruning_selective} prunes all first order and their subsequent
#' branches which have a larger starting diameter or a smaller branch angle than
#' specified.
#'
#' @param qsm An object of class \code{QSM}.
#' @param diameter_m \code{numeric}, branch diameter in meters above which all
#' first order branches should be removed.
#' @param angle_deg \code{numeric}, angle between stem and branches in degree
#' below which all first order branches should be removed.
#' @param remove \code{boolean}, whether the to be pruned cylinders should be
#' removed (\code{TRUE}) or labelled (\code{FALSE}).
#'
#' @return
#' \code{QSM}, with removed or labelled cylinders (column \code{cylinder$prune}).
#'
#' @seealso \code{\link{pruning_conventional}}, \code{\link{pruning_whorlwise}}
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # prune specified branches, without removing
#' not_removed <- pruning_selective(qsm, diameter_m = 0.03, angle_deg = 40, remove = FALSE)
#'
#' # plot qsm
#' plot(not_removed, col_var = "prune")
#'
#' # prune specified branches, with removing
#' removed <- pruning_selective(qsm, diameter_m = 0.03, angle_deg = 40, remove = TRUE)
#'
#' # plot qsm
#' plot(removed)
#'
#' @references
#' Springmann, S., Rogers, R., & Spiecker, H. (2011). Impact of artificial
#' pruning on growth and secondary shoot development of wild cherry
#' (\emph{Prunus avium} L.). Forest Ecology and Management, 261(3), 764-769.
#' doi: \href{https://doi.org/10.1016/j.foreco.2010.12.007}{10.1016/j.foreco.2010.12.007}
#'
#' @export
pruning_selective <- function(qsm, diameter_m = 0.03, angle_deg = 40, remove = FALSE) {

  # make subsets
  cylinder <- qsm@cylinder
  branch <- qsm@branch
  sub <- branch[branch$order == 1,]

  # get branches above certain diameter & angles
  branch_id <- sub$bra_id[sub$angle < angle_deg | sub$diameter > diameter_m]

  # get children of branches
  branch_id <- find_childs_recursive_branch(cylinder, branch_id, TRUE)

  # remove branches
  qsm <- pruning_branches(qsm, branch_id = branch_id, remove = remove)

  # return modified qsm
  return(qsm)
}

################################################################################

#' Whorlwise pruning to a specified number of whorls
#'
#' @description
#' \code{pruning_whorlwise} prunes all first order and their subsequent
#' branches belonging to whorls below a specified number of whorls. Whorls are
#' determined automatically. Single branches without other branches nearby are
#' counted as a single whorl as well.
#'
#' @param qsm An object of class \code{QSM}.
#' @param whorl_m \code{numeric}, maximum whorl stem interval in meters.
#' @param remaining_whorls \code{integer}, number of whorls that should remain
#' at the top of the crown.
#' @param remove \code{boolean}, whether the to be pruned cylinders should be
#' removed (\code{TRUE}) or labelled (\code{FALSE}).
#'
#' @return
#' \code{QSM}, with removed or labelled cylinders (column \code{cylinder$prune}).
#'
#' @seealso \code{\link{pruning_conventional}}, \code{\link{pruning_selective}}
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # prune specified branches, without removing
#' not_removed <- pruning_whorlwise(qsm, whorl_m = 0.3, remaining_whorls = 5, remove = FALSE)
#'
#' # plot qsm
#' plot(not_removed, col_var = "prune")
#'
#' # prune specified branches, with removing
#' removed <- pruning_whorlwise(qsm, whorl_m = 0.3, remaining_whorls = 5, remove = TRUE)
#'
#' # plot qsm
#' plot(removed)
#'
#' @references
#' Springmann, S., Rogers, R., & Spiecker, H. (2011). Impact of artificial
#' pruning on growth and secondary shoot development of wild cherry
#' (\emph{Prunus avium} L.). Forest Ecology and Management, 261(3), 764-769.
#' doi: \href{https://doi.org/10.1016/j.foreco.2010.12.007}{10.1016/j.foreco.2010.12.007}
#'
#' @export
pruning_whorlwise <- function(qsm, whorl_m = 0.3, remaining_whorls = 5, remove = FALSE) {

  # first cylinders of main branches
  cylinder <- qsm@cylinder
  sub <- cylinder[cylinder$BranchOrder == 1 & cylinder$PositionInBranch == 1,]

  # ignore very small branches
  sub <- sub[sub$radius * 2 > 0.01,]

  # get stem length
  stem_length_m  <- apply(sub, 1, function(cyl) {
    parent_pos <- cylinder$PositionInBranch[cylinder$cyl_id == cyl["parent"]]
    stem_below <- cylinder[cylinder$BranchOrder == 0 & cylinder$PositionInBranch <= parent_pos,]
    return(sum(stem_below$length))
  })

  # clustering branches to get whorls
  whorls <- value_cluster(values = stem_length_m, threshold = whorl_m)

  # check if pruning is necessary
  prune <- ifelse(length(whorls) > remaining_whorls, TRUE, FALSE)
  if (prune) {

    # get pruning length
    lowest_whorl <- length(whorls) - remaining_whorls + 1
    pruning_length_m <- min(whorls[[lowest_whorl]])

    # remove branches
    qsm <- pruning_conventional(qsm, threshold_m = pruning_length_m, method = "length", remove = remove)
  } else {
    if (!remove) {
      qsm@cylinder$prune <- FALSE
    }
  }

  # return old / modified qsm
  return(qsm)
}

################################################################################
