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

#' todo
#'
#' @description
#' \code{todo} todo.
#'
#' @param todo \code{todo}, todo.
#'
#' @return
#' \code{todo}, todo.
#'
#' @seealso \code{\link{todo}}
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

#' todo
#'
#' @description
#' \code{todo} todo.
#'
#' @param todo \code{todo}, todo.
#'
#' @return
#' \code{todo}, todo.
#'
#' @seealso \code{\link{todo}}
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # prune specified branches, without removing
#' not_removed <- pruning_branches(qsm, branch_id = c(5,33), remove = FALSE)
#'
#' # plot qsm
#' plot(not_removed, col_var = "prune")
#'
#' # prune specified branches, with removing
#' removed <- pruning_branches(qsm, branch_id = c(3,111,333), remove = TRUE)
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

#' todo
#'
#' @description
#' \code{todo} todo.
#'
#' @param todo \code{todo}, todo.
#'
#' @return
#' \code{todo}, todo.
#'
#' @seealso \code{\link{todo}}
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

#' todo
#'
#' @description
#' \code{todo} todo.
#'
#' @param todo \code{todo}, todo.
#'
#' @return
#' \code{todo}, todo.
#'
#' @seealso \code{\link{todo}}
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

#' todo
#'
#' @description
#' \code{todo} todo.
#'
#' @param todo \code{todo}, todo.
#'
#' @return
#' \code{todo}, todo.
#'
#' @seealso \code{\link{todo}}
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
