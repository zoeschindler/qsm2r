################################################################################
# MAIN FUNCTIONS
################################################################################

#' Read QSM from Matlab file
#'
#' @description
#' \code{readQSM} reads in one QSM from a \code{*.mat} file created with
#' \href{https://github.com/InverseTampere/TreeQSM}{TreeQSM} in Matlab.
#'
#' @param file_path Path to the \code{*.mat} file.
#' @param qsm_var Name of the variable in the file in which the QSM is stored.
#' Uses per default the first variable in the file.
#' @param qsm_idx  Index of the QSM in the variable is stored in.
#' Uses per default the first QSM in the variable.
#'
#' @return
#' A new object of class \code{QSM}.
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#' # or: qsm <- readQSM(file_path, qsm_var = "OptQSM", qsm_idx = 1)
#'
#' # inspect qsm
#' print(qsm)
#' @import data.table
#' @import R.matlab
#' @export
readQSM <- function(file_path, qsm_var = 1, qsm_idx = 1) {

  # data_in: path to a matlab file containing the qsm or the read in matlab file
  # qsm_var: name of the qsm in the matlab file (if there are multiple objects)
  # qsm_idx: which QSM to take, if there are multiple

  # check datatype
  if (is(file_path, "character") & endsWith(file_path, ".mat")) { # TreeQSM

    # read file
    data_mat <- R.matlab::readMat(file_path)

    # read in data
    data_mat <- data_mat[[qsm_var]][,,qsm_idx]

    # extract input parameters
    input_pars_all <- (data_mat$rundata[,,1])$inputs[,,1]
    if ("filter" %in% names(input_pars_all)) {
      filter_parameters <- get_as_df(input_pars_all$filter[,,1])
      input_parameters  <- get_as_df(input_pars_all[names(input_pars_all) != "filter"])
    } else {
      filter_parameters <- list()
      input_parameters <- get_as_df(input_pars_all)
    }

    # extract rundata
    rundata_time <- c(data_mat$rundata[,,1]$time)
    rundata_date_start <- as.POSIXct(
      paste(as.integer(data_mat$rundata[,,1]$date[1,]), collapse = "-"),
      format = "%Y-%m-%d-%H-%M-%S")
    rundata_date_end   <- as.POSIXct(
      paste(as.integer(data_mat$rundata[,,1]$date[2,]), collapse = "-"),
      format = "%Y-%m-%d-%H-%M-%S")
    rundata_version <- as.character(data_mat$rundata[,,1]$version)
    rundata <- list(
      "time" = rundata_time,
      "date_start" = rundata_date_start,
      "date_end" = rundata_date_end,
      "version" = rundata_version
    )

    # prepare cylinder data
    cylinder_names_old <- names(data_mat$cylinder[,,1])
    cylinder_names_new <- c()
    for (idx in 1:length(cylinder_names_old)) {
      cylinder_names_new <- c(cylinder_names_new, ifelse(
        cylinder_names_old[idx] %in% c("start", "axis"),
        list(paste0(cylinder_names_old[idx],c("_X","_Y","_Z"))),
        cylinder_names_old[idx])[[1]])
    }

    # extract cylinder data
    cylinder <- get_as_df(data_mat$cylinder[,,1], dim = 1)
    colnames(cylinder) <- cylinder_names_new
    cylinder <- cbind("cyl_id" = 1:nrow(cylinder), cylinder)

    # extract branch data
    branch <- get_as_df(data_mat$branch[,,1])
    branch <- cbind("bra_id" = 1:nrow(branch), branch)

    # prepare treedata
    tree_mat <- data_mat$treedata[,,1]
    tree_names <- names(tree_mat)
    idx_loc <- which(tree_names == "location")
    tree_overview_mat <- tree_mat[1:(idx_loc - 1)]
    tree_other_mat <- tree_mat[(idx_loc):length(tree_names)]

    # extract treedata - overview
    treedata_overview <- get_as_df(tree_overview_mat)

    # prepare pmdistance data
    pmdist_mat <- data_mat$pmdistance[,,1]
    pmdist_names <- names(pmdist_mat)

    # extract pmdistance data
    pmdist_overview <- get_as_df(pmdist_mat[pmdist_names != "CylDist"])
    pmdist_distance <- data.frame("CylDist" = pmdist_mat[["CylDist"]])
    pmdist <- append(as.list(pmdist_overview), pmdist_distance)

    # create QSM object
    qsm <- new(
      "QSM",
      name = input_parameters[["name"]],
      cylinder = data.table::as.data.table(cylinder),
      branch  = data.table::as.data.table(branch),
      overview = as.list(treedata_overview),
      input_parameters = as.list(input_parameters),
      filter_parameters = as.list(filter_parameters),
      rundata = rundata,
      pmdistance = pmdist
    )

    # return results
    return(qsm)

  } else if (is(file_path, "character") & endsWith(file_path, ".txt")) { # raycloudtool

    # read data
    txt_dat <- read.delim(file_path, skip = 2, header = FALSE)[,1]
    txt_dat <- do.call(rbind, strsplit(strsplit(txt_dat, ", ")[[1]], ","))

    # convert format
    class(txt_dat) <- "numeric"
    txt_dat <- data.frame(txt_dat)
    names(txt_dat) <- c("x","y","z","radius","parent_id","section_id")

    # fix naming starting with 0
    txt_dat$parent_id <- txt_dat$parent_id + 1
    txt_dat$id <- (1:nrow(txt_dat))

    # radius
    radius <- txt_dat$radius[-1]

    # get end coordinates (skip 1st because not a cylinder)
    end_X <- txt_dat$x[-1]
    end_Y <- txt_dat$y[-1]
    end_Z <- txt_dat$z[-1]

    # get start coordinates (+1 because in R, index starts at 1)
    start_X <- txt_dat$x[txt_dat$parent_id]
    start_Y <- txt_dat$y[txt_dat$parent_id]
    start_Z <- txt_dat$z[txt_dat$parent_id]

    # get length
    length <- sqrt((end_X - start_X)^2+(end_Y - start_Y)^2+(end_Z - start_Z)^2)

    # get axis vector
    axis_X <- (end_X - start_X) / length
    axis_Y <- (end_Y - start_Y) / length
    axis_Z <- (end_Z - start_Z) / length

    # correct ids (first entry is no cylinder)
    topology <- data.frame(
      "cyl_id" = txt_dat$id[-1] - 1,
      "parent" = txt_dat$parent_id[-1] - 1,
      "branch" = NA,
      "PositionInBranch" = NA,
      "BranchOrder" = NA)

    # determine branch_id
    for (curr_parent in unique(topology$parent)) {

      # first cylinder is first in stem
      if (curr_parent == 0) {

        # update data
        topology$branch[topology$parent == curr_parent] <- 1
        topology$PositionInBranch[topology$parent == curr_parent] <- 1
        topology$BranchOrder[topology$parent == curr_parent] <- 0
      } else {

        # get ids of cylinders with that parent
        curr_ids <- topology$cyl_id[topology$parent == curr_parent]

        # update data
        topology$branch[topology$cyl_id == curr_ids[1]] <- topology$branch[topology$cyl_id == curr_parent]
        topology$BranchOrder[topology$cyl_id == curr_ids[1]] <- topology$BranchOrder[topology$cyl_id == curr_parent]
        topology$PositionInBranch[topology$cyl_id == curr_ids[1]] <- topology$PositionInBranch[topology$cyl_id == curr_parent] + 1

        # check if more ids
        if (length(curr_ids) > 1) {

          # loop through ids
          for (curr_id in curr_ids[2:length(curr_ids)]) {

            # update data
            topology$branch[topology$cyl_id == curr_id] <- max(topology$branch, na.rm = TRUE) + 1
            topology$BranchOrder[topology$cyl_id == curr_id] <- topology$BranchOrder[topology$cyl_id == curr_parent] + 1
            topology$PositionInBranch[topology$cyl_id == curr_id] <- 1
          }
        }
      }
    }

    # fill cylinder data frame
    cylinder <- data.table::data.table(
      "cyl_id" = topology$cyl_id,
      "radius" = radius,
      "length" = length,
      "start_X" = start_X,
      "start_Y" = start_Y,
      "start_Z" = start_Z,
      "axis_X" = axis_X,
      "axis_Y" = axis_Y,
      "axis_Z" = axis_Z,
      "end_X" = end_X,
      "end_Y" = end_Y,
      "end_Z" = end_Z,
      "parent" = topology$parent,
      "branch" = topology$branch,
      "BranchOrder" = topology$BranchOrder,
      "PositionInBranch" = topology$PositionInBranch,
      "added" = 0)

    # estimate dbh
    dbh <- cylinder$radius[min(which(cumsum(cylinder$length) > 1.3))] * 2

    # create QSM object
    qsm <- new(
      "QSM",
      name = substr(basename(file_path), 1, (nchar(basename(file_path))-4)),
      cylinder = cylinder,
      branch  = data.table::data.table(),
      overview = list("DBHqsm" = dbh),
      input_parameters = list(),
      filter_parameters = list(),
      rundata = list(),
      pmdistance = list()
    )

    # fill missing data
    qsm <- qsm2r::updateQSM(qsm)
  } else {
    stop("input must be from class 'character'")
  }


}

################################################################################
