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
  if (is(file_path, "character")) {
    data_mat <- R.matlab::readMat(file_path)  # read matlab file from path
  } else {
    stop("input must be from class 'character'")
  }

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
}

################################################################################
