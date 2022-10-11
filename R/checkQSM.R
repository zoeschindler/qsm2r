################################################################################
# MAIN FUNCTIONS
################################################################################

#' Check if QSM must be updated
#'
#' @description
#' \code{checkQSM} conducts some superficial tests whether the \code{overview}
#' and \code{branch} data are up-to-date with the \code{cylinder} data of a
#' \code{QSM} object.
#'
#' @param qsm An object of class \code{QSM}.
#' @param precision Number of decimal points to be compared.
#'
#' @return
#' A logical. \code{TRUE} if the object is likely up-to-date, \code{FALSE} if
#' the object needs to be updated.
#'
#' @seealso \code{\link{updateQSM}}
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # check qsm
#' checkQSM(qsm)
checkQSM <- function(qsm, precision = 2L) {

  if (!is(qsm, "QSM")) {
    stop("input must be from class 'QSM'")
  }

  if (!is(precision, "integer")) {
    stop("precision must be an integer")
  } else if (precision < 0 | precision > 10) {
    stop("precision must be between 1 and 10")
  }

  # check dimensions
  empty_branch <- nrow(qsm@branch) == 0
  empty_cylinder <- nrow(qsm@cylinder) == 0

  # check branch id and branch order
  same_branch_id <- all(sort(unique(qsm@cylinder$branch)) == (1:nrow(qsm@branch)))
  same_branch_order <- all(sort(unique(qsm@cylinder$BranchOrder)) == sort(unique(qsm@branch$order)))

  # check basic stats
  same_length <- round(sum(qsm@cylinder$length), precision) == round(qsm@overview$TotalLength, precision)
  same_area   <- round(sum(2 * pi * qsm@cylinder$radius * qsm@cylinder$length), precision) == round(qsm@overview$TotalArea, precision)
  same_volume <- round(sum(pi * qsm@cylinder$length * qsm@cylinder$radius ** 2), precision) == round(qsm@overview$TotalVolume / 1000, precision)

  # print output
  message(paste("\nchecking dimensions ..."))
  message(paste(" - branch not empty -", ifelse(!empty_branch, "yes", "no")))
  message(paste(" - cylinder not empty -", ifelse(!empty_cylinder, "yes", "no")))
  #
  message(paste("checking branch id & order ..."))
  message(paste(" - branch & cylinder have same branch ids -", ifelse(same_branch_id, "yes", "no")))
  message(paste(" - branch & cylinder have same branch orders -", ifelse(same_branch_order, "yes", "no")))
  #
  message(paste("checking basic stats ..."))
  message(paste(" - cylinder & overview have same length -", ifelse(same_length, "yes", "no")))
  message(paste(" - cylinder & overview have same area -", ifelse(same_area, "yes", "no")))
  message(paste(" - cylinder & overview have same volume -", ifelse(same_volume, "yes", "no"), "\n"))

  # all checks fine?
  fine <- all(c(
    !empty_branch, !empty_cylinder,
    same_branch_id, same_branch_order,
    same_length, same_area, same_volume))

  # print result
  if (fine) {
    message("everything seems to be fine\n")
  } else {
    message("please use 'updateQSM()' to update the QSM")
  }

  # return result
  return(invisible(fine))
}

################################################################################
