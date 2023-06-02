################################################################################
# CLASS DEFINITION
################################################################################

#' @import data.table
#' @export
setClass(
  Class = "QSM",
  slots = list(
    name = "character",
    cylinder = "data.table",
    branch  = "data.table",
    overview = "list",
    input_parameters = "list",
    filter_parameters = "list",
    rundata = "list",
    pmdistance = "list"
  ),
  prototype = list(
    name = "tree",
    cylinder = data.table::data.table(),
    branch = data.table::data.table(),
    overview = list(),
    input_parameters = list(),
    filter_parameters = list(),
    rundata = list(),
    pmdistance = list()
  )
)

################################################################################
# CLASS METHODS
################################################################################

setMethod(
  "show",
  "QSM",
  function(object) {
    cat("class:       ", class(object), "\n")
    cat("name:        ", object@name, "\n")
    cat("DBH:         ", round(object@overview$DBHcyl * 100), "cm\n")
    cat("height:      ", round(object@overview$TreeHeight, 2), "m\n")
    cat("cylinder     ", "\n")
    cat(" - count:    ", nrow(object@cylinder), "\n")
    cat(" - length:   ", round(object@overview$TotalLength, 2), "m\n")
    cat(" - area:     ", round(object@overview$TotalArea, 2), "m^2\n")
    cat(" - volume:   ", round(object@overview$TotalVolume / 1000, 2), "m^3\n")
    return(invisible(object))
  }
)

################################################################################

#' Plot a QSM object
#'
#' @description
#' Displays a 3D \code{rgl} plot of a \code{QSM} object. Can be manipulated
#' using the package \code{rgl}.
#'
#' @param x An object of class \code{QSM}.
#' @param y Not used.
#' @param col \code{character}, color used for coloring the whole
#' \code{QSM} in a single color, if trying to color according to a variable,
#' leave it empty.
#' @param col_var \code{character}, column name of the cylinder data that should
#' be used for coloring.
#' @param pal \code{function}, color palette used for coloring.
#' @param sides \code{integer}, number of sides in the polygon cross section.
#' @param lit \code{boolean}, whether the scene should be with or without light.
#' @param center  \code{boolean}, whether all cylinders should be centered at
#' the stem base.
#'
#' @return
#' \code{rgl} plot of a QSM.
#'
#' @seealso \code{\link{readQSM}}
#'
#' @examples
#' # load qsm
#' file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
#' qsm <- readQSM(file_path)
#'
#' # plot qsm
#' plot(qsm, col = "salmon4")
#' plot(qsm, col_var = "PositionInBranch")
#' @export
#' @import rgl
#' @method plot QSM
setGeneric("plot", function(x, y, ...)
  standardGeneric("plot"))

#' @rdname plot
setMethod(
  "plot",
  signature(x = "QSM", y = "missing"),
  function(x, y, col = NULL, col_var = "BranchOrder", pal = grDevices::rainbow,
           sides = 6, lit = FALSE, center = TRUE) {

    # extract relevant data
    qsm <- x

    # extract cylinders
    if (is(qsm, "QSM")) {
      cylinder <- qsm@cylinder
    } else {
      stop("input must be from class 'QSM'")
    }

    # center cylinders
    if (center) {
      cylinder$start_X = cylinder$start_X - cylinder$start_X[1]
      cylinder$start_Y = cylinder$start_Y - cylinder$start_Y[1]
      cylinder$start_Z = cylinder$start_Z - cylinder$start_Z[1]
    }

    # calculate end points of cylinders
    cylinder$end_X = cylinder$start_X + cylinder$axis_X * cylinder$length
    cylinder$end_Y = cylinder$start_Y + cylinder$axis_Y * cylinder$length
    cylinder$end_Z = cylinder$start_Z + cylinder$axis_Z * cylinder$length

    # create color ramp
    cyl_vals <- unique(cylinder[[col_var]])
    col_n <- length(cyl_vals)
    col_vec <- pal(col_n)

    # if single color should be used
    if (!is.null(col)) {
      col_vec <- rep(col, col_n)
    }

    # assign the colors to the cylinders
    cylinder$color <- NA
    for (idx in 1:col_n) {
      cylinder$color[cylinder[[col_var]] == cyl_vals[idx]] <- col_vec[idx]
    }

    # print progress
    message("preparing cylinders ...")

    # prepare cylinders
    # https://stackoverflow.com/a/70684628/13427882
    cylinder_list <- lapply(1:nrow(cylinder), function(i) {
      cyl <- rgl::cylinder3d(
        center = cbind(
          c(cylinder$start_X[i], cylinder$end_X[i]),
          c(cylinder$start_Y[i], cylinder$end_Y[i]),
          c(cylinder$start_Z[i], cylinder$end_Z[i])),
        radius = cylinder$radius[i],
        closed = -2,
        sides = sides)
      cyl$material$color <- cylinder$color[i]
      cyl
    })

    # print progress
    message("plotting cylinders ...")

    # plot cylinders
    rgl::open3d()
    rgl::shade3d(shapelist3d(cylinder_list, plot = FALSE), lit = lit)
  }
)

################################################################################
