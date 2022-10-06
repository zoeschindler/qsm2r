################################################################################
# DEPENDENCIES
################################################################################

library(data.table)
library(rgl)
library(R.matlab)

# library(terra)

################################################################################
# CLASS DEFINITION
################################################################################

setClass(
  "QSM",
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
    cylinder = data.table(),
    branch = data.table(),
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
    cat("DBH:         ", round(obj@overview$DBHcyl * 100), "cm\n")
    cat("height:      ", round(obj@overview$TreeHeight, 2), "m\n")
    cat("cylinders    ", "\n")
    cat("  - count:   ", nrow(object@cylinder), "\n")
    cat("  - length:  ", round(obj@overview$TotalLength, 2), "m\n")
    cat("  - area:    ", round(obj@overview$TotalArea, 2), "m²\n")
    cat("  - volume:  ", round(obj@overview$TotalVolume / 1000, 2), "m³\n")
    return(invisible(object))
  }
)

setMethod(
  "plot",
  "QSM",
  function(x, y = NULL, col = NULL, col_var = "BranchOrder", pal = rainbow,
           bg = "grey20", window = c(500,500), sides = 6) {

    # col:        single color to use for all cylinders
    # col_var:    which variable to use for coloring (e.g. branch, BranchOrder)
    # pal:        color palette to use (function)
    # bg:         background color
    # window:     initial window size
    # sides:      number of sides of each cylinder

    # extract relevant data
    qsm <- x

    # extract cylinders
    if (is(qsm, "QSM")) {
      cylinder <- qsm@cylinder
    } else {
      stop("input must be from class 'QSM'")
    }

    # calculate end points of cylinders
    cylinder$end_X = cylinder$start_X + cylinder$axis_X * cylinder$length
    cylinder$end_Y = cylinder$start_Y + cylinder$axis_Y * cylinder$length
    cylinder$end_Z = cylinder$start_Z + cylinder$axis_Z * cylinder$length

    # create color ramp
    cyl_vals <- unique(cylinder[,get(col_var)])
    col_n <- length(cyl_vals)
    col_vec <- pal(col_n)

    # if single color should be used
    if (!is.null(col)) {
      col_vec <- rep(col, col_n)
    }

    # assign the colors to the cylinders
    cylinder$color <- NA
    for (idx in 1:col_n) {
      cylinder$color[cylinder[,get(col_var)] == cyl_vals[idx]] <- col_vec[idx]
    }

    # print progress
    message("preparing cylinders ...")

    # prepare cylinders
    # https://stackoverflow.com/a/70684628/13427882
    cylinder_list <- lapply(1:nrow(cylinder), function(i) {
      cyl <- cylinder3d(
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
    open3d()
    par3d(windowRect = c(50, 50, window[1] + 50, window[2] + 50))
    bg3d(bg)
    shade3d(shapelist3d(cylinder_list, plot = FALSE), lit = FALSE)
  }
)

# TERRA
# setMethod("ext", "qsm", function(x, ...) { .ext(x) })
# setMethod("crs", "qsm", function(x, asText = FALSE) {})
# setMethod("area", "qsm", function(object) {})

################################################################################
