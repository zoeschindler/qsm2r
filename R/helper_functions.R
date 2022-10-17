################################################################################
# HELPER FUNCTIONS
################################################################################

find_childs_recursive_branch <- function(cylinder, bra_IDs, include_self = TRUE) {
  # get cylinders to start with
  cyl_sub <- cylinder[cylinder$branch %in% bra_IDs,]

  # get all cylinders which are children of the branches
  cyl_childs <- cylinder[cylinder$parent %in% cyl_sub$cyl_id & !(cylinder$branch %in% bra_IDs),]

  # return the branch IDs of the children
  if (nrow(cyl_childs) == 0) {
    return(NULL)
  } else {
    id_childs <- unique(cyl_childs$branch)
    id_childs_childs <- find_childs_recursive_branch(cylinder, id_childs)
    if (include_self) {
      return(c(bra_IDs, id_childs, id_childs_childs))
    } else {
      return(c(id_childs, id_childs_childs))
    }
  }
}

################################################################################

find_childs_recursive_cylinder <- function(cylinder, cyl_IDs, include_self = TRUE) {
  # get cylinders to start with
  cyl_sub <- cylinder[cylinder$cyl_id %in% cyl_IDs,]

  # get all cylinders which are children of the branches
  cyl_childs <- cylinder[cylinder$parent %in% cyl_sub$cyl_id,]

  # return the cylinder IDs of the children
  if (nrow(cyl_childs) == 0) {
    return(NULL)
  } else {
    id_childs <- unique(cyl_childs$cyl_id)
    id_childs_childs <- find_childs_recursive_cylinder(cylinder, id_childs)
    if (include_self) {
      return(c(cyl_IDs, id_childs, id_childs_childs))
    } else {
      return(c(id_childs, id_childs_childs))
    }
  }
}

################################################################################

# turns lists read in from matlab into data frames
# (assumes that all selected list entries are scalars of the same length)
# (transforms each list entry to a data frame column)
get_as_df <- function(target_list, dim = 1) {

  # target_list: list to be converted into a data frame
  # dim:         only used if first entry of list has multiple dimensions,
  #              dim = 1 -> stored as one variable per column
  #              dim = 2 -> stored as one variable per row

  # get list names
  list_names <- names(target_list)

  # prepare empty data frame
  first_entry = target_list[[list_names[1]]]
  if (length(first_entry) %in% dim(first_entry)) {
    list_nrow <- length(first_entry)
  } else {
    list_nrow <- dim(first_entry)[dim]
  }
  new_df <- data.frame(matrix(NA, nrow = list_nrow, ncol = 0))

  # add data column-wise to data frame
  for (list_name in list_names) {
    current_var <- data.frame(matrix(target_list[[list_name]], nrow = list_nrow)) # get data
    if (ncol(current_var) == 0) {
      current_var <- data.frame(matrix(NA, ncol = 1, nrow = 0))
    }
    colnames(current_var) <- list_name
    new_df <- cbind(new_df, current_var)
  }

  # rename columns
  if (!is.null(col) & length(col) == ncol(new_df)) {
    colnames(new_df) <- col
  }

  # return data frame
  return(new_df)
}

################################################################################

# calculates area of a polygon with n - 1 unique points with coordinates x and y
polygon_area <- function(x, y, n) {
  # first and last point of x and y must be the same
  area <- 0.5 * sum(x[1:(n - 1)]*y[2:n] - x[2:n]*y[1:(n - 1)])
  return(area)
}

################################################################################

# calculates centroid of a polygon with n - 1 unique vertices with coordinates x and y
polygon_centroid <- function(x, y, n, area) {
  # first and last point of x and y must be the same
  centroid <- list(
    "X" = sum((x[1:(n - 1)] + x[2:n]) * (x[1:(n - 1)]*y[2:n] - x[2:n]*y[1:(n - 1)]))/6/area,
    "Y" = sum((y[1:(n - 1)] + y[2:n]) * (x[1:(n - 1)]*y[2:n] - x[2:n]*y[1:(n - 1)]))/6/area)
  return(centroid)
}

################################################################################

# create random unit vector orthogonal to another unit vector
random_orth_norm <- function(x, y, z) {

  # get random vector
  orth_x <- rnorm(length(x))
  orth_y <- rnorm(length(y))
  orth_z <- (-x * orth_x - y * orth_y) / z

  # normalize vector
  orth_len <- sqrt(orth_x**2 + orth_y**2 + orth_z**2) # 3D length
  orth_x <- orth_x / orth_len
  orth_y <- orth_y / orth_len
  orth_z <- orth_z / orth_len

  # combine and rename
  orth_vec <- data.table::data.table(
    "orth_X" = orth_x,
    "orth_Y" = orth_y,
    "orth_Z" = orth_z)

  # return result
  return(orth_vec)
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

# create unit vector orthogonal to two other unit vectors
two_vector_orth_norm <- function(x1, y1, z1, x2, y2, z2) {

  # get orthogonal vector
  orth_x <- y1 * z2 - z1 * y2
  orth_y <- z1 * x2 - x1 * z2
  orth_z <- x1 * y2 - y1 * x2

  # normalize vector
  orth_len <- sqrt(orth_x**2 + orth_y**2 + orth_z**2) # 3D length
  orth_x <- orth_x / orth_len
  orth_y <- orth_y / orth_len
  orth_z <- orth_z / orth_len

  # combine and rename
  orth_vec <- data.table::data.table(
    "orth_X" = orth_x,
    "orth_Y" = orth_y,
    "orth_Z" = orth_z)

  # return result
  return(orth_vec)
}

################################################################################
