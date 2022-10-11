################################################################################
# HELPER FUNCTIONS
################################################################################

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
# MAIN FUNCTIONS
################################################################################

updateQSM_basics <- function(qsm) {

  # prepare storage
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
    qsm@cylinder$start_Z[1]
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

  # check method argument
  if(!method %in% c("TreeQSM", "qsm2r")) {
    stop("method must be either 'TreeQSM' or 'qsm2r'")
  }

  # prepare storage
  overview <- qsm@overview
  cylinder <- qsm@cylinder

  ### PART 1: create point cloud from cylinders

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
  trunk_cyl$angle_mod  <- rep(seq(0, 2*pi, length.out = 12), length.out = nrow(trunk_cyl))
  branch_cyl$angle_mod <- rep(seq(0, 2*pi, length.out = 4), length.out = nrow(branch_cyl))

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

  # delete duplicates
  points <- unique(points)

  ### PART 2: vertical profiles

  # I don't care about spreads

  ### PART 3: crown area & diameter

  # get coordinates of convex hull
  chull_idx <- grDevices::chull(points[,c("X", "Y")])
  chull_idx <- c(chull_idx, chull_idx[1])
  chull_n <- length(chull_idx)
  chull_x <- points$X[chull_idx]
  chull_y <- points$Y[chull_idx]

  # get centroid of the polygon
  conv_area <- 0.5 * sum(chull_x[1:(chull_n - 1)]*chull_y[2:chull_n] - chull_x[2:chull_n]*chull_y[1:(chull_n - 1)])
  center_X <- sum((chull_x[1:(chull_n - 1)] + chull_x[2:chull_n]) * (chull_x[1:(chull_n - 1)]*chull_y[2:chull_n] - chull_x[2:chull_n]*chull_y[1:(chull_n - 1)]))/6/conv_area
  center_Y <- sum((chull_y[1:(chull_n - 1)] + chull_y[2:chull_n]) * (chull_x[1:(chull_n - 1)]*chull_y[2:chull_n] - chull_x[2:chull_n]*chull_y[1:(chull_n - 1)]))/6/conv_area

  # choosing method
  if (method == "TreeQSM") {

    # calculate vector from center to cylinder tips
    vec_to_tips <- data.table::data.table(
      "vec_X" = (qsm@cylinder$start_X + qsm@cylinder$length * qsm@cylinder$axis_X) - center_X,
      "vec_Y" = (qsm@cylinder$start_Y + qsm@cylinder$length * qsm@cylinder$axis_Y) - center_Y)

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
      idx_1 <- which(vec_to_tips$angle >= (section-1)*pi/18 & vec_to_tips$angle < section*pi/18)
      len_1 <- ifelse(length(idx_1) > 0, max(vec_to_tips$angle[idx_1]), 0)

      # get length in opposite direction
      idx_2 <- which(vec_to_tips$angle >= ((section-1)*pi/18 + pi) & vec_to_tips$angle < (section*pi/18 + pi))
      len_2 <- ifelse(length(idx_2) > 0, max(vec_to_tips$angle[idx_2]), 0)

      # add crown width
      crown_widths <- c(crown_widths, len_1 + len_2)
    }

    # derive average crown width from section crown widths
    avg_crown_width <- mean(crown_widths)

    # choosing method
  } else if (method == "qsm2r") {

    # derive average crown width from crown projection area
    avg_crown_width <- 2 * sqrt(abs(conv_area) / pi)
  }

  # calculate maximum spread in the data
  max_crown_width <- max(dist(cbind(chull_x, chull_y)))

  # # derive alphahull
  # alpha_value <- max(0.5, avg_crown_width/10)
  # alpha_points <- unique(points[,c("X","Y")])
  # alpha_shape <- alphahull::ashape(x = alpha_points$X, y = alpha_points$Y, alpha = alpha_value)

  #
  #   %% Crown areas from convex hull and alpha shape:
  #   treedata.CrownAreaConv = A;                                               # HERE
  #   alp = max(0.5,treedata.CrownDiamAve/10);
  #   shp = alphaShape(X(:,1),X(:,2),alp);
  #   treedata.CrownAreaAlpha = shp.area;                                       # HERE

  ### PART 4: crown base, length, ratio & volume

  #   %% Crown base
  #   % Define first major branch as the branch whose diameter > min(0.05*dbh,5cm)
  #   % and whose horizontal relative reach is more than the median reach of 1st-ord.
  #   % branches (or at maximum 10). The reach is defined as the horizontal
  #   % distance from the base to the tip divided by the dbh.
  #   dbh = treedata.DBHcyl;
  #   nb = length(branch.order);
  #   HL = zeros(nb,1); % horizontal reach
  #   branches1 = (1:1:nb)';
  #   branches1 = branches1(branch.order == 1); % 1st-order branches
  #   nb = length(branches1);
  #   nc = size(Sta,1);
  #   ind = (1:1:nc)';
  #   for i = 1:nb
  #     C = ind(cylinder.branch == branches1(i));
  #     if ~isempty(C)
  #       base = Sta(C(1),:);
  #       C = C(end);
  #       tip = Sta(C,:)+Len(C)*Axe(C);
  #       V = tip(1:2)-base(1:2);
  #       HL(branches1(i)) = sqrt(V*V')/dbh*2;
  #     end
  #   end
  #   M = min(10,median(HL));
  #
  #   % Sort the branches according to the their heights
  #   Hei = branch.height(branches1);
  #   [Hei,SortOrd] = sort(Hei);
  #   branches1 = branches1(SortOrd);
  #
  #   % Search the first/lowest branch:
  #   d = min(0.05,0.05*dbh);
  #   b = 0;
  #   if nb > 1
  #     i = 1;
  #     while i < nb
  #       i = i+1;
  #       if branch.diameter(branches1(i)) > d && HL(branches1(i)) > M
  #         b = branches1(i);
  #         i = nb+2;
  #       end
  #     end
  #     if i == nb+1 && nb > 1
  #       b = branches1(1);
  #     end
  #   end
  #
  #   if b > 0
  #     % search all the children of the first major branch:
  #     nb = size(branch.parent,1);
  #     Ind = (1:1:nb)';
  #     chi = Ind(branch.parent == b);
  #     B = b;
  #     while ~isempty(chi)
  #       B = [B; chi];
  #       n = length(chi);
  #       C = cell(n,1);
  #       for i = 1:n
  #         C{i} = Ind(branch.parent == chi(i));
  #       end
  #       chi = vertcat(C{:});
  #     end
  #
  #     % define crown base height from the ground:
  #     BaseHeight = max(Sta(:,3)); % Height of the crown base
  #     for i = 1:length(B)
  #       C = ind(cylinder.branch == B(i));
  #       ht = min(Tip(C,3));
  #       hb = min(Sta(C,3));
  #       h = min(hb,ht);
  #       if h < BaseHeight
  #         BaseHeight = h;
  #       end
  #     end
  #     treedata.CrownBaseHeight = BaseHeight-Sta(1,3);                         # HERE
  #
  #     %% Crown length and ratio
  #     treedata.CrownLength = treedata.TreeHeight-treedata.CrownBaseHeight;    # HERE
  #     treedata.CrownRatio = treedata.CrownLength/treedata.TreeHeight;         # HERE
  #
  #     %% Crown volume from convex hull and alpha shape:
  #     I = P(:,3) >= BaseHeight;
  #     X = P(I,:);
  #     [K,V] = convhull(X(:,1),X(:,2),X(:,3));
  #     treedata.CrownVolumeConv = V;                                           # HERE
  #     alp = max(0.5,treedata.CrownDiamAve/5);
  #     shp = alphaShape(X(:,1),X(:,2),X(:,3),alp,'HoleThreshold',10000);
  #     treedata.CrownVolumeAlpha = shp.volume;                                 # HERE
  #
  #     else
  #       % No branches
  #       treedata.CrownBaseHeight = treedata.TreeHeight;                       # HERE
  #       treedata.CrownLength = 0;                                             # HERE
  #       treedata.CrownRatio = 0;                                              # HERE
  #       treedata.CrownVolumeConv = 0;                                         # HERE
  #       treedata.CrownVolumeAlpha = 0;                                        # HERE
  #     end
  #   end

  ### PART 5: gathering new data

  # overwrite old data
  overview$CrownAreaConv    <- conv_area
  overview$CrownDiamAve     <- avg_crown_width
  overview$CrownDiamMax     <- max_crown_width
  overview$CrownAreaAlpha   <- NA
  overview$CrownBaseHeight  <- NA
  overview$CrownLength      <- NA
  overview$CrownRatio       <- NA
  overview$CrownVolumeConv  <- NA
  overview$CrownVolumeAlpha <- NA

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

  # update crown statistics based on cylinder
  message("updating crown statistics ...")
  qsm <- updateQSM_crown(qsm, method)

  # update branch based on cylinder
  message("updating branch statistics ...")
  qsm <- updateQSM_branch(qsm)

  # return results
  return(qsm)
}

################################################################################
