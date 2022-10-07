################################################################################
# HELPER FUNCTIONS
################################################################################

# loops through list entries
# assumes that all selected list entries are scalars of the same length
# transforms each list entry to a data frame column

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
# FUNCTIONS
################################################################################

readQSM <- function(file_path, qsm_var = 1, qsm_idx = 1) {

  # data_in: path to a matlab file containing the qsm or the read in matlab file
  # qsm_var: name of the qsm in the matlab file (if there are multiple objects)
  # qsm_idx: which QSM to take, if there are multiple

  # check datatype
  if (is(file_path, "character")) {
    data_mat <- readMat(file_path)  # read matlab file from path
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
      list(paste0(cylinder_names_old[idx], c("_X","_Y","_Z"))),
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
    cylinder = as.data.table(cylinder),
    branch  = as.data.table(branch),
    overview = as.list(treedata_overview),
    input_parameters = as.list(input_parameters),
    filter_parameters = as.list(filter_parameters),
    rundata = rundata,
    pmdistance = pmdist
  )

  # return results
  return(qsm)
}

checkQSM <- function(qsm, precision = 2L) {

  if (!class(qsm) == "QSM") {
    stop("input must be from class 'QSM'")
  }

  if (!class(precision) == "integer") {
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

updateQSM_crown <- function(qsm) {

  # prepare storage
  overview <- qsm@overview

  # TODO :(

  # $CrownDiamAve
  # $CrownDiamMax
  # $CrownAreaConv
  # $CrownAreaAlpha
  # $CrownBaseHeight
  # $CrownLength
  # $CrownRatio
  # $CrownVolumeConv
  # $CrownVolumeAlpha

  #   function [treedata,spreads] = crown_measures(treedata,cylinder,branch)
  #
  #   %% Generate point clouds from the cylinder model
  #   Axe = cylinder.axis;
  #   Len = cylinder.length;
  #   Sta = cylinder.start;
  #   Tip = Sta+[Len.*Axe(:,1) Len.*Axe(:,2) Len.*Axe(:,3)]; % tips of the cylinders
  #   nc = length(Len);
  #   P = zeros(5*nc,3); % four mid points on the cylinder surface
  #   t = 0;
  #   for i = 1:nc
  #     [U,V] = orthonormal_vectors(Axe(i,:)); # vectors U & V are unit vectors orthogonal to themselves and to the input vector U
  #     U = cylinder.radius(i)*U;
  #     if cylinder.branch(i) == 1
  #       % For stem cylinders generate more points
  #       R = rotation_matrix(Axe(i,:),pi/12);
  #       for k = 1:4
  #         M = Sta(i,:)+k/4*Len(i)*Axe(i,:);
  #         for j = 1:12
  #           if j > 1
  #             U = R*U;
  #           end
  #           t = t+1;
  #           P(t,:) = M+U';
  #         end
  #       end
  #     else
  #       M = Sta(i,:)+0.5*Len(i)*Axe(i,:);
  #       R = rotation_matrix(Axe(i,:),pi/4);
  #       for j = 1:4
  #         if j > 1
  #           U = R*U;
  #         end
  #         t = t+1;
  #         P(t,:) = M+U';
  #       end
  #     end
  #   end
  #
  #   P = P(1:t,:);
  #   I = ~isnan(P(:,1));
  #   P = P(I,:);
  #   P = double([P; Sta; Tip]);
  #   P = unique(P,'rows');
  #
  ############################################################################
  #
  #   %% Vertical profiles (layer diameters/spreads), mean:
  #   bot = min(P(:,3));
  #   top = max(P(:,3));
  #   Hei = top-bot;
  #   if Hei > 10
  #     m = 20;
  #   elseif Hei > 2
  #     m = 10;
  #   else
  #     m = 5;
  #   end
  #   spreads = zeros(m,18);
  #   for j = 1:m
  #     I = P(:,3) >= bot+(j-1)*Hei/m & P(:,3) < bot+j*Hei/m;
  #     X = unique(P(I,:),'rows');
  #     if size(X,1) > 5
  #       [K,A] = convhull(X(:,1),X(:,2));
  #       % compute center of gravity for the convex hull and use it as
  #       % center for computing average diameters
  #       n = length(K);
  #       x = X(K,1);
  #       y = X(K,2);
  #       CX = sum((x(1:n-1)+x(2:n)).*(x(1:n-1).*y(2:n)-x(2:n).*y(1:n-1)))/6/A;
  #       CY = sum((y(1:n-1)+y(2:n)).*(x(1:n-1).*y(2:n)-x(2:n).*y(1:n-1)))/6/A;
  #
  #       V = mat_vec_subtraction(X(:,1:2),[CX CY]);
  #       ang = atan2(V(:,2),V(:,1))+pi;
  #       [ang,I] = sort(ang);
  #       L = sqrt(sum(V.*V,2));
  #       L = L(I);
  #       for i = 1:18
  #         I = ang >= (i-1)*pi/18 & ang < i*pi/18;
  #         if any(I)
  #           L1 = max(L(I));
  #         else
  #           L1 = 0;
  #         end
  #         J = ang >= (i-1)*pi/18+pi & ang < i*pi/18+pi;
  #         if any(J)
  #           L2 = max(L(J));
  #         else
  #           L2 = 0;
  #         end
  #         spreads(j,i) = L1+L2;                                             # SPREADS
  #       end
  #     end
  #   end
  #
  ############################################################################
  #
  #   %% Crown diameters (spreads), mean and maximum:
  #   X = unique(P(:,1:2),'rows');
  #   [K,A] = convhull(X(:,1),X(:,2));
  #   % compute center of gravity for the convex hull and use it as center for
  #   % computing average diameters
  #   n = length(K);
  #   x = X(K,1);
  #   y = X(K,2);
  #   CX = sum((x(1:n-1)+x(2:n)).*(x(1:n-1).*y(2:n)-x(2:n).*y(1:n-1)))/6/A;
  #   CY = sum((y(1:n-1)+y(2:n)).*(x(1:n-1).*y(2:n)-x(2:n).*y(1:n-1)))/6/A;
  #   V = Tip(:,1:2)-[CX CY];
  #   ang = atan2(V(:,2),V(:,1))+pi;
  #   [ang,I] = sort(ang);
  #   L = sqrt(sum(V.*V,2));
  #   L = L(I);
  #   S = zeros(18,1);
  #   for i = 1:18
  #     I = ang >= (i-1)*pi/18 & ang < i*pi/18;
  #     if any(I)
  #       L1 = max(L(I));
  #     else
  #       L1 = 0;
  #     end
  #     J = ang >= (i-1)*pi/18+pi & ang < i*pi/18+pi;
  #     if any(J)
  #       L2 = max(L(J));
  #     else
  #       L2 = 0;
  #     end
  #     S(i) = L1+L2;
  #   end
  #   treedata.CrownDiamAve = mean(S);                                        # HERE
  #   MaxDiam = 0;
  #   for i = 1:n
  #     V = mat_vec_subtraction([x y],[x(i) y(i)]);
  #     L = max(sqrt(sum(V.*V,2)));
  #     if L > MaxDiam
  #       MaxDiam = L;
  #     end
  #   end
  #   treedata.CrownDiamMax = L;                                              # HERE
  #
  ############################################################################
  #
  #   %% Crown areas from convex hull and alpha shape:
  #   treedata.CrownAreaConv = A;                                             # HERE
  #   alp = max(0.5,treedata.CrownDiamAve/10);
  #   shp = alphaShape(X(:,1),X(:,2),alp);
  #   treedata.CrownAreaAlpha = shp.area;                                     # HERE
  #
  ############################################################################
  #
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
  #     treedata.CrownBaseHeight = BaseHeight-Sta(1,3);                       # HERE
  #
  #     %% Crown length and ratio
  #     treedata.CrownLength = treedata.TreeHeight-treedata.CrownBaseHeight;  # HERE
  #     treedata.CrownRatio = treedata.CrownLength/treedata.TreeHeight;       # HERE
  #
  #     %% Crown volume from convex hull and alpha shape:
  #     I = P(:,3) >= BaseHeight;
  #     X = P(I,:);
  #     [K,V] = convhull(X(:,1),X(:,2),X(:,3));
  #     treedata.CrownVolumeConv = V;                                         # HERE
  #     alp = max(0.5,treedata.CrownDiamAve/5);
  #     shp = alphaShape(X(:,1),X(:,2),X(:,3),alp,'HoleThreshold',10000);
  #     treedata.CrownVolumeAlpha = shp.volume;                               # HERE
  #
  #     else
  #       % No branches
  #       treedata.CrownBaseHeight = treedata.TreeHeight;                     # HERE
  #       treedata.CrownLength = 0;                                           # HERE
  #       treedata.CrownRatio = 0;                                            # HERE
  #       treedata.CrownVolumeConv = 0;                                       # HERE
  #       treedata.CrownVolumeAlpha = 0;                                      # HERE
  #     end
  #   end

  # overwrite old data
  qsm@overview <- overview

  # return results
  return(qsm)
}

updateQSM_branch <- function(qsm) {

  # prepare storage
  branch <- data.table()

  # create subsets
  radius <- data.table("cyl_id" = qsm@cylinder$cyl_id, "radius" = qsm@cylinder$radius)
  length <- data.table("cyl_id" = qsm@cylinder$cyl_id, "length" = qsm@cylinder$length)
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

    # if the first cylinder is added to fill a gap,
    # use the second cylinder to compute the angle
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

updateQSM <- function(qsm) {

  # update basic statistics based on cylinder
  message("updating basic statistics ...")
  qsm <- updateQSM_basics(qsm)

  # update crown statistics based on cylinder
  message("updating crown statistics ...")
  qsm <- updateQSM_crown(qsm)

  # update branch based on cylinder
  message("updating branch statistics ...")
  qsm <- updateQSM_branch(qsm)

  # return results
  return(qsm)
}

writeQSM <- function(qsm, file_path) {
  # TODO: which datatype? rds?
  # TODO: change readQSM depending on output datatype
  return(invisible(file_path))
}

################################################################################

get_location <- function(qsm) {

  # get starting coordinates of stem base
  base_id <- qsm@cylinder$cyl_id[qsm@cylinder$BranchOrder == 0 & qsm@cylinder$PositionInBranch == 1]
  location <- list(
    "x" = qsm@cylinder$start_X[qsm@cylinder$cyl_id == base_id],
    "y" = qsm@cylinder$start_Y[qsm@cylinder$cyl_id == base_id],
    "z" = qsm@cylinder$start_Z[qsm@cylinder$cyl_id == base_id])

  # show results
  cat("x:", round(location[["x"]], 2), "\n")
  cat("y:", round(location[["y"]], 2), "\n")
  cat("z:", round(location[["z"]], 2), "\n")

  # return results
  return(invisible(location))
}

set_location <- function(qsm, coordinates) {

  # recenter cylinders from coordinates to (0|0|0)
  qsm@cylinder$start_X <- qsm@cylinder$start_X - coordinates$x
  qsm@cylinder$start_Y <- qsm@cylinder$start_Y - coordinates$y
  qsm@cylinder$start_Z <- qsm@cylinder$start_Z - coordinates$z

  # return results
  return(qsm)
}

get_stemtaper <- function(qsm) {
  # TODO
  print("todo")
}

summary_branchorder <- function(qsm) {
  # TODO
  print("todo")
}

summary_cylinder_diameter <- function(qsm) {
  # TODO
  print("todo")
}

summary_cylinder_height <- function(qsm) {
  # TODO
  print("todo")
}

summary_cylinder_zenith <- function(qsm) {
  # TODO
  print("todo")
}

summary_cylinder_azimuth <- function(qsm) {
  # TODO
  print("todo")
}

summary_branch_diameter <- function(qsm) {
  # TODO
  print("todo")
}

summary_branch_height <- function(qsm) {
  # TODO
  print("todo")
}

summary_branch_zenith <- function(qsm) {
  # TODO
  print("todo")
}

summary_branch_azimuth <- function(qsm) {
  # TODO
  print("todo")
}

################################################################################
