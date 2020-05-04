# ----------------------------- MOVE FUNCTIONS ---------------------------------
# Movement process functions for individual-based model simulation

#' AngleToPoint
#'
#' Calculates the angle (based on a unit circle, where due East is 0 and due
#' West is pi) between an origin (x,y) and a target (x,y).
#'
#' @usage CenterXYInCell(x, y, xmin, ymin, cellsize)
#'
#' @param origin_x starting location x value
#' @param origin_y starting location y value
#' @param target_x finishing location x value
#' @param target_y finishing location y value
#'
#' @return numeric
#' @export
#'
AngleToPoint <- function(origin_x,
                         origin_y,
                         target_x,
                         target_y){
  dx <- c(target_x - origin_x)
  dy <- c(target_y - origin_y)
  abs_angle <- atan2(dy, dx)
  abs_angle <- ifelse(abs_angle < 0, (2*pi) + abs_angle, abs_angle)
}


#' CenterXYInCell
#'
#' Centers x and y values into the center of a raster cell based on the xmin,
#' ymin, and cellsize parameters
#'
#' @usage CenterXYInCell(x, y, xmin, ymin, cellsize)
#'
#' @param x x value
#' @param y y value
#' @param xmin minimum value of x that establishes the grid arrangement
#' @param ymin minimum value of y that establishes the grid arrangement
#' @param cellsize cell size in units of the x and y values
#'
#' @return vector of centered x and y (i.e., (x,y))
#' @export
#'
CenterXYInCell <- function(x,
                           y,
                           xmin,
                           ymin,
                           cellsize) {
    x <- xmin + (floor(((x-xmin)/cellsize))*cellsize) + (cellsize/2)
    y <- ymin + (floor(((y-ymin)/cellsize))*cellsize) + (cellsize/2)
    output <- c(x, y)
    return(output)
}


#' CalculateAngleToPoint
#'
#' Calculates absolute angle in radians [0, 2*pi) of origin location to target
#' point
#'
#' @param origin_long column name of origin point longitude
#' @param origin_lat column name of origin point latitude
#' @param target_long column name of target point longitude
#' @param target_lat column name of target point latitude
#'
#' @return vector of angles
#' @export
#'
#' @details  Coordinates must be in identically-scaled units (e.g. UTM meters).
#' If location data is in lat/long, project coordinates to UTM with 'rgdal'
#' package before running this function.
#'
CalculateAngleToPoint <- function(origin_long,
                                  origin_lat,
                                  target_long,
                                  target_lat){
  dx <- c(target_long - origin_long)
  dy <- c(target_lat - origin_lat)
  abs_angle <- atan2(dy, dx)
  output <- ifelse(abs_angle < 0, (2*pi) + abs_angle, abs_angle)
  return(output)
}

#' ConvertAngle
#'
#' Converts a radian angle that is outside of the Unit Circle range (i.e.,
#'   [0, 2pi]) to the equivalent value that is within the Unit Circle range.
#' @param x radian value
#'
#' @return vector of one value
#' @export
#'
#' @details Used so that the angles never go outside of the Unit Circle range
#'   when angles are being added or subtracted from each other over and over
#'   during a simulation.
ConvertAngle <- function(x) {
  if (x > 2*pi) x <- x-(2*pi)
  if (x < 0) x <- x+(2*pi)
  return(x)
}

#' CoordinatesFromAngleLength
#'
#' Calculates a new x,y location based on an absolute angle in radians [0, 2*pi)
#'   and distance from origin x,y
#'
#' @usage CoordinatesFromAngleLength(xy, angle, length)
#' @param xy numeric, vector of two coordinates c(x,y)
#' @param angle numeric, radians [0, 2*pi) based on unit circle orientation.
#'   Therefore, "due east" is 0, "due west" is pi.
#' @param length numeric, same units as x,y coordinates
#'
#' @return vector of two coordinates c(x,y)
#' @export
#'
#' @details Coordinates must be in identically-scaled units (e.g. UTM meters).
#' If location data is in lat & long, project coordinates to UTM with 'rgdal'
#' package before running this function.
CoordinatesFromAngleLength <- function (xy,
                                        angle,
                                        length){
  new_xy <- c(xy[1] + cos(angle) * length, xy[2] + sin(angle) * length)
  return(new_xy)
}

#' CreateConNestProb
#'
#' Create ConNest Probability Raster
#'
#' CreateConNestProb(con_nest_raster, gamma_shape, gamma_rate, x, y, max_r,
#'   cellsize, base)
#'
#' @param con_nest_raster Raster,
#' @param gamma_shape Numeric
#' @param gamma_rate Numeric
#' @param x Numeric,
#' @param y Numeric,
#' @param max_r Numeric,
#' @param cellsize Numeric
#' @param base Raster,
#'
#' @return Raster
#' @export
#'
CreateConNestProb <- function(con_nest_raster,
                              gamma_shape,
                              gamma_rate,
                              x,
                              y,
                              max_r,
                              cellsize,
                              base){
  max_r_cells <- ceiling(max_r/cellsize)
  xmin <- xmin(base)
  ymin <- ymin(base)
  xy <- CenterXYInCell(x, y, xmin, ymin, cellsize)  # May be unnecessary
  cell_extent <- extent(xy[1]-(cellsize/2), xy[1]+(cellsize/2), xy[2]-
    (cellsize/2), xy[2]+(cellsize/2))
  cell <- setValues(raster(cell_extent, crs=projection(base), res=cellsize),1)
  movement_kernel <- extend(cell, c(max_r_cells, max_r_cells), value=NA)
  con_nest_crop <- crop(con_nest_raster, movement_kernel, snap='in')
#plot(con_nest_crop, col=terrain.colors(255))
  xy_pt <- data.frame(x = xy[1], y = xy[2])
  xy_con_nest <- extract(con_nest_crop, xy_pt)
  con_nest_adjust <- calc(con_nest_crop, fun=function(x){(x - xy_con_nest)/1000})
#  plot(con_nest_adjust, col=terrain.colors(255))
#  points(xy[1], xy[2], pch=20, col="blue")
#  (extract(con_nest_adjust, xy_pt)) #should be zero
  xy_log_scale <- NonlinearRangeRescaleGamma(x=(xy_con_nest/1000),
    shape=gamma_shape, rate=gamma_rate,  min=NULL, max=NULL, lowBound=1,
    upBound=NULL, movement_kernel=movement_kernel, negative=TRUE)
#  curve(LogisticByInflection(x, inflection=0, scale=xy_log_scale), -15, 15)
  LogisticByInflection2 <- function(x){
    x <- LogisticByInflection(x, inflection=0, scale=xy_log_scale)
  }
  con_nest_rescale <- calc(con_nest_adjust, fun=LogisticByInflection2)
#  plot(con_nest_rescale)
#  points(xy[1], xy[2], pch=20, col="blue")
  return(con_nest_rescale)
}

#' CreateParetoKernel
#'
#' Creates a probability matrix based on a Pareto distribution.
#'
#' @usage CreateParetoKernel(scale, shape, max_r, cellsize)
#'
#' @param scale scale parameter of Pareto distribution
#' @param shape shape parameter of Pareto distribution
#' @param max_r maximum radius of kernel in meters, default = 100
#' @param cellsize cell size in meters, default = 1
#'
#' @return RasterLayer
#' @export
#'
CreateParetoKernel <- function(scale,
                               shape,
                               max_r = 100,
                               cellsize = 1) {
  max_r_cells <- max_r/cellsize
  size = ceiling(max_r_cells) * 2 + 1
  center = ceiling(max_r_cells) + 1
  kernel <- new("matrix", 0, size, size)
  for (i in 1:size) for (j in 1:size) {
    r = sqrt((i - center)^2 + (j - center)^2) * cellsize
    if (r <= max_r)
      kernel[i, j] <- VGAM::dgpd(r,scale=scale,shape=shape,log=FALSE)
  }
  kernel[center, center] <- 1/scale
  kernel <- kernel / sum(kernel)
  # This last part deletes the cells at the edge if they are all zero
  if (all(kernel[1, ] == 0, kernel[, 1] == 0,
          kernel[nrow(kernel),] == 0, kernel[, ncol(kernel)] == 0)) {
    kernel <- kernel[2:(nrow(kernel) - 1), 2:(ncol(kernel) - 1)]
  }
  return(kernel)
}

#' CreateMoveKernel
#'
#' Create a movement kernel matrix based on a wrapped Cauchy distribution
#'   for direction and a Pareto distribution for distance.
#'
#' @usage CreateMoveKernel(max_r, cellsize, mu, rho, shape, scale,
#'    ignore_cauchy, ignore_pareto)
#'
#' @param max_r maximum radius of kernel in meters, default = 300
#' @param cellsize cell size in meters, default = 30
#' @param mu mu parameter of wrapped Cauchy distribution, 0 radians is due east
#'   because everything is based on the Unit Circle
#' @param rho rho parameter of wrapped Cauchy distribution
#' @param shape shape parameter of Pareto distribution
#' @param scale scale parameter of Pareto distribution
#' @param ignore_cauchy logical, removes cauchy kernel's contribution to output
#'   raster. Default is FALSE.
#' @param ignore_pareto logical, removes pareto kernel's contribution to output
#'    raster. Default is FALSE.
#'
#' @return matrix
#' @export
CreateMoveKernel <- function(max_r = 300,
                               cellsize = 30,
                               mu,
                               rho,
                               shape,
                               scale,
                               ignore_cauchy = FALSE,
                               ignore_pareto = FALSE) {
  # Create the empty kernel objects
  max_r_cells <- ceiling(max_r/cellsize)
  size <- max_r_cells * 2 + 1
  center <- max_r_cells + 1
  wrpc_kernel <- new("matrix", 0, size, size)
  gpd_kernel <- new("matrix", 0, size, size)
  for (i in 1:size) {
    for (j in 1:size) {
      r = sqrt((i - center)^2 + (j - center)^2) * cellsize
      b = AngleToPoint(center, center, j, i)
      if(r <= max_r){
        wrpc_kernel[i, j] <- round(suppressWarnings(circular::dwrappedcauchy(b,
          mu=mu, rho=rho)), 5)
        gpd_kernel[i, j] <- texmex::dgpd(r, sigma=scale, xi=shape, log=FALSE)
      }
    }
  }
  wrpc_kernel <- apply(wrpc_kernel, 2, rev)
  gpd_kernel[center, center] <- 1/scale
  # This last part deletes the cells at the edge if they are all zero
  if (all(wrpc_kernel[1, ] == 0, wrpc_kernel[, 1] == 0,
    wrpc_kernel[nrow(wrpc_kernel),] == 0, wrpc_kernel[, ncol(wrpc_kernel)] ==0))
    wrpc_kernel <- wrpc_kernel[2:(nrow(wrpc_kernel) - 1), 2:(ncol(wrpc_kernel)
      - 1)]
  if (all(gpd_kernel[1, ] == 0, gpd_kernel[, 1] == 0,
    gpd_kernel[nrow(gpd_kernel),] == 0, gpd_kernel[, ncol(gpd_kernel)] == 0))
    gpd_kernel <- gpd_kernel[2:(nrow(gpd_kernel) - 1), 2:(ncol(gpd_kernel) - 1)]
  # Multiply the two kernels together and re-normalize
  if (ignore_cauchy) wrpc_kernel <- 1
  if (ignore_pareto) gpd_kernel <- 1
  move_kernel <- gpd_kernel*wrpc_kernel
  move_kernel <- move_kernel/sum(move_kernel)
  return(move_kernel)
}

#' CreateMoveKernelWeibull
#'
#' Create a movement kernel matrix based on a wrapped Cauchy distribution
#'   for direction and a Weibull distribution for distance.
#'
#' @usage CreateMoveKernelWeibull(max_r, cellsize, mu, rho, shape, scale,
#'    ignore_cauchy, ignore_weibull)
#'
#' @param max_r maximum radius of kernel in meters, default = 300
#' @param cellsize cell size in meters, default = 30
#' @param mu mu parameter of wrapped Cauchy distribution, 0 radians is due east
#'   because everything is based on the Unit Circle
#' @param rho rho parameter of wrapped Cauchy distribution
#' @param shape shape parameter of Pareto distribution
#' @param scale scale parameter of Pareto distribution
#' @param ignore_cauchy logical, removes cauchy kernel's contribution to output
#'   raster. Default is FALSE.
#' @param ignore_weibull logical, removes Weibull kernel's contribution to
#'   output raster. Default is FALSE.
#'
#' @return matrix
#' @details the Weibull parameters need to be based on distance in meters
#'
#' @export
#'
CreateMoveKernelWeibull <- function(max_r = 300,
                                      cellsize = 30,
                                      mu,
                                      rho,
                                      shape,
                                      scale,
                                      ignore_cauchy = FALSE,
                                      ignore_weibull = FALSE) {
  if (is.null(max_r)) max_r <- qweibull(.99, shape, scale) #* 1000
  max_r_cells <- ceiling(max_r/cellsize)
  size <- max_r_cells * 2 + 1
  center <- max_r_cells + 1
  angle_matrix <- new("matrix", 0, size, size)
  row_matrix <- row(angle_matrix)
  col_matrix <- col(angle_matrix)
  distance_matrix <- new("matrix", 0, size, size)
  weibull_kernel <- new("matrix", 0, size, size)
  dx <- row_matrix - center
  dy <- col_matrix - center
  abs_angle <- atan2(dx, dy)
  angle_matrix <- ifelse(abs_angle < 0, (2*pi) + abs_angle, abs_angle)
  wrpc_kernel <- suppressWarnings(circular::dwrappedcauchy(angle_matrix,
    mu=mu, rho=rho))
  wrpc_kernel <- apply(wrpc_kernel, 2, rev)
  distance_matrix <- (sqrt((row_matrix - center)^2 + (col_matrix - center)^2) *
      cellsize) #/ 1000
  weibull_kernel[] <- dweibull(as.vector(distance_matrix), shape=shape,
    scale=scale)
  weibull_kernel[center, center] <- 0 # probabilty at cell cell = Inf
  distance_matrix[distance_matrix > (max_r)] <- NA # (max_r/1000)] <- NA
  distance_matrix[!is.na(distance_matrix)] <- 1
  distance_matrix[is.na(distance_matrix)] <- 0
  distance_matrix[center, center] <- 0 # forces agent to move from center cell
  wrpc_kernel <- distance_matrix*wrpc_kernel
  weibull_kernel <- distance_matrix*weibull_kernel
  # This last part deletes the cells at the edge if they are all zero
  if (all(wrpc_kernel[1, ] == 0, wrpc_kernel[, 1] == 0,
    wrpc_kernel[nrow(wrpc_kernel),] == 0, wrpc_kernel[, ncol(wrpc_kernel)]==0)){
    wrpc_kernel <- wrpc_kernel[2:(nrow(wrpc_kernel) - 1), 2:(ncol(wrpc_kernel)
      - 1)]
  }
  if (all(weibull_kernel[1, ] == 0, weibull_kernel[, 1] == 0,
    weibull_kernel[nrow(weibull_kernel),] == 0, weibull_kernel[,
      ncol(weibull_kernel)] == 0)){
    weibull_kernel <- weibull_kernel[2:(nrow(weibull_kernel) - 1),
      2:(ncol(weibull_kernel) - 1)]
  }
  # Multiply the two kernels together and re-normalize
  if (ignore_cauchy) wrpc_kernel <- 1
  if (ignore_weibull) weibull_kernel <- 1
  move_kernel <- weibull_kernel*wrpc_kernel
  move_kernel <- move_kernel/sum(move_kernel)
  return(move_kernel)
}

#' CreateMoveKernelWeibullVonMises
#'
#' Create a movement kernel matrix based on a mixed von Mises distribution
#'   for direction and a Weibull distribution for distance.
#'
#' @usage CreateMoveKernelWeibullVonMises(max_r, cellsize, mu1, mu2, kappa1,
#'    kappa2, mix, rho, shape, scale, ignore_cauchy, ignore_weibull)
#'
#' @param max_r maximum radius of kernel in meters, default = 300
#' @param cellsize cell size in meters, default = 30
#' @param pars dataframe with columns for parameters, default is NULL. Only uses
#'   first row. Overrides the paramaters below (mu1:scale)
#' @param mu1 mu1 parameter of mixed von Mises distribution, 0 radians is due
#'   East because everything is based on the Unit Circle
#' @param mu2 mu2 parameter of mixed von Mises distribution, 0 radians is due
#'   East because everything is based on the Unit Circle
#' @param kappa1 kappa1 parameter of mixed von Mises distribution
#' @param kappa2 kappa2 parameter of mixed von Mises distribution
#' @param mix mixture (p) parameter of mixed von Mises distribution
#' @param shape shape parameter of Weibull distribution
#' @param scale scale parameter of Weibull distribution
#' @param ignore_von_mises logical, removes mixed von Mises kernel's
#'   contribution to output raster. Default is FALSE.
#' @param ignore_weibull logical, removes Weibull kernel's contribution to
#'   output raster. Default is FALSE.
#'
#' @return matrix
#' @details the Weibull parameters need to be based on distance in meters
#'
#' @export
#'

CreateMoveKernelWeibullVonMises <- function(max_r = 300,
                                            cellsize = 30,
                                            pars = NULL,
                                            mu1,
                                            mu2,
                                            kappa1,
                                            kappa2,
                                            mix,
                                            shape,
                                            scale,
                                            ignore_von_mises = FALSE,
                                            ignore_weibull = FALSE) {
  if(!is.null(pars)){
    mu1 = pars$mvm_mu1[1]
    mu2 = pars$mvm_mu2[1]
    kappa1 = pars$mvm_kappa1[1]
    kappa2 = pars$mvm_kappa2[1]
    mix = pars$mvm_prop[1]
    shape = pars$weibull_shape[1]
    scale = pars$weibull_scale[1]
  }
  if(is.null(max_r)) max_r <- qweibull(.99, shape, scale) #* 1000
  max_r_cells <- ceiling(max_r/cellsize)
  size <- max_r_cells * 2 + 1
  center <- max_r_cells + 1
  angle_matrix <- new("matrix", 0, size, size)
  row_matrix <- row(angle_matrix)
  col_matrix <- col(angle_matrix)
  distance_matrix <- new("matrix", 0, size, size)
  weibull_kernel <- new("matrix", 0, size, size)
  dx <- row_matrix - center
  dy <- col_matrix - center
  abs_angle <- atan2(dx, dy)
  angle_matrix <- ifelse(abs_angle < 0, (2*pi) + abs_angle, abs_angle)
  mvm_kernel <- suppressWarnings(CircStats::dmixedvm(angle_matrix,
    mu1=mu1, mu2=mu2, kappa1=kappa1, kappa2=kappa2, p=mix))
  mvm_kernel <- apply(mvm_kernel, 2, rev)
  distance_matrix <- (sqrt((row_matrix - center)^2 + (col_matrix - center)^2) *
      cellsize) #/ 1000
  weibull_kernel[] <- dweibull(as.vector(distance_matrix), shape=shape,
    scale=scale)
  weibull_kernel[center, center] <- 0 # probabilty at cell cell = Inf
  absolute_distance_matrix <- distance_matrix
  distance_matrix[distance_matrix > (max_r)] <- NA # (max_r/1000)] <- NA
  distance_matrix[!is.na(distance_matrix)] <- 1
  distance_matrix[is.na(distance_matrix)] <- 0
  distance_matrix[center, center] <- 0 # forces agent to move from center cell
  mvm_kernel <- distance_matrix*mvm_kernel
  weibull_kernel <- distance_matrix*weibull_kernel
  # Multiply the two kernels together and re-normalize
  if (ignore_von_mises) mvm_kernel <- 1
  if (ignore_weibull) weibull_kernel <- 1
  move_kernel <- weibull_kernel*mvm_kernel
  move_kernel <- move_kernel/sum(move_kernel)
  move_kernel[absolute_distance_matrix > (max_r)] <- NA  # only NA in last step
  return(move_kernel)
}

#' MovementSubModel
#'
#' Movement submodel that predicts movements based on current location,
#'   homerange_kernel, nest_return probability, and others.
#'
#' @usage MovementSubModel(sim, agent_states, step_data, step)
#'
#' @param sim 'sim' object
#' @param agent_states 'agent_states' list object
#' @param step_data 'step_data' dataframe
#' @param step step interval
#'
#'
#' @return 'step_data' object
#' @export
#'
MovementSubModel <- function(sim = sim,
                             agent_states = agent_states,
                             step_data = step_data,
                             step = step) {
  sim <- sim
  base <- sim$spatial$base
  cellsize <- raster::res(sim$spatial$base)[1]
  step_start <- lubridate::int_start(step)
  step_end <- lubridate::int_end(step)
  sex <- agent_states$sex
  season <- FindSeasonFromDatetime(datetime = step_start,
    seasons = sim$pars$global$sim_seasons)
  home_return <-
    sim$pars$classes[[sex]]$julian[lubridate::yday(lubridate::int_start(step)),
    2]
  step_max_r <- sim$pars$classes[[sex]]$season[[season]]$step_max_r
  step_cauchy_mu <- sim$pars$classes[[sex]]$constant$fixed$step_cauchy_mu
  step_cauchy_rho <- sim$pars$classes[[sex]]$constant$fixed$step_cauchy_rho
  step_weibull_shape <-
    sim$pars$classes[[sex]]$season[[season]]$step_weibull_shape
  step_weibull_scale <-
    sim$pars$classes[[sex]]$season[[season]]$step_weibull_scale

  connest_gamma_shape <-
    sim$pars$classes[[sex]]$constant$fixed$nestcon_gamma_shape
  connest_gamma_rate <-
    sim$pars$classes[[sex]]$constant$fixed$nestcon_gamma_rate

  con_nest_raster <- sim$spatial$con_nest_raster[[agent_states$nest_id]]

#  homerange_kernel <- sim$spatial$homerange_kernel[[agent_states$nest_id]]
#  landcover <- sim$spatial$landcover
#  hydro_dist <- sim$spatial$hydro_dist

  if (nrow(step_data) == 1) {
    i <- 1
    step_data[i+1, "datetime"] <- lubridate::int_end(step)
    step_data[1, "exp_angle"] <- sample(x=seq(from=0, to=(2*pi), by=(2*pi/360)),
      size=1)
    go_home <- FALSE
  } else {
    i <- nrow(step_data)
    step_data[i+1, "datetime"] <- lubridate::int_end(step)
    step_data$exp_angle[i] <- step_data$abs_angle[i-1]
  #  go_home <- rbinom(1, 1, home_return)
    go_home <- FALSE
  }
  if (go_home == TRUE) {
    home_xy <- c(agent_states$start_x, agent_states$start_y)
    step_data$x[i+1] <- home_xy[[1]]
    step_data$y[i+1] <- home_xy[[2]]
    step_data$abs_angle[i] <- CalculateAngleToPoint(step_data$x[i],
      step_data$y[i], step_data$x[i+1], step_data$y[i+1])
    step_data$step_length[i] <- as.integer(sqrt((step_data[i, "x"] -
      step_data[i+1, "x"])^2 + (step_data[i, "y"]-step_data[i+1, "y"])^2))
  } else {
    move_kernel <- CreateMoveKernelWeibull(max_r=step_max_r, cellsize=cellsize,
      mu=step_data$exp_angle[i], rho=step_cauchy_rho, shape=step_weibull_shape,
      scale=step_weibull_scale)
    r <- (cellsize*((nrow(move_kernel)-1)/2))+(cellsize/2)
    move_raster <- raster::raster(move_kernel, xmn=-r, xmx=r, ymn=-r, ymx=r)
    move_shift <- raster::shift(move_raster, x=step_data$x[i],
      y=step_data$y[i])
    ### PLACE TO ADD IN OTHER PROBABILITY LAYERS

    print(paste("x:", step_data$x[i], " y:", step_data$y[i]))

    move_shift <- raster::crop(move_shift, base, snap="in")
    ### NEW ConNestProb Raster
    con_nest <- CreateConNestProb(con_nest_raster,
      gamma_shape=connest_gamma_shape, gamma_rate=connest_gamma_rate,
      x=step_data$x[i], y=step_data$y[i], max_r=step_max_r, cellsize=cellsize,
      base=base)
    print(paste0("con_nest:", as.vector(raster::extent(con_nest)),
      "move_shift:", as.vector(raster::extent(move_shift))))
    con_nest_crop <- raster::crop(con_nest, move_shift, snap="out")
#    landcover_crop <- crop(landcover, move_shift, snap="out")
#    hydro_dist_crop <- crop(hydro_dist, move_shift, snap="out")
#    homerange_crop <- crop(homerange_kernel, move_shift, snap="out")
#    prob_raster <- overlay(move_shift, landcover_crop, hydro_dist_crop,
#      homerange_crop, fun=function(a,b,c,d) {return(a*b*c*d)}, recycle=FALSE)
#    prob_raster <- prob_raster/cellStats(prob_raster, stat="sum")
#    prob_raster <- move_shift
#    prob_raster <- prob_raster/cellStats(prob_raster, stat="sum")

    # USE GEOMETRIC MEAN for final probability layer?

    prob_raster <- raster::overlay(move_shift, con_nest_crop,
      fun=function(a,b){return(a*b)}, recycle=FALSE)
    prob_raster <- prob_raster/raster::cellStats(prob_raster, stat="sum")
    print("prob_min:", raster::minValue(prob_raster))
  #  prob_raster[prob_raster <= .000001] <- 0

    raster::crs(prob_raster) <- raster::crs(sim$spatial$base)
#    ExportKMLRasterOverlayWithTime(raster = prob_raster, time = step,
#      alpha = .8, color_pal= jet2.col(20),
#      outfile = paste0(agent_states$id, "_", i),
#      output_dir= file.path(getwd(), "Prob_Rasters"))

    ### END OF OTHER PROBABILITY LAYERS

#    plot(prob_raster)

    destination_cell <- suppressWarnings(sampling::strata(data=data.frame(cell=
      1:ncell(prob_raster)), stratanames=NULL, size=1, method="systematic",
      pik=prob_raster@data@values))

    while(is.na(destination_cell[1,1])) {
      destination_cell <- suppressWarnings(sampling::strata(data=data.frame(
        cell=1:ncell(prob_raster)), stratanames=NULL, size=1,
        method="systematic", pik=prob_raster@data@values))
    }

    print(paste("destination_cell:", destination_cell[1,1]))

    destination_xy <- raster::xyFromCell(prob_raster, destination_cell[1,1])

    print(paste("x:", destination_xy[1], " y:", destination_xy[2]))

    step_data[i+1, "x"] <- destination_xy[1]
    step_data[i+1, "y"] <- destination_xy[2]
    step_data$abs_angle[i] <- CalculateAngleToPoint(step_data$x[i],
      step_data$y[i], step_data$x[i+1], step_data$y[i+1])
    step_data$step_length[i] <- as.integer(sqrt((step_data[i, "x"] -
      step_data[i+1, "x"])^2 + (step_data[i, "y"] - step_data[i+1, "y"])^2))
  }
  step_data[i+1, "id"] <- step_data[i, "id"]
  return(step_data)
}

#' MovementSubModelBAEA
#'
#' Movement submodel that predicts movements based on current location,
#'   homerange_kernel, nest_return probability, and others.
#'
#' @usage MovementSubModelBAEA(sim, agent_states, step_data, step)
#'
#' @param sim 'sim' object
#' @param agent_states 'agent_states' list object
#' @param step_data 'step_data' dataframe
#' @param step step interval
#'
#'
#' @return 'step_data' object
#' @export
#'
MovementSubModelBAEA <- function(sim = sim,
                             agent_states = agent_states,
                             step_data = step_data,
                             step = step) {
  sim <- sim
  base <- sim$spatial$base
  cellsize <- raster::res(sim$spatial$base)[1]
  #step_start <- lubridate::int_start(step)
  #step_end <- lubridate::int_end(step)
  i <- which(step_data$datetime == step)
  current_behavior <- as.numeric(step_data[i, "behavior"])
  next_behavior <- as.numeric(step_data[i + 1, "behavior"])
  behavior_trans <- paste0(current_behavior, "_", next_behavior)
  sex <- agent_states$sex
  season <- FindSeasonFromDatetime(datetime = step,
    seasons = sim$pars$global$sim_seasons)
#  homerange_kernel <- sim$spatial$homerange_kernel[[agent_states$nest_id]]
#  landcover <- sim$spatial$landcover
#  hydro_dist <- sim$spatial$hydro_dist

  if (i == 1) {
    step_data[1, "exp_angle"] <- sample(x=seq(from=0, to=(2*pi), by=(2*pi/360)),
      size=1)
  } else {
    step_data$exp_angle[i] <- step_data$abs_angle[i-1]
  }
  #print(paste("behavior_trans:", behavior_trans))

  if (i == nrow(step_data)){
    print(paste("Last step for :", agent_states$id))
  } else if (behavior_trans %in% c("3_3", "5_5")){
    step_type <- "None"
  } else if (next_behavior == 3){
    step_type <- "To Nest"
  } else if (behavior_trans == "4_4"){
    perch_perch_pars <- sim$pars$classes$male$constant$fixed$move_pars %>%
      filter(behavior_behavior == "Perch -> Perch")
    bern_p <- perch_perch_pars$bern_p[1]
    step_type <- ifelse(rbinom(1, 1, bern_p), "Move", "None")
  } else {
    step_type <- "Move"
  }
  #print(paste("step_type:", step_type))
  if (i == nrow(step_data)){
    step_data$abs_angle[i] <- NA
    step_data$step_length[i] <- NA
  } else if (step_type == "None") { # no movement
    step_data$x[i+1] <- step_data$x[i]
    step_data$y[i+1] <- step_data$y[i]
    step_data$abs_angle[i] <- 0
    step_data$step_length[i] <- 0
  } else if (step_type == "To Nest" & i != nrow(step_data)) {
    home_xy <- c(agent_states$start_x, agent_states$start_y)
    step_data$x[i+1] <- home_xy[[1]]
    step_data$y[i+1] <- home_xy[[2]]
    step_data$abs_angle[i] <- CalculateAngleToPoint(step_data$x[i],
      step_data$y[i], step_data$x[i+1], step_data$y[i+1])
    step_data$step_length[i] <- as.integer(sqrt((step_data[i, "x"] -
      step_data[i+1, "x"])^2 + (step_data[i, "y"]-step_data[i+1, "y"])^2))
  } else if (step_type == "Move"){
    move_kernel <-sim$spatial$classes[[sex]][["move_kernels"]][[behavior_trans]]
    move_shift <- raster::shift(move_kernel, x=step_data$x[i],
      y=step_data$y[i])
    ### PLACE TO ADD IN OTHER PROBABILITY LAYERS
    #print(paste("x:", step_data$x[i], " y:", step_data$y[i]))
    move_shift <- raster::crop(move_shift, base, snap="in")
    ### NEW ConNestProb Raster
#    con_nest <- CreateConNestProb(con_nest_raster,
#      gamma_shape=connest_gamma_shape, gamma_rate=connest_gamma_rate,
#      x=step_data$x[i], y=step_data$y[i], max_r=step_max_r, cellsize=cellsize,
#      base=base)
#    print(paste0("con_nest:", as.vector(raster::extent(con_nest)),
#      "move_shift:", as.vector(raster::extent(move_shift))))
#    con_nest_crop <- raster::crop(con_nest, move_shift, snap="out")
#    landcover_crop <- crop(landcover, move_shift, snap="out")
#    hydro_dist_crop <- crop(hydro_dist, move_shift, snap="out")
#    homerange_crop <- crop(homerange_kernel, move_shift, snap="out")
#    prob_raster <- overlay(move_shift, landcover_crop, hydro_dist_crop,
#      homerange_crop, fun=function(a,b,c,d) {return(a*b*c*d)}, recycle=FALSE)
#    prob_raster <- prob_raster/cellStats(prob_raster, stat="sum")
    prob_raster <- move_shift
#    prob_raster <- prob_raster/cellStats(prob_raster, stat="sum")
#    prob_raster <- raster::overlay(move_shift, con_nest_crop,
#      fun=function(a,b){return(a*b)}, recycle=FALSE)
    prob_raster <- prob_raster/raster::cellStats(prob_raster, stat = "sum")
    prob_raster[prob_raster <= .000001] <- 0

      # USE GEOMETRIC MEAN for final probability layer?

    raster::crs(prob_raster) <- raster::crs(sim$spatial$base)
#    ExportKMLRasterOverlayWithTime(raster = prob_raster, time = step,
#      alpha = .8, color_pal= jet2.col(20),
#      outfile = paste0(agent_states$id, "_", i),
#      output_dir= file.path(getwd(), "Prob_Rasters"))

    ### END OF OTHER PROBABILITY LAYERS

    destination_cell <- suppressWarnings(sampling::strata(data=data.frame(cell=
      1:raster::ncell(prob_raster)), stratanames=NULL, size=1,
      method="systematic", pik=prob_raster@data@values))

    while(is.na(destination_cell[1,1])) {
      destination_cell <- suppressWarnings(sampling::strata(data=data.frame(
        cell=1:ncell(prob_raster)), stratanames=NULL, size=1,
        method="systematic", pik=prob_raster@data@values))
    }
    destination_xy <- raster::xyFromCell(prob_raster, destination_cell[1,1])
    step_data[i+1, "x"] <- destination_xy[1]
    step_data[i+1, "y"] <- destination_xy[2]
    step_data$abs_angle[i] <- CalculateAngleToPoint(step_data$x[i],
      step_data$y[i], step_data$x[i+1], step_data$y[i+1])
    step_data$step_length[i] <- as.integer(sqrt((step_data[i, "x"] -
      step_data[i+1, "x"])^2 + (step_data[i, "y"] - step_data[i+1, "y"])^2))

  }
  return(step_data)
}

#' NonlinearRangeRescaleGamma
#'
#' Non-linear Range Rescale Gamma Distribution
#'
#' @param x numeric, x value
#' @param shape numeric, Gamma shape
#' @param rate numeric
#' @param min numeric, value
#' @param max numeric
#' @param lowBound numeric
#' @param upBound numeric
#' @param movement_kernel kernel
#' @param negative logical
#'
#' @return vector
#' @export
#'
NonlinearRangeRescaleGamma <- function(x,
                                       shape = shape,
                                       rate = rate,
                                       min = NULL, #e.g., .001; pgamma
                                       max = NULL, #e.g., .99; pgamma
                                       lowBound = 1,
                                       upBound = NULL,
                                       movement_kernel = movement_kernel,
                                       negative = TRUE){
  if(is.null(min)){
    max_distance <- qgamma(0.999, shape=shape, rate=rate)
    min <- pgamma(max_distance, shape=shape, rate=rate, lower.tail=FALSE)
  }
  if(is.null(max)){
    max <- pgamma(.075, shape=shape, rate=rate, lower.tail=FALSE)
  }
  if(is.null(upBound)){
    upBound <- sqrt((xmin(movement_kernel)-xmax(movement_kernel))^2+
                      (ymin(movement_kernel)-ymax(movement_kernel))^2)/1000
  }
  # Get predicted y (gamma pdf) for x
  y_pred <- ppareto(x, shape=shape, rate=rate, lower.tail=FALSE)
  rescale <- lowBound + (((y_pred-min)/(max-min)) * (upBound-lowBound))
  if(negative==TRUE){
    rescale <- rescale*-1
  }
  return(rescale)
}

#' NonlinearRangeRescalePareto
#'
#' Non-linear Range Rescale Pareto Distribution
#'
#' @param x numeric, x (independent variable) value
#' @param scale numeric, scale (sigma) parameter for texmex::gpd
#' @param shape numeric, shape (xi) parameter for texmex::gpd
#' @param min_prob numeric, minimum Pareto probability value
#' @param max_prob numeric, maximum Pareto probability value
#' @param lower_bound numeric, lower limit for rescaled y values
#' @param upper_bound numeric, upper limit for rescaled y values
#' @param movement_kernel kernel, may be used to get the upper_bound
#' @param negative logical, whether or not to multiple outcome by -1
#'
#' @return vector
#' @export
#'
NonlinearRangeRescalePareto <- function(x,
                                        scale = shape,
                                        shape = shape,
                                        min_prob = NULL,
                                        max_prob = NULL,
                                        lower_bound = 1,
                                        upper_bound = NULL,
                                        movement_kernel = movement_kernel,
                                        negative = TRUE){
  if(is.null(min_prob)){
    max_distance <- texmex::qgpd(0.999, sigma=scale, xi=shape)
    min_prob <- texmex::pgpd(max_distance, sigma=scale, xi=shape, lower.tail=FALSE)
  }
  if(is.null(max_prob)){
    max_prob <- texmex::pgpd(.075, sigma=shape, xi=scale, lower.tail=FALSE)
  }
  if(is.null(upper_bound)){
    upper_bound <- sqrt((xmin(movement_kernel)-xmax(movement_kernel))^2+
                      (ymin(movement_kernel)-ymax(movement_kernel))^2)/1000
  }
  # Get predicted y (gamma pdf) for x
  #x = 5
  y_pred <- texmex::pgpd(x, sigma=scale, xi=shape, lower.tail=FALSE)
  rescale <- lower_bound + (((y_pred-min)/(max-min)) * (upper_bound-lower_bound))
  if(negative==TRUE){
    rescale <- rescale*-1
  }
  return(rescale)
}


