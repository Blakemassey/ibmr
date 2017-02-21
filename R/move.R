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
#' If location data is in lat & long, project coordinates to UTM with 'rgdal'
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

#' CreateRedistKernel
#'
#' Create a redistribution kernel matrix based on a wrapped Cauchy distribution
#'   for direction and a Pareto distribution for distance.
#'
#' @usage CreateRedistKernel(max_r, cellsize, mu, rho, shape, scale,
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
CreateRedistKernel <- function(max_r = 300,
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
  redist_kernel <- gpd_kernel*wrpc_kernel
  redist_kernel <- redist_kernel/sum(redist_kernel)
  return(redist_kernel)
}

#' ExportKMLRasterOverlayWithTime
#'
#' Export KML Raster function
#'
#' @usage ExportKMLRasterOverlayWithTime(raster, time, color_pal, alpha,
#'     maxpixels, blur, colNA, outfile, output_dir)
#'
#' @param raster a Raster* object
#' @param time Interval* object, default = NULL
#' @param color_pal color palette, can be a color ramp (e.g., c("white", "red")
#'    or a specific palette (e.g., "SAGA_pal[[1]]")
#' @param alpha numeric (0-1), transparency level of the .kml. Default is 1.
#' @param method method used to compute values for the new RasterLayer. Either
#'    'ngb' (nearest neighbor), which is useful for categorical variables, or
#'    'bilinear' (bilinear interpolation; the default value), which is
#'    appropriate for continuous variables.
#' @param overwrite logical, overwrite output file
#' @param maxpixels maximum number of pixels. If ncell(raster) > maxpixels,
#'    sample is used to reduce the number of pixels.
#' @param blur integer (default=10). Higher values help avoid blurring of
#'    isolated pixels (at the expense of a png file that is blur^2 times
#'    larger).
#' @param colNA color to use for the background (default is transparent)
#' @param outfile name of KML, default is to use name of raster
#' @param output_dir output folder location, default is getwd()
#' @param zip logical, whether or not to convert .kml to .kmz
#'
#' @return KML of a Raster
#' @export
#'
#' @details Modified from functions in the 'kml' and 'raster' packages
#'
ExportKMLRasterOverlayWithTime <- function(raster = raster,
                                           time = NULL,
                                           color_pal = rev(terrain.colors(255)),
                                           alpha = 1,
                                           method = "ngb",
                                           overwrite = TRUE,
                                           maxpixels = 500000,
                                           blur = 10,
                                           colNA = "transparent",
                                           outfile = NULL,
                                           output_dir= getwd(),
                                           zip = TRUE) {
  suppressPackageStartupMessages(require(raster))
  x <- raster
  if (nlayers(x) > 1) {
    x <- x[[1]]
  }
  if(!is.null(outfile)){
    name <- outfile
    outfile <- paste(output_dir, "/", name, ".kml", sep="")
  } else {
    name <- names(x)
    if (name == "layer") {
      name <- deparse(substitute(raster))
    }
    outfile <- paste(output_dir, "/", name, ".kml", sep="")
  }
  stopifnot(hasValues(x))
  x <- projectRaster(x, crs="+proj=longlat +datum=WGS84", method=method)
  unique_x <- length(unique(getValues(x)))
  col <- colorRampPalette(color_pal, alpha=TRUE)(unique_x)
  cols <- adjustcolor(col, alpha)
  #  if (unique_x > 250) unique_x <- 250
  filename <- extension(outfile, ".kml")
  x <- sampleRegular(x, size=maxpixels, asRaster=TRUE, useGDAL=TRUE)
  imagefile <- filename
  extension(imagefile) <- ".png"
  kmlfile <- kmzfile <- filename
  extension(kmlfile) <- ".kml"
  if (file.exists(kmlfile)) {
    if (overwrite) {
      file.remove(kmlfile)
    } else {
      stop("kml file exists, use \"overwrite=TRUE\" to overwrite it")
    }
  }
  png(filename=imagefile, width=max(480, blur * ncol(x)), height=max(480,
    blur * nrow(x)), bg="transparent", type="cairo-png")
  if (!is.na(colNA)) {
    par(mar = c(0, 0, 0, 0), bg = colNA)
  } else {
    par(mar = c(0, 0, 0, 0))
  }
  x[x== 0] <- NA
  image(x, col=cols, axes=FALSE, useRaster=TRUE, maxpixels=maxpixels)
  #plot(x, colNA="transparent")
  dev.off()
  kml <- c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
    "<kml xmlns=\"http://www.opengis.net/kml/2.2\">", "<GroundOverlay>")
  kmlname <- paste("<name>", name, "</name>", sep="")
  if(!is.null(time)){
    start_time <- paste0(strftime(int_start(time), "%Y-%m-%d", tz = time@tzone),
      "T", strftime(int_start(time), "%H:%M", tz=time@tzone), "Z")
    end_time <- paste0(strftime(int_end(time), "%Y-%m-%d", tz = time@tzone),
      "T", strftime(int_end(time), "%H:%M", tz=time@tzone), "Z")
    timespan <- paste0("<TimeSpan>", "<begin>", start_time, "</begin>", "<end>",
      end_time, "</end>", "</TimeSpan>")
  } else {
    timespan <- ""
  }
  icon <- paste("<Icon><href>", basename(imagefile),"</href><viewBoundScale>",
    "0.75</viewBoundScale></Icon>", sep="")
  e <- extent(x)
  latlonbox <- c("\t<LatLonBox>", paste("\t\t<north>", e@ymax,"</north><south>",
    e@ymin, "</south><east>", e@xmax, "</east><west>", e@xmin, "</west>",
    sep = ""), "\t</LatLonBox>")
  footer <- "</GroundOverlay></kml>"
  kml <- c(kml, kmlname, timespan, icon, latlonbox, footer)
  cat(paste(kml, sep = "", collapse = "\n"), file = kmlfile, sep = "")
  if(zip) ZipKML(kmlfile, imagefile)
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
#' @return 'step_data' object
#' @export
#'
MovementSubModel <- function(sim = sim,
                             agent_states = agent_states,
                             step_data = step_data,
                             step = step) {
  sim <- sim
  base <- sim$spatial$base
  cellsize <- res(sim$spatial$base)[1]
  step_start <- int_start(step)
  step_end <- int_end(step)
  sex <- agent_states$sex
  season <- FindSeasonFromDatetime(datetime = step_start,
    seasons = sim$pars$global$sim_seasons)
  home_return <- sim$pars$classes[[sex]]$julian[yday(int_start(step)),
    "home_return"]
  step_max_r <- sim$pars$classes[[sex]]$season[[season]]$step_max_r
  step_cauchy_mu <- sim$pars$classes[[sex]]$constant$fixed$step_cauchy_mu
  step_cauchy_rho <- sim$pars$classes[[sex]]$constant$fixed$step_cauchy_rho
  step_pareto_shape <-sim$pars$classes[[sex]]$season[[season]]$step_pareto_shape
  step_pareto_scale <-sim$pars$classes[[sex]]$season[[season]]$step_pareto_scale

#  homerange_kernel <- sim$spatial$homerange_kernel[[agent_states$nest_id]]
#  landcover <- sim$spatial$landcover
#  hydro_dist <- sim$spatial$hydro_dist

  if (nrow(step_data) == 1) {
    i <- 1
    step_data[i+1, "datetime"] <- int_end(step)
    step_data[1, "exp_angle"] <- sample(x=seq(from=0, to=(2*pi), by=(2*pi/360)),
      size=1)
    go_home <- FALSE
  } else {
    i <- nrow(step_data)
    step_data[i+1, "datetime"] <- int_end(step)
    step_data$exp_angle[i] <- step_data$abs_angle[i-1]
    go_home <- rbinom(1, 1, home_return)
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
    redist <- CreateRedistKernel(max_r=step_max_r, cellsize=cellsize,
      mu=step_data$exp_angle[i], rho=step_cauchy_rho, shape=step_pareto_shape,
      scale=step_pareto_scale)
    r <- (cellsize*((nrow(redist)-1)/2))+(cellsize/2)
    redist_raster <- raster(redist, xmn=-r, xmx=r, ymn=-r, ymx=r)
    redist_shift <- shift(redist_raster, x=step_data$x[i], y=step_data$y[i])

    ### PLACE TO ADD IN OTHER PROBABILITY LAYERS
    redist_shift <- crop(redist_shift, base, snap="in")

#    landcover_crop <- crop(landcover, redist_shift, snap="out")
#    hydro_dist_crop <- crop(hydro_dist, redist_shift, snap="out")
#    homerange_crop <- crop(homerange_kernel, redist_shift, snap="out")
#    prob_raster <- overlay(redist_shift, landcover_crop, hydro_dist_crop,
#      homerange_crop, fun=function(a,b,c,d) {return(a*b*c*d)}, recycle=FALSE)
#    prob_raster <- prob_raster/cellStats(prob_raster, stat="sum")

    prob_raster <- redist_shift
    prob_raster <- prob_raster/cellStats(prob_raster, stat="sum")

    crs(prob_raster) <- crs(sim$spatial$base)
#    ExportKMLRasterOverlayWithTime(raster = prob_raster, time = step,
#      alpha = .8, color_pal= jet2.col(20),
#      outfile = paste0(agent_states$id, "_", i),
#      output_dir= file.path(getwd(), "Prob_Rasters"))

    ### END OF OTHER PROBABILITY LAYERS

    destination_cell <- suppressWarnings(sampling::strata(data=data.frame(cell=
      1:ncell(prob_raster)), stratanames=NULL, size=1, method="systematic",
      pik=prob_raster@data@values))
    destination_xy <- xyFromCell(prob_raster, destination_cell[1,1])
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


#' ZipKML
#'
#' Zip a kml file
#'
#' @usage ZipKML(kml, image)
#'
#' @param kml filename for RasterStack or RasterBrick
#' @param image filename for image file
#'
#' @return a .zip file
#' @export
#'
#' @details Based on function from 'raster' package
#'
ZipKML <- function(kml,
                   image) {
	wd <- getwd()
	on.exit(setwd(wd))
 	setwd(dirname(kml))
	kml <- basename(kml)
	kmz <- extension(kml, '.kmz')
	image <- basename(image)
	if (file.exists(kmz)) {
		x <- file.remove(kmz)
	}
	kmzzip <- extension(kmz, '.zip')
	cmd <- paste('7z', 'a', kmzzip, kml, image, collapse=" ")
  sss <- try(system(cmd, intern=TRUE), silent=TRUE)
  file.rename(kmzzip, kmz)
  if (file.exists(kmz)) {
		file.remove(kml, image)
		return(invisible(kmz))
	} else {
		return(invisible(kml))
	}
}
