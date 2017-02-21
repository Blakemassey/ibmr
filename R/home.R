# -------------------- HOMERANGE FUNCTIONS -------------------------------------
# General functions for creation, analysis, and simulation of home ranges


#' CreateHomeRangeKernels
#'
#' Creates RasterLayers of homerange kernels based on home centroids
#'
#' @usage CreateHomeRangeKernels(df_all, df_home, base, max_r, home_inflection,
#' home_scale, avoid_inflection, avoid_scale, output_dir, write_distance,
#' write_homerange)
#'
#' @param df_all dataframe of all homerange centroids
#' @param df_home dataframe of homerange centroids to calculate homerange
#'  kernels (these must be a subset of the df_all dataframe), Default is to use
#'  'df_all' dataframe
#' @param base base Raster that sets the projection, extent, and dimensions of
#'  the study area
#' @param max_r maximum radius to calculate the homerange raster from
#' @param home_inflection inflection point of the Logistic function that governs
#'   the home kernel
#' @param home_scale scale parameter of the Logistic function that governs the
#'   scale parameter
#' @param avoid_inflection inflection point of the Logistic function that
#'   governs the conspecific avoidance kernel
#' @param avoid_scale scale parameter of the Logistic function that governs the
#'   conspecific avoidance kernel
#' @param output_dir directory for output files (distance, homerange)
#' @param write_distance logical, write distance raster to file. Default is
#'   FALSE.
#' @param write_homerange logical, write home range raster to file. Default is
#'   FALSE.
#'
#' @return List, contains homerange kernel Rasters for all the df_home
#'   centroids
#'
CreateHomeRangeKernels <- function(df_all,
                                   df_home = df_all,
                                   base,
                                   max_r,
                                   home_inflection,
                                   home_scale,
                                   avoid_inflection,
                                   avoid_scale,
                                   output_dir,
                                   write_distance = FALSE,
                                   write_homerange = FALSE) {
  source('C:/Work/R/Functions/sim/move.R')
  cellsize <- res(base)[1]
  max_r_cells <- ceiling(max_r/cellsize)
  xmin <- xmin(base)
  ymin <- ymin(base)
  df_sp <- SpatialPointsDataFrame(df_all[,c("x","y")], df_all,
    proj4string=crs(base))
  writeLines(noquote(paste("Calculating global distance")))
  if (write_distance == TRUE) {
    ifelse(exists("i"), i <- i , i <- 1)
    filename <- paste0(output_dir,"/global_dist_", sprintf("%03d", i), ".tif")
    global_dist <- distanceFromPoints(base, df_sp, filename=filename,
        overwrite=TRUE)
    writeLines(noquote(paste("Writing:", filename)))
  } else {
    global_dist <- distanceFromPoints(base, df_sp)
  }
  homerange_ids <- df_home$nest_id
  total <- length(homerange_ids)
  homerange_kernels <- as.list(setNames(rep(NA,length(homerange_ids)),
    homerange_ids), homerange_ids)
  for (j in 1:nrow(df_home)) {
    writeLines(noquote(paste("Calculating homerange", j, "of", total)))
    home <- df_home[j,]
    home_sp <- SpatialPointsDataFrame(home[,c("x","y")], home,
      proj4string=crs(base))
    xy <- CenterXYInCell(home_sp@coords[,"x"], home_sp@coords[,"y"],
      xmin, ymin, cellsize)
    cell_extent <- extent(xy[1]-(cellsize/2), xy[1]+(cellsize/2), xy[2]-
      (cellsize/2), xy[2]+(cellsize/2))
    cell <- setValues(raster(cell_extent, crs=projection(base), res=cellsize),j)
    home_ext <- extend(cell, c(max_r_cells, max_r_cells), value=NA)
    home_dist <- distanceFromPoints(home_ext, home[,c("x","y")])
    home_kern <- calc(home_dist, fun = function(x){(1/(exp((-(x -
      home_inflection)) / home_scale) + 1))})
    global_dist_crop <- crop(global_dist, home_dist)
    cent_dist <- overlay(home_dist, global_dist_crop, fun=function (x,y)
      {ifelse(x != y, NA, x)})
    cent_bounds <- boundaries(cent_dist)
    cent_bounds <- subs(cent_bounds, data.frame(from=c(0,1), to=c(NA,1)))
    edge_dist_abs <- distance(cent_bounds)
    edge_dist <- overlay(cent_dist, edge_dist_abs, fun=function(x,y)
      {ifelse(!is.na(x), y*-1, y*1)})
    avoid_kern <- calc(edge_dist, fun = function(x){(1/(exp((-(x -
      avoid_inflection)) / avoid_scale) + 1))})
    homerange_kern <- overlay(avoid_kern, home_kern, fun=function(x,y){
      p <- y*x #    y+x-1
    #       ifelse(p<=0, 0, p)
      })
    if (write_distance == TRUE) {
      k <- home[,"id"]
      filename <- paste0(output_dir,"/homerange_", k, "_",sprintf("%03d",
        i), ".tif")
      writeLines(noquote(paste0("Writing: ", filename)))
      writeRaster(homerange_kern, filename=filename, overwrite=TRUE)
    }
    homerange_kernels[[j]] <- homerange_kern
    names(homerange_kernels[[j]]) <- df_home[j,"nest_id"]
  }
  return(homerange_kernels)
}

#' PlotHomeRange
#'
#' Plot a homerange raster
#'
#' @usage PlotHomeRange(homerange, col, extent)
#'
#' @param homerange home range probability raster
#' @param col color palette for base raster, default is jet2.col()
#' @param extent an extent object for the plot, default is homerange extent
#'
#' @return plot of locations and pathways
#' @export
#'
PlotHomeRange <- function (homerange,
                           col = NULL,
                           extent = NULL){
  if(is.null(col)) col <- jet2.col(length(unique(homerange)))
  homerange_df <- data.frame(rasterToPoints(homerange))
  colnames(homerange_df)[3] <- "Probability"
  g <- ggplot(homerange_df, aes(x=x, y=y))+ geom_raster(aes(fill=Probability)                                          ,
    interpolate=TRUE) + scale_fill_gradientn(limits=c(0,1), colours=col)
  g <- g + xlab("X") + ylab("Y") +
  coord_fixed(ratio = 1) +
  theme(text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(colour="black")) +
  theme(axis.title.x = element_text(angle = 0, vjust = 0, hjust=0.5)) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  if (!is.null(extent)) {
    g <- g + scale_x_continuous(limits=extent[c(1,2)]) + scale_y_continuous(
      limits=extent[c(3,4)])
  }
  return(g)
}

#' PlotSimHomeRange
#'
#' Plot locations, pathways, and bias values of SimulateHomeRange() output
#' dataframes
#'
#' @usage PlotSimHomeRange(sim, home_df, homerange, col, label, extent)
#'
#' @param sim dataframe with x, y columns
#' @param df_home dataframe with x, y columns
#' @param homerange home range probability raster
#' @param col color palette for base raster, default is jet2.col()
#' @param label column name for point labels
#' @param extent an extent object for the plot, default is homerange extent
#'
#' @return plot of locations and pathways
#' @export
#'
PlotSimHomeRange <- function (sim,
                              df_home,
                              homerange = NULL,
                              col = NULL,
                              label = NULL,
                              extent = NULL){
  sim <- sim
  if (!is.null(label)) sim$label <- sim[,label]
  if (!is.null(homerange)) {
    homerange_df <- data.frame(rasterToPoints(homerange))
    colnames(homerange_df)[3] <- "Probability"
    if(is.null(col)) col <- jet2.col(length(unique(homerange)))
    g <- ggplot(homerange_df, aes(x=x, y=y))+ geom_raster(aes(fill=Probability)                                          ,
      interpolate=TRUE) + scale_fill_gradientn(limits=c(0,1),colours=col)
  } else {
    g <- ggplot(data=sim) + theme(legend.position="none")
  }
  g <- g + xlab("X") + ylab("Y") +
  geom_point(data=sim[1,], aes(x=x, y=y),size=3, color="green",
    fill="white", shape=18) +
  geom_point(data=sim[nrow(sim), ], aes(x=x, y=y),size=3, color="red",
    fill="black", shape=18) +
  geom_path(data=sim, aes(x=x, y=y), color= "gray90", arrow = arrow(angle=20,
    length=unit(.01, "npc"))) +
  geom_point(data=sim, aes(x=x, y=y), size=2, color="gray20", shape=20) +
  coord_fixed(ratio = 1) +
  theme(text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(colour="black")) +
  theme(axis.title.x = element_text(angle = 0, vjust = 0, hjust=0.5)) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5))
#  theme(panel.background = element_rect(fill = "white"))
  if (!is.null(label)) {
    g <- g + geom_text(aes(x=x, y=y, label=signif(label,3)), vjust=1.5,
      hjust=-.5, size=3)
  }
  if (!is.null(extent)) {
    g <- g + scale_x_continuous(limits=extent[c(1,2)]) + scale_y_continuous(
      limits=ext[c(3,4)])
  }
  return(g)
}

#' SimulateHomeRange
#'
#' Runs a simulation of home range behavior based on nest and conspecific
#'   locations and a logistic function for angle bias
#'
#' @usage SimulateHomeRange(steps, start_xy, homerange, mu, rho, extent, extend,
#'   redist_plots, arrow_plots, alpha)
#'
#' @param steps number of locations
#' @param start_xy animal's starting longitude and latitude.
#' @param homerange homerange kernel raster
#' @param mu wrapped Cauchy mu parameter, default is 0
#' @param rho wrapped Cauchy rho parameter
#' @param col color palette for base raster, default is jet2.col()
#' @param extend integer, changes the extent of the plots to a certain number of
#'   cells around the prob_raster. Default is 0.
#' @param extent sets a fixed extent for plots, default is NULL which uses the
#'   extent of the homerange raster
#' @param redist_plots save all of the redistribution kernel probability plots
#'   in the working directory. Default is FALSE.
#' @param arrow_plots save all of the redistribution kernel probability plots
#'   with movement arrows in the working directory. Default is FALSE.
#' @param alpha alpha (i.e., transparency) value for plotting the redistribution
#'   kernels. Default is .9
#'
#' @return Dataframe with simulation data
#' @export
#'
#' @details Still a work in progress. May modify code so all of the
#'   redistribution kernels are saved as a RasterStack. This would allow the
#'   legend scale to be set at fixed values and to potentially create ArcGIS
#'   rasters and KMLs after the simultion is run.
#'
SimulateHomeRange <- function(steps,
                              start_xy,
                              homerange,
                              mu = 0,
                              rho,
                              col = NULL,
                              extend = 0,
                              extent = NULL,
                              redist_plots = FALSE,
                              arrow_plots = FALSE,
                              alpha = .9) {
  cellsize <- res(homerange)[1]
  n <- steps
  df <- data.frame(step = seq(1:n))
  df$x <- NA
  df$y <- NA
  df$exp_angle <- NA
  df$abs_angle <- NA
  for (i in 1:n) {
    if (i == 1) {
      df$x[i] <- CenterXYInCell(start_xy[1], start_xy[2], xmin(homerange),
        ymin(homerange), cellsize)[1]
      df$y[i] <- CenterXYInCell(start_xy[1], start_xy[2], xmin(homerange),
        ymin(homerange), cellsize)[2]
      df$exp_angle[i] <- sample(x=seq(from=0, to=(2*pi), by=(2*pi/360)), size=1)
    }
    if (i > 1) df$exp_angle[i] <- df$abs_angle[i-1]
    if (i < n) {
      writeLines(noquote(paste0("Exp angle: ", round(df$exp_angle[i], 2))))
      # Create redistribution kernel raster
      redist <- CreateRedistKernel(max_r=45, cellsize=cellsize,
        mu=df$exp_angle[i], rho=rho, center_zero = TRUE)
      r <- (cellsize*((nrow(redist)-1)/2))+(cellsize/2)
      redist_raster <- raster(redist, xmn=-r, xmx=r, ymn=-r, ymx=r)
      redist_shift <- shift(redist_raster, x=df$x[i], y=df$y[i])
      homerange_crop <- crop(homerange, redist_shift)
      # Redistribution raster X homerange raster = probability raster
      writeLines(noquote(paste0("Calculating prob raster for step ", i)))
      prob_raster <- overlay(redist_shift, homerange_crop, fun=function(x,y)
        {return(x*y)}, recycle=FALSE)
      if (redist_plots == TRUE) {
        redist_shift2 <- extend(redist_shift, extend, value = 1)
        homerange_crop2 <- crop(homerange, redist_shift2)
        prob_raster2 <- overlay(redist_shift2, homerange_crop2,fun=function(x,y)
          {return(x*y)}, recycle=FALSE)
  #      prob_raster_NA[prob_raster_NA <= 0] <- NA
        if(is.null(col)) col <- jet2.col(length(unique(homerange)))
        prob_raster_df <- data.frame(rasterToPoints(prob_raster2))
        colnames(prob_raster_df)[3] <- "Probability"
        g <- ggplot(prob_raster_df, aes(x=x, y=y)) + geom_raster(aes(fill=
          Probability), interpolate=TRUE) + scale_fill_gradientn(limits=c(0,1),
          colours=col)
        g <- g + xlab("X") + ylab("Y") + coord_fixed(ratio = 1) +
        theme(plot.title=element_text(size=22)) +
        theme(text=element_text(size=20, colour="black")) +
        theme(axis.text=element_text(colour="black")) +
        xlab("X") + ylab("Y") + labs(title =  paste0("Step ", i))
        if (!is.null(extent)) {
          g <- g + scale_x_continuous(limits=extent[c(1,2)]) +
            scale_y_continuous(limits= extent[c(3,4)])
        }
        leading_zeros <- paste0("%0", nchar(steps), "d")
        SaveGGPlot(filename = paste0("Step", sprintf(leading_zeros, i),".jpeg"))
      }
      destination_cell <- suppressWarnings(strata(data=data.frame(cell=
        1:ncell(prob_raster)), stratanames=NULL, size=1, method="systematic",
        pik=prob_raster@data@values))
      destination_xy <- as.vector(xyFromCell(prob_raster,destination_cell[1,1]))
      df$x[i+1] <- destination_xy[1]
      df$y[i+1] <- destination_xy[2]
      df$abs_angle[i] <- CalculateAngleToPoint(df$x[i], df$y[i], df$x[i+1],
        df$y[i+1])
      if (arrow_plots == TRUE) {
        df2 <- cbind.data.frame(x=df$x[i], y=df$y[i], xend=df$x[i+1],
          yend=df$y[i+1])
        g2 <- g + geom_segment(data=df2, aes(x=x, y=y, xend=xend, yend=yend),
          arrow = arrow(length = unit(0.01, "npc")))
        SaveGGPlot(filename = paste0("Step",sprintf(leading_zeros, i),"a.jpeg"))
      }
    }
    if (i == n) df$abs_angle[i] <- NA
  }
  return(df)
}
