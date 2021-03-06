% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/home.R
\name{SimulateHomeRange}
\alias{SimulateHomeRange}
\title{SimulateHomeRange}
\usage{
SimulateHomeRange(steps, start_xy, homerange, mu, rho, extent, extend,
  redist_plots, arrow_plots, alpha)
}
\arguments{
\item{steps}{number of locations}

\item{start_xy}{animal's starting longitude and latitude.}

\item{homerange}{homerange kernel raster}

\item{mu}{wrapped Cauchy mu parameter, default is 0}

\item{rho}{wrapped Cauchy rho parameter}

\item{col}{color palette for base raster, default is jet2.col()}

\item{extend}{integer, changes the extent of the plots to a certain number of
cells around the prob_raster. Default is 0.}

\item{extent}{sets a fixed extent for plots, default is NULL which uses the
extent of the homerange raster}

\item{redist_plots}{save all of the redistribution kernel probability plots
in the working directory. Default is FALSE.}

\item{arrow_plots}{save all of the redistribution kernel probability plots
with movement arrows in the working directory. Default is FALSE.}

\item{alpha}{alpha (i.e., transparency) value for plotting the redistribution
kernels. Default is .9}
}
\value{
Dataframe with simulation data
}
\description{
Runs a simulation of home range behavior based on nest and conspecific
  locations and a logistic function for angle bias
}
\details{
Still a work in progress. May modify code so all of the
  redistribution kernels are saved as a RasterStack. This would allow the
  legend scale to be set at fixed values and to potentially create ArcGIS
  rasters and KMLs after the simultion is run.
}
