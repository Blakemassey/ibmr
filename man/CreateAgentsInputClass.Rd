% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\name{CreateAgentsInputClass}
\alias{CreateAgentsInputClass}
\title{CreateAgentsInputClass}
\usage{
CreateAgentsInputClass(n, class, sex_ratio, age_min, age_max)
}
\arguments{
\item{n}{number of individuals}

\item{class}{"class" of individuals}

\item{sex_ratio}{numeric [1-0], male/female sex ratio}

\item{age_min}{integer, minimum age (lubridate 'period' is determined later)}

\item{age_max}{integer, maximum age (lubridate 'period' is determined later)}

\item{start}{"random" or "regular"}

\item{base}{RasterLayer, base}
}
\value{
The input data frame with a column "birth_date" containing each
agent's birth_date as a POSIXct object.
}
\description{
Creates a data frame of input agents based on class, sex_ration, and min/max
age
}
