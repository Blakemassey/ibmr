# ----------------------- SIMULATION FUNCTIONS ---------------------------------
# General functions for simulation construction and analysis of simulation.
#

#' AgingSubModel
#'
#' @usage AgingSubModel(agent_states, dateime)
#'
#' @param agent_states agent_states object
#' @param step_data step_data object
#' @param step step object
#'
#' @return agent_states
#' @export
#'
AgingSubModel <- function(agent_states = agent_states,
                          step_data = step_data,
                          step = step){
  agent_states <- agent_states
  step_data <- step_data
  step <- step
  if (nrow(step_data) == 1){
    agent_states$current_age <- CalculateCurrentAge(datetime=int_start(step),
      birth_date=agent_states$birth_date)
  } else {
    current_step_day <- day(step_data[nrow(step_data), "datetime"])
    previous_step_day <- day(step_data[nrow(step_data)-1, "datetime"])
    if(current_step_day != previous_step_day) {
      agent_states$current_age <- CalculateCurrentAge(datetime=int_start(step),
        birth_date=agent_states$birth_date)
    }
  }
  if (is.na(agent_states$start_datetime)) {
    agent_states[["start_datetime"]] <- as.character(int_start(step))
  }
  return(agent_states)
}

#' CalculateCurrentAge
#'
#' Calculates current age of agent at a given step
#'
#' @usage CalculateCurrentAge(datetime, birthdate)
#'
#' @param datetime datetime of step
#' @param birth_date birthdate of the agent
#'
#' @return current_age in days
#' @export
#'
CalculateCurrentAge <- function(datetime,
                                birth_date) {
  current_age <- as.numeric(difftime(datetime, birth_date,
    tz=tz(sim$pars$global$sim_start)),
    units = c("days"))
  return(current_age)
}

#' CalculateReportAge
#'
#' Creates and updates population interval dataframe
#'
#' @usage CalculateReportAge(agent_states, step_data)
#' @param agent_states = agent_states
#' @param step_data = step data dataframe
#'
#' @return report_date in 'report_age_period' time Period
#' @export
#'
#' @examples
#'
CalculateReportAge <- function(agent_states,
                               step_data) {
  suppressPackageStartupMessages(require(zoo))
  birth_date <- ymd(agent_states$birth_date)
  datetime <- step_data[nrow(step_data), "datetime"]
  age_period <- sim$pars$global$report_age_period
  if (age_period == "day" || age_period == "days") {
    age <- floor(as.numeric(difftime(datetime, birth_date, units="days")))
  }
  if (age_period == "week" || age_period == "weeks") {
    age <- floor(as.numeric(difftime(datetime, birth_date, units="weeks")))
  }
  if (age_period == "month" || age_period == "months") {
    age <- floor((as.yearmon(datetime) - as.yearmon(birth_date))*12)
  }
  if (age_period == "year" || age_period == "years") {
    age <- floor(as.yearmon(datetime) - as.yearmon(birth_date))
  }
  return(age)
}

#' CompileAllAgentsStepData
#'
#' Creates a compiled dataframe of all the agent$step_data objects within
#' sim$agents$all
#'
#' @usage CompileAllAgentsStepData(sim)
#'
#' @param sim  a 'sim' list object
#'
#' @return a dataframe of all the step_data objects
#' @export
#'
#' @examples
#'
CompileAllAgentsStepData <- function(sim = sim) {
  all_step_data <- data.frame()
  all <- sim$agents$all
  for (i in 1:length(all)) {
    if (exists("step_data", where=sim$agents$all[[i]])){
      all_step_data <- rbind(all_step_data, sim$agents$all[[i]]$step_data)
    }
  }
  return(all_step_data)
}

#' ConvertStepDataCoordinates
#'
#' Converts "x", "y" columns to coordinates to lat/long WGS84 (for Google)
#'
#' @usage ConvertStepDataCoordinates(df)
#'
#' @param df input dataframe with "x", "y" columns
#' @param crs string of projection for "x", "y" columns, default is for
#'   Maine BAEA GPS data (UTM Zone 19N).
#'
#' @return A dataframe with added "long", "lat" columns
#' @export
#'

ConvertStepDataCoordinates <- function(df,
                                       crs = "+proj=utm +zone=19 ellps=WGS84"){
  suppressPackageStartupMessages(library(rgdal))
  suppressPackageStartupMessages(library(sp))
  df <- df
  coordinates(df) <- c("x", "y")
  proj4string(df) <- CRS(crs)
  res <- spTransform(df, CRS("+proj=longlat +datum=WGS84"))
  long_lat <- coordinates(res)
  colnames(long_lat) <- c("long", "lat")
  output <- cbind.data.frame(df, long_lat)
  output$optional <- NULL
  return(output)
}


#' CreateAgentsInputClass
#'
#' Creates a data frame of input agents based on class, sex_ration, and min/max
#' age
#'
#' @usage CreateAgentsInputClass(n, class, sex_ratio, age_min, age_max)
#'
#' @param n number of individuals
#' @param class "class" of individuals
#' @param sex_ratio  numeric [1-0], male/female sex ratio
#' @param age_min integer, minimum age (lubridate 'period' is determined later)
#' @param age_max integer, maximum age (lubridate 'period' is determined later)
#' @param start "random" or "regular"
#' @param base RasterLayer, base
#'
#' @return The input data frame with a column "birth_date" containing each
#' agent's birth_date as a POSIXct object.
#' @export
#'
CreateAgentsInputClass <- function(n = 100,
                                   class = "turtle",
                                   sex_ratio = .5,
                                   age_min = 1,
                                   age_max = 10,
                                   start = "random",
                                   base = base){
  id <- 1:n
  sex <- sample(c("male", "female"), n, replace=TRUE, prob = c(sex_ratio,
    1-sex_ratio))
  age <- sample(age_min:age_max, n, replace=TRUE)
  input <- cbind.data.frame(class, sex, age)
  if (start == "random"){
    base_sample <- raster::sampleRandom(x=base, size=n, na.rm=TRUE, xy=TRUE)
  }
  if (start == "regular"){
     base_sample <- raster::sampleRandom(x=base, size=n, na.rm=TRUE, xy=TRUE,
       sp=TRUE)
  }
  start_x <- base_sample[,1]
  start_y <- base_sample[,2]
  input <- cbind.data.frame(id, class, sex, age, start_x, start_y,
    stringsAsFactors = FALSE)
  return(input)
}

#' CreateBase
#'
#' Creates a 'base' RasterLayer
#'
#' @usage CreateAgentsInputClass(n, class, sex_ratio, age_min, age_max)
#'
#' @param x_min raster x minimum
#' @param y_min raster y minimum
#' @param x_max raster x maximum
#' @param y_max raster y maximum
#' @param cell_size integer, cell size (meters)
#'
#' @return RasterLayer
#' @export
#'
CreateBase <- function(x_min,
                       y_min,
                       x_max,
                       y_max,
                       cell_size){
  base <- raster::raster(xmn=x_min, xmx=x_max, ymn=y_min, ymx=y_max,
    resolution=cell_size, vals=NULL)
  base[] <- 1
#  base[] <- sample(1:10, ncell(base), replace=TRUE)
  return(base)
}

#' CreateBirthDate
#'
#' Creates a column in the input data frame specifying each agent's birth date
#' baed on its specified age
#'
#' @usage CreateBirthDate(input)
#'
#' @param input a data frame containing the agents used to start the model.
#' Must containing a column specifying "age" as a numeric value.
#'
#' @return The input data frame with a column "birth_date" containing each
#' agent's birth_date as a POSIXct object.
#' @export
#'
CreateBirthDate <- function(input = sim$agents$input){
  # Only proceed if there is no birth_date column
  if(is.null(input$birth_date)){
    # Loop through each row in the input
    input_age_period <- sim$pars$global$input_age_period
    birth_day <- sim$pars$global$birth_day
    for(a in 1:nrow(input)){
      # Is the age_period a year?
      if(input_age_period == "year" || input_age_period == "years") {
        # Determine the first sim_start date after the birth_day
        one_year <- as.period(1, "year")
        s0 <- as.Date(sim$pars$global$sim_start - (one_year*input$age[a]))
        # Set the format of the birth_day
        birth_day_format <- tail(guess_formats(birth_day, orders ="dm"), 1)
        birth_day_format <- paste(birth_day_format,"%Y",sep="")
        # Determine the first birth_day after s0
        s1 <- as.Date(paste(birth_day,year(s0),sep=""), format=birth_day_format)
        if(s0 >= s1) {
          input$birth_date[a] <- as.character(s1)
        } else {
          input$birth_date[a] <- as.character(s1-one_year)
        }
      } else {
        # If age period is not a year
        age_period_unit <- as.period(1, input_age_period)
        input$birth_date[a] <- as.character(sim$pars$global$sim_start -
          (age_period_unit*input$age[a]))
      }
    }
  }
  return(input)
}
#' CreateParsClassConstant
#'
#' Creates a parameter class 'constant'
#'
#' @usage CreateBirthDate(step_cauchy_mu, step_cauchy_rho)
#'
#' @param step_cauchy_mu step cauchy mu
#' @param step_cauchy_rho step cauchy rho
#'
#' @return list
#' @export
#'
CreateParsClassConstant <- function(step_cauchy_mu = 0,
                                    step_cauchy_rho = .5){
  step_cauchy_mu <- step_cauchy_mu
  step_cauchy_rho <- step_cauchy_rho
  fixed <- NamedList(step_cauchy_mu, step_cauchy_rho)
  constant <- NamedList(fixed)
  return(constant)
}

#' CreateParsClassSeason
#'
#' Creates a parameter class 'season'
#'
#' @usage CreateParsClassSeason(step_pareto_scale,
#'   step_pareto_shape, step_max_r)
#'
#' @param step_pareto_scale pareto scale
#' @param step_pareto_shape pareto shape
#' @param step_max_r max radius
#'
#' @return list
#' @export
#'
CreateParsClassSeason <- function(step_pareto_scale,
                                  step_pareto_shape,
                                  step_max_r){
  step_pareto_scale = step_pareto_scale
  step_pareto_shape = step_pareto_shape
  step_max_r = step_max_r
  season <- NamedList(step_pareto_scale, step_pareto_shape, step_max_r)
  return(season)
}

#' CreateParsClassSeason
#'
#' Creates a parameter class 'season'
#'
#' @usage CreateParsClassSeason(input)
#'
#' @param sim_start POSIXct default = as.POSIXct("2017-01-01", tz = "UTC")
#' @param sim_period 'period', default = period(100, "days")
#' @param sim_end POSIXct, default = as.POSIXct("2017-01-01", tz = "UTC")
#' @param rep_period 'period', default = period(1, "days")
#' @param rep_interval character string of dates, default = c("01Jan", "01Jul")
#' @param rep_interval_custom = character string of dates, dafault = NULL
#' @param step_period 'period', default = period(1, "day")
#' @param time_step_period 'period', default = period(6, "hour")
#' @param birth_day = character string of date, default = "01Jan"
#' @param input_age_period = character string of period, default = "years"
#' @param report_age_period = character string of period, default = "months"
#' @param sim_seasons = dataframe, columns = c(start, season), default are
#'   typical sesasons
#'
#' @import lubridate
#'
#' @return list
#' @export
#'
CreateParsGlobal <- function(sim_start = as.POSIXct("2015-01-01", tz = "UTC"),
                             sim_period = period(100, "days"),
                             sim_end = NULL,
                             rep_period = period(1, "days"),
                             rep_interval = c("01Jan", "01Jul"),
                             rep_interval_custom = NULL,
                             step_period = period(1, "day"),
                             time_step_period = period(6, "hour"),
                             birth_day = "01Jan",
                             input_age_period = "years",
                             report_age_period = "months",
                             sim_seasons = NULL){
  sim_start <- sim_start
  sim_period <- sim_period
  sim_end <- sim_end
  rep_period <- rep_period
  rep_interval <- rep_interval
  rep_interval_custom <- rep_interval_custom
  step_period <- step_period
  time_step_period <- time_step_period
  birth_day <- birth_day
  input_age_period = input_age_period
  report_age_period = report_age_period
  if (is.null(sim_seasons)){
    sim_seasons <- data.frame(start = c("20Mar", "21Jun", "23Sep", "22Dec"),
      season = c("spring", "summer", "winter", "fall"), stringsAsFactors=FALSE)
  }
  global <- NamedList(sim_start, sim_period, sim_end, rep_period, rep_interval,
    rep_interval_custom, step_period, time_step_period, sim_seasons, birth_day,
    input_age_period, report_age_period)
  return(global)
}

#' CreateRunsList
#'
#' Helper function for RunSimulation(). Creates a list of empty objects, each
#' one named individually as "run_xx" where the number of leading zeroes matches
#' the maximum digit width of the number of runs.
#'
#' @param runs = numeric, number of runs in RunSimulation()
#'
#' @return a list
#' @export
#'
#' @note Where there are < 10 runs, there will have no leading zeros, when
#' there are between 10 and 99 runs, there will be a leading zero for runs
#' numbered < 10.

CreateRunsList <- function(runs = runs) {
  format <- paste0("%0", nchar(runs), "d")
  runs_list <- list()
  for (i in 1:runs) {
    runs_list[[i]] <- NA
    names(runs_list)[[i]] <- paste0("run_", sprintf(format, i))
  }
  return(runs_list)
}



#' CreateReportIntervals
#'
#' Creates summary intervals for the simulation input
#'
#' @usage CreateReportIntervals(sim)
#'
#' @param sim a list that contains a "pars" list which contains a "global" list
#'   that sets the arguments for the Report Interval. See Description section
#'   below.
#'
#' @return list of intervals
#' @export
#'
#' @description
#' sim list that contains a "pars" list which contains a "global" list that
#' contains the following objects:
#'   sim_start = a POSIXct object for simulation's start, required
#'   sim_period = a Period object (from 'lubridate' package), e.g. period(10,
#'     "year"), period(10, "week"), period(10, "week") that sets the length of
#'     time the simulation runs
#'   sim_end = a POSIXct object for simulation's end. Setting this will override
#'     the 'sim_period' parameter. Default = NULL.
#'   rep_period = a Period object that sets the summary interval
#'   rep_interval = vector or dateframe (with date in first column) of the
#'     format %B%d, %b%d, %d%B, or %B%d (see ?strptime for more details on
#'     formatting). Example: c("April01", "May15", "Sep1", "Oct15"). Overrides
#'     'rep_period' parameter.
#'   rep_interval_custom = a vector of POSIXct objects that set the start of the
#'     summary periods. If sim_start precedes the first value, the first value
#'     will go until from start_sim to rep_interval_custom. Last interval will
#'     go from the last date in rep_interval_custom to end of sim_period or
#'     end_sim. Overrides 'rep_period' or 'rep_interval' parameters.
#'
#' @note Either 'sim_period' or 'sim_end' must be set, and either 'rep_period',
#' 'rep_interval', or 'rep_interval_custom' must be set. For 'rep_interval' the
#' dates can be set as numeric for both day and month, but there is a
#' possibility that the day and month will be inverted because dates with both
#' values <12 can be confused (e.g. 02-10 can be either February 10th or October
#' 2nd. Date format: "10Feb" or "02Oct" is preferred.

CreateReportIntervals <- function(sim = sim) {
  suppressPackageStartupMessages(require(lubridate))
  sim_start <- sim$pars$global$sim_start
  sim_period <- sim$pars$global$sim_period
  sim_end <- sim$pars$global$sim_end
  rep_period <- sim$pars$global$rep_period
  rep_interval <- sim$pars$global$rep_interval
  rep_interval_custom <- sim$pars$global$rep_interval_custom
  if (is.null(sim_period) && is.null(sim_end)) {
    stop("must provide either 'sim_period' or 'sim_end' parameter value")
  }
  if (!is.null(sim_period))  sim_end <- sim_start + sim_period
  if (!is.null(rep_period) && !is.null(rep_interval)){
    print("'rep_period' used instead of 'rep_interval'")
  }
  if (!is.null(rep_period)) {
    rep_period_unit <- ExtractUnitFromPeriod(rep_period)
    first_int <- as.interval(sim_start, ceiling_date(sim_start+period(1,
      "second"), unit=rep_period_unit))
    rep_intervals <- sim_start
    end_int <- int_end(first_int)
    rep_intervals <- with_tz(append(rep_intervals, end_int), tz(sim_start))
    while (end_int < sim_end) {
      end_int <- end_int + rep_period
      rep_intervals <- with_tz(append(rep_intervals, end_int), tz(sim_start))
    }
    if (tail(rep_intervals, 1) > sim_end) {
      rep_intervals[length(rep_intervals)] <- sim_end
    }
    rep_interval <- NULL
  }
  if (!is.null(rep_interval)){
      if (is.data.frame(rep_interval)) rep_interval <- rep_interval[, 1]
      if (max(names(guess_formats(rep_interval, c("md", "dm")))) == "md") {
        md_format <- max(guess_formats(rep_interval, c("md", "dm")))
        rep_interval <- as.character(as.Date(rep_interval, format = md_format),
          "%d%b")  # date must be in day month order
      }
      dates <- dmy(paste0(rep_interval, 2000)) # creates POSIXct
      rep_interval  <- rep_interval[order(dates)]  # orders rep_interval dates
      first_rep_int <- FindFirstReportInterval(sim_start, rep_interval)
      end_int <- int_end(first_rep_int)
      rep_intervals <- with_tz(append(sim_start, end_int), tz(sim_start))
      while (end_int < sim_end) {
        end_int_jul <- yday(end_int)
        rep_int_jul <- yday(dmy(paste0(rep_interval, year(end_int))))
        rep_int_position <- match(end_int_jul, rep_int_jul)
        if (rep_int_position == length(rep_interval)) {
          next_int <- interval(end_int, dmy(paste0(rep_interval[1],
            (year(end_int)))) + period(1, "year"))
        } else {
          next_int <- interval(end_int, dmy(paste0(rep_interval
            [rep_int_position+1], year(end_int))))
        }
        end_int <- int_end(next_int)
        rep_intervals <- with_tz(append(rep_intervals, end_int), tz(sim_start))
      }
      if (tail(rep_intervals, 1) > sim_end) {
        rep_intervals[length(rep_intervals)] <- sim_end
      }
      rep_interval_custom <- NULL
  }
  if(!is.null(rep_interval_custom)){
    if (is.data.frame(rep_interval)) rep_interval <- rep_interval_custom[,1]
    rep_interval <- rep_interval_custom # needed if rep_interval_custom != df
    rep_intervals <- as.POSIXct(rep_interval, tz=tz(sim_start))
    if (tail(rep_intervals, 1) > sim_end && !is.null(rep_period)) {
      rep_intervals[length(rep_intervals)] <- sim_end
      warning("'sum_inverval_custom' exceeds 'sim_period' time length")
    }
    if (tail(rep_intervals, 1) > sim_end && is.null(rep_period)) {
      rep_intervals[length(rep_intervals)] <- sim_end
      warning("'sum_inverval_custom' exceeds 'sim_end' time length")
    }
    if (tail(rep_intervals, 1) < sim_end) {
      rep_intervals[length(rep_intervals)+1] <- sim_end
    }
  }
  rep_intervals_list <- list()
  for (i in 1:(length(rep_intervals)-1)) {
    interval <- as.interval(rep_intervals[i], rep_intervals[i+1])
    rep_intervals_list[[length(rep_intervals_list)+1]] <- interval
  }
  return(rep_intervals_list)
}

#' CreateStepIntervals
#'
#' Create step intervals within a summary interval
#'
#' @usage CreateStepIntervals(rep_interval, step_period)
#'
#' @param rep_interval = an Interval object from lubridate that represents the
#' current summary interval
#' @param step_period = a Period object specifying the length of the step
#' period. The step_period should be of equal or greater length than the
#' time_step_period.
#'
#' @return a list of intervals
#' @export
#'
CreateStepIntervals <- function(rep_interval = rep_interval,
                                step_period = sim$pars$global$step_period) {
  suppressPackageStartupMessages(require(lubridate))
  step_period <- step_period
  step_intervals <- list()
  interval_counter <- 1
  current_start <- int_start(rep_interval)
  current_end <- (current_start + step_period)
  stop_point <- int_end(rep_interval)
  while(current_start < (stop_point)) {
    current_end <- (current_start+step_period)
    step_intervals[[interval_counter]] <- interval(current_start,
      current_end)
    interval_counter <- interval_counter + 1
    current_start <- current_start + step_period
    }
  if (int_end(step_intervals[[length(step_intervals)]]) > stop_point) {
    step_intervals[[length(step_intervals)]] <-
      interval(int_start(step_intervals[[length(step_intervals)]]), stop_point)
    }
  return(step_intervals)
}

#' CreateTimeSteps
#'
#' Create the time steps within a step interval
#'
#' @usage CreateTimeSteps(step_interval)
#'
#' @param step_interval = an Interval object from lubridate that represents the
#' current step interval
#' @param time_step_period = = a Period object specifying the length of the time
#' step period. Note that if the length of the current_step_interval is equal to
#' the time_step_period the function will return a POSIXct object with the same
#' date as the end of the current_step_interval.
#'
#' @return a list of intervals
#' @export
#'
CreateTimeSteps <- function(step_interval = step_interval,
                            time_step_period =
                              sim$pars$global$time_step_period) {
  options(lubridate.verbose=FALSE)
  step_interval_end <- int_end(step_interval)
  steps <- int_start(step_interval) # creates output list object w/start time
  end_int <- int_end(as.interval(time_step_period, int_start(step_interval)))
  steps <- with_tz(append(steps, end_int), tz(int_end(step_interval))) # 2nd t
  while (end_int < step_interval_end) {
    end_int <- end_int + time_step_period
    steps <- with_tz(append(steps, end_int), tz(int_end(step_interval)))
  }
  if (tail(steps, 1) > step_interval_end) {
    steps[length(steps)] <- step_interval_end
  }
  time_step_list <- list()
  for (i in 1:(length(steps)-1)) {
    interval <- as.interval(steps[i], steps[i+1])
    time_step_list[[length(time_step_list)+1]] <- interval
  }
  return(time_step_list)
}

#' CustomAgentsReportData
#'
#' Creates custom report data
#'
#' @param agent_states = agent_states
#' @param report_data = report_data
#' @param step_data = step_data
#' @param q = counter for reporting interval
#'
#' @return 'report_data' object
#' @export
#'
CustomAgentsReportData <- function(agent_states,
                                   report_data,
                                   step_data,
                                   q) {
    report_data <- report_data
#    report_data[q, "sex"] <- agent_states$sex
#    report_data[q, "move_length"] <- sum(step_data$step_length, na.rm=TRUE)
    return(report_data)
}

#' CustomPopReportData
#'
#' Creates and updates custom pop_report data.
#'
#' @usage CustomPopReportData(report_data, q, agent_states, step_data)
#'
#' @param pop_report = population report
#' @param report_data =  report_data
#' @param q = report_data row for new data
#'
#' @return a population report
#' @export
#'

CustomPopReportData <- function(pop_report,
                                report_data,
                                q) {
  pop_report <- pop_report
#  pop_report[q, "mf_ratio"] <- nrow(subset(report_data, sex=="male")) /
#    nrow(subset(report_data, sex=="female"))
  return(pop_report)
}

#' ExtractUnitFromPeriod
#'
#' Helper function for CreateSumInterval(). Finds unit of Period object.
#'
#' @usage ExtractUnitFromPeriod(period)
#'
#' @param period = period object
#'
#' @return period object's unit as a chararcter
#' @export
#'
ExtractUnitFromPeriod <- function(period) {
  period <- period
  if (period$year > 0) unit <- "year"
  if (period$month > 0) unit <- "month"
  if (period$day == 7 ) unit <- "week"
  if (period$day > 0 && period@day < 7) unit <- "day"
  if (period$day > 7) unit <- "day"
  if (period$hour > 0) unit <- "hour"
  if (period$minute > 0) unit <- "minute"
  if (period$.Data > 0) unit <- "second"
  return(unit)
}

#' FindFirstReportInterval
#'
#' Helper function for CreateSumInterval. Finds first interval based on
#' sim_start datetime and rep_interval dates
#'
#' @usage FindFirstReportInterval(sim_start, rep_interval)
#'
#' @param sim_start  = a POSIXct object for simulation's start, required
#' @param rep_interval = vector or dateframe (with date in first column) of the
#'   format %B%d, %b%d, %d%B, or %B%d (see ?strptime for more details on
#'   formatting). Example: c("April01", "May15", "Sep1", "Oct15"). Required.
#'
#' @return an Interval object
#' @export

FindFirstReportInterval <- function(sim_start,
                                    rep_interval){
  require(lubridate)
  start_year <- as.Date(0, origin=as.Date(floor_date(sim_start, "year")))
  if (length(rep_interval) == 1) {
    first_int <- as.interval(period(1, "year"), dmy(paste0(rep_interval[1],
      year(start_year))))
    if (sim_start %within% first_int) {
      first_rep_int <- interval(sim_start, dmy(paste0(rep_interval[1],
        year(start_year + period(1, "year")))))
    } else {
      first_rep_int <- interval(sim_start, dmy(paste0(rep_interval[1],
        year(start_year))))
    }
  }
  if (length(rep_interval) == 2) {
    first_int <- as.interval(dmy(paste0("0101", year(start_year))),
      dmy(paste0(rep_interval[1], year(start_year))))
    if ((sim_start + period(1, "second")) %within% first_int) {
      first_rep_int <- interval(sim_start, dmy(paste0(rep_interval[1],
        year(start_year))))
    }
    if (exists("first_rep_int") == FALSE) {
      second_int <- as.interval(dmy(paste0(rep_interval[1], year(start_year))),
        dmy(paste0(rep_interval[2], year(start_year))))
      if ((sim_start + period(1, "second")) %within% second_int) {
        first_rep_int <- interval(sim_start, dmy(paste0(rep_interval[2],
          year(start_year))))
      } else {
        first_rep_int <- interval(sim_start, dmy(paste0(rep_interval[1],
          year(start_year+period(1, "year")))))
      }
    }
  }
  if (length(rep_interval) > 2) {
    first_int <- as.interval(dmy(paste0("0101", year(start_year))),
        dmy(paste0(rep_interval[1], year(start_year))))
    if ((sim_start + period(1, "second")) %within% first_int) {
      first_rep_int <- interval(sim_start, dmy(paste0(rep_interval[1],
        year(start_year))))
    }
    for (i in 2:(length(rep_interval)-1)) {
      mid_int <- as.interval(dmy(paste0(rep_interval[i], year(start_year))),
        dmy(paste0(rep_interval[i+1], year(start_year))))
      if ((sim_start + period(1, "second")) %within% mid_int) {
        first_rep_int <- interval(sim_start, dmy(paste0(rep_interval[i+1],
          year(start_year))))
      }
    }
    for (i in length(rep_interval)) {
      last_int <- as.interval(dmy(paste0(rep_interval[i], year(start_year))),
        dmy(paste0("0101", (year(start_year)+1) )))
      if ((sim_start + period(1, "second")) %within% last_int) {
        first_rep_int <- interval(sim_start, dmy(paste0(rep_interval[1],
          (year(start_year)+1))))
      }
    }
  }
  return(first_rep_int)
}

#' FindSeasonFromDatetime
#'
#' Determines which season (from sim$pars$global$sim_seasons) a datetime is
#' within
#'
#' @usage FindSeasonFromDatetime(datetime, season)
#'
#' @param datetime = a POSIXct object
#' @param seasons = a dataframe of starting dates ("start" column) and season
#' names ("season" column), usually located at: sim$pars$global$sim_seasons
#'
#' @return  an atomic character object
#' @export
#'
#'
FindSeasonFromDatetime <- function(datetime = datetime,
                                   seasons = sim$pars$global$sim_seasons) {
  require(lubridate)
  start_year <- as.Date(0, origin=as.Date(floor_date(datetime, "year")))
  if (max(names(guess_formats(seasons[, "start"], c("md", "dm")))) == "md") {
    md_format <- max(guess_formats(seasons[, "start"], c("md", "dm")))
    seasons[, "start"] <- as.character(as.Date(seasons[, "start"], format =
      md_format), "%d%b")  # date must be in day month order
  }
  dates <- dmy(paste0(seasons[,"start"], 2000)) # creates POSIXct
  seasons  <- seasons[with(seasons, order(dates)), ]  # orders seasons by dates
  if (nrow(seasons) == 2) {
    first_season <- interval(dmy(paste0(paste0(seasons[1, "start"],
      year(start_year)))), dmy(paste0(seasons[2, "start"], year(start_year)))-1)
    if (datetime %within% first_season) {
      season <- seasons[1, "season"]
    } else {
      season <- seasons[2, "season"]
    }
  }
  if (nrow(seasons) > 2) {
    last_season <- TRUE
    for (i in 1:(nrow(seasons)-1)) {
      mid_season <- interval(dmy(paste0(seasons[i, "start"], year(
        start_year))), dmy(paste0(seasons[i+1, "start"], year(start_year)))-1)
      if (datetime %within% mid_season) {
        season <- seasons[i, "season"]
        last_season = FALSE
      }
    }
    if(last_season == TRUE) season <- seasons[nrow(seasons), "season"]
  }
  return(season)
}

#' NamedList
#'
#' Creates a list with named objects, and it tries not to replace any
#' already-named arguments.
#'
#' @usage NamedList(...)
#'
#' @param ... = objects to add to list
#'
#' @return  a list
#' @export
#'
#' @note Original function came from Ben Bolker's answer on StackOverflow:
#' "http://stackoverflow.com/questions/16951080"
#'
NamedList <- function(...) {
  list <- list(...)
  str_name <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(name <- names(list))) name <- str_name
  if (any(no_names <- name == "")) name[no_names] <- str_name[no_names]
    setNames(list, name)
}

#' RemoveExcept
#'
#' Helper function that removes all the objects in the environment except ones
#' listed in argument
#'
#' @usage RemoveExcept(object)
#'
#' @param object = objects to keep
#'
#' @return  environment with only the kept objects
#' @export
#'
RemoveExcept <- function(object = object){
  if (length(setdiff(ls(pos = .GlobalEnv), object)) > 0) {
    rm(list=setdiff(ls(pos = .GlobalEnv), object), pos = .GlobalEnv)
  }
}

#' ReproductionSubModel
#'
#' Reproduction submodel
#'
#' @param agent_states = a 'agent_states' list object
#' @param step_data = list object
#'
#' @return 'agent_states' list object
#' @export
#'
#' @note Just a placeholder for now.
#'
ReproductionSubModel <- function(agent_states = agent_states,
                                 step_data = step_data) {
  if(dummy){

  } else {


  }
}

#' ReturnActiveSeq
#'
#' Create vector of positions of agents that had a step_data$datetime within
#' the current reporting interval
#'
#' @param all = 'all' list
#' @param rep_interval = reporting interval
#'
#' @return 'active_seq' vector
#' @export
#'
ReturnActiveSeq <- function(all = sim$agents$all,
                            rep_interval = rep_interval){
  all <- all
  active_seq <- vector()
  for (i in 1:length(all)) {
    step_data <- all[[i]]$step_data
#    step_data <- sim$agents$all[[i]]$step_data
#    print(paste("The step_data:", step_data))
#    print(paste("The active_seq:", rep_interval))
    if (any(step_data$datetime %within% rep_interval)) {
      active_seq <- alive_seq <- append(active_seq, i)
    }
  }
  return(active_seq)
}

#' ReturnAliveIds
#'
#' Create vector of ids of alive agents
#'
#' @param all = 'all' list
#'
#' @return a vector
#' @export
#'
ReturnAliveIds <- function(all = sim$agents$all){
  all <- all
  alive_ids <- vector()
  for (i in 1:length(all)) {
    agent <- all[[i]]
    if (is.na(agent$states$died)) alive_ids <- append(alive_ids,agent$states$id)
  }
  return(alive_ids)
}

#' ReturnAliveSeq
#'
#' Create vector of alive agents positions
#'
#' @param sim = 'sim' list
#'
#' @return a vector
#' @export
#'

ReturnAliveSeq <- function(sim = sim){
  all <- sim$agents$all
  alive_seq <- vector()
  for (i in 1:length(all)) {
    agent <- all[[i]]
    if (is.na(agent$states$died)) alive_seq <- append(alive_seq, i)
  }
  return(alive_seq)
}

#' ReturnPopReportData
#'
#' Helper function for UpdatePopReport(). Finds intervals with start dates
#' within the report interval.
#'
#' @usage ReturnPopReportData(sim, agents, rep_interval)
#'
#' @param sim = 'sim' list
#' @param agents = agents object
#' @param rep_interval = counter for reporting interval (q)
#'
#' @return period object's unit as a chararcter
#' @export
#'

ReturnPopReportData <- function(sim=sim,
                                agents = agents,
                                rep_interval=rep_interval) {
  sim = sim #
  rep_interval = rep_interval #
  report_data = agents$agents_report
#  print(paste0("Inside PopReportData ", report_data))
#  print(paste0("RepInterval ", rep_interval))
  pop_report_data <- data.frame()
  for (s in 1:nrow(report_data)) {
    if (as.POSIXct(report_data[s, "start"], tz = tz(sim$pars$global$sim_start))
      %within% rep_interval) {
        pop_report_data <- rbind(pop_report_data, report_data[s, ])
    }
  }
  return(pop_report_data)
}

#' ReturnReportStepData
#'
#' Returns agent's step data within the current reporting interval
#'
#' @param agent agent = 'agent' list
#' @param rep_interval = counter for reporting interval (q)
#'
#' @return a step dataframe
#' @export
#'
ReturnReportStepData <- function(agent=agent,
                                 rep_interval = rep_interval) {
  step_data = agent$step_data
  report_step_data <- data.frame()
  for (s in 1:nrow(step_data)) {
    if (step_data[s, "datetime"] %within% rep_interval) report_step_data<-
      rbind(report_step_data, step_data[s, ])
  }
  return(report_step_data)
}

#' SurvivalSubModel
#'
#' Survival submodel
#'
#' @param agent_states = 'agent_states' list object
#' @param step_data = list object
#'
#' @return agent_states
#' @export
#'
#' @note Just a placeholder for now.

SurvivalSubModel <- function(agent_states = agent_states,
                             step_data = step_data) {
  agent_states <- agent_states
  return(agent_states)
}

#' UpdateAgentParsData
#'
#' Creates and updates agents parameters dataframe
#'
#' @usage CreateAgentParsData(sim, init)
#'
#' @param sim = sim list, MUST be included in function's parameter argument
#' needed when: when init == TRUE
#' @param init = logical, whether or not this is the initation step
#'
#' @return 'sim' list
#' @export
#'
UpdateAgentParsData <- function(sim = sim,
                                init = FALSE){
  sim=sim
  if (init == TRUE) {
    input <- sim$agents$input
    all <- sim$agents$all
    for (i in 1:nrow(input)){
      agent <- all[[i]]
      pars_data <- data.frame(id=agent$states$id, rep_int = NA)
      all[[i]] <- append(agent, NamedList(pars_data))
    }
    sim$agents$all <- all
    return(sim)
  } else {
    pars_data <- sim$agents$all
    # placeholder for updates
    sim$agents$all <- pars_data
    return(sim)
  }
}

#' UpdateAgentsReport
#'
#' Creates and updates agents interval data dataframe
#'
#' @usage CreateAgentInData(agents, init)
#'
#' @param sim = 'sim' list
#' @param rep_interval = rep_interval
#' @param step_intervals = step_intervals
#'
#' @return an 'agents' list
#' @export
#'
#' @note rep_data, init are internal to RunSimualtion()
UpdateAgentsReport <- function(sim = sim,
                               rep_interval = rep_interval,
                               step_intervals = step_intervals,
                               custom_agent_report = FALSE) {
  agents <- sim$agents
  alive_seq  <- ReturnActiveSeq(all=sim$agents$all, rep_interval=rep_interval)
  for (p in alive_seq) {
    if(!("report_data" %in% names(agents$all[[p]]))) {
      report_data <- data.frame()
    } else {
      report_data <- agents$all[[p]]$report_data
    }
    agent <- agents$all[[p]]
    q <- nrow(report_data) + 1
    agent_states <- agent$states
    step_data <- ReturnReportStepData(agent, rep_interval)
    pars_data <- agent$pars
    report_data[q, "start"] <- as.character(int_start(step_intervals[[1]]))
    report_data[q, "end"] <- as.character(int_end(step_intervals[[length(
      step_intervals)]]))
    report_data[q, "id"] <- agent_states$id
    report_data[q, "age"] <- CalculateReportAge(agent_states, step_data)
    if (custom_agent_report == TRUE) {
      report_data <- CustomAgentsReportData(agent_states, report_data,
        step_data, q)
    }
    agents$all[[p]][["report_data"]] <- report_data
  }
  if(!("agents_report" %in% names(agents))) {
    agents_report <- data.frame()
    print("Created agents_report")
  } else {
    agents_report <- sim$agents$agents_report
  }
  active_seq <- ReturnActiveSeq(all=sim$agents$all, rep_interval=rep_interval)
  for (r in active_seq) {
    report_data <- agents$all[[r]]$report_data
    agents_report  <- rbind(agents_report, report_data[nrow(report_data), ])
  }
  print("Updated agents_report")
  agents[["agents_report"]] <- agents_report
   return(agents)
}

#' UpdateAgentStates
#'
#' Creates and updates agents list from the sim object
#'
#' @param agent_states a list from sim$agents$all$agent$states
#' @param sim sim list, MUST be included in function's parameter arguments
#' needed when: when init == TRUE
#' @param init logical, whether or not this is the initation step
#'
#' @return an 'agents' list
#' @export
#'
UpdateAgentStates <- function(agent_states = NULL,
                              sim = sim,
                              init = FALSE) {
  if (init == TRUE) {
    input <- sim$agents$input
    input <- CreateBirthDate(input)
    input_columns <- colnames(input)
    na_columns <- c("start_datetime", "died")
    all <- list()
    for (i in 1:nrow(input)) {
      states <- list()
      for (j in input_columns) states <- append(states, input[i, j])
      for (k in 1:length(na_columns)) states <- append(states, NA)
      states <- setNames(states, c(input_columns, na_columns))
      agent <- NamedList(states)
      all <- append(all, NamedList(agent))
    }
    sim$agents <- append(sim$agents, NamedList(all))
    return(sim)
  } else {
    agent_states <- agent_states
    return(agent_states)
  }
}

#' UpdateAgentStepData
#'
#' Creates and updates agents step data database
#'
#' @usage CreateAgentStepData(agents, sim, init)
#'
#' @param step_data = 'step data' dataframe
#' @param sim = 'sim' object
#' @param init = logical, whether or not this is the initation step
#'
#' @return an 'agents' list
#' @export
#'
UpdateAgentStepData <- function(step_data = NULL,
                                sim = sim,
                                init = FALSE) {
  if (init == TRUE) {
    sim_start <- sim$pars$global$sim_start
    all <- sim$agents$all
    for (i in 1:length(all)) {
      agent <- all[[i]]
      step_data <- data.frame(id=agent$states$id, datetime=sim_start,
        x=agent$states$start_x, y=agent$states$start_y, exp_angle=NA,
        abs_angle=NA)
      agent  <-  append(agent, NamedList(step_data))
      all[[i]] <- agent
    }
    sim$agents$all <- all
    return(sim)
  } else {
    step_data <- step_data
    return(step_data)
  }
}

#' UpdatePopReport
#'
#' Creates and updates population interval dataframe
#'
#' @usage UpdatePopReportData(sim, rep_interval, step_intervals)
#'
#' @param sim  'sim' list
#' @param rep_interval rep_interval
#' @param step_intervals  step_intervals
#'
#' @return a 'pop_report' dataframe
#' @export
#'
UpdatePopReport <- function(sim = sim,
                            rep_interval=rep_interval,
                            step_intervals=step_intervals,
                            custom_pop_report = FALSE) {
  agents <- sim$agents
  rep_interval = rep_interval
  if(!("pop_report" %in% names(agents))) {
    pop_report <- data.frame()
    print("Created pop_report")
  } else {
    pop_report <- sim$agents$pop_report
  }
  q <- nrow(pop_report) + 1
  report_data <- ReturnPopReportData(sim=sim, agents=agents, rep_interval=rep_interval)
  pop_report[q, "start"] <- as.character(int_start(step_intervals[[1]]))
  pop_report[q, "end"] <- as.character(int_end(step_intervals[[length(
    step_intervals)]]))
  pop_report[q, "total_n"] <- length(unique(report_data$id))
  if (custom_pop_report == TRUE) {
    pop_report <- CustomPopReportData(pop_report, report_data, q)
  }
  agents[["pop_report"]] <- pop_report
  return(agents)
}

#' UpdateSpatial
#'
#' Creates and updates population interval dataframe
#'
#' @usage UpdateSpatial(spatial, init)
#'
#' @param sim = sim list, MUST be included in function's parameter argument
#' needed when: when init == TRUE
#' @param init = logical, whether or not this is the initation step
#'
#' @return an 'agents' list
#' @export
#'
UpdateSpatial <- function(sim = sim,
                          init = TRUE) {
  sim <- sim
  if (init == TRUE) {
    spatial <- sim$spatial
    sim$spatial <- spatial
    return(sim)
  } else {
#   spatial_timer <- UpdateSpatialTimer(spatial_timer)
#   if (spatial_timer == timer_number) UpdateSpatial(); rm(spatial_timer)
#   spatial <- append()
    spatial <- sim$spatial
    sim$spatial <- spatial
    return(sim)
  }
}

#' WriteSimList
#'
#' Writes "sim" list to a directory,
#'
#' @usage WriteSimList(write, run, sim, output_dir, components)
#'
#' @param write  = logical, whether or not to write the sim list to a file.
#' Default is TRUE.
#' @param run = name of run. Default uses the name from the runs object to
#' ensure that the proper number of leading zeros is used. If default is not
#' use, the name will simply be: "sim_(run).RData"
#' @param sim = a 'sim' list object
#' @param output_dir = output directory, default is working directory
#' @param components = components of 'sim' to write. Options include: "agents",
#' "pars", and "spatial". Default is all.
#'
#' @return Writes a file to the output_dir
#' @export
#'
#' @note NOT FINISHED - Code not implemented for writing subset of components
WriteSimList <- function(write = TRUE,
                         run = names(runs[j]),
                         sim = sim,
                         output_dir = getwd(),
                         components = "all") {
  sim = sim
  if (write == TRUE) {
    if (components == "agents") {
      file_path = file.path(output_dir, paste0("sim_",run,".RData"))
      print(paste0("Writing: ", file_path, " - agents only"))
      save(sim$agents, file = file_path)
    }
    if (components == "all") {
      file_path = file.path(output_dir, paste0("sim_",run,".RData"))
      print(paste0("Writing: ", file_path))
      save(sim, file = file_path)
    }
  }
}

#' TestMe
#' This is to test the order things load
#'
#' @return print statement
#' @export

TestMe <- function(){
  print("Testing from sim.R. Double-check")
}
