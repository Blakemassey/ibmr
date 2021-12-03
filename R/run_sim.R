###  Creates and updates agents list from the sim object
###  Usage: CreateAgentStates(agent_states, sim, init)
###  Arguments: agent_states = a list from sim$agents$all$agent$states
###             sim = sim list, MUST be included in function's parameter
###               arguments needed when: when init == TRUE
###             init = logical, whether or not this is the initation step
###             current_agent = numeric, the id (a numeric value) of the
###               current agent
###             add_agents = logical, whether or not new agents should be added
###               to 'agents'
###             new_intput = not implemented yet, currently set to NULL
###  Returns: an 'agents' list
###  Notes: New columns can be via the "add_columns" vector w. The arguments
###     'current_agent' and 'add_agents' are currently designed for use in the
###     indigo snake IBM where new agents (i.e., new births) are being
###     continuously added to the model.
###  Blake Massey and Javan Bauder
###  2015.03.24
UpdateAgentStates <- function(agent_states = NULL,
                              sim = sim,
                              init = FALSE,
                              current_agent = NULL,
                              add_agents = FALSE,
                              new_input = NULL) {
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
    clutch_data <- data.frame(nest_id = numeric(), parent_id = character(),
      date = character(), x = numeric(), y = numeric(), offspring = numeric())
    sim$agents <- append(sim$agents, NamedList(all))
    sim$agents <- append(sim$agents, NamedList(clutch_data))
    return(sim)
  } else {
    if(add_agents == TRUE){
      new_input$id <- seq((max(sim$agents$input$id)+1),
        (max(sim$agents$input$id) + nrow(new_input)))
      input_columns <- colnames(new_input)
      na_columns <- c("died")
      all <- list()
      for (i in 1:nrow(new_input)) {
        states <- list()
        for (j in input_columns) states <- append(states, new_input[i, j])
        for (k in 1:length(na_columns)) states <- append(states, NA)
        states <- setNames(states, c(input_columns, na_columns))
        #states$start_datetime <- as.character(new_input[i,"start_date"])
        agent <- NamedList(states)
        all <- append(all, NamedList(agent))
      }
      sim$agents$all <- append(sim$agents$all, all)
      if(ncol(sim$agents$input)<ncol(new_input)){
        new_cols <- colnames(new_input)[-(which(colnames(new_input) %in%
            colnames(sim$agents$input)))]
        new_data <- data.frame(matrix(nrow = nrow(sim$agents$input),
          ncol = length(new_cols)))
        colnames(new_data) <- new_cols
        sim$agents$input <- cbind(sim$agents$input, new_data)
      }
      sim$agents$input <- rbind(sim$agents$input, new_input)
      return(sim)
    } else {
      current_step_data <- sim$agents$all[[current_agent]]$step_data
      if(nrow(current_step_data) == sim$pars$global$burnin){

        total_center <- colMeans(sim$agents$all[[current_agent]]$step_data[,
          c("x","y")],na.rm=T)
        agent_states <- append(agent_states, NamedList(total_center))
        return(agent_states)
      } else {
        agent_states <- agent_states
        return(agent_states)
      }
    }
  }
}

###  Creates and updates agents step data database
###  Usage: CreateAgentStepData(agents, sim, init)
###  Arguments: agents = agents object
###             sim = sim list, MUST be included in function's parameter
###               arguments needed when: when init == TRUE
###             init = logical, whether or not this is the initation step
###             step_attributes = character vector of column names of step-level
###               attributes to be recorded (e.g., c("exp_anlge","abs_angle"))
###  Returns: an 'agents' list
###  Notes:
###  Blake Massey
###  2015.03.27
###  Modified by Javan Bauder, 2016.09.05
UpdateAgentStepData <- function(step_data = NULL,
                                sim = sim,
                                init = FALSE,
                                step_attributes = NULL,
                                add_roads = FALSE,
                                add_centers = FALSE) {
  if (init == TRUE) {
    #sim_start <- sim$pars$global$sim_start
    all <- sim$agents$all
    for (i in 1:length(all)) {
      agent <- all[[i]]
      if(any(names(agent) == "step_data")){
        next
      } else {
        if(is.na(agent$states$start_datetime)){
          start_datetime <- sim$pars$global$sim_start
        } else {
          start_datetime <- as.POSIXct(agent$states$start_datetime, tz = "UTC")
        }
        step_data <- data.frame(id = agent$states$id,
          class = agent$states$class, datetime = start_datetime,
          x = agent$states$start_x, y = agent$states$start_y)
        step_data_col <- ncol(step_data)
        if(!is.null(step_attributes)){
          step_data <- cbind(step_data, matrix(ncol = length(step_attributes),
            nrow = nrow(step_data), data = NA))
          colnames(step_data)[-(1:step_data_col)] <- step_attributes
        }

        agent  <-  append(agent, NamedList(step_data))
        if(add_roads==TRUE){
          Roads <- list()
          agent  <-  append(agent, NamedList(Roads))
        }
        if(add_centers==TRUE){
          Home_range_centers <- data.frame(datetime = sim$pars$global$sim_start)
          agent <- append(agent, NamedList(Home_range_centers))
        }

        all[[i]] <- agent
      }
    }
    sim$agents$all <- all
    return(sim)
  } else {
    step_data <- step_data
    return(step_data)
  }
}

### Writes "sim" list to a directory,
### Usage: WriteSimList(write, sim, output_dir, components)
### Arguments: write = logical, whether or not to write the sim list to a file.
###              Default is TRUE.
###            run = name of run. Default uses the name from the runs object to
###              ensure that the proper number of leading zeros is used. If
###              default is not use, the name will simply be: "sim_(run).RData"
###            sim = a 'sim' list object
###            output_dir = output directory, default is working directory
###            components = components of 'sim' to write. Options include:
###              "agents", "pars", and "spatial". Default is all.
###            file.type = character, use ".rds" to write sim as a RDS file
###              otherwise sim is written as a .Rdata file.
### Returns: Writes a file to the output_dir
### Notes: NOT FINISHED - Code not implemented for writing subset of components
### Blake Massey and Javan Bauder
### 2015.03.29
WriteSimList <- function(write = TRUE,
                         run = names(runs[j]),
                         sim1 = sim,
                         output_dir = getwd(),
                         components = "all",
                         file.type = "rds") {
  sim = sim1
  if (write == TRUE) {
    if (components == "agents") {

      if(file.type=="rds"){
        file_path = file.path(output_dir, paste0("sim_",run,".rds"))
        print(paste0("Writing: ", file_path, " - agents only"))
        saveRDS(sim$agents, file = file_path)
      } else {
        file_path = file.path(output_dir, paste0("sim_",run,".RData"))
        print(paste0("Writing: ", file_path, " - agents only"))
        save(sim$agents, file = file_path)
      }
    }
    if (components == "all") {
      if(file.type=="rds"){
        file_path = file.path(output_dir, paste0("sim_",run,".rds"))
        print(paste0("Writing: ", file_path))
        saveRDS(sim, file = file_path)
      } else {
        file_path = file.path(output_dir, paste0("sim_",run,".RData"))
        print(paste0("Writing: ", file_path))
        save(sim, file = file_path)
      }
    }
  }
}

###  Combines multiple raster objects into a single, normalized redistribution
###    kernel.
###  Usage: CreateRedistributionKernel(raster_stack, combine_method,
###    adjust_distance, current_location, base)
###  Arguments: raster_stack = a RasterStack object containing all raster
###               objects to be combined
###             combine method = character. "product" takes the product of all
###               rasters within the RasterStack. "geometric" takes the
###               geometric mean of all rasters within the raster stack.
###             adjust_distance = logical. indicates whether or not the final
###               redistribution kernel should be adjusted based on expected
###               cell counts.
###             current_location = numeric. the x and y coordinates of current
###               location.
###             base = base Raster that sets the projection, extent, and
###               dimensions of the study area
###  Returns: a raster object
###  Notes:
###  Javan Bauder
###  2017.04.05
CreateRedistributionKernel <- function(raster_stack,
                                       combine_method,
                                       adjust_distance = TRUE,
                                       current_location,
                                       base){
  cellsize <- res(base)[1]
  if(sim$pars$global$combine_method == "product"){
    prob_raster <- raster::overlay(raster_stack, fun=function(a,b,c,d,e)
      {return(a*b*c*d*e)}, recycle = FALSE, unstack = TRUE)
  } else {
    prob_raster <- raster::overlay(raster_stack, fun = function(a,b,c,d,e)
      {return((a*b*c*d*e)^(1/5))}, recycle = FALSE, unstack = TRUE)
  }
  if(adjust_distance==TRUE){
    center <- cellFromXY(prob_raster, matrix(c(current_location[1],
      current_location[2]), nrow = 1))
    distance_raster <- prob_raster
    distance_raster[center] <- 999
    distance_raster[distance_raster!=999] <- NA
    distance_raster <- distance(distance_raster)
    no_cells <- 2*pi*distance_raster/cellsize
    adjusted_raster <- prob_raster/no_cells
    # adjusted_raster[center] <- NA
    # adjusted_raster[center] <- cellStats(adjusted_raster,"max")
  }
  adjusted_raster[cellFromXY(adjusted_raster, matrix(c(current_location[1],
    current_location[2]), nrow=1))] <- 0
  adjusted_raster <- adjusted_raster/cellStats(adjusted_raster, stat="sum")
  crs(adjusted_raster) <- crs(base)
  adjusted_raster[is.na(adjusted_raster)] <- 0
  return(adjusted_raster)
}

#' RunSimulation
#'
#' The backbone function for an individual-based model simulation.
#'
#' @usage RunSimulation(sim, runs, write, output_dir)
#'
#' @param sim = list of (agents, pars, and spatial)
#' @param runs = number of runs
#' @param write = write
#' @param output_dir = output files directory
#'
#' @return A list object
#' @export
#'
RunSimulation <- function(sim = sim,
                          runs = 1,
                          write = FALSE,
                          output_dir = getwd()) {
  runs <- CreateRunsList(runs)
  for (i in 1:length(runs)) {
    rep_intervals <- CreateReportIntervals(sim)
    sim <- UpdateAgentStates(init=TRUE, sim=sim)
    sim <- UpdateAgentStepData(init=TRUE, sim=sim)
    sim <- UpdateAgentParsData(init=TRUE, sim=sim)
    sim <- UpdateSpatial(init=TRUE, sim=sim)
    for (j in 1:length(rep_intervals)) {
      step_intervals <- CreateStepIntervals(rep_intervals[[j]])
      for (k in 1:length(step_intervals)) {
        time_steps <- CreateTimeSteps(step_intervals[[k]])
        for (m in 1:length(time_steps)) {
          step <- time_steps[[m]]
          alive_seq <- ReturnAliveSeq(sim)
          sim$agents$all <- UpdateAgentParsData(sim$agents$all)
          # potentially add if statement about alive_seq == 0 then don't run
          for (n in alive_seq) {
            print(paste("alive_seq:", n))
            agent_states <- sim$agents$all[[n]][["states"]]
            step_data <- sim$agents$all[[n]][["step_data"]]
            pars_data <- sim$agents$all[[n]][["pars_data"]]
            # START Submodels #
            agent_states <- AgingSubModel(agent_states, step_data, step)
            step_data <- MovementSubModel(sim, agent_states, step_data, step)
            agent_states <- SurvivalSubModel(agent_states, step_data)
            agent_states <- ReproductionSubModel(agent_states, step_data)
            # END Submodels #
            sim$agents$all[[n]][["step_data"]] <- UpdateAgentStepData(step_data)
            sim$agents$all[[n]][["states"]] <- UpdateAgentStates(agent_states)
          } # end of alive_seq[[n]]
          print(paste("end of time_step:", time_steps[[m]]))
          # Hatchling/Dispersal Submodel
          # add alive_seq
          sim$spatial <- UpdateSpatial(sim$spatial)
        } # end of time_steps[[m]]
        print(paste("end of step_interval:", step_intervals[[k]]))
      } # end of step_interval[[k]]
      sim$agents <- UpdateAgentsReport(sim, rep_intervals[[j]], step_intervals)
      sim$agents <- UpdatePopReport(sim, rep_intervals[[j]], step_intervals)
      print(paste("end of rep_interval:", rep_intervals[[j]]))
    } # end of rep_interval[[j]]
    runs[[i]] <- sim
    WriteSimList(write = write, run = names(runs[j]), sim = sim,
      output_dir = getwd(), components = "all")
  }
  return(runs)
}
