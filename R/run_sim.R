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
            sim$agents$all[[n]][["states"]] <- UpdateAgentStates(agent_states)
            sim$agents$all[[n]][["step_data"]] <- UpdateAgentStepData(step_data)
          } # end of alive_seq[[n]]
          print(paste("end of time_step:", time_steps[[m]]))
          # Hatchling/Dispersal Submodel
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
