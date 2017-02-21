############################## SPATIAL #########################################

base <- CreateBase(0, 0, 30000, 30000, 30)
spatial <- NamedList(base)

################################ AGENTS ########################################

input <- CreateAgentsInputClass(n=10, class="turtle", sex_ratio=.5, age_min=1,
  age_max=20, base=base)
agents <- NamedList(input)

################################# PARS #########################################

# Global
library(lubridate)
global <- CreateParsGlobal(
  sim_start = as.POSIXct("2015-01-01", tz = "UTC"),
  sim_period = period(25, "days"),
  sim_end = NULL,
  rep_period = period(1, "days"),
  rep_interval = c("01Jan", "01Jul"),
  rep_interval_custom = NULL,
  step_period = period(1, "day"),
  time_step_period = period(6, "hour"),
  birth_day = "01Jan",
  input_age_period = "years",
  report_age_period = "months",
  sim_seasons = NULL)

# Classes
# constant
constant <- CreateParsClassConstant(step_cauchy_mu = 0, step_cauchy_rho = .5)

# Season
spring <- CreateParsClassSeason(300, .2, 300)
summer <- CreateParsClassSeason(100, .1, 300)
fall <- CreateParsClassSeason(300, .2, 300)
winter <- CreateParsClassSeason(100, .1, 300)
season <- NamedList(spring, summer, fall, winter)

# Julian
x <- 1:365
y <- sin(3 * pi * x / 365) * - .003 + runif(length(x), .005, .010)
julian <- data.frame(day = x)
julian$home_return <- predict(loess(y ~ x, span=.75, data.frame(x=x, y=y)),
  data.frame(x=x))

male <- NamedList(constant, season, julian)
female <- NamedList(constant, season, julian)
classes <- NamedList(male, female)

pars <- NamedList(global, classes)

################################## SIM #########################################

sim <- NamedList(agents, pars, spatial)
RemoveExcept("sim")

############################# RUN SIMULATION ###################################
sim = sim
runs = 1
write = FALSE
output_dir = getwd()

runs <- RunSimulation(sim, runs, write, output_dir)
sim <- runs[[1]]

sim_step_data <- CompileAllAgentsStepData(sim = sim)

library(ggplot2)
library(plot3D)
ggplot(sim_step_data, aes(x, y))

ggplot() + geom_tile(data = as.data.frame(sim$spatial$base, xy=TRUE),
  aes(x=x, y=y, fill=layer)) +
  geom_point(data=sim_step_data, aes(x, y))






################################################################################
############################ OLD CODE ##########################################
################################################################################

sim <- NamedList(agents, pars, spatial)
RemoveExcept("sim")
save(sim, file="C:/Work/R/Data/Simulation/sim.RData")

plot(sim$spatial$homerange_kernels[[1]])
source('C:/Work/R/Functions/gen.R')
SavePlot("homerange_kernel_1.png", height = 5, width=5)

redist <- CreateRedistKernel(max_r = 300, cellsize = 30, mu=0, rho=.5, shape=300,
   scale=100, ignore_cauchy = FALSE, ignore_pareto = FALSE)

ggplot(redist) + geom_raster(aes())

require(reshape2); require(ggplot2)
dataL = melt(redist, id="x")

qplot(redist, x=Var1, y=value, data=dataL, group=Var2)