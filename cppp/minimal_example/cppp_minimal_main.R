library(nimble)

# load compiled model
model <- readRDS("cppp/minimal_example/saved_inputs/compiled_model.rds")

# load original samples
MCMCOutput <- readRDS("cppp/minimal_example/saved_inputs/samples.rds")

# load data
y <- readRDS("cppp/minimal_example/saved_inputs/y.rds")

dataNames <- "y"

# if *not* conditioning on latent state
paramNames <- c("p", "psi")
simNodes <- c("z", "y")

# if conditioning on latent state
paramNames <- c("p", "psi", "z")
simNodes <- "y"

discrepancyFunctions <- list(chisqDiscFunction)
discrepancyFunctionsArgs <- list(list(nVisits = nVisits))
