library(nimble)

# load compiled model
model <- readRDS("cppp/minimal_example/saved_inputs/compiled_model.rds")

# load original samples
samples <- readRDS("cppp/minimal_example/saved_inputs/samples.rds")

dataNames <- "y"
paramNames <- c("p", "psi")
simNodes <- c("z", "y")