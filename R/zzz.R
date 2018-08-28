.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "medshift v", utils::packageDescription("medoutcon")$Version,
    ": Causal Mediation Analysis Under Mediator-Outcome Confounding"
  ))
}
