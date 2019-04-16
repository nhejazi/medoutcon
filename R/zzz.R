.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "medoutcon v", utils::packageDescription("medoutcon")$Version,
    ": Causal Mediation Analysis Under Mediator-Outcome Confounding"
  ))
}
