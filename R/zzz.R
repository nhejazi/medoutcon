.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "medoutcon v", utils::packageDescription("medoutcon")$Version,
    ": Mediation Analysis with Mediator-Outcome Confounding by Exposure"
  ))
}
