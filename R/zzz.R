.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "medoutcon v", utils::packageDescription("medoutcon")$Version,
    ": Efficient Causal Mediation Analysis Under Intermediate Confounding"
  ))
}
