.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "medoutcon v", utils::packageDescription("medoutcon")$Version,
    ": ", utils::packageDescription("medoutcon")$Title
  ))
}
