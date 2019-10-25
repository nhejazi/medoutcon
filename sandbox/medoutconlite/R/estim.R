#' Stochastic Mediation Effect of Contrasts Under Intermediate Confounding
#'
#' @param data ...
#' @param contrast ...
#' @param weights ...
#' @param candidatesg ...
#' @param candidatese ...
#' @param candidatesm ...
#' @param candidatesq ...
#' @param candidatesr ...
#' @param candidatesu ...
#' @param candidatesv ...
#' @param nfolds ...
#' @param family.outcome ...
#'
#' @importFrom stats glm coef offset predict var sd
#' @importFrom stats gaussian binomial plogis qlogis
#'
#' @export
estimatheta <- function(data, contrast, weights,
                        candidatesg,
                        candidatese,
                        candidatesm,
                        candidatesq,
                        candidatesr,
                        candidatesu,
                        candidatesv,
                        nfolds,
                        family.outcome) {
  # contrasts
  aprime <- contrast[1]
  astar <- contrast[2]

  # extract data
  A <- data[, "A"]
  M <- data[, substr(names(data), 1, 1) == "M"]
  Z <- data[, substr(names(data), 1, 1) == "Z"]
  Y <- data[, "Y"]
  W <- data[, substr(names(data), 1, 1) == "W"]

  # observations and cross-validation
  n <- length(A)
  valSets <- split(sample(seq_len(n)), rep(seq_len(nfolds), length = n))

  # fit nuisance functions
  fitg <- mySL(
    Y = A,
    X = W,
    family = stats::binomial(),
    SL.library = candidatesg,
    validRows = valSets
  )
  fite <- mySL(
    Y = A,
    X = data.frame(W, M),
    family = stats::binomial(),
    SL.library = candidatese,
    validRows = valSets
  )
  fitm <- mySL(
    Y = Y,
    X = data.frame(W, M, Z, A),
    family = family.outcome,
    SL.library = candidatesm,
    validRows = valSets
  )
  fitq <- mySL(
    Y = Z,
    X = data.frame(W, A),
    family = stats::binomial(),
    SL.library = candidatesq,
    validRows = valSets
  )
  fitr <- mySL(
    Y = Z,
    X = data.frame(W, M, A),
    family = stats::binomial(),
    SL.library = candidatesr,
    validRows = valSets
  )

  # compute u pseudo outcome and fit u function
  gone <- fitg$SL.predict[, 1]
  eone <- fite$SL.predict[, 1]

  gprime <- gone * aprime + (1 - gone) * (1 - aprime)
  gstar <- gone * astar + (1 - gone) * (1 - astar)
  eprime <- eone * aprime + (1 - eone) * (1 - aprime)
  estar <- eone * astar + (1 - eone) * (1 - astar)

  qoneprime <- stats::predict(fitq,
                              newdata = data.frame(W, A = aprime))$pred[, 1]
  roneprime <- stats::predict(fitr,
                              newdata = data.frame(W, M, A = aprime))$pred[, 1]

  qprime <- Z * qoneprime + (1 - Z) * (1 - qoneprime)
  rprime <- Z * roneprime + (1 - Z) * (1 - roneprime)

  mprime <- stats::predict(fitm,
                           newdata = data.frame(W, M, Z, A = aprime))$pred[, 1]
  hstar <- gprime / gstar * qprime / rprime * estar / eprime
  upseudo <- mprime * hstar

  # in case values are all very close to each other
  if (sd(upseudo) < .Machine$double.eps) {
    candidatesu <- c("SL.mean")
  }

  # predict on estimated u
  fitu <- mySL(
    Y = upseudo,
    X = data.frame(W, Z, A),
    family = stats::gaussian(),
    SL.library = candidatesu,
    validRows = valSets
  )

  uprimeone <-
    stats::predict(fitu, newdata = data.frame(W, Z = 1, A = aprime))$pred[, 1]
  uprimezero <-
    stats::predict(fitu, newdata = data.frame(W, Z = 0, A = aprime))$pred[, 1]
  uprime <- Z * uprimeone + (1 - Z) * uprimezero


  # compute v pseudo outcome and fit v function
  mprimeone <-
    stats::predict(fitm,
                   newdata = data.frame(W, M, Z = 1, A = aprime))$pred[, 1]
  mprimezero <-
    stats::predict(fitm,
                   newdata = data.frame(W, M, Z = 0, A = aprime))$pred[, 1]

  # build v nuisance function
  vpseudo <- mprimeone * qoneprime + mprimezero * (1 - qoneprime)
  if (stats::sd(vpseudo) < .Machine$double.eps) {
    candidatesv <- c("SL.mean")
  }

  # fit nuisance function v
  fitv <- mySL(
    Y = vpseudo,
    X = data.frame(W, A),
    family = stats::gaussian(),
    SL.library = candidatesv,
    validRows = valSets
  )
  vstar <- stats::predict(fitv, newdata = data.frame(W, A = astar))$pred[, 1]

  # compute one step
  ipwprime <- as.numeric(A == aprime) / gprime
  ipwstar <- as.numeric(A == astar) / gstar

  # compute EIF components
  eify <- ipwprime * hstar / mean(ipwprime * hstar) * (Y - mprime)
  eifu <- ipwprime / mean(ipwprime) * (uprimeone - uprimezero) *
    (Z - qoneprime)
  eifv <- ipwstar / mean(ipwstar) * (vpseudo - vstar)
  eif <- weights * (eify + eifu + eifv + vstar)
  os <- mean(eif)
  eifos <- eif

  # now, compute the TMLE
  stopcrit <- FALSE
  iter <- 1

  # iterative TMLE
  while (!stopcrit) {
    # compute contrast difference for u
    uprimediff <- (uprimeone - uprimezero)

    # evaluate h' under various contrasts
    hstar <- gprime / gstar * qprime / rprime * estar / eprime
    hstarone <- gprime / gstar * qoneprime / roneprime * estar / eprime
    hstarzero <- gprime / gstar * (1 - qoneprime) / (1 - roneprime) *
      estar / eprime

    # first fluctuation/tilting
    suppressWarnings(
      tiltm <- stats::glm(
        as.formula("Y ~ -1 + offset(mprime_logit) + hstar"),
        data = data.frame(list(Y = Y, A = A,
                               mprime_logit = stats::qlogis(mprime),
                               hstar = hstar)),
        subset = A == aprime,
        weights = weights / gprime,
        family = stats::binomial()
      )
    )

    # second fluctuation/tilting
    suppressWarnings(
      tiltq <- stats::glm(
        stats::as.formula("Z ~ -1 + offset(qoneprime_logit) + uprimediff"),
        data = data.frame(list(Z = Z, A = A,
                               qoneprime_logit = stats::qlogis(qoneprime),
                               uprimediff = uprimediff)),
        subset = A == aprime,
        weights = weights / gprime,
        family = stats::binomial()
      )
    )

    # extract epsilon and set to zero if failed fluctuation/tilting
    coefq <- stats::coef(tiltq)
    coefm <- stats::coef(tiltm)
    if (is.na(coefq)) coefq <- 0
    if (is.na(coefm)) coefm <- 0

    mprime <- stats::plogis(stats::qlogis(mprime) + coefm * hstar)
    mprimeone <- stats::plogis(stats::qlogis(mprimeone) + coefm * hstarone)
    mprimezero <- stats::plogis(stats::qlogis(mprimezero) + coefm *
      hstarzero)

    qoneprime <- stats::plogis(stats::qlogis(qoneprime) + coefq *
      uprimediff)
    qprime <- Z * qoneprime + (1 - Z) * (1 - qoneprime)

    # iterate iterator
    iter <- iter + 1

    # note: interesting stopping criterion
    stopcrit <- max(abs(c(coefm, coefq))) < 0.001 / n^(0.6) | iter > 6
  }

  vpseudo <- mprimeone * qoneprime + mprimezero * (1 - qoneprime)
  vstar[vstar < 1e-3] <- 1e-3
  vstar[vstar > 1 - 1e-3] <- 1 - 1e-3

  suppressWarnings(
    tiltv <- stats::glm(
      as.formula("vpseudo ~ offset(vstar_logit)"),
      data = data.frame(list(vpseudo = vpseudo, A = A,
                             vstar_logit = stats::qlogis(vstar))),
      subset = A == astar,
      weights = weights / gstar,
      family = stats::binomial()
    )
  )
  vstar <- stats::plogis(stats::qlogis(vstar) + stats::coef(tiltv))

  qprime <- Z * qoneprime + (1 - Z) * (1 - qoneprime)
  hstar <- gprime / gstar * qprime / rprime * estar / eprime
  upseudo <- mprime * hstar
  vpseudo <- mprimeone * qoneprime + mprimezero * (1 - qoneprime)
  ipwprime <- as.numeric(A == aprime) / gprime
  ipwstar <- as.numeric(A == astar) / gstar

  eify <- ipwprime * hstar / mean(ipwprime * hstar) * (Y - mprime)
  eifu <- ipwprime / mean(ipwprime) * (uprimeone - uprimezero) *
    (Z - qoneprime)
  eifv <- ipwstar / mean(ipwstar) * (vpseudo - vstar)

  eif <- weights * (eify + eifu + eifv + vstar)

  tmle <- mean(eif)
  eiftmle <- weights * (eify + eifu + eifv + vstar)

  return(list(
    estimates = c(os = os, tmle = tmle),
    eifs = cbind(os = eifos, tmle = eiftmle)
  ))
}

#' Stochastic Direct/Indirect Effects Under Intermediate Confounding
#'
#' @param data ...
#' @param weights ...
#' @param candidatesg ...
#' @param candidatese ...
#' @param candidatesm ...
#' @param candidatesq ...
#' @param candidatesr ...
#' @param candidatesu ...
#' @param candidatesv ...
#' @param nfolds ...
#' @param family.outcome ...
#'
#' @importFrom stats var
#'
#' @export
mediation <- function(data, weights, candidatesg, candidatese,
                      candidatesm, candidatesq, candidatesr,
                      candidatesu, candidatesv, nfolds,
                      family.outcome) {
  # preliminaries
  n <- nrow(data)

  # compute under different contrasts
  theta11 <- estimatheta(data,
    contrast = c(1, 1), weights,
    candidatesg, candidatese, candidatesm,
    candidatesq, candidatesr, candidatesu,
    candidatesv, nfolds, family.outcome
  )
  theta10 <- estimatheta(data,
    contrast = c(1, 0), weights,
    candidatesg, candidatese, candidatesm,
    candidatesq, candidatesr, candidatesu,
    candidatesv, nfolds, family.outcome
  )
  theta00 <- estimatheta(data,
    contrast = c(0, 0), weights,
    candidatesg, candidatese, candidatesm,
    candidatesq, candidatesr, candidatesu,
    candidatesv, nfolds, family.outcome
  )

  # point estimates
  indirect <- theta11$estimates - theta10$estimates
  direct <- theta10$estimates - theta00$estimates

  # standard error estimates
  seindirect <- sqrt(diag(stats::var(theta11$eifs - theta10$eifs)) / n)
  sedirect <- sqrt(diag(stats::var(theta10$eifs - theta00$eifs)) / n)

  # output
  return(list(
    effects = rbind(indirect, direct),
    ses = rbind(seindirect, sedirect)
  ))
}
