context("Variance estimate is larger when specifying level of independence using 'ids'")

set.seed(4353645)

n <- 1000
j <- 20
k <- n / j
re <- rep(rnorm(j, sd = 2), each = k)

id <- rep(1:j, each = k)
W1 <- rbinom(n, 1, 0.5)
W2 <- rbinom(n, 1, 0.33)

# treatment on the group level
A <- rep(rbinom(j, 1, plogis(-1 + log(2)*mean(W1) + log(1.5)*mean(W2))), each = k)
M <- rbinom(n, 1, plogis(2 - log(3)*A + log(1.5)*W1 - log(2)*W2))

tmp <- data.frame(re, W1, W2, A, M)
# Y is drawn with random effects
Y <- unlist(lapply(split(tmp, id),
                   function(x) with(x, rbinom(k, 1, plogis(re - log(0.5)*A - log(2)*M + log(1.5)*W1 + log(1.2)*W2)))))

tml_direct_no_id <- medoutcon(W = data.frame(W1, W2),
                              A = A,
                              Z = NULL,
                              M = M,
                              Y = Y,
                              effect = "direct")

tml_direct_id <- medoutcon(W = data.frame(W1, W2),
                           A = A,
                           Z = NULL,
                           M = M,
                           Y = Y,
                           effect = "direct",
                           ids = id)

tml_indirect_no_id <- medoutcon(W = data.frame(W1, W2),
                                A = A,
                                Z = NULL,
                                M = M,
                                Y = Y,
                                effect = "indirect")

tml_indirect_id <- medoutcon(W = data.frame(W1, W2),
                             A = A,
                             Z = NULL,
                             M = M,
                             Y = Y,
                             effect = "indirect",
                             ids = id)

os_direct_no_id <- medoutcon(W = data.frame(W1, W2),
                             A = A,
                             Z = NULL,
                             M = M,
                             Y = Y,
                             effect = "direct",
                             estimator = "onestep")

os_direct_id <- medoutcon(W = data.frame(W1, W2),
                          A = A,
                          Z = NULL,
                          M = M,
                          Y = Y,
                          effect = "direct",
                          ids = id,
                          estimator = "onestep")

os_indirect_no_id <- medoutcon(W = data.frame(W1, W2),
                               A = A,
                               Z = NULL,
                               M = M,
                               Y = Y,
                               effect = "indirect",
                               estimator = "onestep")

os_indirect_id <- medoutcon(W = data.frame(W1, W2),
                            A = A,
                            Z = NULL,
                            M = M,
                            Y = Y,
                            effect = "indirect",
                            ids = id,
                            estimator = "onestep")

test_that("NDE: Variance is larger with 'ids', TMLE", {
  expect_gt(tml_direct_id$var, tml_direct_no_id$var)
})

test_that("NDE: Variance is larger with 'ids', one-step", {
  expect_gt(os_direct_id$var, os_direct_no_id$var)
})

test_that("NIE: Variance is larger with 'ids', TMLE", {
  expect_gt(tml_indirect_id$var, tml_indirect_no_id$var)
})

test_that("NIE: Variance is larger with 'ids', one-step", {
  expect_gt(os_indirect_id$var, os_indirect_no_id$var)
})
