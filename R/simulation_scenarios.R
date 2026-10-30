scqte_draw_equal <- function(n, error, seed) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n, 0, 2)
  A <- rbinom(n, 1, plogis(-0.5 + 0.5 * x1 + 0.5 * x2))
  epsilon <- if (error == "normal") rnorm(n) else runif(n)^(-1 / 5)
  data.frame(x1 = x1, x2 = x2, A = A,
             Y = 1 + A + x1 + x2 + (1 + x1 + 1.5 * A) * epsilon)
}

scqte_draw_shift <- function(n, seed) {
  set.seed(seed)
  sites <- 20L
  per_site <- n / sites
  t <- seq(0, 1, length.out = sites)
  x1 <- unlist(lapply(t, function(v) rbeta(per_site, 1 + 4 * v, 5 - 4 * v)),
               use.names = FALSE)
  x2 <- 2 * unlist(lapply(t, function(v) rbeta(per_site, 5 - 4 * v, 1 + 4 * v)),
                   use.names = FALSE)
  A <- rbinom(n, 1, plogis(-0.5 + 0.5 * x1 + 0.5 * x2))
  epsilon <- rnorm(n)
  data.frame(x1 = x1, x2 = x2, A = A,
             Y = 1 + A + x1 + x2 + (1 + x1 + 1.5 * A) * epsilon)
}

scqte_reverse_blocks <- function(data, sizes) {
  ends <- cumsum(sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  order <- unlist(rev(Map(seq.int, starts, ends)), use.names = FALSE)
  data[order, , drop = FALSE]
}

scqte_scenario <- function(design = "equal", error = "normal", sites = 20L,
                           pattern = "equal", replication = 1L, n = 10000L) {
  if (!design %in% c("equal", "shift", "unequal") ||
      !error %in% c("normal", "pareto") ||
      length(replication) != 1L || replication < 1L ||
      replication != as.integer(replication) ||
      length(n) != 1L || n < 1L || n != as.integer(n)) {
    stop("Invalid scenario argument.")
  }
  sites <- as.integer(sites)
  replication <- as.integer(replication)
  if (design == "equal") {
    if (pattern != "equal" || sites < 2L || n %% sites != 0L) {
      stop("Equal-site design requires pattern='equal' and n divisible by sites.")
    }
    seed <- if (error == "normal") 666L + replication else 100666L + replication
    data <- scqte_draw_equal(n, error, seed)
    sizes <- rep(n / sites, sites)
  } else if (design == "shift") {
    if (error != "normal" || pattern != "equal" || sites != 20L || n %% sites != 0L) {
      stop("The published covariate-shift design uses normal error and 20 equal sites.")
    }
    data <- scqte_draw_shift(n, 200666L + replication)
    sizes <- rep(n / sites, sites)
  } else {
    if (error != "normal" || sites != 200L || n != 10000L ||
        !pattern %in% c("50-50", "60-40", "70-30", "40-60", "30-70")) {
      stop("The published unequal-site design uses normal error, n=10000, M=200, and a listed pattern.")
    }
    data <- scqte_draw_equal(n, "normal", 666L + replication)
    base_pattern <- switch(pattern, `40-60` = "60-40", `30-70` = "70-30", pattern)
    first <- as.integer(sub("-.*", "", base_pattern))
    second <- as.integer(sub(".*-", "", base_pattern))
    sizes <- c(rep(first, 100L), rep(second, 100L))
    if (pattern %in% c("40-60", "30-70")) {
      data <- scqte_reverse_blocks(data, sizes)
      sizes <- rev(sizes)
    }
  }
  list(data = data, sizes = as.integer(sizes), design = design,
       error = error, sites = sites, pattern = pattern, replication = replication)
}
