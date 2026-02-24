chain_graph <- function(p) {
  G <- list(conn = list(), weight = list())
  G$conn[1:(p - 1)] <- as.integer(2:p) 
  G$conn[p] <- as.integer(1)
  G$weight[1:p] <- rep(1.0, p)
  class(G) <- "FusionGraph"
  G
}

conn_from_adj <- function(adjacency_matrix) {
  G <- list(conn = list(), weight = list())
  pairs <- which(upper.tri(adjacency_matrix), arr.ind = TRUE)
  for (idx in 1:nrow(pairs)) {
    i = unname(pairs[idx,][1])
    j = unname(pairs[idx,][2])
    if (adjacency_matrix[i, j] == 1) {
      G$conn[[i]] = j
      G$conn[[j]] = i
      G$weight[[i]] = G$weight[[j]] <- as.numeric(1)
    }
  }
  class(G) <- "FusionGraph"
  G
}

model_default <- function()
  list(verbose     = 0, # default control options
       intercept   = TRUE,
       normalize   = TRUE,
       nlambda1    = 100,
       min_ratio   = 1e-2
  )


optim_enet_default <- function(d)
  list(verbose     = 0, # default control options
       timer       = FALSE,
       maxiter     = max(50, d),
       method      = "quadra",
       threshold   = 1e-6,
       monitor     = 0,
       usechol     = TRUE
  )

optim_breg_default <- function(d)
  list(verbose     = 0, # default control options
       timer       = FALSE,
       maxiter     = 10,
       method      = "quadra",
       threshold   = 1e-4,
       monitor     = 0,
       usechol     = FALSE
  )

optim_fused_default <- function(d)
  list(verbose       = 0, # default control options
       timer         = FALSE,
       maxiterout    = 100,
       maxiterin     = 10000,
       maxactivation = 10, 
       accuracy      = 1e-6,
       monitor       = 0,
       fusioncheck   = "all" ## c("all","active","none", "naive", "huber")
  )

status_to_message <- function(status) {
  message <- switch(as.character(status),
                    "0"  = "converged",
                    "1"  = "max # of iterate reached",
                    "2"  = "max # of feature reached",
                    "3"  = "system has become singular",
                    "Return status not recognized"
  )
  message
}

