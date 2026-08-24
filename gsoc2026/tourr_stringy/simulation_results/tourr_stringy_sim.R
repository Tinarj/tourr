library(tourr)
library(cassowaryr)
library(spinebil)
library(ferrn)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)


# --------------------------------------------------
# Match mean and SD
# --------------------------------------------------

match_mean_sd <- function(x, target) {
  (x - mean(x)) / sd(x) * sd(target) + mean(target)
}


# --------------------------------------------------
# Generate polynomial datasets
# --------------------------------------------------

make_poly_data <- function(
    n = 500,
    dimensions = c(4, 6, 8, 10, 12),
    seed = 1050,
    noise_sd = 0.15,
    signal_sd = c(0.01, 0.01)
) {
  
  set.seed(seed)
  
  t <- seq(-2, 2, length.out = n)
  
  poly_signal <- poly(
    t,
    degree = 2,
    raw = TRUE
  )
  
  # Generate signal once
  signal <- cbind(
    poly_signal[, 1] + rnorm(n, sd = signal_sd[1]),
    poly_signal[, 2] + rnorm(n, sd = signal_sd[2])
  )
  
  # Generate noise once for largest dimension
  max_p <- max(dimensions)
  
  noise <- matrix(
    rnorm(n * (max_p - 2), sd = noise_sd),
    nrow = n
  )
  
  # Match signal magnitude to noise
  target <- as.vector(noise)
  
  signal[, 1] <- match_mean_sd(signal[, 1], target)
  signal[, 2] <- match_mean_sd(signal[, 2], target)
  
  # Create nested datasets
  out <- lapply(dimensions, function(p) {
    
    x <- cbind(
      noise[, seq_len(p - 2), drop = FALSE],
      signal
    )
    
    colnames(x) <- paste0("V", seq_len(p))
    
    x
  })
  
  names(out) <- paste0("poly", dimensions)
  
  out
}


poly_data <- make_poly_data()

poly4  <- poly_data$poly4
poly6  <- poly_data$poly6
poly8  <- poly_data$poly8
poly10 <- poly_data$poly10
poly12 <- poly_data$poly12


# --------------------------------------------------
# Stringy05 index
# --------------------------------------------------

rescale_stringy05 <- function(z, n) {
  
  lb <- 0.05 + 3.86 / sqrt(n)
  
  pmax(
    0,
    (z - lb) / (1 - lb)
  )
}


stringy05_index <- function(mat) {
  
  z <- cassowaryr::sc_stringy05(
    mat[, 1],
    mat[, 2]
  )
  
  rescale_stringy05(
    z,
    nrow(mat)
  )
}


# --------------------------------------------------
# Repeated Jellyfish simulations
# --------------------------------------------------

run_jellyfish_simulation <- function(
    data,
    n_jellies,
    n_runs = 50,
    seed_start = 1050,
    max_tries = 40,
    store_full_result = FALSE
) {
  
  data <- as.matrix(data)
  
  p <- ncol(data)
  
  # final two variables contain the signal
  true_basis <- spinebil::basis_matrix(
    p - 1,
    p,
    p
  )
  
  true_index <- stringy05_index(
    data %*% true_basis
  )
  
  
  simulation_results <- purrr::map_dfr(
    seq_len(n_runs),
    
    function(run) {
      
      # Seed increases by 1 for each run
      seed <- seed_start + run - 1
      
      set.seed(seed)
      
      message(
        "Dimension = ", p,
        " | Jellies = ", n_jellies,
        " | Run = ", run, "/", n_runs,
        " | Seed = ", seed
      )
      
      
      jellyfish_res <- animate_xy(
        data,
        guided_tour(
          index_f = stringy05_index,
          search_f = search_jellyfish,
          n_jellies = n_jellies,
          max.tries = max_tries
        )
      )
      
      
      # Find loop containing overall best result
      best_loop <- jellyfish_res |>
        ferrn::get_best() |>
        pull(loop)
      
      
      # Extract bases from best loop
      bases_best_loop <- jellyfish_res |>
        filter(loop == best_loop) |>
        pull(basis) |>
        check_dup(0.1)
      
    
      best_jelly <- jellyfish_res |>
        ferrn::get_best()
      
      
      best_basis <- best_jelly$basis[[1]]
      
      best_projection <- data %*% best_basis
      
      
      tibble(
        dimension = p,
        n_jellies = n_jellies,
        max_tries = max_tries,
        
        run = run,
        seed = seed,
        
        true_2d_index = true_index,
        
        jellyfish_index =
          best_jelly$index_val[[1]],
        
        index_difference =
          abs(
            best_jelly$index_val[[1]] -
              true_index
          ),
        
        best_loop = best_loop,
        
        best_tries =
          best_jelly$tries[[1]],
        
        jellyfish_basis =
          list(best_basis),
        
        best_loop_bases =
          list(bases_best_loop),
        
        jellyfish_projection =
          list(best_projection),
        
        jellyfish_result =
          if (store_full_result) {
            list(jellyfish_res)
          } else {
            list(NULL)
          }
      )
    }
  )
  
  simulation_results
}





results_8d_100 <- run_jellyfish_simulation(
  data = poly8,
  n_jellies = 100,
  n_runs = 20,
  seed_start = 1050,
  max_tries = 40
)


saveRDS(
  results_8d_100,
  "sim_res_8d_100.rds"
)



plot_jelly_runs <- function(results, ncol = 2) {
  
  plots <- purrr::map(seq_len(nrow(results)), function(i) {
    
    proj <- as.data.frame(
      results$jellyfish_projection[[i]]
    )
    
    names(proj) <- c("x", "y")
    
    ggplot(proj, aes(x, y)) +
      geom_point(size = 0.7, alpha = 0.7) +
      coord_equal() +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank()
      ) +
      labs(
        x = NULL,
        y = NULL,
        title = paste0("Run ", results$run[i]),
        subtitle = paste0(
          "seed = ", results$seed[i],
          ", loop = ", results$best_loop[i],
          ", tries = ", results$best_tries[i],
          "\nstringy05 = ",
          round(results$jellyfish_index[i], 3),
          ", true = ",
          round(results$true_2d_index[i], 3)
        )
      )
  })
  
  patchwork::wrap_plots(
    plots,
    ncol = ncol
  )
}


plot_jelly_runs(
  results_test,
  ncol = 2
)
