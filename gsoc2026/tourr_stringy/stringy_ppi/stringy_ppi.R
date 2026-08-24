library(tourr)

set.seed(1050)

n <- 300
t <- seq(-1, 1, length.out = n)

poly_signal <- poly(t, degree = 2, raw = TRUE)

poly_data <- data.frame(
  V1 = rnorm(n),
  V2 = rnorm(n),
  V3 = poly_signal[, 1],
  V4 = rnorm(n),
  V5 = poly_signal[, 2]
)


saveRDS(poly_data, "blog/stringy_ppi/poly_data.rds")

# Raw stringy05
idx_stringy_raw <- stringy05()

set.seed(1050)

poly_history_raw <- save_history(
  poly_data,
  guided_tour(idx_stringy_raw)
)

saveRDS(poly_history_raw, "blog/stringy_ppi/poly_history_raw.rds")

render_gif(
  poly_data,
  planned_tour(poly_history_raw),
  display_xy(),
  gif_file = "blog/stringy_ppi/gifs/stringy05_raw.gif"
)

# Rescaled stringy05

idx_stringy_rescaled <- stringy05(rescale = TRUE)

set.seed(1050)

poly_history_rescaled <- save_history(
  poly_data,
  guided_tour(idx_stringy_rescaled)
)

saveRDS(poly_history_rescaled, "blog/stringy_ppi/poly_history_rescaled.rds")

render_gif(
  poly_data,
  planned_tour(poly_history_rescaled),
  display_xy(),
  gif_file = "blog/stringy_ppi/gifs/stringy05_rescaled.gif"
)




# Raw stringy05 with jellyfish search 

idx_stringy_raw <- stringy05()

set.seed(1050)

poly_history_raw_jellyfish <- save_history(
  poly_data,
  guided_tour(
    idx_stringy_raw,
    search_f = search_jellyfish
  )
)

saveRDS(
  poly_history_raw_jellyfish,
  "blog/stringy_ppi/poly_history_raw_jellyfish.rds"
)

render_gif(
  poly_data,
  planned_tour(poly_history_raw_jellyfish),
  display_xy(),
  gif_file = "blog/stringy_ppi/gifs/stringy05_raw_jellyfish.gif"
)