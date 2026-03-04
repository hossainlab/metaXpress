# Generate hex sticker logo for metaXpress
# Run from package root: Rscript inst/logo/make_logo.R

library(ggplot2)
library(hexSticker)

# --- Central icon: DNA helix with connected network nodes ---
# DNA helix backbone curves
t <- seq(0, 2 * pi, length.out = 100)
helix1 <- data.frame(x = sin(t) * 0.6, y = t / (2 * pi) * 2 - 1,
                      strand = "a")
helix2 <- data.frame(x = -sin(t) * 0.6, y = t / (2 * pi) * 2 - 1,
                      strand = "b")

# Rungs connecting strands
n_rungs <- 6
rung_t <- seq(0, 2 * pi, length.out = n_rungs + 2)[2:(n_rungs + 1)]
rungs <- data.frame(
  x    = sin(rung_t) * 0.6,
  xend = -sin(rung_t) * 0.6,
  y    = rung_t / (2 * pi) * 2 - 1,
  yend = rung_t / (2 * pi) * 2 - 1
)

# Network nodes at rung midpoints + ends
nodes <- data.frame(
  x = c(0, 0, 0, 0.55, -0.55, 0.55, -0.55),
  y = c(-0.7, -0.15, 0.4, 0.7, -0.5, -0.1, 0.3)
)

# Network edges
edges <- data.frame(
  x    = c(0,     0,     0.55,  -0.55, 0.55),
  y    = c(-0.7,  -0.15, 0.7,   -0.5,  -0.1),
  xend = c(0,     0,     0,     0,     0),
  yend = c(-0.15, 0.4,   0.4,   -0.15, -0.15)
)

p <- ggplot() +
  # DNA helix strands
  geom_path(data = helix1, aes(x, y), colour = "#89D2DC", linewidth = 0.8,
            alpha = 0.9) +
  geom_path(data = helix2, aes(x, y), colour = "#89D2DC", linewidth = 0.8,
            alpha = 0.9) +
  # DNA rungs
  geom_segment(data = rungs, aes(x = x, y = y, xend = xend, yend = yend),
               colour = "#5AB8B0", linewidth = 0.4, alpha = 0.6) +
  # Network edges
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
               colour = "#A8E6CF", linewidth = 0.3, alpha = 0.5) +
  # Network nodes
  geom_point(data = nodes, aes(x, y), colour = "#2EC4B6", fill = "#A8E6CF",
             shape = 21, size = 1.8, stroke = 0.4) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1.2, 1.2)) +
  theme_void() +
  theme_transparent()

# --- Build hex sticker ---
sticker(
  p,
  package   = "metaXpress",
  p_size    = 18,
  p_color   = "#FFFFFF",
  p_y       = 1.45,
  p_family  = "sans",
  s_x       = 1,
  s_y       = 0.78,
  s_width   = 1.4,
  s_height  = 1.2,
  h_fill    = "#1B2A4A",
  h_color   = "#2EC4B6",
  h_size    = 1.4,
  url       = "github.com/hossainlab/metaXpress",
  u_size    = 3.5,
  u_color   = "#89D2DC",
  filename  = "man/figures/logo.png",
  dpi       = 300
)

message("Logo saved to man/figures/logo.png")
