# ------------------------------------------------------------------------------
# Plot colors and other things to be used across all plots
# TS O'Leary
# ------------------------------------------------------------------------------

# Acclimation colors
acc_colors <- c(
  "#8698C0", 
  "#D09F7B"
)

# Cell-type colors
cell_type_colors <- c(
  "grey90",
  "#ADD9F4",
  "#57A4B2",
  "#D39C9E",
  "#FEF29A",
  "#F9DCEE",
  "#819FC5",
  "#A7BF9B",
  "#bfa3a4"
)

# Cell-type colors
color_cell_type <- tibble(
  cell_type = c(
    "germ cell", 
    "peripheral nervous system prim.",
    "ectoderm prim.",
    "mesoderm prim.",
    "endoderm prim.",
    "foregut & hindgut prim.",
    "ventral nerve cord prim.",
    "tracheal prim.",
    "amnioserosa"
  ), 
  colors = c(
    "grey90", 
    "#ADD9F4",
    "#57A4B2",
    "#D39C9E",
    "#FEF29A",
    "#F9DCEE",
    "#819FC5",
    "#A7BF9B",
    "#bfa3a4"
  )
)
