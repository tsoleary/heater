# ------------------------------------------------------------------------------
# Pseudobulk analysis plots
# June 14, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)

# Load data
res <- readRDS(here::here("output/degs/pseudobulk_DESeq_res.rds"))

# Quick volcano plot
p <- res |>
  as_tibble(rownames = "gene") |>
  ggplot(aes(label = gene)) +
  geom_hline(yintercept = -log10(0.05),
             color = "grey50",
             linetype = 2) +
  geom_point(aes(y = -log10(padj),
                 x = log2FoldChange,
                 fill = padj < 0.05,
                 size = padj < 0.05),
             color = "grey90",
             stroke = 0.2,
             shape = 21) +
  scale_size_manual(values = c(1, 2)) +
  scale_fill_manual(values = c("grey80", "firebrick")) +
  scale_y_continuous(expand = c(0, 0.5),
                     breaks = c(0, 5, 10, 15),
                     labels = c("0", "5", "10", "15")) +
  cowplot::theme_cowplot() +
  theme(legend.position = "none")

# Save volcano
ggsave(here::here("output/figs/degs/pseudobulk_volcano.pdf"),
       plot = p,
       height = 10,
       width = 10,
       units = "cm")
ggsave(here::here("output/figs/degs/pseudobulk_volcano.png"),
       plot = p,
       height = 10,
       width = 10,
       units = "cm")

# Quick plotly output
plotly::ggplotly(p, tooltip = c("gene"))