# ------------------------------------------------------------------------------
# HEATER Project Shiny App 
# May 31, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries -----
require(tidyverse)
require(Seurat)
require(Signac)
require(cowplot)
library(shiny)
library(shinyWidgets)
library(shinyThings)
library(shinybusy)
library(bslib)
library(bsicons)
library(plotly)
library(DESeq2)


# Set themes and stuff
tso_theme <- bs_theme(primary = "#FFA6F7", 
                      font_scale = NULL, 
                      bootswatch = "united")


# Load data etc ----------------------------------------------------------------
# Latest Seurat Object
dat <- readRDS(
  here::here("data/processed/seurat_object/10_dat_linked.rds")
)
# Pseudobulk analysis from DESeq2
res <- readRDS(
  here::here("output/degs/pseudobulk_DESeq_res.rds")
)


# Cards ------------------------------------------------------------------------
test_card <- card(
  card_header("Test",
              class = "bg-dark"),
  "This is a test body"
)

coverage_plot_card <- card(
  height = 700,
  full_screen = TRUE,
  card_header("Chromatin accessibility and gene expression", class = "bg-dark"),
  plotOutput("Coverage_Plot"),
)

dim_plot_card <- card(
  height = 700,
  full_screen = TRUE,
  card_header("UMAP plot", class = "bg-dark"),
  plotOutput("Dim_Plot"),
)

gene_error_card <- card(
  max_height = 200,
  max_width = 200,
  card_header("Error", class = "error"),
  "Please enter a valid gene name from the dropdown menu."
)

volcano_card <- card(
  height = 400,
  full_screen = TRUE,
  card_header("Pseudobulk DESeq Volcano Plot", class = "bg-dark"),
  plotlyOutput("Volcano_Plot"),
)
  

# Value boxes ------------------------------------------------------------------
cell_number_box <- value_box(
  title = "Number of total cells", 
  value = length(Cells(dat)),
  showcase = bs_icon("record-circle"),
  theme_color = "success",
  max_width = 400,
  max_height = 200
)
qc_card <- card(
  card_header("Quality control"),
  lorem::ipsum(1),
  width = 400,
  max_height = 200
)



################################################################################
# Define UI for application ----------------------------------------------------
################################################################################
ui <- page_navbar(
  
  # Main options ------------
  theme = tso_theme,
  title = "HEATER Project",

  # Panel setups ---------------------------------------
  tabPanel(title = "Home",
           volcano_card),
  
  # Quality Control Tab
  tabPanel(title = "Quality Control",
           layout_column_wrap(
             width = 1/3, height = 300, 
             cell_number_box, qc_card
           )
  ),
  
  # GOI on UMAP Feature Plots
  tabPanel(title = "UMAP Plot",
           sidebarLayout(
             sidebarPanel(
               tags$style(".well {background-color:#ffffff;}"),
               # Split Acclimation Temp
               shinyThings::radioSwitchButtons(
                 inputId = "split",
                 choices = c("Split", "Together"),
                 selected_background = "#00589a",
                 not_selected_background = "grey90"
               ),
               width = 3
             ),
             # Show generated plot
             mainPanel(
               as_fill_carrier(dim_plot_card),
               width = 9
             )
           )
  ),
  
  
  # Coverage Plot Panel
  tabPanel(title = "Coverage Plot",
    sidebarLayout(
        sidebarPanel(
          tags$style(".well {background-color:#ffffff;}"),
            # A select input for gene of interest
            selectizeInput(inputId = "gene",
                           label = "Gene:",
                           choices = NULL),
            # Split Acclimation Temp
            shinyThings::radioSwitchButtons(
              inputId = "Split",
              choices = c("Split", "Together"),
              selected_background = "#00589a",
              not_selected_background = "grey90"
            ),
            width = 3
        ),
        # Show generated plot
        mainPanel(
          as_fill_carrier(coverage_plot_card),
           width = 9
        )
    )
  ),
  nav_spacer(),
  tabPanel(title = "About",
           test_card)
)



################################################################################
# Define server logic ----------------------------------------------------------
################################################################################
server <- function(input, output, session) {
  
  # Coverage plot --------------------------------------------------------------
  
  # Server side selectize for faster rendering
  updateSelectizeInput(session,
                 inputId = "gene",
                 label = "Gene:",
                 choices = sort(dat@assays$SCT@counts@Dimnames[[1]]),
                 server = TRUE)
  
    output$Coverage_Plot <- renderPlot({
      
      # Validate for coverage plot
      validate(
        need(
          input$gene %in% dat@assays$SCT@counts@Dimnames[[1]],
          message = "Please enter a valid gene name from the dropdown menu."
          )
        )
      
      # Modal to show while rendering the plot
      show_modal_spinner(spin = "fingerprint",
                         color = "#782c54",
                         text = "Rendering plot")
      
      # Set Default Assay for the Coverage Plot
      DefaultAssay(dat) <- "ATAC"
      
      # Create plot for just the 18°C
      p18 <- CoveragePlot(
        object = dat |> subset(acc_temp == "18°C"),
        region = input$gene,
        features = input$gene,
        expression.assay = "SCT",
        extend.upstream = 500,
        extend.downstream = 1000
      ) + cowplot::theme_cowplot(font_size = 24)
      
      # Create plot for just the 25°C
      p25 <- CoveragePlot(
        object = dat |> subset(acc_temp == "25°C"),
        region = input$gene,
        features = input$gene,
        expression.assay = "SCT",
        extend.upstream = 500,
        extend.downstream = 1000
      ) + cowplot::theme_cowplot(font_size = 24)
      
      # Combine 18C and 25C plots
      p <- cowplot::plot_grid(
        p18, 
        p25,
        labels = c("A. 18°C", "B. 25°C"),
        label_size = 20
      )
      
      # Create title for the plot
      title_gene <- input$gene
      title <- cowplot::ggdraw() + 
        cowplot::draw_label(
          input$gene,
          size = 24,
          fontface = "bold.italic",
          x = 0,
          hjust = 0
        ) +
        theme(
          plot.margin = margin(0, 0, 0, 7)
        )
      
      # Put plot together with the title
      p <- cowplot::plot_grid(
        title, p,
        ncol = 1,
        rel_heights = c(0.1, 1)
      )
      remove_modal_spinner()
      # Print out plot
      p
    })
    # End of coverage plot -----------------------------------------------------
    
    # Dim Plot or Features -----------------------------------------------------
    output$Dim_Plot <- renderPlot({
      # Create plot
      ifelse(
        # Conditional test
        input$split == "Split",
        
        # Split plot
        p <- DimPlot(dat, 
                reduction = "umap",
                split.by = "acc_temp",
                pt.size = 1)  +
          labs(title = element_blank()) +
          theme_cowplot() +
          theme(legend.position = "bottom",
                title = element_blank(),
                axis.title = element_blank(),
                axis.line = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank()) +
          guides(color = guide_legend(byrow = TRUE,
                                      nrow = 1,
                                      override.aes = list(size = 2))),
        # Together plot
        p <- DimPlot(dat, 
                group.by = "acc_temp",
                pt.size = 2) +
          scale_color_manual(name = element_blank(),
                             values = c("#43aa8b", "#f3722c")) +
          labs(title = element_blank()) +
          theme_nothing() +
          theme(legend.position = "bottom") +
          guides(color = guide_legend(byrow = TRUE,
                                      nrow = 1,
                                      override.aes = list(size = 2)))
        
        
      )
      # Plot Output
      p
      
    })
    # End of Dim Plot ----------------------------------------------------------
    
    
    # QC plots -----------------------------------------------------------------
    output$QC_plots <- renderPlot({
      # Create QC plot
      
      
    })
    # End of QC plots ----------------------------------------------------------
    
    # Volcano plots -----------------------------------------------------------------
    output$Volcano_Plot <- renderPlot({
      # Create Volcano plot
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
                           labels = c("", "5", "10", "15")) +
        ggh4x::coord_axes_inside(labels_inside = TRUE) +
        cowplot::theme_cowplot() +
        theme(legend.position = "none",
              axis.title.x = element_blank())
      
      plotly::ggplotly(p, tooltip = c("gene"))
    })
    # End of QC plots ----------------------------------------------------------
    
}

# Run the application ----------------------------------------------------------
shinyApp(
  ui = ui, 
  server = server
)
