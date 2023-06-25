# ------------------------------------------------------------------------------
# HEATER Project Shiny App 
# May 31, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries ---------------------------------------------------------------
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

# Set themes and stuff ---------------------------------------------------------
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
about_card <- card(
  card_header("About",
              class = "bg-info"),
  "This is an interactive Shiny App for the HEATER project. Instructions and overview to appear on this page.",
  max_height = 400
)

dev_card <- card(
  card_header("Under Development",
              class = "bg-warning"),
  "This app is under development",
  max_height = 400
)

analysis_card <- card(
  card_header("Analysis",
              class = "bg-dark"),
  markdown('Here is a [link](https://tsoleary.github.io/heater/docs/analysis.html) to an overview of all the analysis involved in the project'),
  max_height = 400
)

# Plot cards -----
dim_plot_card <- card(
  height = 700,
  full_screen = TRUE,
  card_header(
    class = "d-flex justify-content-between",
    "UMAP plot",
    checkboxInput(
      "split", 
      "Split by acclimation", 
      TRUE)
  ),
  plotOutput("Dim_Plot"),
)

feature_qc_plot_card <- card(
  height = 700,
  full_screen = TRUE,
  card_header("Quality contro feature plot", class = "bg-dark"),
  plotOutput("Feature_QC_Plot"),
)

coverage_plot_card <- card(
  height = 400,
  full_screen = TRUE,
  card_header("Chromatin accessibility and gene expression", class = "bg-dark"),
  plotOutput("Coverage_Plot"),
)

feat_plot_card <- card(
  height = 400,
  full_screen = TRUE,
  card_header("Feature plot", class = "bg-dark"),
  plotOutput("Feature_Plot"),
)

vln_plot_card <- card(
  height = 400,
  full_screen = TRUE,
  card_header("Violin Plot", class = "bg-dark"),
  plotOutput("Violin_Plot"),
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
  title = "HEATER Project App",

  # Panel setups ---------------------------------------
  
  # Home page
  tabPanel(title = "Home",
           fluidPage(
             fluidRow(
               column(
                 3,
                 about_card,
                 analysis_card
               ),
               column(
                 2,
                 dev_card
               ),
             )
           )
           ),
  
  # Quality Control Tab
  tabPanel(title = "Quality Control",
           layout_column_wrap(
             width = 1/3, height = 300, 
             cell_number_box, qc_card
           ),
           # A select input for gene of interest
           selectizeInput(inputId = "feature",
                          label = "QC Feature:",
                          select = "nCount_RNA",
                          choices = colnames(dat@meta.data)[c(2:3,6:12)]),
           feat_plot_card
  ),
  
  # Cell types and clusters panel
  tabPanel(title = "Cell types",
           as_fill_carrier(dim_plot_card)
    ),
  
  
  # Gene plotting panel
  tabPanel(title = "Gene of Interest",
    sidebarLayout(
        sidebarPanel(
          tags$style(".well {background-color:#ffffff;}"),
            # A select input for gene of interest
            selectizeInput(inputId = "gene",
                           label = "Gene:",
                           select = "lola",
                           choices = sort(dat@assays$SCT@counts@Dimnames[[1]])),
          width = 3
          ),
          
        # Show generated plot
        mainPanel(
          fluidRow(
            column(8,
                   as_fill_carrier(coverage_plot_card)
            ),
            column(4,
                   as_fill_carrier(feat_plot_card)
            )
          ),
          fluidRow(
            as_fill_carrier(vln_plot_card)
          ),
           width = 9
        ),
    ),
  ),

  nav_spacer(),


  tabPanel(title = "About",
           about_card)
)



################################################################################
# Define server logic ----------------------------------------------------------
################################################################################
server <- function(input, output, session) {
  
  # Coverage plot --------------------------------------------------------------
  
  # # Server side selectize for faster rendering
  # updateSelectizeInput(session,
  #                inputId = "gene",
  #                label = "Gene:",
  #                choices = sort(dat@assays$SCT@counts@Dimnames[[1]]),
  #                server = TRUE)
  
    output$Coverage_Plot <- renderPlot({
      
      # Set Default assay to peaks
      DefaultAssay(dat) <- "peaks"
      
      # # Validate for coverage plot
      # validate(
      #   need(
      #     input$gene %in% dat@assays$SCT@counts@Dimnames[[1]],
      #     message = "Please enter a valid gene name from the dropdown menu."
      #     )
      #   )
      # 
      # # Modal to show while rendering the plot
      # show_modal_spinner(spin = "fingerprint",
      #                    color = "#782c54",
      #                    text = "Rendering plot")
      
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
      #remove_modal_spinner()
      # Print out plot
      p
    })
    # End of coverage plot -----------------------------------------------------
    
    # Feature plot -------------------------------------------------------------
    output$Feature_Plot <- renderPlot({
      
      # # Validate for coverage plot
      # validate(
      #   need(
      #     input$gene %in% dat@assays$SCT@counts@Dimnames[[1]],
      #     message = "Please enter a valid gene name from the dropdown menu."
      #   )
      # )
      
      # Create Feature Plot
      DefaultAssay(dat) <- "SCT"
      
      FeaturePlot(dat, 
                  features = input$gene,
                  reduction = "umap",
                  split.by = "acc_temp",
                  pt.size = 1)
    })
    # End Feature plot ---------------------------------------------------------
    
    # Violin plot --------------------------------------------------------------
    output$Violin_Plot <- renderPlot({
      
      # # Validate for coverage plot
      # validate(
      #   need(
      #     input$gene %in% dat@assays$SCT@counts@Dimnames[[1]],
      #     message = "Please enter a valid gene name from the dropdown menu."
      #   )
      # )
      
      # Create Feature Plot
      DefaultAssay(dat) <- "SCT"
      
      VlnPlot(dat, 
              features = input$gene,
              split.by = "acc_temp",
              pt.size = 0)
    })
    # End of Violin plots ------------------------------------------------------
    
    
    
    # Dim Plot -----------------------------------------------------------------
    output$Dim_Plot <- renderPlot({
      # Create plot
      ifelse(
        # Conditional test
        input$split,
        
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
    
    # Volcano plots ------------------------------------------------------------
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
    # End of Volcano plots -----------------------------------------------------
    
    
    # QC plots -----------------------------------------------------------------
    output$Feature_QC_Plot <- renderPlot({
      # QC Plot of Features
      dat |>
        FeaturePlot(
          features = "nCount_RNA",
          split.by = "acc_temp"
        )

    })
    # End of Volcano plots -----------------------------------------------------
}

# Run the application ----------------------------------------------------------
shinyApp(
  ui = ui, 
  server = server
)
