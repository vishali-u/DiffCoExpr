# RStudio Inc. (2013). Tabsets. Shiny Gallery. 
# URL:https://shiny.rstudio.com/gallery/tabsets.html

library(shiny)
library(shinyalert)

# Increase the file upload size limit to 100 MB
options(shiny.maxRequestSize = 100 * 1024^2)

# Define UI for random distribution app ----
ui <- fluidPage(
  
  # Change title
  titlePanel("DiffCoExpr"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      tags$p("This is a simple Shiny App that is part of the DiffCoExpr R
             package. Its purpose to demonstrate how DiffCoExpr can be used
             to construct and compare coexpression networks."),
      
      # Add vertical space
      br(),
      
      tags$p("Description: DiffCoExpr can be used construct cell-type-specific
             coexpression networks, and subsequently compare coexpression
             under different cell types. This Shiny App is part of the
             DiffCoExpr package. It permits the user to upload a 10X dataset,
             generate expression matrices, correlation matrices, and 
             coexpression networks. Subsequently, it also allows the user to
             compare expression levels, under two different cell types (or more
             generally, two conditions), of two coexpressed genes."),
      
      br(),
      
      # Accept input from the user
      tags$p("Instructions: To build coexpression networks for a cell type,
             Below, enter or select values required to perform the analysis
             Default values are shown. Then press 'Generate Coexpression 
             Networks'. Note this part may take up to 10 minutes. Navigate through 
             the first 2 tabs to the right to explore the results."),
      
      br(),
      
      # Input files
      useShinyalert(force = TRUE),  # Set up shinyalert
      uiOutput("tab1"),
      actionButton(inputId = "data1",
                   label = "Input 1 Details"),
      
      textInput("geneBCPath",
                "Gene BC Matrices: Please download a directory containing 10X
                Gene BC matrices. This directory should include
                matrix.mtx, genes.tsv or features.tsv, and barcodes.csv.
                Example data is given above; download the directory, and paste
                the path to the directory here."),
      
      uiOutput("tab2"),
      actionButton(inputId = "data2",
                   label = "Input 2 Details"),
      
      fileInput("cellTypesPath", 
                "Cell Types: Please upload a .csv file mapping cell type labels 
                to marker genes. The first column must be the marker genes and 
                the second column must be cell type labels. Example data is 
                given above; download the .csv file and upload it here.",
                accept = c(".csv")),
      
      selectInput("cellTypeA", 
                  "CellTypeA: Select a cell type from the file you uploaded.", 
                  choices = NULL),
      
      selectInput("cellTypeB", 
                  "CellTypeB: Select another cell type from the file you 
                  uploaded. This cell type will be used for the second part
                  (differential coexpression) below.", 
                  choices = NULL),
      
      textInput(inputId = "minCellCount",
                label = "Enter the minimum number of cells that you require in
                your analysis after filtering. This should be a positive value
                greater than 2:", "5"),
      
      textInput(inputId = "minGeneCount",
                label = "Enter the minimum number of genes that you require in
                your analysis after filtering. This should be a positive value
                greater than 2:", "5"),
      
      sliderInput(inputId = "thresholdVariation",
                  label = "Choose a value between 0 and 1. This value will be the
                  percentile used to filter out genes with low variation. For
                  example, the default is 0.20, so any genes that have a 
                  variation lower than the variation at the 20th percentile will
                  be filtered out:", 
                  min = 0,
                  max = 1,
                  value = 0.2),
      
      sliderInput(inputId = "thresholdCorrelation",
                  label = "Choose a value between 0 and 1. This value will be the
                  percentile used to filter out pairs of genes with low
                  correlations. For example, the default is 0.80, so any pairs of
                  genes that have a correlation lower than the correlation at the 
                  80th percentile will be filtered out:", 
                  min = 0, 
                  max = 1,
                  value = 0.8),
      
      textInput(inputId = "thresholdFC",
                label = "Enter a threshold for correlation. This value will be
                the minimum log2 fold change of correlation values calculated
                from two coexpression networks. For example, the default
                value is 0.50, so any pairs of genes that have a log2 fold 
                change less than 0.50 will be filtered out.", "0.5"),
      
      # actionButton
      actionButton(inputId = "button1",
                   label = "Generate Coexpression Networks"),
      
      # Add vertical spacing
      br(),
      br(),
      br(),
      
      tags$p("Instructions: To see pairs of differentially coexpressed genes 
             between the two chosen cell types, navigate to tab 3. Then, 
             to plot expression levels across cells for two genes that are
             differentially coexpressed, select two genes below. See tab 3
             for possible options. Press 'Run Differential Coexpression'."),
      
      selectizeInput("gene1", 
                     "Gene1: Select a gene.", 
                  choices = NULL),
      
      selectizeInput("gene2", 
                     "Gene2: Select a different gene.", 
                     choices = NULL),
      
      actionButton(inputId = "button2",
                   label = "Generate DiffCoExpr Plot"),
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      
      tags$style(type = "text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),
      
      # Output: only show the coexpression network dataframe, the network
      # graph, the differential coexpression dataframe, and the differential
      # coexpression scatterplot
      tabsetPanel(type = "tabs",
                  tabPanel("Coexpressed Genes",
                           h3("Coexpression Network for Cell Type A:"),
                           br(),
                           verbatimTextOutput("CoexprNetwork")),
                  
                  tabPanel("Coexpression Network Plot",
                           h3("Coexpression Network Plot for Cell Type A:"),
                           br(),
                           plotOutput("CoexprNetworkPlot")),
                  
                  tabPanel("Differentially Coexpressed Genes",
                           h3("Genes that are differentially coexpressed between 
                              Cell Type A and Cell Type B:"),
                           br(),
                           verbatimTextOutput("DiffCoexpTable")),
                  
                  tabPanel("Differential Coexpression Plot",
                           h3("Expression levels for two genes that are 
                              differentially coexpressed between Cell Type A 
                              and Cell Type B:"),
                           br(),
                           plotOutput("DiffCoExprPlot"))
                  
      )
    )
  )
)

# Define server logic for random distribution app
server <- function(input, output, session) {
  
  # URLs for downloading data
  url1 <- a("Example Dataset: Gene BC Matrices", 
            href = paste("https://github.com/vishali-u/DiffCoExpr/tree/main",
                         "/inst/extdata/filtered_gene_bc_matrices", sep = ''))
  output$tab1 <- renderUI({
    tagList("Download:", url1)
  })
  
  url2 <- a("Example Dataset: Cell Types", 
            href = paste("https://github.com/vishali-u/DiffCoExpr/tree/main",
                         "/inst/extdata/cellTypes.csv", sep = ''))
  output$tab2 <- renderUI({
    tagList("Download:", url2)
  })
  
  observeEvent(input$data1, {
    # Show a modal when the button is pressed
    shinyalert(title = "Example Input 1: Gene-BC Matrices",
               text = "A single-cell RNA-seq dataset from peripheral blood 
               mononuclear cells (PBMCs) from a healthy donor that was 
               published May 26th, 2016. 
               Citation: 10X Genomics. 2016. “3k PBMCs from a Healthy Donor.”. 
               URL https://www.10xgenomics.com/resources/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0.",
               type = "info")
  })
  
  observeEvent(input$data2, {
    # Show a modal when the button is pressed
    shinyalert(title = "Example Input 2: Cell Types",
               text = "A .csv file mapping cell type labels to marker genes that
               were used to identify the cell type label. The data was generated
               using information gather from a Seurat tutorial. 
               Citation: Butler, A. (2015). Seurat: Tools for Single Cell Genomics.",
               type = "info")
  })
  
  # Add to cellType drop-down menu
  
  cellTypes <- reactiveVal(NULL)
  
  # First dropdown menu
  observe({
    req(input$cellTypesPath)
    
    cellTypeFile <- input$cellTypesPath
    df <- read.csv(cellTypeFile$datapath)
    
    # Update cellTypes
    cellTypes(unique(df[[2]]))  # assuming cell types are in the second column
    updateSelectInput(session, 'cellTypeA', choices = cellTypes())
    updateSelectInput(session, 'cellTypeB', choices = cellTypes())
  })
  
  # Update second dropdown menu based on the first selection
  observe({
    req(input$cellTypeA)
    availableChoices <- setdiff(cellTypes(), input$cellTypeA)
    updateSelectInput(session, 'cellTypeB', choices = availableChoices)
  })
  
  # Prepare data separately so it is only done once
  # This step is the most time consuming
  startPrepareData <- eventReactive(eventExpr = input$button1, {
    
    DiffCoExpr::prepareData(geneMatrixPath = input$geneBCPath,
                            cellTypesPath = input$cellTypesPath$datapath)
  })
  
  # Create expression matrix separately because it will be needed again for
  # button2
  startExprMatrixA <- eventReactive(eventExpr = input$button1, {
    
    if (! is.null(startPrepareData)) {
      srat <- startPrepareData()
      DiffCoExpr::getExpressionMatrix(srat = srat,
                                      cellType = input$cellTypeA)
    }
  
  })
  
  # Construct coexpression network for cell type A
  startCoExprNetA <- eventReactive(eventExpr = input$button1, {
    
    if ((! is.null(startPrepareData)) && (! is.null(startExprMatrixA))) {
      
      exprMatrixA <- startExprMatrixA()
      
      corrMatrixA <-
        DiffCoExpr::getCorrelationMatrix(expressionMatrix = exprMatrixA,
                                         minCellCount = as.numeric(input$minCellCount),
                                         minGeneCount = as.numeric(input$minGeneCount),
                                         minPt = as.numeric(input$thresholdVariation))
      
      DiffCoExpr::getCoexpressionNetwork(correlationMatrix = corrMatrixA,
                                         thresholdCorrelation = as.numeric(input$thresholdCorrelation))
    
    }
  })
  
  # Repeat for the second cell type
  startExprMatrixB <- eventReactive(eventExpr = input$button1, {
    
    if (! is.null(startPrepareData)) {
      DiffCoExpr::getExpressionMatrix(srat = startPrepareData(),
                                      cellType = input$cellTypeB)
    }
  })
  
  startCoExprNetB <- eventReactive(eventExpr = input$button1, {
    
    if ((! is.null(startPrepareData)) && (! is.null(startExprMatrixB))) {
      
      exprMatrixB <- startExprMatrixB()
      
      corrMatrixB <-
        DiffCoExpr::getCorrelationMatrix(expressionMatrix = exprMatrixB,
                                         minCellCount = as.numeric(input$minCellCount),
                                         minGeneCount = as.numeric(input$minGeneCount),
                                         minPt = as.numeric(input$thresholdVariation))
      
      DiffCoExpr::getCoexpressionNetwork(correlationMatrix = corrMatrixB,
                                         thresholdCorrelation = as.numeric(input$thresholdCorrelation))
      
    }
  })
  
  # Add gene choices to the gene dropdown menus in the second part
  gene1Choices <- reactiveVal(NULL)
  gene2Choices <- reactiveVal(NULL)
  
  # Get differentially coexpressed genes
  startDiffCoExpr <- eventReactive(eventExpr = input$button1, {
    if ((! is.null(startCoExprNetA)) && (! is.null(startCoExprNetA))) {
      
      coexprNetA <- startCoExprNetA()
      coexprNetB <- startCoExprNetB()
      
      diffCoexpTable <- 
        DiffCoExpr::getDifferentialCoexpression(networkA = coexprNetA,
                                                networkB = coexprNetB,
                                                thresholdLogFC = as.numeric(input$thresholdFC))
      
      
      gene1Choices(unique(diffCoexpTable$Gene1))
      gene2Choices(unique(diffCoexpTable$Gene2))
      
      updateSelectInput(session, "gene1", choices = gene1Choices())
      updateSelectInput(session, "gene2", choices = gene2Choices())
      
      diffCoexpTable
    }
  })
  
  # observe({
  #   updateSelectInput(session, "gene1", choices = gene1Choices())
  #   updateSelectInput(session, "gene2", choices = gene2Choices())
  # })
  
  # Plot for differentially coexpressed genes
  startDiffCoExprPlot <- eventReactive(eventExpr = input$button2, {
    if ((! is.null(startDiffCoExpr)) && (! is.null(startExprMatrixA)) &&
        (! is.null(startExprMatrixB))) {
      
      diffCoexpTable <- startDiffCoExpr()
      
      gene1 <- input$gene1
      gene2 <- input$gene2
      
      exprMatrixB <- startExprMatrixB()
      exprMatrixA <- startExprMatrixA()
      
      DiffCoExpr::plotDifferentialCoexpression(diffCoexpTable = diffCoexpTable,
                                               expressionMatrixA = exprMatrixA,
                                               expressionMatrixB = exprMatrixB,
                                               gene1 = gene1,
                                               gene2 = gene2,
                                               conditionA = input$cellTypeA,
                                               conditionB = input$cellTypeB)
      
    }
  })
  
  # Only output one of the coexpression networks (the first one)
  output$CoexprNetwork <- renderPrint({
    if (! is.null(startCoExprNetA)) {
      coexprNetwork <- startCoExprNetA()
      coexprNetwork
    }
  })
  
  # Coexpression plotting
  output$CoexprNetworkPlot <- renderPlot({
    if (! is.null(startCoExprNetA))
      edgeList <- startCoExprNetA()
      DiffCoExpr::plotCoexpressionNetwork(edgeList = edgeList)
  })
  
  # Table of differentially coexpressed genes
  output$DiffCoexpTable <- renderPrint({
    if (! is.null(startDiffCoExpr)) {
      diffCoexpTable <- startDiffCoExpr()
      diffCoexpTable
    }
  })
  
  output$DiffCoExprPlot <- renderPlot({
    if (! is.null(startDiffCoExprPlot)) {
      startDiffCoExprPlot()
    }
  })
}

# Create Shiny app ----
shinyApp(ui, server)
# [END]