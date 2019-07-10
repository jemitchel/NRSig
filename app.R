library(shiny)
library(DT)
library(shinydashboard)
library(affy)
library(ggplot2)

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(
    width = 350,
    fileInput(inputId = "files",
              label = "Upload .CEL files for one sample or a group of replicates (replicates will be averaged)",
              multiple = TRUE),
    actionButton(inputId = "example",
                 label = "Use Example Input"),
    verbatimTextOutput("fileName"),
    textOutput("addLine1"),
    fileInput(inputId = "crossProbes",
              label = "For xenograft samples only: Upload .csv of cross-hybridized probes to be removed (single column)",
              multiple = FALSE),
    actionButton(inputId = "exampleCrossProbes",
                 label = "Load probes crosshybridized with athymic nude mice - Hollingshead et al. study",
                 style="white-space: normal;
                 height:45px;
                 width:250px;
                 font-size: 12px;"),
    verbatimTextOutput("hybridFileName"),
    textOutput("addLine2"),
    verbatimTextOutput("errorBox"),
    actionButton(inputId = "clear",
                 label = "Clear All Input"),
    actionButton(inputId = "compute",
                 label = "Compute",
                 style="color: #fff; background-color: #008000; border-color: #008000"),
    uiOutput("download1"),
    uiOutput("download2"),
    uiOutput("download3"),
    uiOutput("download4")
    ),
  dashboardBody(
    fluidRow(
      box(title = "NRSig", width = 12, textOutput("intro")),
      column(6,
             textOutput("tableTitle"),
             textOutput("placeHolder"),
             tags$head(tags$style("#tableTitle{color: blue;
                                  font-size: 20px;
                                  }")),
             DT::dataTableOutput("resTable")),
      column(6,div(style = 'overflow-x: auto; height:850px;', plotOutput("bplots"))),
    div(style = "height:1200px;"))
    )
)


# ui <- fluidPage(
#   sidebarLayout(
#     sidebarPanel(width = 4,
#       fileInput(inputId = "files",
#                 label = "Upload .CEL files for one sample or a group of replicates (replicates will be averaged)",
#                 multiple = TRUE),
#       actionButton(inputId = "example",
#                    label = "Use Example Input"),
#       verbatimTextOutput("fileName"),
#       textOutput("addLine1"),
#       fileInput(inputId = "crossProbes",
#                 label = "For xenograft samples only: Upload .csv of cross-hybridized probes to be removed (single column)",
#                 multiple = FALSE),
#       actionButton(inputId = "exampleCrossProbes",
#                    label = "Load probes crosshybridized with athymic nude mice - Hollingshead et al. study",
#                    style="white-space: normal;
#                    height:45px;
#                    width:250px;
#                    font-size: 12px;"),
#       verbatimTextOutput("hybridFileName"),
#       textOutput("addLine2"),
#       verbatimTextOutput("errorBox"),
#       actionButton(inputId = "clear",
#                    label = "Clear All Input"),
#       actionButton(inputId = "compute",
#                    label = "Compute",
#                    style="color: #fff; background-color: #008000; border-color: #008000"),
#       uiOutput("download1"),
#       uiOutput("download2"),
#       uiOutput("download3"),
#       uiOutput("download4")
#     ),
# 
#     mainPanel(
#       titlePanel("NRSig"),
#       textOutput("intro"),
#       # box(title = "NRSig", width = 12, textOutput("intro")),
#       column(width = 6,
#              textOutput("tableTitle"),
#              textOutput("placeHolder"),
#              tags$head(tags$style("#tableTitle{color: blue;
#                                   font-size: 20px;
#                                   }")),
#              DT::dataTableOutput("resTable",width = "100%")),
#       column(width = 4, plotOutput("bplots")),
#       div(style = "height:5000px;")
#     )
#   )
# )


# ui <- fixedPage(
#   fixedRow(
#     column(width=4,
#       fileInput(inputId = "files",
#                 label = "Upload .CEL files for one sample or a group of replicates (replicates will be averaged)",
#                 multiple = TRUE),
#       actionButton(inputId = "example",
#                    label = "Use Example Input"),
#       verbatimTextOutput("fileName"),
#       textOutput("addLine1"),
#       fileInput(inputId = "crossProbes",
#                 label = "For xenograft samples only: Upload .csv of cross-hybridized probes to be removed (single column)",
#                 multiple = FALSE),
#       actionButton(inputId = "exampleCrossProbes",
#                    label = "Load probes crosshybridized with athymic nude mice - Hollingshead et al. study",
#                    style="white-space: normal;
#                    height:45px;
#                    width:250px;
#                    font-size: 12px;"),
#       verbatimTextOutput("hybridFileName"),
#       textOutput("addLine2"),
#       verbatimTextOutput("errorBox"),
#       actionButton(inputId = "clear",
#                    label = "Clear All Input"),
#       actionButton(inputId = "compute",
#                    label = "Compute",
#                    style="color: #fff; background-color: #008000; border-color: #008000"),
#       uiOutput("download1"),
#       uiOutput("download2"),
#       uiOutput("download3"),
#       uiOutput("download4")),
#     column(width=8, title = "NRSig", 
#            textOutput("intro"),
#            fixedRow(
#              column(5,DT::dataTableOutput("resTable",width="70%")),
#              column(3,plotOutput("bplots"))
#            )
#     )
#   )
# )
#   # column(width = 9,
#   #   fixedRow(box(title = "NRSig", width = 12, textOutput("intro")),
#   #     column(width = 10,box(title = "NRSig", width = 12, textOutput("intro"))),
#   #   fixedRow(
#   #     column(width = 5,
#   #            textOutput("tableTitle"),
#   #            textOutput("placeHolder"),
#   #            tags$head(tags$style("#tableTitle{color: blue;
#   #                                 font-size: 20px;
#   #                                 }")),
#   #            DT::dataTableOutput("resTable")),
#   #     column(width = 6, plotOutput("bplots"))
#   #   )
#   # )))


server <- function(input, output, session) {
  # increases allowable size of file uploads
  options(shiny.maxRequestSize=35*1024^2)
  
  rv <- reactiveValues() # container for selected files
  rv2 <- reactiveValues() # container for crosshybrid probes file
  nclicks <- reactiveVal(0) # reactive number times compute is clicked
  status <- reactiveValues("finished"="Results will be displayed here","title"="") 
  
  
  output$intro <- renderText({
    "This program computes differentially expressed genes between a query MCF7 microarray sample
    and prior distributions generated from serum-starved MCF7 samples. This allows for
    identification of transcriptionally active nuclear receptors in the query sample by
    calculating enriched NR-target gene sets among the differentially expressed genes. Currently,
    NRSig only accepts samples generated on the HG-U133 Plus 2.0 platform."
  })
  
  observe({
    rv$data <- input$files$datapath
    rv$name <- input$files$name
    rv$errorMessage <- NULL
    nclicks(0)
  })
  
  observe({
    rv2$data <- input$crossProbes$datapath
    rv2$name <- input$crossProbes$name
  })
  
  observeEvent(input$example, {
    rv$data <- c("./data/example.CEL")
    rv$name <- c("example.CEL")
    rv$errorMessage <- NULL
    nclicks(0)
  })
  
  observeEvent(input$exampleCrossProbes, {
    rv2$data <- c("./data/crosshybrid.csv")
    rv2$name <- c("./data/crosshybrid.csv")
  })
  
  observeEvent(input$clear, {
    rv$data <- NULL
    rv$name <- NULL
    rv2$data <- NULL
    rv2$name <- NULL
  })
  
  output$fileName <- renderText({
    # tmp <- "No Files Uploaded"
    tmp <- NULL
    if (!is.null(rv$name)) {
      tmp <- "Uploaded Files:\n"
      for (i in basename(rv$name)) {
        tmp <- paste(tmp,i,"\n")
      }
    }
    return(tmp)
  })
  
  output$errorBox <- renderText({
    rv$errorMessage
  })
  
  output$addLine1 <- renderText({
    "--------------------------------------------------------------------------------"
  })
  
  output$hybridFileName <- renderText({
    # tmp <- "No Files Uploaded"
    tmp <- NULL
    if (!is.null(rv2$name)) {
      tmp <- paste("Uploaded File:\n",basename(rv2$name),"")
    }
    return(tmp)
  })
  
  output$addLine2 <- renderText({
    "--------------------------------------------------------------------------------"
  })
  
  output$placeHolder <- renderText({
    status$finished
  })
  
  output$tableTitle <- renderText({
    status$title
  })
  
  
  # determines if there is a problem with uploaded files
  celFileError <- function(flist) {
    if (is.null(flist)) {
      return("Error: Please upload files to process")
    } else {
      for (i in 1:length(flist)) {
        len <- nchar(flist[[i]])
        last4 <- substr(flist[[i]],len-3,len)
        if (!identical(last4,".CEL")) {
          return("Error: Upload includes file(s) not ending in .CEL ")
        }
        cdfName <- whatcdf(flist[[i]])
        if (!identical(cdfName,"HG-U133_Plus_2")) {
          return("Error: Input file(s) are for wrong platform. Use CEL files for HG-U133_Plus_2 only.")
        }
      }
    }
    return(NULL)
  }
  
  
  csvFileError <- function(fl) {
    if (is.null(fl)) {
      return(NULL)
    } else {
      len <- nchar(fl)
      last4 <- substr(fl,len-3,len)
      if (!identical(last4,".csv")) {
        return("Error: Uploaded probe list not a .csv")
      }
      dimension <- ncol(read.csv(fl))
      if (dimension > 1) {
        return("Error: Uploaded probe list has more than 1 column. Ensure list is a single column of Affymetrix probe symbols.")
      }
    }
    return(NULL)
  }
  
  
  res <- reactiveVal() # holds all calculation results
  chosenPlt <- reactiveVal() # holds output z-score figure
  
  observeEvent(input$compute, {
    err1 <- celFileError(rv$data)
    err2 <- csvFileError(rv2$data)
    
    if (!is.null(err1)) {
      rv$errorMessage <- err1
      return(NULL)
    } else if (!is.null(err2)) {
      rv$errorMessage <- err2
      return(NULL)
    }
    
    if(nclicks() != 0){
      return(NULL)
    }
    
    # Increment clicks and prevent concurrent analyses
    nclicks(nclicks() + 1)
    
    if (identical(rv$name,"example.CEL")) {
      # uses precomputed results
      tmpres <- readRDS("./data/exres2.rds")
      res(tmpres)
    } else {
      progress <- Progress$new(session)
      progress$set(message = 'Preprocessing',
                   detail = 'This may take a few minutes...')
      
      source("preprocess.R")
      samples_matrix <- pre_proc(rv$data,rv$name,rv2$data)
      progress$set(message = 'Preprocessing Completed', detail = "")
      on.exit(progress$close())
      
      withProgress(message = 'Computing Enrichment of NR-Target Gene Sets...', value = 0,{
        source("compute_enriched_NRs.R")
        results <- CalcEnrich(samples_matrix,rv2$data)
        res(results)
      })
    }
    
    # finished("")
    # finished("Enriched NR-Target Gene Sets")
    status$title <- "Enriched NR-Target Gene Sets"
    status$finished <- ""
  })
  
  # makes inputs for buttons on output results table
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  # produces the results dataframe with column of buttons
  tmptbl <- reactive({
    if (is.null(res())){
      return(NULL)
    }
    tbl <- res()[[2]]
    tbl[,2:3] <- round(as.numeric(tbl[,2:3]),4) 
    data.frame(tbl,
               Actions = shinyInput(actionButton, nrow(res()[[2]]), 'button_', label = "See Targets", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' )
    )
  })
  
  # renders the results dataframe with column of buttons
  output$resTable <- DT::renderDataTable({
    if (is.null(res())){
      return(NULL)
    }
    tmptbl()
  },server = FALSE, escape = FALSE, selection = 'none', rownames = FALSE, 
  options = list(pageLength = 15, dom = 't'),
  colnames = c('NR','P Value','Adjusted P Value','Target Genes'))
  
  
  # chooses the correct plot depending which button is pushed
  observeEvent(input$select_button, {
    if (is.null(res())){
      return(NULL)
    }
    selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    selectedNR <- rownames(tmptbl())[selectedRow]
    chosenPlt(res()[[1]][[selectedNR]])
  })
  
  # shows the selected table
  observe({
    if (is.null(res())){
      return(NULL)
    }
    if (is.null(chosenPlt())) {
      chosenPlt(res()[[1]][[1]])
    } 
    
    output$bplots <- renderPlot({
      chosenPlt()[[1]]
    },height = (3000/74)*chosenPlt()[[2]],width = 500)
    
    # output$bplots <- renderPlot({
    #   chosenPlt()[[1]]
    # },width = "auto",
    # height = function() {
    #   session$clientData$output_bplots_width*(3000/74)*chosenPlt()[[2]]/500
    # })
    
  })
  
  
  # download buttons and handlers
  output$download1 <- renderUI({
    if(!is.null(res())) {
      downloadButton('downloadfRMA', label = "fRMA Preprocessed Input",
                     style="color: #fff; background-color: #A9A9A9; border-color: #A9A9A9")
    }
  })
  
  output$download2 <- renderUI({
    if(!is.null(res())) {
      downloadButton('downloadDiffex', label = "Differential Expression Results",
                     style="color: #fff; background-color: #A9A9A9; border-color: #A9A9A9")
    }
  })
  
  output$download3 <- renderUI({
    if(!is.null(res())) {
      downloadButton('downloadEnrich', label = "Enrichment Results",
                     style="color: #fff; background-color: #A9A9A9; border-color: #A9A9A9")
    }
  })
  
  output$download4 <- renderUI({
    if(!is.null(chosenPlt())) {
      downloadButton('downloadPlot', label = "Boxplot Figure",
                     style="color: #fff; background-color: #A9A9A9; border-color: #A9A9A9")
    }
  })
  
  output$downloadfRMA <- downloadHandler(
    filename = function() { paste('input_preprocessed', '.csv', sep='') },
    content = function(file) {
      write.csv(res()[[4]], file, row.names = FALSE)
    }
  )
  
  output$downloadDiffex <- downloadHandler(
    filename = function() { paste('diff_exprs', '.csv', sep='') },
    content = function(file) {
      write.csv(res()[[3]], file, row.names = FALSE)
    }
  )
  
  output$downloadEnrich <- downloadHandler(
    filename = function() { paste('enriched', '.csv', sep='') },
    content = function(file) {
      write.csv(res()[[2]], file, row.names = FALSE)
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(chosenPlt()[[1]]$labels$title, '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = chosenPlt()[[1]], device = "png",
             height = (45/74)*chosenPlt()[[2]], width = 6.75,limitsize = FALSE)
    }
  )
  
  }

shinyApp(ui, server)




























































































