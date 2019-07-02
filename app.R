library(shiny)
library(DT)
library(shinydashboard)
library(affy)

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(
    width = 350,
    fileInput(inputId = "files",
              label = "Upload .CEL files for one sample or a group of replicates (replicates will be averaged)",
              multiple = TRUE),
    verbatimTextOutput("fileName"),
    actionButton(inputId = "example",
                 label = "Use Example Input"),
    textOutput("addLine1"),
    fileInput(inputId = "crossProbes",
              label = "Xenograft samples only: Upload .csv of cross-hybridized probes to be removed (single column)",
              multiple = FALSE),
    verbatimTextOutput("hybridFileName"),
    textOutput("addLine2"),
    verbatimTextOutput("errorBox"),
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
      column(width = 6, DT::dataTableOutput("resTable"),textOutput("placeHolder")),
      column(width = 6, plotOutput("bplots"))
    ),
    div(style = "height:5000px;")
  )
)


server <- function(input, output, session) {
  # increases allowable size of file uploads
  options(shiny.maxRequestSize=30*1024^2)

  rv <- reactiveValues() # container for selected files
  nclicks <- reactiveVal(0) # reactive number times compute is clicked
  finished <- reactiveVal("Results will be displayed here") 

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

  observeEvent(input$example, {
    rv$data <- c("./data/example.CEL")
    rv$name <- c("example.CEL")
    rv$errorMessage <- NULL
    nclicks(0)
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
    if (!is.null(input$crossProbes$name)) {
      tmp <- paste("Uploaded File:\n",basename(input$crossProbes$name),"")
    }
    return(tmp)
  })
  
  output$addLine2 <- renderText({
    "--------------------------------------------------------------------------------"
  })
  
  output$placeHolder <- renderText({
    finished()
  })
  
  
  # determines if there is a problem with uploaded files
  FileError <- function(flist) {
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
  

  res <- reactiveVal() # holds all calculation results

  observeEvent(input$compute, {
    rv$errorMessage <- FileError(rv$data)
    if (!is.null(rv$errorMessage)) {
      return(NULL)
    }

    if(nclicks() != 0){
      return(NULL)
    }
    
    # Increment clicks and prevent concurrent analyses
    nclicks(nclicks() + 1)
    
    if (identical(rv$name,"example.CEL")) {
      # uses precomputed results
      tmpres <- readRDS("./data/exres.rds")
      res(tmpres)
    } else {
      progress <- Progress$new(session)
      progress$set(message = 'Preprocessing',
                   detail = 'This may take a few minutes...')
      
      source("preprocess.R")
      samples_matrix <- pre_proc(rv$data,rv$name,input$crossProbes$datapath)
      progress$set(message = 'Preprocessing Completed', detail = "")
      on.exit(progress$close())
      
      withProgress(message = 'Computing Enrichment of NR-Target Gene Sets...', value = 0,{
        source("compute_enriched_NRs.R")
        results <- CalcEnrich(samples_matrix,input$crossProbes$datapath)
        res(results)
      })
    }
    
    finished("")
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
    data.frame(res()[[2]],
                      Actions = shinyInput(actionButton, 15, 'button_', label = "See Targets", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' )
               )
  })

  # renders the results dataframe with column of buttons
  output$resTable <- DT::renderDataTable({
    if (is.null(res())){
      return(NULL)
    }
    tmptbl()
  },server = FALSE, escape = FALSE, selection = 'none', rownames = FALSE)
  

  # chooses the correct plot depending which button is pushed
  chosenPlt <- eventReactive(input$select_button, {
    if (is.null(res())){
      return(NULL)
    }
    selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    selectedNR <- rownames(tmptbl())[selectedRow]
    return(res()[[1]][[selectedNR]])
  })

  # shows the selected table
  observe({
    if (is.null(res())){
      return(NULL)
    }
    output$bplots <- renderPlot({
      chosenPlt()[[1]]
    },height = (3000/74)*chosenPlt()[[2]],width = 500)
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
             height = (60/74)*chosenPlt()[[2]], width = 5,limitsize = FALSE)
    }
  )
  
}

shinyApp(ui, server)































































































