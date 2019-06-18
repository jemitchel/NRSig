library(shiny)
library(DT)
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(
    fileInput(inputId = "files",
              label = "Upload .CEL files for one sample or multiple replicates (replicates will be averaged)",
              multiple = TRUE),
    actionButton(inputId = "example",
                 label = "Use Example Input"),
    verbatimTextOutput("fileName"),
    actionButton(inputId = "compute",
                 label = "Compute",
                 style="color: #fff; background-color: #008000; border-color: #008000"),
    uiOutput("download1"),
    uiOutput("download2"),
    uiOutput("download3")
  ),
  dashboardBody(
    fluidRow(
      column(width = 6,
             DT::dataTableOutput("resTable")),
      column(width = 6,
             plotOutput("bplots"))
    )
  )
)

# ui <- fluidPage(
#   
#   sidebarLayout(
#     sidebarPanel(
#       fileInput(inputId = "files",
#                 label = "Upload .CEL files for one sample or multiple replicates (replicates will be averaged)",
#                 multiple = TRUE),
#       actionButton(inputId = "example",
#                    label = "Use Example Input"),
#       verbatimTextOutput("fileName"),
#       actionButton(inputId = "compute",
#                    label = "Compute",
#                    style="color: #fff; background-color: #008000; border-color: #008000"),
#       uiOutput("download1"),
#       uiOutput("download2"),
#       uiOutput("download3")
#     ),
#     mainPanel(
#       column(width = 6,
#              DT::dataTableOutput("resTable")),
#       column(width = 6,
#              plotOutput("bplots"))
#     )
#   )
# )


server <- function(input, output, session) {
  # increases allowable size of file uploads
  options(shiny.maxRequestSize=30*1024^2)
  
  rv <- reactiveValues() # container for selected files
  
  observe({
    rv$data <- input$files$datapath
    rv$name <- input$files$name
  })
  
  observeEvent(input$example, {
    rv$data <- c("C:/Users/jonat/Documents/R/NRSig-app/data/example.CEL")
    rv$name <- c("C:/Users/jonat/Documents/R/NRSig-app/data/example.CEL")
  })
  
  output$fileName <- renderText({ 
    tmp <- "No Files Uploaded"
    if (!is.null(rv$name)) {
      tmp <- "Uploaded Files:\n"
      for (i in basename(rv$name)) {
        tmp <- paste(tmp,i,"\n")
      }
    } 
    return(tmp)
  })
  
  # errorMessage <- reactive({ 
  #   if (res() == "no file error") {
  #     return("No input files uploaded")
  #   } else if (res() == "wrong file error") {
  #     return("Wrong type of file uploaded. Need .CEL file.")
  #   } else if (res() == "wrong platform error") {
  #     return(".CEL file is for wrong platform. Only upload files for HG-U133 Plus 2.0")
  #   } else if (is.null(res())) {
  #     return(NULL)
  #   }
  # })
  
  FileTypeError <- function(flist) {
    if (is.null(flist)) {
      return("")
    }
    for (i in 1:length(flist)) {
      len <- nchar(flist[[i]])
      last4 <- substr(flist[[i]],len-3,len)
      if (!identical(last4,".CEL")) {
        return("Upload includes file(s) not ending in .CEL ")
      }
      cdfName <- whatcdf(flist[[i]])
      if (!identical(cdfName,"HG-U133_Plus_2")) {
        return("Input file(s) are for wrong platform. Use CEL files for HG-U133_Plus_2 only.")
      }
    }
  }
  
  res <- eventReactive(input$compute, {
    print(rv$data)
    shiny::validate(
      need(!is.null(rv$name), message = "Please upload files to process"),
      FileTypeError(rv$name)
    )
    cool <- 4
    progress <- Progress$new(session)
    progress$set(message = 'Preprocessing',
                 detail = 'This may take a few minutes...')
    source("C:/Users/jonat/Documents/R/NRSig-app/preprocess.R")
    samples_matrix <- pre_proc(rv$data,rv$name)
    progress$set(message = 'Preprocessing Completed', detail = "")
    on.exit(progress$close())

    withProgress(message = 'Computing Enrichment of NR-Target Gene Sets...', value = 0,{
      source("C:/Users/jonat/Documents/R/NRSig-app/compute_enriched_NRs.R")
      results <- CalcEnrich(samples_matrix)
      return(results)
    })

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
    data.frame(res()[[2]],
                      Actions = shinyInput(actionButton, 15, 'button_', label = "See Targets", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' )
               )
  })
  
  # renders the results dataframe with column of buttons
  output$resTable <- DT::renderDataTable({
    tmptbl()
  },server = FALSE, escape = FALSE, selection = 'none', rownames = FALSE)

  # chooses the correct plot depending which button is pushed
  chosenPlt <- eventReactive(input$select_button, {
    selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    selectedNR <- rownames(tmptbl())[selectedRow]
    return(res()[[1]][[selectedNR]])
  })

  # shows the selected table
  observe({
    output$bplots <- renderPlot({
      chosenPlt()[[1]]
    },height = (5000/74)*chosenPlt()[[2]],width = 500)
  })
  
  
  output$download1 <- renderUI({
    if(!is.null(res())) {
      downloadButton('downloadfRMA', label = "fRMA Preprocessed Input")
    }
  })

  output$download2 <- renderUI({
    if(!is.null(res())) {
      downloadButton('downloadDiffex', label = "Differential Expression Results")
    }
  })
    
  output$download3 <- renderUI({
    if(!is.null(chosenPlt())) {
      downloadButton('downloadPlot', label = "Boxplot Figure")
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
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(chosenPlt()[[1]]$labels$title, '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = chosenPlt()[[1]], device = "png",
             height = (60/74)*chosenPlt()[[2]], width = 5,limitsize = FALSE)
    }
  )
  
}

shinyApp(ui, server)






