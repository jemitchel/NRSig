library(shiny)
library(DT)
library(shinydashboard)
library(affy)
library(ggplot2)

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(
    width = 350,
    fileInput(
      inputId = "cel_files",
      label = "Upload .CEL files for one sample or a group of
              replicates",
      multiple = TRUE
    ),
    actionButton(
      inputId = "example",
      label = "Use Example Input"
    ),
    verbatimTextOutput("cel_file_name"),
    textOutput("add_line_1"),
    fileInput(
      inputId = "cross_probes",
      label = "For xenograft samples only: Upload .csv of
              cross-hybridized probes to be removed (single column)",
      multiple = FALSE
    ),
    actionButton(
      inputId = "example_cross_probes",
      label = "Load probes crosshybridized with athymic nude mice
                (from Hollingshead et al. study)",
      style = "white-space: normal;
                 height:45px;
                 width:250px;
                 font-size: 12px;"
    ),
    verbatimTextOutput("hybrid_file_name"),
    textOutput("add_line_2"),
    verbatimTextOutput("error_box"),
    actionButton(
      inputId = "clear",
      label = "Clear All Input"
    ),
    actionButton(
      inputId = "compute",
      label = "Compute",
      style = "color: #fff; background-color: #008000;
                 border-color: #008000"
    ),
    uiOutput("download1"),
    uiOutput("download2"),
    uiOutput("download3"),
    uiOutput("download4")
  ),
  dashboardBody(
    fluidRow(
      box(title = "NRSig", width = 12, textOutput("intro")),
      column(
        5,
        textOutput("table_title"),
        textOutput("placeholder"),
        tags$head(tags$style("#table_title{color: blue;
                                  font-size: 20px;
                                  }")),
        DT::dataTableOutput("res_table")
      ),
      column(7, div(
        style = "overflow-x: auto; height:850px;",
        plotOutput("bplots")
      )),
      div(style = "height:1200px;")
    )
  )
)


server <- function(input, output, session) {

  # increases allowable upload file size to accomodate large .cel files
  options(shiny.maxRequestSize = 35 * 1024^2)

  rv_cel <- reactiveValues() # container for selected files
  rv_cross <- reactiveValues() # container for crosshybrid probes file
  nclicks <- reactiveVal(0) # reactive number times compute is clicked
  status <- reactiveValues(
    "finished" = "Results will be displayed here",
    "title" = "",
    "error_message" = NULL
  )


  output$intro <- renderText({
    "This program computes differentially expressed genes between a query MCF7
    microarray sample and prior distributions generated from serum-starved MCF7
    samples. This allows for identification of transcriptionally active nuclear
    receptors in the query sample by calculating enriched NR-target gene sets
    among the differentially expressed genes. NRSig currently only accepts
    samples generated on the HG-U133 Plus 2.0 Affymetrix platform. Runtime is
    approximately 10-15 minutes."
  })

  observe({
    rv_cel$data <- input$cel_files$datapath
    rv_cel$name <- input$cel_files$name
    status$error_message <- NULL
    nclicks(0)
  })

  observe({
    rv_cross$data <- input$cross_probes$datapath
    rv_cross$name <- input$cross_probes$name
  })

  observeEvent(input$example, {
    rv_cel$data <- "precomputed"
    rv_cel$name <- c("example.CEL")
    status$error_message <- NULL
    nclicks(0)
  })

  observeEvent(input$example_cross_probes, {
    rv_cross$data <- c("./data/crosshybrid.csv")
    rv_cross$name <- c("./data/crosshybrid.csv")
  })

  observeEvent(input$clear, {
    rv_cel$data <- NULL
    rv_cel$name <- NULL
    rv_cross$data <- NULL
    rv_cross$name <- NULL
  })

  output$cel_file_name <- renderText({
    # tmp <- "No Files Uploaded"
    tmp <- NULL
    if (!is.null(rv_cel$name)) {
      tmp <- "Uploaded Files:\n"
      for (i in basename(rv_cel$name)) {
        tmp <- paste(tmp, i, "\n")
      }
    }
    return(tmp)
  })

  output$error_box <- renderText({
    status$error_message
  })

  output$add_line_1 <- renderText({
    "--------------------------------------------------------------------------------"
  })

  output$hybrid_file_name <- renderText({
    tmp <- NULL
    if (!is.null(rv_cross$name)) {
      tmp <- paste("Uploaded File:\n", basename(rv_cross$name), "")
    }
    return(tmp)
  })

  output$add_line_2 <- renderText({
    "--------------------------------------------------------------------------------"
  })

  output$placeholder <- renderText({
    status$finished
  })

  output$table_title <- renderText({
    status$title
  })


  # determines if there is a problem with uploaded files
  celFileError <- function(flist) {
    if (is.null(flist)) {
      return("Error: Please upload files to process")
    } else {
      for (i in 1:length(flist)) {
        len <- nchar(flist[[i]])
        last4 <- substr(flist[[i]], len - 3, len)
        if (!identical(last4, ".CEL")) {
          return("Error: Upload includes file(s) not ending in .CEL ")
        }
        cdfName <- whatcdf(flist[[i]])
        if (!identical(cdfName, "HG-U133_Plus_2")) {
          return("Error: Input file(s) are for wrong platform. Use CEL files
                 for HG-U133_Plus_2 only.")
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
      last4 <- substr(fl, len - 3, len)
      if (!identical(last4, ".csv")) {
        return("Error: Uploaded probe list not a .csv")
      }
      dimension <- ncol(read.csv(fl))
      if (dimension > 1) {
        return("Error: Uploaded probe list has more than 1 column. Ensure list
               is a single column of Affymetrix probe symbols.")
      }
    }
    return(NULL)
  }


  res <- reactiveVal() # holds all calculation results
  chosen_plt <- reactiveVal() # holds output z-score figure

  observeEvent(input$compute, {
    if (nclicks() != 0) {
      return(NULL)
    }

    # Increment clicks and prevent concurrent analyses
    nclicks(nclicks() + 1)

    if (identical(rv_cel$data, "precomputed")) {
      # uses precomputed results
      rv_cel$data <- readRDS("./data/exres.rds")
      res(rv_cel$data)
    } else {

      # checks for file errors
      err1 <- celFileError(rv_cel$data)
      err2 <- csvFileError(rv_cross$data)

      if (!is.null(err1)) {
        status$error_message <- err1
        return(NULL)
      } else if (!is.null(err2)) {
        status$error_message <- err2
        return(NULL)
      }

      progress <- Progress$new(session)
      progress$set(
        message = "Preprocessing",
        detail = "This may take a few minutes..."
      )

      source("preprocess.R")
      samples_matrix <- pre_proc(rv_cel$data, rv_cel$name, rv_cross$data)
      progress$set(message = "Preprocessing Completed", detail = "")
      on.exit(progress$close())

      withProgress(
        message = "Computing Enrichment of NR-Target Gene Sets...",
        value = 0, {
          source("compute_enriched_NRs.R")
          results <- CalcEnrich(samples_matrix, rv_cross$data)
          res(results)
        }
      )
    }

    status$title <- "Enriched NR-Target Gene Sets"
    status$finished <- ""
  })

  # makes inputs for buttons on output results table
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    return(inputs)
  }

  # produces the results dataframe with column of buttons
  tmptbl <- reactive({
    if (is.null(res())) {
      return(NULL)
    }
    tbl <- res()[[2]]
    tbl[, 2:3] <- round(as.numeric(tbl[, 2:3]), 4)
    data.frame(tbl,
      Actions = shinyInput(actionButton, nrow(res()[[2]]), "button_",
        label = "See Targets",
        onclick = 'Shiny.onInputChange(
                                    \"select_button\",  this.id)'
      )
    )
  })

  # renders the results dataframe with column of buttons
  output$res_table <- DT::renderDataTable({
    if (is.null(res())) {
      return(NULL)
    }
    tmptbl()
  },
  server = FALSE, escape = FALSE, selection = "none", rownames = FALSE,
  options = list(pageLength = 15, dom = "t"),
  colnames = c("NR", "P Value", "Adjusted P Value", "Target Genes")
  )


  # chooses the correct plot depending which button is pushed
  observeEvent(input$select_button, {
    if (is.null(res())) {
      return(NULL)
    }
    selected_row <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    selected_NR <- rownames(tmptbl())[selected_row]
    chosen_plt(res()[[1]][[selected_NR]])
  })

  # shows the selected table
  observe({
    if (is.null(res())) {
      return(NULL)
    }
    if (is.null(chosen_plt())) {
      chosen_plt(res()[[1]][[1]])
    }

    output$bplots <- renderPlot({
      chosen_plt()[[1]]
    }, height = (3000 / 74) * chosen_plt()[[2]], width = 500)
  })


  # download buttons and handlers
  output$download1 <- renderUI({
    if (!is.null(res())) {
      downloadButton("downloadfRMA",
        label = "fRMA Preprocessed Input",
        style = "color: #fff; background-color: #A9A9A9; 
                     border-color: #A9A9A9"
      )
    }
  })

  output$download2 <- renderUI({
    if (!is.null(res())) {
      downloadButton("downloadDiffex",
        label = "Differential Expression Results",
        style = "color: #fff; background-color: #A9A9A9;
                     border-color: #A9A9A9"
      )
    }
  })

  output$download3 <- renderUI({
    if (!is.null(res())) {
      downloadButton("downloadEnrich",
        label = "Enrichment Results",
        style = "color: #fff; background-color: #A9A9A9; 
                     border-color: #A9A9A9"
      )
    }
  })

  output$download4 <- renderUI({
    if (!is.null(chosen_plt())) {
      downloadButton("download_plot",
        label = "Boxplot Figure",
        style = "color: #fff; background-color: #A9A9A9; 
                     border-color: #A9A9A9"
      )
    }
  })

  output$downloadfRMA <- downloadHandler(
    filename = function() {
      paste("input_preprocessed", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(res()[[4]], file, row.names = FALSE)
    }
  )

  output$downloadDiffex <- downloadHandler(
    filename = function() {
      paste("diff_exprs", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(res()[[3]], file, row.names = FALSE)
    }
  )

  output$downloadEnrich <- downloadHandler(
    filename = function() {
      paste("enriched", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(res()[[2]], file, row.names = FALSE)
    }
  )

  output$download_plot <- downloadHandler(
    filename = function() {
      paste(chosen_plt()[[1]]$labels$title, ".png",
        sep = ""
      )
    },
    content = function(file) {
      ggsave(file,
        plot = chosen_plt()[[1]], device = "png",
        height = (45 / 74) * chosen_plt()[[2]], width = 6.75,
        limitsize = FALSE
      )
    }
  )
}

shinyApp(ui, server)
