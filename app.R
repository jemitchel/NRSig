library(shiny)
library(DT)

ui <- fluidPage(

  sidebarLayout(

    sidebarPanel(
      fileInput(inputId = "files",
                label = "Upload .CEL files for one sample or multiple replicates (replicates will be averaged)",
                multiple = TRUE),
      verbatimTextOutput("fileName"),
      actionButton(inputId = "compute",
                   label = "Compute")
    ),

    mainPanel(
      DTOutput("table")
    )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)

  output$fileName <- renderText({
    tmp <- "Uploaded Files:\n"
    if (!is.null(input$files$datapath)) {
      for (i in basename(input$files$name)) {
        tmp <- paste(tmp,i,"\n")
      }
      return(tmp)
    } else {
      return(NULL)
    }
  })


  data <- eventReactive(input$compute, {
    progress <- Progress$new(session)
    progress$set(message = 'Preprocessing',
                 detail = 'This may take a few minutes...')
    source("C:/Users/jonat/Documents/R/NRSig-app/preprocess.R")
    samples_matrix <- pre_proc(input$files$datapath)
    progress$set(message = 'Preprocessing Completed', detail = "")
    on.exit(progress$close())

    withProgress(message = 'Computing Enrichment of NR-Target Gene Sets...', value = 0,{
      source("C:/Users/jonat/Documents/R/NRSig-app/compute_enriched_NRs.R")
      results <- CalcEnrich(samples_matrix)
      return(results[[2]])
    })

  })

  output$table <- renderDT(data())

}

shinyApp(ui, server)



