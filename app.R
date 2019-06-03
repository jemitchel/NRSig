library(shiny)

ui <- fluidPage(
  fileInput(inputId = "files",
            label = "Select all .CEL files to process",
            multiple = TRUE),
  verbatimTextOutput("fileName"),
  actionButton(inputId = "compute",
               label = "Compute"),
  tableOutput("table")
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  
  output$fileName <- renderText({
    tmp <- "Uploaded Files\n"
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
      return(results[1])
    })
    
  })
  
  
  output$table <- renderTable(data())
  
}

shinyApp(ui, server)

