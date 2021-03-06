options(stringsAsFactors=FALSE)
options(shiny.maxRequestSize=100*1024^2)
## NOTE: library folder is "./rlib/4.0" (check .Rprofile)
library(shiny)
library(shinyjs)
library(shinydashboard)
library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)

# define available chains
chains <- data.frame(rbind(
  c(org  ="human",
    old  ="hg19",
    new  ="hg38",
    chain="hg19ToHg38.over.chain",
    bsgenome.old="BSgenome.Hsapiens.UCSC.hg19",
    bsgenome.new="BSgenome.Hsapiens.UCSC.hg38"),
  c(org  ="human",
    old  ="hg38",
    new  ="hg19",
    chain="hg38ToHg19.over.chain",
    bsgenome.old="BSgenome.Hsapiens.UCSC.hg38",
    bsgenome.new="BSgenome.Hsapiens.UCSC.hg19")
))

# define available formats
formats <- list("bed"     =c(ext=".bed"     , fun="export.bed"),
                "bedGraph"=c(ext=".bedGraph", fun="export.bedGraph"),
                "bigWig"  =c(ext=".bw"      , fun="export.bw"),
                "gff"     =c(ext=".gff"     , fun="export.gff"),
                "wig"     =c(ext=".wig"     , fun="export.wig"))

# Define UI for application that draws a histogram
ui <- dashboardPage(
  dashboardHeader(title="liftOver coordinates", titleWidth=300),
  dashboardSidebar(disable=TRUE),
  dashboardBody(
    useShinyjs(),
    fluidRow(
      column(6,
             box(width=NULL, title="", status="warning",
                 fileInput("file1", "Choose track file",
                           multiple=FALSE,
                           accept=c(".gff", ".bed", "bedpe", "bedGraph", ".bw", ".wig")),
                 hr(),
                 checkboxInput("liftover", "Lift over", value=TRUE),
                 fluidRow(
                   column(4, selectInput("org", "Organism:", "")),
                   column(4, selectInput("old", "Original assembly:", "")),
                   column(4, selectInput("new", "New assembly:", ""))
                 ),
                 hr(),
                 fluidRow(
                   column(6, selectInput("format", "Output format:", names(formats))),
                   column(3, actionButton("clear", label="Clear")),
                   column(3, actionButton("convert", label="Convert!"))
                 )
             )
      ),
      column(6,
        box(width=NULL, title="", status="warning",
          fluidRow(
            column(4, downloadButton("downloadData", "Download converted file")),
            column(8, verbatimTextOutput("downloadMsg"))
          ),
          hr(),
          tags$b("Input coordinates:"),
          verbatimTextOutput("previewInput"),
          hr(),
          tags$b("Output coordinates:"),
          verbatimTextOutput("previewOutput")
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  session$onSessionEnded(stopApp)
  
  # enable chains if user wants to liftover
  observe({
    toggleState("org", input$liftover)
    toggleState("old", input$liftover)
    toggleState("new", input$liftover)
  })
  
  # update available chains depending on the organism selected
  observe({
    x <- unique(chains$org)
    updateSelectInput(session, "org", choices=x, selected=x[1])
  })
  observe({
    x <- unique(chains$old[chains$org == input$org])
    updateSelectInput(session, "old", choices=x, selected=x[1])
  })
  observe({
    x <- unique(chains$new[chains$org == input$org & chains$old == input$old])
    updateSelectInput(session, "new", choices=x, selected=x[1])
  })
  
  # convert track to new assembly
  track.new <- reactiveVal(NULL)
  observeEvent(input$convert, {
    isolate({
      req(track.old())
      if(input$liftover) {
        withProgress(message="Lifting over", value=0, {
          i <- chains$org == input$org & chains$old == input$old & chains$new == input$new
  
          # load bsgenome info about the to/from genomes
          genome.old <- eval(parse(text=paste0(chains$bsgenome.old[i], "::", chains$bsgenome.old[i])))
          genome.new <- eval(parse(text=paste0(chains$bsgenome.new[i], "::", chains$bsgenome.new[i])))
          
          # add genome info to/from track
          x.old <- track.old()
          seqlevelsStyle(x.old) <- "UCSC"
          seqlevels(x.old) <- seqlevels(genome.old)
          seqinfo(x.old)   <- seqinfo(genome.old)
          
          # convert
          x.new <- unlist(liftOver(x.old, rtracklayer::import.chain(as.character(chains$chain[i]))))
          seqlevels(x.new) <- seqlevels(genome.new)
          seqinfo(x.new)   <- seqinfo(genome.new)
          
          # drop the overlapping ranges
          hits <- findOverlaps(x.new, drop.self=TRUE)
          if(length(hits) > 0) {
            x.new <- x.new[-queryHits(hits)]
          }
          if("score" %in% colnames(x.new)) {
            x.new <- x.new[!is.na(x.new$score)] 
          }
        })
      } else {
        withProgress(message="Converting", value=0, {
          x.new <- track.old()
          
          # drop the overlapping ranges
          hits <- findOverlaps(x.new, drop.self=TRUE)
          if(length(hits) > 0) {
            x.new <- x.new[-queryHits(hits)]
          }
          if("score" %in% colnames(x.new)) {
            x.new <- x.new[!is.na(x.new$score)] 
          }
        })
      }
      track.new(x.new)
    })
  })
  
  # upload input
  track.old <- reactiveVal(NULL)
  observe({
    req(input$file1)
    tryCatch( {
      track.old(rtracklayer::import(input$file1$datapath))  # in principle, rtracklayer will detect the format automatically
    },
    error=function(e) {
      stop(safeError(e))
    })
  })
  
  # clear control
  observeEvent(input$clear, {
    track.old(NULL)
    track.new(NULL)
  })
  
  # download
  output$downloadData <- downloadHandler(
    filename=function() {
      if(input$liftover) {
        paste0(sub("\\..+", "", input$file1$name), "_", input$new, formats[[input$format]]["ext"])
      } else {
        paste0(sub("\\..+", "", input$file1$name), formats[[input$format]]["ext"])
      }
    },
    content=function(file) {
      req(track.new())
      do.call(formats[[input$format]]["fun"], list(track.new(), file))
    }
  )
  
  # previews
  output$previewInput <- renderPrint({
    req(track.old())
    print(track.old())
  })
  
  output$previewOutput <- renderPrint({
    req(track.new())
    print(track.new())
  })
}

# Run the application 
shinyApp(ui=ui, server=server)
