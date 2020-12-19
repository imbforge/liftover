options(stringsAsFactors = FALSE)
library(shiny)
library(shinydashboard)
## NOTE: this will go away once the Rstudio-Shiny container is fixed
library(rtracklayer, lib.loc=paste0(getwd(), "/rlib"))
library(BSgenome, lib.loc=paste0(getwd(), "/rlib"))
library(BSgenome.Hsapiens.UCSC.hg19, lib.loc=paste0(getwd(), "/rlib"))
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc=paste0(getwd(), "/rlib"))
#library(rtracklayer)
## END NOTE

# define available chains
chains <- data.frame(rbind(
  c(org  ="human",
    from ="hg19",
    to   ="hg38",
    chain="hg19ToHg38.over.chain",
    bsgenomefrom="BSgenome.Hsapiens.UCSC.hg19",
    bsgenometo  ="BSgenome.Hsapiens.UCSC.hg38"),
  c(org  ="human",
    from ="hg38",
    to   ="hg19",
    chain="hg38ToHg19.over.chain",
    bsgenomefrom="BSgenome.Hsapiens.UCSC.hg38",
    bsgenometo  ="BSgenome.Hsapiens.UCSC.hg19")
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
    fluidRow(
      column(6,
             box(width=NULL, title="", status="warning",
                 fileInput("file1", "Choose track file",
                           multiple = FALSE,
                           accept = c(".gff", ".bed", "bedpe", "bedGraph", ".bw", ".wig")),
                 hr(),
                 fluidRow(
                   column(4, selectInput("org", "Organism:", "")),
                   column(4, selectInput("from", "From:", "")),
                   column(4, selectInput("to", "To:", ""))
                 ),
                 hr(),
                 fluidRow(
                   column(8, selectInput("format", "Output format:", names(formats))),
                   column(4, actionButton("convert", label = "Convert!"))
                 )
             )
      ),
      column(6,
        box(width=NULL, title="", status="warning",
          downloadButton("downloadData", "Download converted file"),
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
  
  # update available chains depending on the organism selected
  observe({
    x <- unique(chains$org)
    updateSelectInput(session, "org", choices=x, selected=x[1])
  })
  observe({
    x <- unique(chains$from[chains$org == input$org])
    updateSelectInput(session, "from", choices=x, selected=x[1])
  })
  observe({
    x <- unique(chains$to[chains$org == input$org & chains$from == input$from])
    updateSelectInput(session, "to", choices=x, selected=x[1])
  })
  
  trackConverted <- reactive({
    input$convert
    
    isolate({
      req(track())
      withProgress(message = "Lifting over", value = 0, {
        i <- chains$org == input$org & chains$from == input$from & chains$to == input$to

        # load bsgenome info about the to/from genomes
        genomefrom <- eval(parse(text=paste0(chains$bsgenomefrom[i], "::", chains$bsgenomefrom[i])))
        genometo   <- eval(parse(text=paste0(chains$bsgenometo[i], "::", chains$bsgenometo[i])))
        
        # add genome info tofrom track
        xfrom <- track()
        seqlevelsStyle(xfrom) <- "UCSC"
        seqlevels(xfrom) <- seqlevels(genomefrom)
        seqinfo(xfrom)   <- seqinfo(genomefrom)
        
        # convert
        xto <- unlist(liftOver(xfrom, rtracklayer::import.chain(as.character(chains$chain[i]))))
        seqlevels(xto) <- seqlevels(genometo)
        seqinfo(xto)   <- seqinfo(genometo)
        
        # drop the overlapping ranges
        hits <- findOverlaps(xto, drop.self=TRUE)
        xto <- xto[-queryHits(hits)]
        if("score" %in% colnames(xto)) {
          xto <- xto[!is.na(xto$score)] 
        }
        xto
      })
    })
  })
  
  # upload input
  track <- reactive({
    req(input$file1)
    tryCatch( {
      rtracklayer::import(input$file1$datapath)  # in principle, rtracklayer will detect the format automatically
    },
    error = function(e) {
      stop(safeError(e))
    })
  })
  
  # download
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(sub("\\..+", "", input$file1$name), "_", input$to, formats[[input$format]]["ext"])
    },
    content = function(file) {
      do.call(formats[[input$format]]["fun"], list(trackConverted(), file))
    }
  )
  
  # previews
  output$previewInput <- renderPrint({
    req(track())
    print(track())
  })
  
  output$previewOutput <- renderPrint({
    req(trackConverted())
    print(trackConverted())
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
