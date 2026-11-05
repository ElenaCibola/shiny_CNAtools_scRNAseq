# Shiny app: guided tool selection for CNA inference from scRNA-seq
# Launch: shiny::runApp("~/Downloads/files")
library(shiny)
library(DT)

# Load helper functions from R/engine.R (needed for QUESTIONS, run_engine, etc.)
source("R/engine.R", local = TRUE)

# Definition of the safe markdown function (must be placed before ui)
includeMarkdown_safe <- function() HTML("
<p><b>Exclusion</b>: the tool cannot satisfy one of your constraints (required input you don't have, output not produced, unsupported species...).</p>
<p><b>Warning</b>: the tool is usable, but your scenario touches a known limit for methodological reasons (e.g. sparse allelic signal with few cells).</p>
<p><b>Score</b> = (no. of bonuses) &minus; 0.5 &times; (no. of warnings). Bonuses reward tools that exploit information you have. <u>Not</u> a measure of accuracy: benchmarks do not enter the ranking.</p>
<p><b>Missing data</b> (<code>unknown</code>) = information that cannot be derived from the review, to be completed in <code>data/tools_kb.csv</code>.</p>")

kb    <- read.csv("data/tools_kb.csv", na.strings = c("", "NA"), stringsAsFactors = FALSE)
bench <- read.csv("data/benchmark_evidence.csv", stringsAsFactors = FALSE)

stages <- unique(vapply(QUESTIONS, `[[`, "", "stage"))

make_input <- function(q) {
  w <- switch(q$type,
              radio  = radioButtons(q$id, q$label, choices = q$choices, selected = q$choices[[1]]),
              multi  = checkboxGroupInput(q$id, q$label, choices = q$choices),
              select = selectInput(q$id, q$label, choices = q$choices),
              num    = numericInput(q$id, q$label, value = q$value, min = 0, step = 100),
              chk    = checkboxInput(q$id, q$label, value = isTRUE(q$value)))
  tagList(w, if (!is.null(q$help)) div(class = "help-block", style = "margin-top:-8px", q$help))
}

ui <- fluidPage(
  titlePanel("Which tool to infer CNAs from scRNA-seq?"),
  sidebarLayout(
    sidebarPanel(width = 4, style = "max-height:92vh; overflow-y:auto;",
                 lapply(stages, function(s) tagList(h4(s), lapply(Filter(function(q) q$stage == s, QUESTIONS), make_input), hr())),
                 checkboxInput("show_bench", "Show benchmark evidence (informative only, does not affect the ranking)", FALSE),
                 actionButton("reset", "Reset answers")),
    mainPanel(width = 8, tabsetPanel(
      tabPanel("Recommendation", br(), uiOutput("top"), h4("All compatible tools"),
               DTOutput("ok_tbl"), br(), downloadButton("dl", "Download (CSV)")),
      tabPanel("Decision path", br(),
               p("Number of tools still compatible after each question."),
               plotOutput("funnel", height = "520px"), DTOutput("trace_tbl")),
      tabPanel("Excluded tools", br(), DTOutput("excl_tbl")),
      tabPanel("Benchmark evidence", br(), uiOutput("bench_msg"), DTOutput("bench_tbl")),
      tabPanel("Knowledge base", br(),
               p("Tool attributes used by the rules. 'unknown' = data to be verified: the tool is NOT excluded, but a warning appears."),
               DTOutput("kb_tbl")),
      tabPanel("How to read the results", br(), includeMarkdown_safe())
    ))
  )
)

server <- function(input, output, session) {
  ids <- vapply(QUESTIONS, `[[`, "", "id")
  answers <- reactive(setNames(lapply(ids, function(i) input[[i]]), ids))
  res  <- reactive(run_engine(kb, answers()))
  tabs <- reactive(summarise_results(kb, res()))
  
  output$top <- renderUI({
    ok <- tabs()$ok
    if (!nrow(ok)) return(div(class = "alert alert-danger",
                              "No tool satisfies all the constraints. See 'Excluded tools' for the reasons and consider which constraint to relax."))
    top <- ok[ok$Score == max(ok$Score), ]
    div(h4("Recommended tools (best match with your constraints)"),
        lapply(seq_len(nrow(top)), function(i) div(class = "panel panel-success",
                                                   div(class = "panel-heading", strong(top$Tool[i]), span(style = "float:right", top$Classe[i])),
                                                   div(class = "panel-body",
                                                       if (nzchar(top$Perche[i])) p(strong("Pros: "), top$Perche[i]),
                                                       if (nzchar(top$Avvertenze[i])) tagList(strong("Warnings:"),
                                                                                              tags$ul(lapply(strsplit(top$Avvertenze[i], "; ", fixed = TRUE)[[1]], tags$li)))
                                                       else p("No warnings for this scenario.")))))
  })
  output$ok_tbl   <- renderDT(datatable(tabs()$ok,   rownames = FALSE, options = list(pageLength = 15, dom = "tip")))
  output$excl_tbl <- renderDT(datatable(tabs()$excl, rownames = FALSE, options = list(pageLength = 15, dom = "tip")))
  output$trace_tbl <- renderDT(datatable(res()$trace, rownames = FALSE, options = list(pageLength = 25, dom = "tp")))
  output$kb_tbl   <- renderDT(datatable(kb, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15)))
  
  output$funnel <- renderPlot({
    tr <- res()$trace
    op <- par(mar = c(4, 14, 2, 1)); on.exit(par(op))
    bp <- barplot(rev(tr$remaining), names.arg = rev(tr$step), horiz = TRUE, las = 1, cex.names = .85,
                  xlim = c(0, nrow(kb) + 1), xlab = "Compatible tools", col = "#4C78A8", border = NA)
    text(rev(tr$remaining), bp, rev(tr$remaining), pos = 4, cex = .85)
  })
  
  output$bench_msg <- renderUI(if (!isTRUE(input$show_bench))
    p("Enable 'Show benchmark evidence' in the sidebar. The evidence is shown only for tools that are still compatible and does NOT change exclusions or ranking."))
  output$bench_tbl <- renderDT({
    req(input$show_bench)
    b <- bench[bench$tool %in% c(kb$tool[res()$alive], "all"), ]
    datatable(b, rownames = FALSE, options = list(pageLength = 20, dom = "tip"))
  })
  
  output$dl <- downloadHandler("cna_tool_ranking.csv", function(file) write.csv(tabs()$ok, file, row.names = FALSE))
  observeEvent(input$reset, session$reload())
}

shinyApp(ui, server)
