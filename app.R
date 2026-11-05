# Shiny app: selezione guidata di tool per inferenza CNA da scRNA-seq
# Avvio: shiny::runApp("~/Downloads/files")
library(shiny)
library(DT)

# Carica le funzioni di supporto da R/engine.R (necessario per QUESTIONS, run_engine, ecc.)
source("R/engine.R", local = TRUE)

# Definizione della funzione safe per il markdown (deve stare prima di ui)
includeMarkdown_safe <- function() HTML("
<p><b>Esclusione</b>: il tool non può soddisfare un tuo vincolo (input richiesto che non hai, output non prodotto, specie non supportata...).</p>
<p><b>Avvertenza</b>: il tool è utilizzabile, ma il tuo scenario tocca un limite noto per ragioni metodologiche (es. segnale allelico sparso con poche cellule).</p>
<p><b>Score</b> = (n. bonus) &minus; 0.5 &times; (n. avvertenze). I bonus premiano i tool che sfruttano informazioni che hai. <u>Non</u> è una misura di accuratezza: i benchmark non entrano nel ranking.</p>
<p><b>Dato mancante</b> (<code>unknown</code>) = informazione non ricavabile dalla review, da completare in <code>data/tools_kb.csv</code>.</p>")

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
  titlePanel("Quale tool per inferire CNA da scRNA-seq?"),
  sidebarLayout(
    sidebarPanel(width = 4, style = "max-height:92vh; overflow-y:auto;",
                 lapply(stages, function(s) tagList(h4(s), lapply(Filter(function(q) q$stage == s, QUESTIONS), make_input), hr())),
                 checkboxInput("show_bench", "Mostra evidenze da benchmark (solo informative, non influenzano il ranking)", FALSE),
                 actionButton("reset", "Reset risposte")),
    mainPanel(width = 8, tabsetPanel(
      tabPanel("Raccomandazione", br(), uiOutput("top"), h4("Tutti i tool compatibili"),
               DTOutput("ok_tbl"), br(), downloadButton("dl", "Scarica (CSV)")),
      tabPanel("Percorso decisionale", br(),
               p("Numero di tool ancora compatibili dopo ogni domanda."),
               plotOutput("funnel", height = "520px"), DTOutput("trace_tbl")),
      tabPanel("Tool esclusi", br(), DTOutput("excl_tbl")),
      tabPanel("Evidenze benchmark", br(), uiOutput("bench_msg"), DTOutput("bench_tbl")),
      tabPanel("Knowledge base", br(),
               p("Attributi dei tool usati dalle regole. 'unknown' = dato da verificare: il tool NON viene escluso, ma compare un'avvertenza."),
               DTOutput("kb_tbl")),
      tabPanel("Come leggere i risultati", br(), includeMarkdown_safe())
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
                              "Nessun tool soddisfa tutti i vincoli. Guarda 'Tool esclusi' per i motivi e valuta quale vincolo rilassare."))
    top <- ok[ok$Score == max(ok$Score), ]
    div(h4("Tool consigliati (miglior corrispondenza con i tuoi vincoli)"),
        lapply(seq_len(nrow(top)), function(i) div(class = "panel panel-success",
                                                   div(class = "panel-heading", strong(top$Tool[i]), span(style = "float:right", top$Classe[i])),
                                                   div(class = "panel-body",
                                                       if (nzchar(top$Perche[i])) p(strong("Punti a favore: "), top$Perche[i]),
                                                       if (nzchar(top$Avvertenze[i])) tagList(strong("Avvertenze:"),
                                                                                              tags$ul(lapply(strsplit(top$Avvertenze[i], "; ", fixed = TRUE)[[1]], tags$li)))
                                                       else p("Nessuna avvertenza per questo scenario.")))))
  })
  output$ok_tbl   <- renderDT(datatable(tabs()$ok,   rownames = FALSE, options = list(pageLength = 15, dom = "tip")))
  output$excl_tbl <- renderDT(datatable(tabs()$excl, rownames = FALSE, options = list(pageLength = 15, dom = "tip")))
  output$trace_tbl <- renderDT(datatable(res()$trace, rownames = FALSE, options = list(pageLength = 25, dom = "tp")))
  output$kb_tbl   <- renderDT(datatable(kb, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15)))
  
  output$funnel <- renderPlot({
    tr <- res()$trace
    op <- par(mar = c(4, 14, 2, 1)); on.exit(par(op))
    bp <- barplot(rev(tr$remaining), names.arg = rev(tr$step), horiz = TRUE, las = 1, cex.names = .85,
                  xlim = c(0, nrow(kb) + 1), xlab = "Tool compatibili", col = "#4C78A8", border = NA)
    text(rev(tr$remaining), bp, rev(tr$remaining), pos = 4, cex = .85)
  })
  
  output$bench_msg <- renderUI(if (!isTRUE(input$show_bench))
    p("Attiva 'Mostra evidenze da benchmark' nella barra laterale. Le evidenze sono mostrate solo per i tool ancora compatibili e NON modificano esclusioni né ranking."))
  output$bench_tbl <- renderDT({
    req(input$show_bench)
    b <- bench[bench$tool %in% c(kb$tool[res()$alive], "all"), ]
    datatable(b, rownames = FALSE, options = list(pageLength = 20, dom = "tip"))
  })
  
  output$dl <- downloadHandler("cna_tool_ranking.csv", function(file) write.csv(tabs()$ok, file, row.names = FALSE))
  observeEvent(input$reset, session$reload())
}

shinyApp(ui, server)