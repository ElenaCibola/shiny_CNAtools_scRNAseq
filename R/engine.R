# =====================================================================
# Motore di regole per la scelta dei tool CNA da scRNA-seq
# - NON usa i risultati dei benchmark per escludere/ordinare (gap noti)
# - Esclusioni = vincoli "di progetto" del tool (input richiesti, output, specie...)
# - Avvertenze = caveat meccanicistici, non punteggi di performance
# - Bonus = il tool sfrutta un'informazione che l'utente ha a disposizione
# Per aggiungere una domanda: aggiungere un elemento a QUESTIONS.
# =====================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
is_ <- function(x, v) !is.na(x) & x %in% v           # NA-safe
unk <- function(x) is.na(x) | x == "unknown"          # dato mancante nella KB
fmt <- function(a) if (is.null(a) || !length(a) || all(is.na(a))) "-" else paste(a, collapse = ", ")

has_out <- function(kb, vals) vapply(kb$output, function(o)
  !is.na(o) && any(trimws(strsplit(o, ";")[[1]]) %in% vals), logical(1), USE.NAMES = FALSE)

# regola: ex = esclusione, wn = avvertenza, bo = bonus (vettori di lunghezza n tool)
rule <- function(kb, ex = FALSE, why = "", wn = FALSE, wmsg = "", bo = FALSE, bmsg = "") {
  n <- nrow(kb)
  list(ex = rep_len(ex, n), why = rep_len(why, n), wn = rep_len(wn, n),
       wmsg = rep_len(wmsg, n), bo = rep_len(bo, n), bmsg = rep_len(bmsg, n))
}
miss <- function(lbl) paste0(lbl, ": dato mancante nella KB, verificare")

QUESTIONS <- list(

  # ------------------------------------------------------------ 1. OBIETTIVO
  list(id = "goals", short = "Obiettivi", stage = "1. Obiettivo dell'analisi", type = "multi",
    label = "Cosa ti serve in output? (seleziona tutto ciò che serve)",
    choices = c("Classificazione cellule tumorali vs normali" = "classification",
                "Sottocloni definiti da CNA" = "subclones",
                "Albero filogenetico / storia evolutiva" = "phylogeny",
                "Stati allele-specifici (cnLOH, LOH)" = "cnloh",
                "Distinguere amplificazioni bilanciate/sbilanciate" = "balanced"),
    help = "Ogni voce esclude i tool che non la implementano.",
    apply = function(kb, a, all) {
      a <- a %||% character(0); out <- list()
      if ("classification" %in% a)
        out[[length(out) + 1]] <- rule(kb,
          ex = is_(kb$classification, "none"), why = "non produce classificazione tumore/normale",
          wn = is_(kb$classification, "threshold") | unk(kb$classification),
          wmsg = ifelse(is_(kb$classification, "threshold"),
                        "fornisce solo uno score per cellula: serve scegliere una soglia",
                        miss("classificazione")))
      feat <- c(subclones = "sottocloni", phylogeny = "filogenesi",
                cnloh = "stati allele-specifici (cnLOH/LOH)", balanced = "amplificazioni bilanciate/sbilanciate")
      col  <- c(subclones = "subclones", phylogeny = "phylogeny", cnloh = "cnloh", balanced = "balanced_amp")
      for (k in intersect(names(feat), a))
        out[[length(out) + 1]] <- rule(kb, ex = is_(kb[[col[[k]]]], "no"),
          why = paste("non implementa:", feat[[k]]),
          wn = unk(kb[[col[[k]]]]), wmsg = miss(feat[[k]]))
      if (length(out)) out else rule(kb)
    }),

  list(id = "output", short = "Tipo di output", stage = "1. Obiettivo dell'analisi", type = "radio",
    label = "Tipo di output richiesto",
    choices = c("Nessuna preferenza" = "any",
                "Stati discreti (anche solo gain/neutral/loss)" = "discrete",
                "Stati discreti a più livelli (deletion/loss/neutral/gain/amp)" = "multi"),
    help = "Gli score continui richiedono una soglia scelta dall'utente.",
    apply = function(kb, a, all) {
      if (a == "any") return(rule(kb))
      vals <- if (a == "discrete") c("3-state", "multi-state", "allele-specific") else c("multi-state", "allele-specific")
      rule(kb, ex = !unk(kb$output) & !has_out(kb, vals), why = "output solo continuo / livelli insufficienti",
           wn = unk(kb$output), wmsg = miss("tipo di output"))
    }),

  list(id = "level", short = "Livello della chiamata", stage = "1. Obiettivo dell'analisi", type = "radio",
    label = "Livello a cui vuoi la chiamata CNA",
    choices = c("Indifferente" = "any", "Per singola cellula" = "cell", "Per clone/sottoclone" = "clone"),
    apply = function(kb, a, all) {
      if (a == "any") return(rule(kb))
      bad <- if (a == "cell") "clone" else "cell"
      rule(kb, ex = is_(kb$call_level, bad), why = paste("la chiamata CNA è solo a livello", bad),
           wn = unk(kb$call_level), wmsg = miss("livello della chiamata"))
    }),

  list(id = "res", short = "Risoluzione", stage = "1. Obiettivo dell'analisi", type = "radio",
    label = "Risoluzione genomica necessaria",
    choices = c("Braccio cromosomico sufficiente" = "arm", "Segmenti / sub-braccio" = "segment",
                "Eventi focali (pochi Mb)" = "focal"),
    help = "Tutti i tool applicano smoothing: più smoothing = più robustezza ma meno risoluzione.",
    apply = function(kb, a, all) {
      if (a == "arm") return(rule(kb))
      out <- list(rule(kb, ex = is_(kb$resolution, "arm"), why = "chiama CNA solo per braccio cromosomico"))
      if (a == "focal") out[[2]] <- rule(kb, wn = TRUE,
        wmsg = "eventi focali: lo smoothing su finestre di geni ne riduce la rilevabilità")
      out
    }),

  list(id = "k", short = "N. sottocloni noto?", stage = "1. Obiettivo dell'analisi", type = "radio",
    label = "(Se vuoi sottocloni) conosci a priori il numero di sottocloni?",
    choices = c("No / non lo so" = "unknown", "Sì, posso specificarlo" = "known"),
    apply = function(kb, a, all) {
      if (a == "known" || !("subclones" %in% (all$goals %||% character(0)))) return(rule(kb))
      rule(kb, ex = is_(kb$k_mode, "apriori"), why = "richiede il numero di sottocloni a priori",
           wn = unk(kb$k_mode), wmsg = miss("gestione del numero di sottocloni"))
    }),

  # ------------------------------------------------------------ 2. DATI
  list(id = "format", short = "Formato dati", stage = "2. Dati disponibili", type = "radio",
    label = "Che dati hai?",
    choices = c("BAM/FASTQ disponibili (conteggi allelici possibili)" = "bam",
                "Solo matrice gene x cellula" = "matrix"),
    help = "I tool allele-aware estraggono REF/ALT per SNP e cellula dal BAM.",
    apply = function(kb, a, all) {
      if (a == "bam") return(rule(kb))
      rule(kb, ex = is_(kb$allele_input, "required"), why = "richiede BAM per i conteggi allelici")
    }),

  list(id = "phasing", short = "Pannello phasing", stage = "2. Dati disponibili", type = "radio",
    label = "Hai un pannello di riferimento per il phasing (popolazione) adatto alla specie/ancestry?",
    choices = c("Sì" = "yes", "No" = "no", "Non so" = "unknown"),
    help = "Alcuni tool importano aplotipi da pannelli di popolazione invece di fasare dai dati scRNA-seq.",
    apply = function(kb, a, all) {
      if (a == "yes") return(rule(kb))
      if (a == "no") return(rule(kb, ex = is_(kb$phasing_panel, "yes"), why = "richiede un pannello di phasing esterno"))
      rule(kb, wn = is_(kb$phasing_panel, "yes"), wmsg = "richiede un pannello di phasing esterno: verificarne la disponibilità")
    }),

  list(id = "bulk", short = "DNA bulk abbinato", stage = "2. Dati disponibili", type = "radio",
    label = "Hai DNA bulk abbinato (WGS/WES)?",
    choices = c("No" = "no", "Sì" = "yes"),
    apply = function(kb, a, all) {
      if (a == "yes") return(rule(kb, bo = is_(kb$bulk_dna, c("optional", "recommended")),
                                  bmsg = "può integrare breakpoint/profili da DNA bulk"))
      rule(kb, wn = is_(kb$bulk_dna, "recommended"), wmsg = "prestazioni ottimali con WGS abbinato, che non hai")
    }),

  list(id = "ref", short = "Riferimento normale", stage = "2. Dati disponibili", type = "radio",
    label = "Riferimento copy-neutral a priori",
    choices = c("No, nessuna annotazione a priori" = "none",
                "Sì, cellule normali annotate nello stesso campione" = "insample",
                "Sì, riferimento esterno (altro dataset/atlante)" = "external"),
    help = "Senza annotazioni restano solo i tool che identificano il riferimento in automatico.",
    apply = function(kb, a, all) {
      if (a != "none") return(rule(kb))
      rule(kb, ex = is_(kb$ref_mode, "required"), why = "richiede un riferimento normale fornito dall'utente",
           wn = unk(kb$ref_mode), wmsg = miss("gestione del riferimento"))
    }),

  list(id = "normal", short = "Cellule normali nel campione", stage = "2. Dati disponibili", type = "radio",
    label = "Il campione contiene cellule non tumorali (immuni/stromali)?",
    choices = c("Sì / probabile" = "yes", "No (es. linea cellulare pura)" = "no", "Non so" = "unknown"),
    help = "Conta solo se non hai un riferimento esterno.",
    apply = function(kb, a, all) {
      if (a == "yes" || (all$ref %||% "none") != "none") return(rule(kb))
      rule(kb, ex = is_(kb$runs_without_normal, "no") & a == "no",
           why = "senza cellule normali né riferimento esterno non può stimare la baseline",
           wn = unk(kb$runs_without_normal),
           wmsg = miss("funzionamento senza cellule normali"))
    }),

  # ------------------------------------------------------------ 3. CAMPIONE
  list(id = "species", short = "Specie", stage = "3. Campione e sequenziamento", type = "radio",
    label = "Specie", choices = c("Umano" = "human", "Topo" = "mouse", "Altra" = "other"),
    help = "Molti tool usano annotazioni geniche predefinite (umano/topo, poche versioni di genoma).",
    apply = function(kb, a, all) {
      v <- kb[[paste0("sp_", a)]]
      rule(kb, ex = is_(v, "no"), why = paste("specie non supportata:", a),
           wn = unk(v) | is_(v, "limited"),
           wmsg = ifelse(is_(v, "limited"), "supporto limitato per questa specie (SNP eterozigoti + phasing)", miss("supporto specie")))
    }),

  list(id = "het", short = "Eterozigosità", stage = "3. Campione e sequenziamento", type = "radio",
    label = "L'organismo/ceppo ha sufficiente variabilità eterozigote germinale?",
    choices = c("Sì (es. umano, topo outbred)" = "yes", "No (es. ceppo inbred)" = "no"),
    apply = function(kb, a, all) {
      if (a == "yes") return(rule(kb))
      rule(kb, ex = is_(kb$allele_input, "required"), why = "servono SNP eterozigoti informativi")
    }),

  list(id = "sample", short = "Tipo di campione", stage = "3. Campione e sequenziamento", type = "select",
    label = "Tipo di campione",
    choices = c("Tumore solido primario" = "primary", "Metastasi" = "met", "Ematologico" = "heme",
                "Linea cellulare" = "line", "Organoide" = "organoid", "PDX / xenotrapianto" = "pdx"),
    apply = function(kb, a, all) {
      if (a == "heme")
        return(rule(kb, wn = is_(kb$ref_mode, "required") & (all$ref %||% "none") != "external",
          wmsg = "ematologico: il riferimento tipico (cellule immuni) non è copy-neutral; meglio un riferimento esterno"))
      if (a == "pdx")
        return(rule(kb, wn = TRUE, wmsg = "PDX: separare le letture umane da quelle murine prima dell'analisi"))
      rule(kb)
    }),

  list(id = "purity", short = "Purezza tumorale", stage = "3. Campione e sequenziamento", type = "radio",
    label = "Frazione di cellule tumorali",
    choices = c("Bilanciata / non so" = "mid", "Alta (tumore dominante)" = "high", "Bassa (tumore minoritario)" = "low"),
    apply = function(kb, a, all) {
      if (a == "high" && (all$ref %||% "none") == "none")
        return(rule(kb, wn = is_(kb$ref_mode, c("auto", "optional")) & kb$signal == "expression",
          wmsg = "poche cellule normali: il riferimento automatico basato su espressione può risultare mal centrato"))
      if (a == "low") return(rule(kb, wn = TRUE, wmsg = "frazione tumorale bassa: segnale CNA diluito, sottocloni rari poco rilevabili"))
      rule(kb)
    }),

  list(id = "ploidy", short = "Ploidia / eventi copy-neutral", stage = "3. Campione e sequenziamento", type = "chk",
    label = "Atteso WGD / aneuploidia estesa / eventi copy-neutral",
    value = FALSE,
    apply = function(kb, a, all) {
      if (!isTRUE(a)) return(rule(kb))
      rule(kb, wn = kb$signal == "expression",
           wmsg = "solo espressione: stesso CN totale con diversa composizione allelica è indistinguibile",
           bo = kb$signal %in% c("allele", "both"), bmsg = "il segnale allelico risolve configurazioni con CN totale invariato")
    }),

  list(id = "batch", short = "Batch / multi-piattaforma", stage = "3. Campione e sequenziamento", type = "chk",
    label = "Più batch o più piattaforme da analizzare insieme", value = FALSE,
    apply = function(kb, a, all) {
      if (!isTRUE(a)) return(rule(kb))
      rule(kb, wn = TRUE, wmsg = "batch: correggere sui conteggi (es. ComBat); i metodi di integrazione scRNA restituiscono embedding, non conteggi")
    }),

  list(id = "platform", short = "Piattaforma", stage = "3. Campione e sequenziamento", type = "select",
    label = "Piattaforma di sequenziamento",
    choices = c("Droplet 3' (es. 10x 3')" = "droplet3", "Droplet 5' / altro tag-based" = "droplet5",
                "Plate-based full-length (es. Smart-seq)" = "plate", "Altro / misto" = "other"),
    apply = function(kb, a, all) {
      if (a == "droplet3")
        return(rule(kb, wn = kb$signal %in% c("allele", "both"),
          wmsg = "protocollo 3': letture concentrate verso il 3', meno SNP coperti"))
      rule(kb)
    }),

  list(id = "cells", short = "N. cellule", stage = "3. Campione e sequenziamento", type = "num",
    label = "Numero di cellule", value = 3000,
    apply = function(kb, a, all) {
      if (is.null(a) || is.na(a)) return(rule(kb))
      rule(kb, wn = kb$signal %in% c("allele", "both") & a < 500,
           wmsg = "poche cellule: il segnale allelico si aggrega tra cellule/SNP, pochi SNP informativi")
    }),

  list(id = "depth", short = "Profondità", stage = "3. Campione e sequenziamento", type = "select",
    label = "Profondità per cellula",
    choices = c("Media / non so" = "mid", "Bassa" = "low", "Alta" = "high"),
    apply = function(kb, a, all) {
      if (a != "low") return(rule(kb))
      list(rule(kb, wn = TRUE, wmsg = "bassa profondità: dropout elevato, serve smoothing più forte (meno risoluzione)"),
           rule(kb, wn = kb$signal %in% c("allele", "both"), wmsg = "bassa profondità: segnale allelico molto sparso"))
    }),

  # ------------------------------------------------------------ 4. RISORSE
  list(id = "lang", short = "Linguaggio", stage = "4. Risorse e vincoli", type = "radio",
    label = "Linguaggio/ambiente", choices = c("Indifferente" = "any", "R" = "R", "Python" = "Python"),
    apply = function(kb, a, all) {
      if (a == "any") return(rule(kb))
      rule(kb, ex = !unk(kb$language) & kb$language != a, why = paste("non disponibile in", a),
           wn = unk(kb$language), wmsg = miss("linguaggio"))
    }),

  list(id = "threads", short = "Multi-threading", stage = "4. Risorse e vincoli", type = "chk",
    label = "Mi serve il multi-threading", value = FALSE,
    apply = function(kb, a, all) {
      if (!isTRUE(a)) return(rule(kb))
      rule(kb, ex = is_(kb$multithread, "no"), why = "nessun supporto al multi-threading",
           wn = unk(kb$multithread), wmsg = miss("multi-threading"))
    })
)

# ---------------------------------------------------------------------
run_engine <- function(kb, answers) {
  n <- nrow(kb); alive <- rep(TRUE, n)
  why <- warn <- bnote <- vector("list", n); bonus <- numeric(n); trace <- list()
  for (q in QUESTIONS) {
    a <- answers[[q$id]]
    res <- q$apply(kb, a, answers)
    if (!is.null(res$ex)) res <- list(res)
    before <- alive
    for (r in res) {
      for (i in which(r$ex)) why[[i]]  <- union(why[[i]],  r$why[i])
      for (i in which(r$wn)) warn[[i]] <- union(warn[[i]], r$wmsg[i])
      for (i in which(r$bo)) if (!(r$bmsg[i] %in% bnote[[i]])) {
        bonus[i] <- bonus[i] + 1; bnote[[i]] <- c(bnote[[i]], r$bmsg[i]) }
      alive <- alive & !r$ex
    }
    trace[[length(trace) + 1]] <- data.frame(step = q$short, answer = fmt(a), remaining = sum(alive),
      removed = paste(kb$tool[before & !alive], collapse = ", "), stringsAsFactors = FALSE)
  }
  nwarn <- vapply(warn, length, integer(1))
  list(alive = alive, why = why, warn = warn, bnote = bnote, bonus = bonus,
       score = bonus - 0.5 * nwarn, trace = do.call(rbind, trace))
}

summarise_results <- function(kb, r) {
  cls <- c(expression = "expression-centric", allele = "allele-aware", both = "expression + variante allelica")
  idx <- which(r$alive); ex <- which(!r$alive)
  ok <- data.frame(Tool = kb$tool[idx], Score = r$score[idx], Classe = unname(cls[kb$signal[idx]]),
    Perche = vapply(idx, function(i) paste(r$bnote[[i]], collapse = "; "), character(1)),
    Avvertenze = vapply(idx, function(i) paste(r$warn[[i]], collapse = "; "), character(1)),
    stringsAsFactors = FALSE)
  ok <- ok[order(-ok$Score, ok$Tool), , drop = FALSE]
  excl <- data.frame(Tool = kb$tool[ex],
    Motivi = vapply(ex, function(i) paste(r$why[[i]], collapse = "; "), character(1)), stringsAsFactors = FALSE)
  list(ok = ok, excl = excl)
}
