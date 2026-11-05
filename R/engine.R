# =====================================================================
# Rule engine for choosing CNA tools from scRNA-seq
# - Does NOT use benchmark results to exclude/rank (known gaps)
# - Exclusions = tool "design" constraints (required inputs, outputs, species...)
# - Warnings = mechanistic caveats, not performance scores
# - Bonus = the tool exploits information the user has available
# To add a question: add an element to QUESTIONS.
# =====================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
is_ <- function(x, v) !is.na(x) & x %in% v           # NA-safe
unk <- function(x) is.na(x) | x == "unknown"          # missing data in the KB
fmt <- function(a) if (is.null(a) || !length(a) || all(is.na(a))) "-" else paste(a, collapse = ", ")

has_out <- function(kb, vals) vapply(kb$output, function(o)
  !is.na(o) && any(trimws(strsplit(o, ";")[[1]]) %in% vals), logical(1), USE.NAMES = FALSE)

# rule: ex = exclusion, wn = warning, bo = bonus (vectors of length n tools)
rule <- function(kb, ex = FALSE, why = "", wn = FALSE, wmsg = "", bo = FALSE, bmsg = "") {
  n <- nrow(kb)
  list(ex = rep_len(ex, n), why = rep_len(why, n), wn = rep_len(wn, n),
       wmsg = rep_len(wmsg, n), bo = rep_len(bo, n), bmsg = rep_len(bmsg, n))
}
miss <- function(lbl) paste0(lbl, ": missing data in the KB, please verify")

QUESTIONS <- list(

  # ------------------------------------------------------------ 1. GOAL
  list(id = "goals", short = "Goals", stage = "1. Analysis goal", type = "multi",
    label = "What do you need as output? (select everything you need)",
    choices = c("Classification of tumor vs normal cells" = "classification",
                "Subclones defined by CNAs" = "subclones",
                "Phylogenetic tree / evolutionary history" = "phylogeny",
                "Allele-specific states (cnLOH, LOH)" = "cnloh",
                "Distinguish balanced/unbalanced amplifications" = "balanced"),
    help = "Each item excludes the tools that do not implement it.",
    apply = function(kb, a, all) {
      a <- a %||% character(0); out <- list()
      if ("classification" %in% a)
        out[[length(out) + 1]] <- rule(kb,
          ex = is_(kb$classification, "none"), why = "does not produce tumor/normal classification",
          wn = is_(kb$classification, "threshold") | unk(kb$classification),
          wmsg = ifelse(is_(kb$classification, "threshold"),
                        "provides only a per-cell score: a threshold must be chosen",
                        miss("classification")))
      feat <- c(subclones = "subclones", phylogeny = "phylogeny",
                cnloh = "allele-specific states (cnLOH/LOH)", balanced = "balanced/unbalanced amplifications")
      col  <- c(subclones = "subclones", phylogeny = "phylogeny", cnloh = "cnloh", balanced = "balanced_amp")
      for (k in intersect(names(feat), a))
        out[[length(out) + 1]] <- rule(kb, ex = is_(kb[[col[[k]]]], "no"),
          why = paste("does not implement:", feat[[k]]),
          wn = unk(kb[[col[[k]]]]), wmsg = miss(feat[[k]]))
      if (length(out)) out else rule(kb)
    }),

  list(id = "output", short = "Output type", stage = "1. Analysis goal", type = "radio",
    label = "Required output type",
    choices = c("No preference" = "any",
                "Discrete states (even just gain/neutral/loss)" = "discrete",
                "Multi-level discrete states (deletion/loss/neutral/gain/amp)" = "multi"),
    help = "Continuous scores require a user-chosen threshold.",
    apply = function(kb, a, all) {
      if (a == "any") return(rule(kb))
      vals <- if (a == "discrete") c("3-state", "multi-state", "allele-specific") else c("multi-state", "allele-specific")
      rule(kb, ex = !unk(kb$output) & !has_out(kb, vals), why = "continuous-only output / insufficient levels",
           wn = unk(kb$output), wmsg = miss("output type"))
    }),

  list(id = "level", short = "Call level", stage = "1. Analysis goal", type = "radio",
    label = "Level at which you want the CNA call",
    choices = c("Doesn't matter" = "any", "Per single cell" = "cell", "Per clone/subclone" = "clone"),
    apply = function(kb, a, all) {
      if (a == "any") return(rule(kb))
      bad <- if (a == "cell") "clone" else "cell"
      rule(kb, ex = is_(kb$call_level, bad), why = paste("the CNA call is only at", bad, "level"),
           wn = unk(kb$call_level), wmsg = miss("call level"))
    }),

  list(id = "res", short = "Resolution", stage = "1. Analysis goal", type = "radio",
    label = "Required genomic resolution",
    choices = c("Chromosome arm is sufficient" = "arm", "Segments / sub-arm" = "segment",
                "Focal events (a few Mb)" = "focal"),
    help = "All tools apply smoothing: more smoothing = more robustness but less resolution.",
    apply = function(kb, a, all) {
      if (a == "arm") return(rule(kb))
      out <- list(rule(kb, ex = is_(kb$resolution, "arm"), why = "calls CNAs only per chromosome arm"))
      if (a == "focal") out[[2]] <- rule(kb, wn = TRUE,
        wmsg = "focal events: smoothing over gene windows reduces their detectability")
      out
    }),

  list(id = "k", short = "Known no. of subclones?", stage = "1. Analysis goal", type = "radio",
    label = "(If you want subclones) do you know the number of subclones a priori?",
    choices = c("No / I don't know" = "unknown", "Yes, I can specify it" = "known"),
    apply = function(kb, a, all) {
      if (a == "known" || !("subclones" %in% (all$goals %||% character(0)))) return(rule(kb))
      rule(kb, ex = is_(kb$k_mode, "apriori"), why = "requires the number of subclones a priori",
           wn = unk(kb$k_mode), wmsg = miss("handling of the number of subclones"))
    }),

  # ------------------------------------------------------------ 2. DATA
  list(id = "format", short = "Data format", stage = "2. Available data", type = "radio",
    label = "What data do you have?",
    choices = c("BAM/FASTQ available (allelic counts possible)" = "bam",
                "Gene x cell matrix only" = "matrix"),
    help = "Allele-aware tools extract REF/ALT per SNP and cell from the BAM.",
    apply = function(kb, a, all) {
      if (a == "bam") return(rule(kb))
      rule(kb, ex = is_(kb$allele_input, "required"), why = "requires BAM for allelic counts")
    }),

  list(id = "phasing", short = "Phasing panel", stage = "2. Available data", type = "radio",
    label = "Do you have a (population) reference panel for phasing suited to the species/ancestry?",
    choices = c("Yes" = "yes", "No" = "no", "Don't know" = "unknown"),
    help = "Some tools import haplotypes from population panels instead of phasing from scRNA-seq data.",
    apply = function(kb, a, all) {
      if (a == "yes") return(rule(kb))
      if (a == "no") return(rule(kb, ex = is_(kb$phasing_panel, "yes"), why = "requires an external phasing panel"))
      rule(kb, wn = is_(kb$phasing_panel, "yes"), wmsg = "requires an external phasing panel: verify its availability")
    }),

  list(id = "bulk", short = "Matched bulk DNA", stage = "2. Available data", type = "radio",
    label = "Do you have matched bulk DNA (WGS/WES)?",
    choices = c("No" = "no", "Yes" = "yes"),
    apply = function(kb, a, all) {
      if (a == "yes") return(rule(kb, bo = is_(kb$bulk_dna, c("optional", "recommended")),
                                  bmsg = "can integrate breakpoints/profiles from bulk DNA"))
      rule(kb, wn = is_(kb$bulk_dna, "recommended"), wmsg = "optimal performance with matched WGS, which you do not have")
    }),

  list(id = "ref", short = "Normal reference", stage = "2. Available data", type = "radio",
    label = "A priori copy-neutral reference",
    choices = c("No, no a priori annotation" = "none",
                "Yes, annotated normal cells in the same sample" = "insample",
                "Yes, external reference (other dataset/atlas)" = "external"),
    help = "Without annotations, only tools that identify the reference automatically remain.",
    apply = function(kb, a, all) {
      if (a != "none") return(rule(kb))
      rule(kb, ex = is_(kb$ref_mode, "required"), why = "requires a user-provided normal reference",
           wn = unk(kb$ref_mode), wmsg = miss("reference handling"))
    }),

  list(id = "normal", short = "Normal cells in sample", stage = "2. Available data", type = "radio",
    label = "Does the sample contain non-tumor cells (immune/stromal)?",
    choices = c("Yes / likely" = "yes", "No (e.g. pure cell line)" = "no", "Don't know" = "unknown"),
    help = "Only matters if you have no external reference.",
    apply = function(kb, a, all) {
      if (a == "yes" || (all$ref %||% "none") != "none") return(rule(kb))
      rule(kb, ex = is_(kb$runs_without_normal, "no") & a == "no",
           why = "without normal cells or an external reference it cannot estimate the baseline",
           wn = unk(kb$runs_without_normal),
           wmsg = miss("operation without normal cells"))
    }),

  # ------------------------------------------------------------ 3. SAMPLE
  list(id = "species", short = "Species", stage = "3. Sample and sequencing", type = "radio",
    label = "Species", choices = c("Human" = "human", "Mouse" = "mouse", "Other" = "other"),
    help = "Many tools use predefined gene annotations (human/mouse, few genome versions).",
    apply = function(kb, a, all) {
      v <- kb[[paste0("sp_", a)]]
      rule(kb, ex = is_(v, "no"), why = paste("unsupported species:", a),
           wn = unk(v) | is_(v, "limited"),
           wmsg = ifelse(is_(v, "limited"), "limited support for this species (heterozygous SNPs + phasing)", miss("species support")))
    }),

  list(id = "het", short = "Heterozygosity", stage = "3. Sample and sequencing", type = "radio",
    label = "Does the organism/strain have sufficient germline heterozygous variability?",
    choices = c("Yes (e.g. human, outbred mouse)" = "yes", "No (e.g. inbred strain)" = "no"),
    apply = function(kb, a, all) {
      if (a == "yes") return(rule(kb))
      rule(kb, ex = is_(kb$allele_input, "required"), why = "informative heterozygous SNPs are needed")
    }),

  list(id = "sample", short = "Sample type", stage = "3. Sample and sequencing", type = "select",
    label = "Sample type",
    choices = c("Primary solid tumor" = "primary", "Metastasis" = "met", "Hematologic" = "heme",
                "Cell line" = "line", "Organoid" = "organoid", "PDX / xenograft" = "pdx"),
    apply = function(kb, a, all) {
      if (a == "heme")
        return(rule(kb, wn = is_(kb$ref_mode, "required") & (all$ref %||% "none") != "external",
          wmsg = "hematologic: the typical reference (immune cells) is not copy-neutral; an external reference is better"))
      if (a == "pdx")
        return(rule(kb, wn = TRUE, wmsg = "PDX: separate human reads from murine reads before the analysis"))
      rule(kb)
    }),

  list(id = "purity", short = "Tumor purity", stage = "3. Sample and sequencing", type = "radio",
    label = "Fraction of tumor cells",
    choices = c("Balanced / don't know" = "mid", "High (dominant tumor)" = "high", "Low (minority tumor)" = "low"),
    apply = function(kb, a, all) {
      if (a == "high" && (all$ref %||% "none") == "none")
        return(rule(kb, wn = is_(kb$ref_mode, c("auto", "optional")) & kb$signal == "expression",
          wmsg = "few normal cells: the automatic expression-based reference may be poorly centered"))
      if (a == "low") return(rule(kb, wn = TRUE, wmsg = "low tumor fraction: diluted CNA signal, rare subclones hard to detect"))
      rule(kb)
    }),

  list(id = "ploidy", short = "Ploidy / copy-neutral events", stage = "3. Sample and sequencing", type = "chk",
    label = "Expected WGD / extensive aneuploidy / copy-neutral events",
    value = FALSE,
    apply = function(kb, a, all) {
      if (!isTRUE(a)) return(rule(kb))
      rule(kb, wn = kb$signal == "expression",
           wmsg = "expression only: same total CN with different allelic composition is indistinguishable",
           bo = kb$signal %in% c("allele", "both"), bmsg = "the allelic signal resolves configurations with unchanged total CN")
    }),

  list(id = "batch", short = "Batch / multi-platform", stage = "3. Sample and sequencing", type = "chk",
    label = "Multiple batches or platforms to analyze together", value = FALSE,
    apply = function(kb, a, all) {
      if (!isTRUE(a)) return(rule(kb))
      rule(kb, wn = TRUE, wmsg = "batch: correct on the counts (e.g. ComBat); scRNA integration methods return embeddings, not counts")
    }),

  list(id = "platform", short = "Platform", stage = "3. Sample and sequencing", type = "select",
    label = "Sequencing platform",
    choices = c("Droplet 3' (e.g. 10x 3')" = "droplet3", "Droplet 5' / other tag-based" = "droplet5",
                "Plate-based full-length (e.g. Smart-seq)" = "plate", "Other / mixed" = "other"),
    apply = function(kb, a, all) {
      if (a == "droplet3")
        return(rule(kb, wn = kb$signal %in% c("allele", "both"),
          wmsg = "3' protocol: reads concentrated toward the 3' end, fewer SNPs covered"))
      rule(kb)
    }),

  list(id = "cells", short = "No. of cells", stage = "3. Sample and sequencing", type = "num",
    label = "Number of cells", value = 3000,
    apply = function(kb, a, all) {
      if (is.null(a) || is.na(a)) return(rule(kb))
      rule(kb, wn = kb$signal %in% c("allele", "both") & a < 500,
           wmsg = "few cells: the allelic signal is aggregated across cells/SNPs, few informative SNPs")
    }),

  list(id = "depth", short = "Depth", stage = "3. Sample and sequencing", type = "select",
    label = "Depth per cell",
    choices = c("Medium / don't know" = "mid", "Low" = "low", "High" = "high"),
    apply = function(kb, a, all) {
      if (a != "low") return(rule(kb))
      list(rule(kb, wn = TRUE, wmsg = "low depth: high dropout, stronger smoothing needed (less resolution)"),
           rule(kb, wn = kb$signal %in% c("allele", "both"), wmsg = "low depth: very sparse allelic signal"))
    }),

  # ------------------------------------------------------------ 4. RESOURCES
  list(id = "lang", short = "Language", stage = "4. Resources and constraints", type = "radio",
    label = "Language/environment", choices = c("Doesn't matter" = "any", "R" = "R", "Python" = "Python"),
    apply = function(kb, a, all) {
      if (a == "any") return(rule(kb))
      rule(kb, ex = !unk(kb$language) & kb$language != a, why = paste("not available in", a),
           wn = unk(kb$language), wmsg = miss("language"))
    }),

  list(id = "threads", short = "Multi-threading", stage = "4. Resources and constraints", type = "chk",
    label = "I need multi-threading", value = FALSE,
    apply = function(kb, a, all) {
      if (!isTRUE(a)) return(rule(kb))
      rule(kb, ex = is_(kb$multithread, "no"), why = "no multi-threading support",
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
  cls <- c(expression = "expression-centric", allele = "allele-aware", both = "expression + allelic variant")
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
