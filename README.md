# CNA tool selector (Shiny)

Selezione guidata dei tool per l'inferenza di CNA da scRNA-seq, basata su un albero decisionale
su (i) scopo dell'analisi, (ii) dati/input disponibili, (iii) tipo di campione e sequenziamento, (iv) risorse.
I benchmark NON entrano in esclusioni né ranking (restano consultabili in una tab separata).

## Avvio
```r
install.packages(c("shiny", "DT"))
shiny::runApp("cna_tool_selector")
```

## Struttura
- `app.R` – interfaccia e server (le domande sono generate da `QUESTIONS`).
- `R/engine.R` – domande + regole (esclusione / avvertenza / bonus) e funzione `run_engine()`.
- `data/tools_kb.csv` – una riga per tool; `unknown` = da verificare (nessuna esclusione, solo avvertenza).
- `data/benchmark_evidence.csv` – evidenze dai benchmark, solo informative.

## Percorso decisionale (ordine delle domande)
1. Obiettivo: classificazione, sottocloni, filogenesi, cnLOH, amplificazioni bilanciate; tipo di output; livello cellula/clone; risoluzione; N sottocloni noto.
2. Dati: BAM vs matrice; pannello di phasing; DNA bulk; riferimento normale a priori; presenza di cellule normali.
3. Campione: specie; eterozigosità; tipo di campione; purezza; ploidia/eventi copy-neutral; batch; piattaforma; N cellule; profondità.
4. Risorse: linguaggio; multi-threading.

## Estendere
- Nuovo tool: aggiungere una riga a `tools_kb.csv`.
- Nuova condizione: aggiungere una colonna alla KB e un elemento a `QUESTIONS` con `apply(kb, a, all)` che restituisce `rule(...)`.
- Soglie euristiche (es. 500 cellule per gli allele-aware) sono in `engine.R`: da validare.
