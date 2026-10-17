# ProfileAlignments.jl

Pacchetto indipendente per allineamenti di sequenze e profili, scritto in Julia.
Julia ≥1.10; usa BioSequences.jl per l'interoperabilità dei tipi biologici,
senza libreria FFT esterna, eseguibile MAFFT, R o servizio remoto. AASV lo usa
come motore di allineamento.

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/ProfileAlignments.jl")
using ProfileAlignments, BioSequences

records = [Sequence("reference", aa"ACDEFGHIKLMNPQ"),
           Sequence("fragment", aa"DEFGHIK"),
           Sequence("insertion", aa"ACDEFGWHIKLMNPQ")]
result = align_sequences(records)
write_fasta(result, "aligned.fasta")
```

Da questa cartella, il pacchetto si può usare e testare autonomamente:

```powershell
julia --project=. examples/protein_alignment.jl
julia --project=. -e 'using Pkg; Pkg.test()'
julia -t auto --project=. benchmark/threading.jl
```

## Algoritmo

Il disegno segue la famiglia di strategie progressive e iterative descritta
negli [algoritmi MAFFT](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html).
È una nuova implementazione, con scelte e parametri propri.

1. Collassa le sequenze identiche, conservando i nomi originali.
2. Stima le distanze con conteggi k-mer nell'alfabeto proteico ridotto a sei
   gruppi Dayhoff. Usa la correzione per lunghezza documentata da MAFFT.
3. Costruisce un albero guida UPGMA deterministico. È disponibile anche
   `linkage=:mixed`, con miscela di distanza minima e media (`mix=0.1`).
4. Allinea progressivamente i profili seguendo l'albero. I punteggi sono
   medie pesate BLOSUM62 e i costi dei gap sono affini.
5. Ricostruisce opzionalmente l'albero dalle distanze nell'MSA e ripete
   l'allineamento. Mantiene il nuovo MSA se migliora il punteggio.
6. Raffina il risultato separando i gruppi sui rami dell'albero e riallineandoli.
   Accetta soltanto miglioramenti del punteggio esatto sum-of-pairs pesato.
7. Ripristina duplicati, identificativi e ordine iniziale. I duplicati dividono
   il peso del loro rappresentante, evitando di sovrappesare copie identiche.

`sequence_weighting=:tree` usa pesi derivati dalle lunghezze dei rami divise
tra le foglie discendenti; `:uniform` assegna lo stesso peso ai rappresentanti.
I pesi rimangono fissi durante il confronto dei candidati e il raffinamento.

## API

| Funzione / tipo | Uso |
|---|---|
| `Sequence(id, sequence)` | sequenza con identificativo |
| `Alignment(records; weights=...)` | profilo rettangolare con pesi positivi normalizzati |
| `ScoringModel(...)` | matrice, penalità e trattamento delle estremità |
| `pairalign(a, b; model=...)` | allineamento affine di due sequenze, risultato `.alignment` e `.score` |
| `merge_profiles(a, b; model=...)` | fusione di due profili con conservazione delle colonne |
| `align_sequences(records; ...)` | MSA progressivo e raffinato |
| `align_sequence_sets(groups; ...)` | allineamento parallelo di insiemi indipendenti |
| `refine_alignment(alignment; ...)` | raffinamento di un MSA esistente |
| `alignment_score(alignment; model=...)` | punteggio affine sum-of-pairs pesato |
| `kmer_distances(records; ...)` | distanze iniziali approssimate |
| `guide_tree(distances; ...)` | albero guida |
| `fft_lags(a, b; ...)`, `fft_anchors(a, b; ...)` | diagonali candidate e catena di ancore |
| `read_fasta`, `write_fasta`, `align_fasta` | input/output e flusso FASTA completo |

`Sequence`, `pairalign` e `align_sequences` accettano sia stringhe sia
`BioSequence`, compresi `LongAA`, `LongDNA{2}` e `LongDNA{4}`. Al confine del
motore i simboli vengono convertiti una sola volta nella rappresentazione ASCII
usata dai cicli DP ottimizzati; i risultati conservano quindi l'API storica a
stringhe. L'alfabeto BioSequences riconosce anche `O`, `U` e `J`, ma il modello
BLOSUM62 predefinito continua a rifiutarli perché non possiede punteggi definiti.

```julia
# Sequenza corta inserita in un riferimento: solo i gap terminali della query
# sono gratuiti. La sequenza query è sempre il primo argomento.
model = ScoringModel(mode=:query)
pair = pairalign("DEFGHIK", "ACDEFGHIKLMNPQ"; model)

# DNA con match/mismatch scalari, senza reverse complement automatico.
dna = align_sequences(["ACGTACGT", "CGTAC", "ACGTTACGT"];
                      model=ScoringModel(alphabet=:dna))

# Traccia dei punteggi dei soli candidati accettati.
trace = align_sequences(records; retree=2, max_iterations=5,
                        return_trace=true)
trace.scores
trace.alignment
```

### Strategie di calcolo

`align_sequences` accetta preset espliciti. Il default resta `:accurate`,
quindi l'aggiornamento non cambia automaticamente l'algoritmo:

```julia
accurate = align_sequences(records; strategy=:accurate)
balanced = align_sequences(records; strategy=:balanced)
large    = align_sequences(records; strategy=:large)
automatic = align_sequences(records; strategy=:auto)
```

- `:accurate`: due passaggi dell'albero, due cicli di raffinamento, DP completa;
- `:balanced`: un albero, un ciclo di raffinamento, ancore FFT;
- `:large`: un albero, nessun raffinamento globale, ancore FFT.
- `:auto`: sceglie i tre preset a 100 e 500 sequenze distinte; le soglie si
  modificano con `strategy_thresholds=(accurate_max, balanced_max)`.

`retree`, `max_iterations` e `anchors`, se specificati, prevalgono sul preset.
Le strategie FFT sono euristiche e possono produrre un allineamento differente.
I tipi `AccurateStrategy`, `BalancedStrategy` e `LargeStrategy` selezionano lo
stesso comportamento tramite dispatch.

Per molti insiemi indipendenti, `align_sequence_sets` distribuisce gli insiemi
tra task Julia preservando l'ordine dei risultati. Con più insiemi assegna un
thread a ciascun allineamento, evitando parallelismo annidato; con un solo
insieme usa l'intero budget di thread.

Il calcolo usa i thread Julia disponibili per default (`threads=:auto`). Avvia
Julia con `julia -t auto` oppure imposta `JULIA_NUM_THREADS` per abilitarne più
di uno. `threads=1` forza l'esecuzione seriale; un intero positivo limita il
numero di thread usati dal pacchetto. Le distanze k-mer, le distanze ricavate
dall'MSA e i punteggi pairwise vengono distribuiti dinamicamente. I due rami
indipendenti dell'albero guida vengono allineati in parallelo con un budget di
task limitato. L'ordine delle riduzioni e delle decisioni di raffinamento resta
fisso, quindi lo stesso input produce lo stesso allineamento e gli stessi
punteggi con `threads=1` e `threads=:auto`.

Il guadagno dipende dalla forma dell'albero, dalla lunghezza e dal numero delle
sequenze. Gli ultimi merge del progressivo hanno dipendenze inevitabilmente
seriali; per questo non ci si deve attendere uno speedup proporzionale ai core.
`benchmark/threading.jl [numero_sequenze] [lunghezza]` confronta le due modalità
sulla macchina corrente senza dipendenze di benchmark aggiuntive.
`benchmark/scaling.jl [dimensioni...]` confronta i tre preset e riporta tempo,
allocazioni e larghezza finale.

## Rappresentazione e memoria

L'API pubblica conserva `Sequence` e FASTA come stringhe ASCII. I kernel interni
leggono direttamente i `codeunits` `UInt8`, senza copiare i residui in matrici
`Char`. Output e colonne filtrate sono costruiti in buffer di byte di dimensione
esatta.

I k-mer comuni sono codificati senza collisioni in `UInt64` (base 7 per i sei
gruppi Dayhoff, base 5 per DNA); valori di `k` troppo grandi usano il percorso
generico a tuple. La costruzione UPGMA e la cache dei punteggi usano matrici
simmetriche compattate con diagonale implicita. Le modalità terminali del DP
sono parametri di tipo di `ScoringModel`, così Julia compila kernel distinti
per allineamento globale, overlap e query. La FFT radix-2 resta Julia puro e
riusa i buffer tra i canali.

`trace.tree` indicizza i rappresentanti unici elencati in `trace.representatives`;
il risultato conserva tutte le sequenze originali. `trace.scores` riporta
i punteggi del risultato completo, incluso il termine costante dei duplicati.

Per `ScoringModel`, le modalità sono `:global` (penalizza tutte le estremità),
`:overlap` (estremità gratuite da entrambi i lati, default) e `:query`.
Quest'ultima vale per `pairalign` e `merge_profiles`, non per MSA/raffinamento.
Un gap di lunghezza L costa **gap_open + (L−1) × gap_extend**. Default proteici:
BLOSUM62, `gap_open=-10`, `gap_extend=-1`, `terminal_gap=0`. DNA:
`match=2`, `mismatch=-3`; `N` ha punteggio zero. Le altre ambiguità DNA
seguono i punteggi scalari; non vengono espanse in probabilità IUPAC.
È possibile passare `matrix` e `symbols` personalizzati; la matrice deve essere
finita e simmetrica. Proteine: 20 residui standard, `B`, `Z`, `X` e `*`.

## FFT opzionale

```julia
fast = align_sequences(records; anchors=:fft, anchor_window=12)
```

Una FFT radix-2 Julia calcola le correlazioni delle proprietà volume/polarità
normalizzate per le proteine, secondo il principio del
[lavoro MAFFT del 2002](https://pmc.ncbi.nlm.nih.gov/articles/PMC135756/).
Per DNA si usano quattro canali indicatori. Tra i picchi di correlazione vengono
cercate finestre ad alta identità; una catena monotona dei loro centri divide
la matrice in blocchi più piccoli. In assenza di ancore viene usato il DP completo.

Questa accelerazione **forza le corrispondenze delle ancore** e può perdere
allineamenti migliori. È disabilitata per default. Non riproduce il metodo
di selezione dei segmenti o il sistema di punteggio esatto di MAFFT.
I parametri non sono intercambiabili con `--op`, `--ep` e `--maxiterate`.

## Correttezza e limiti

- L'allineamento affine tra due sequenze è verificato contro l'enumerazione
  indipendente di tutti gli allineamenti di piccoli input, in tutte e tre
  le modalità terminali. La suite verifica anche FFT contro trasformate e
  correlazioni dirette, colonne dei profili, duplicati, determinismo e
  monotonicità del raffinamento.
- Il DP tra profili usa il prodotto delle frequenze pesate per le sostituzioni
  e l'occupazione delle colonne per i gap. È un'approssimazione del problema
  affine sum-of-pairs: non garantisce il migliore allineamento globale di profili.
  Il raffinamento usa invece il punteggio affine effettivo delle coppie, ignorando
  le colonne in cui entrambe le sequenze hanno un gap.
- Con DP completo, tempo O(L₁×L₂) e memoria di traceback O(L₁×L₂);
  i punteggi intermedi usano vettori a memoria O(L₁), oltre ai profili.
  `max_cells=100_000_000` limita le celle del confronto, o del blocco più grande
  quando si usano ancore. Albero e distanze richiedono memoria O(N²);
  la coda delle distanze rende la costruzione UPGMA O(N² log N). Il costo del
  raffinamento dipende dalle bipartizioni e dalle coppie che le attraversano;
  `max_iterations=0` lo disabilita.
- Non implementa PartTree, punteggi di consistenza L/G/E-INS-i, SIMD o
  allineamento strutturale dell'RNA. Non è un clone completo
  di MAFFT e non ne è dimostrata l'equivalenza biologica o prestazionale.
- I controlli attuali sono algoritmici e su dati sintetici. Prima dell'uso
  scientifico definitivo sono necessari benchmark di allineamenti reali
  e un confronto sull'effetto delle scelte di allineamento in AASV.

Licenza GPL ≥2, testo in `COPYING`. Implementazione indipendente basata
su descrizioni pubblicate; nessun sorgente MAFFT copiato. La matrice BLOSUM62
è il dato di sostituzione aminoacidica della versione Biostrings del confronto
R, in ordine ARNDCQEGHILKMFPSTWYVBZX*. Per `B`, `Z` e `X` si usa quindi la
convenzione Biostrings; `X` contro un residuo ha punteggio −1. Vedere
`THIRD_PARTY_NOTICES.md` per provenienza e riferimenti.
