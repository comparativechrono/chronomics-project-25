# Analysis — how the project data and answer key were made

These scripts generated everything in `data/examples/` and the "checkpoint" numbers quoted in the guide. They are the reproducible record of the supervisor-side analysis: the student reproduces the same results from the bundle, step by step, through the guide. They are **not** needed to run or host the project.

They were run in a Linux environment with Python 3.13 and `pysam`, `pybigtools`, `numpy`, `pandas`, `scipy`, `matplotlib`. Intermediate files were written to a working directory (`raw/`, `out/`, `fig/`); the paths at the top of each script point there. Network access to Ensembl, PanelApp, the 1000 Genomes FTP, the UCSC phyloP mirror and JASPAR is required.

## Order

| Script | Does | Key outputs |
|---|---|---|
| `01_genes_and_promoters.py` | PanelApp green genes → canonical TSS from Ensembl → promoter windows (TSS −1000/+500, strand-aware, GRCh38) | `raw/promoters.json` |
| `02_sequences_and_eboxes.py` | Fetch each promoter's forward-strand sequence from Ensembl; locate canonical `CACGTG` | `raw/promoter_seqs.json`, `out/eboxes.json` |
| `03_ebox_tables.py` | Canonical and PWM (JASPAR ARNTL MA0603.2, rel-score ≥ 0.85) E-box tables + BEDs | `out/eboxes_canonical.*`, `out/eboxes_pwm.*`, `out/pwm_model.json` |
| `03b_enrichment_counts.py` | Observed vs CpG-aware (Markov-1) expected E-box counts | `out/motif_counts.json` |
| `04_variants.py` (+ `04_variants_driver.sh`) | Pull every promoter SNV from the 1000 Genomes 30× phased panel by remote tabix, with per-super-population allele frequencies. Resumable; the driver reruns it until all genes are done (remote reads are flaky). | `out/promoter_snvs.csv` |
| `05_phylop.py` | phyloP (Zoonomia 470-way and 100-way) for every promoter base, from the UCSC S3 mirror via `pybigtools` | `raw/phylop.npz` |
| `06_overlap_selection_shortlist.py` | Variant × E-box overlap; per-variant PWM score change (motifbreakR-style); the selection test (E-box vs CpG-matched control); the shortlist | `out/catalogue.csv`, `out/selection_stats.json`, `out/shortlist.csv`, `out/summary.json` |
| `07_figures.py` | The figures (enrichment, position, conservation, variation) | `fig/*.png` |
| `08_build_bundle.py` | Writes the student-facing bundle into `data/examples/` and the starter notebooks | the shipped data |
| `09_verify_bundle.py` | Recomputes every guide checkpoint from the **shipped bundle** using the guide's own method — the consistency check | prints the numbers |

Scope note: genes on non-nuclear contigs are dropped (one mitochondrial gene, `MT-TK`), leaving **641** nuclear promoters.

## The answer key (what the student should get)

- **641** promoters, 1500 bp each, GC ≈ 59%.
- Canonical `CACGTG` E-boxes: **205** (0 in the 21 channelopathy genes); enriched **1.66×** over a CpG-aware expectation (205 vs ~124).
- PWM sites (ARNTL, rel-score ≥ 0.85, left-to-right non-overlapping): **≈ 1240**, of which **≈ 28** in **13** of the 21 channelopathy genes.
- Promoter SNVs (1000 Genomes 30×): **22,836**; **27** fall inside a canonical E-box; PWM-break **212**, PWM-create **210**.
- **Selection:** mean phyloP470 at E-box positions **2.68** vs **−0.18** at control CpGs (vs 0.81 over all positions); Mann–Whitney p ≈ 3 × 10⁻¹³⁴. Variant density **22.0** vs **42.9** per 1000 bp at control CpGs (odds ratio ≈ 0.50). E-boxes are conserved and depleted of variation — the signature of purifying selection.
- **Shortlist:** disrupting variants at highly conserved positions in epilepsy/channelopathy genes (e.g. PRRT2, TRPM7, CHD4, CUX2, PURA).

The control is CpG-matched because `CACGTG` contains a CpG; the guide tags the C of each CpG as the control position, and `09_verify_bundle.py` uses that same definition so the shipped numbers and the guide agree.

## Data sources

1000 Genomes Project 30× high-coverage phased panel (3,202 genomes); Ensembl/GENCODE GRCh38 gene models and sequence; Genomics England PanelApp (panels "Brain channelopathy" and "Early-onset or syndromic epilepsy", green genes); UCSC/Zoonomia phyloP (470-way, 100-way) on GRCh38; JASPAR 2022 matrix MA0603.2 (Arntl). All public; cited in the guide's reading list.
