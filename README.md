# Clock E-boxes in brain genes — a guided research project

**An undergraduate genomics research project that runs entirely in a web browser.** A student with no computational background works through a left-hand guide, using a terminal, Python and R notebooks, gene databases and a genome browser on the right, to answer a real, novel question from the Chronomic Medicine Group:

> Do common human variants disrupt the **CLOCK:BMAL1 E-box** motifs (`CACGTG`) in the promoters of **brain channelopathy and epilepsy genes**, and do those sites show the hallmarks of **purifying selection** — conservation across mammals and a shortage of human variation?

It is the neurobiology counterpart of the group's work on E-box variants in the cardiac gene *SCN5A*. Everything the student needs is bundled, nothing is installed, and every step produces a real result.

This project is built on the **[bioinformatics workbench](https://github.com/comparativechrono/bioinformatics-workbench)**: the same in-browser terminal, notebooks, databases and genome browser, with the tabs that this project does not need removed (3D viewer, networks, Galaxy, the RStudio/Spyder IDE), and the help pane replaced by a sequential project guide.

---

## What the student does

Six parts, roughly one per week of a part-time term:

| Part | What they do | Skills |
|---|---|---|
| **Start here** | The clock, E-boxes, the question; a tour of the tools | orientation |
| **1 · Genes and promoters** | Read the PanelApp gene lists; promoter windows; pull sequence with `samtools faidx` | tables, FASTA/BED, coordinates |
| **2 · Finding E-boxes** | Scan for `CACGTG`; then a JASPAR **PWM** to find real (degenerate) sites; test enrichment against a **CpG-aware** baseline | motif scanning, Python, statistical thinking |
| **3 · Human variation** | 1000 Genomes variants; `bcftools`; `bedtools intersect`; score each variant's effect on the motif (motifbreakR-style) | VCF, allele frequencies, populations |
| **4 · Conservation and selection** | phyloP (Zoonomia 470-way); compare E-boxes with matched controls; a Mann–Whitney test; a `ggplot2` figure | conservation, R, hypothesis testing |
| **5 · Shortlist and write-up** | Rank candidate variants for the 100,000 Genomes Project; a reproducible report | synthesis, reproducibility |

The headline findings the guide leads them to (the "answer key" is built into the checkpoints on each page):

- Canonical E-boxes are **enriched ~1.7×** in these promoters over a CpG-aware expectation.
- E-box positions are **markedly more conserved** than matched control CpGs and carry **fewer human variants** — the signature of purifying selection.
- A **catalogue** of E-box-disrupting/creating variants with population frequencies, and a ranked **shortlist** for patient follow-up.

## The tabs that are kept

Notebook (Python + R), Terminal (`samtools`, `bcftools`, `bedtools`, `seqtk`, all real programs compiled to WebAssembly), Files, Databases (live Ensembl / NCBI), Genome browser (IGV), and an optional AI assistant. Removed from the base workbench: 3D viewer, Networks, Galaxy, IDE.

Some settings in `assets/js/config.js` adapt the workbench to this project (the matching code changes are marked `chronomics-project-25` in `igv-panel.js`, `workbench.js`, `page-wb.js` and `assistant.js`):

- `igvGenome` and `igvTracks` — the workbench's example genome (`data/examples/reads/`) is not shipped here, so the Genome browser starts on the human genome hg38 hosted by the IGV team (it needs the internet), with the project's promoters, variants and phyloP 470-way conservation as tracks. The guide's genome-browser buttons go there too.
- `helpAliases` — the help buttons of the tabs, the *My files* menu and the Databases cards point at workbench help pages that this guide does not have; these name the guide pages they open instead.
- `aiChatStarters`, `aiAgentTasks` and `aiAgentExample` — the example questions and tasks of the AI tab, which in the workbench use its own example files.

## Publish it on GitHub Pages

1. This repository → **Settings → Pages → Build and deployment**: Source **Deploy from a branch**, branch **main**, folder **/ (root)**. Save.
2. After a minute the site is at `https://comparativechrono.github.io/chronomics-project-25/` (or your custom domain). Share that link with the student.

`index.html` is at the top level with `assets/` and `data/` beside it, and the empty `.nojekyll` file stops GitHub from processing the site. The first visit loads once more by itself to install the service worker that makes the page cross-origin isolated (needed by a couple of the threaded programs) — see the workbench README for the details, which all still apply.

## Edit the guide

The left-hand guide is built from the pages in `build/help/` by `build/make_page.py`, exactly as in the workbench:

```bash
python3 build/make_page.py        # writes index.html
```

Each page is one `build/help/NN-name.html` file; its first line is a comment with its id, title, group, summary and keywords. The parts (groups) and their order are the `GROUPS` list in `make_page.py`. Buttons in the guide run things in the workbench: `data-run="term|py|r"` on a `<pre class="code">` adds **Run in the terminal / Add to notebook**; `data-goto`, `data-igv` (a gene or `chr:start-end` on hg38), `data-igv-track` (a file of the student's, e.g. `~/eboxes_canonical.bed`, as a genome-browser track), `data-api` (e.g. `ensembl:gene:SCN1A`), `data-bench` and `data-open` drive the other tabs. References use `<!--CITE: key …-->` (a References section), `<!--CITE-LIST: key …-->` (the list alone, as in the reading list) and `<!--CITE-SHORT: key-->` (author and year) against `build/citations.json`. `make_page.py` warns about links to pages that do not exist, links with no address, and `/data/examples` files that are not in the bundle.

## The data, and how it was made

Everything the student uses is in `data/examples/`, prepared from public sources so nothing large is downloaded in the browser (`data/examples` is described in the guide's **Reference → The project data**). It was generated by the analysis scripts kept in [`analysis/`](analysis/) (see `analysis/README.md`): gene lists from **Genomics England PanelApp**; canonical TSS and promoter sequence from **Ensembl/GENCODE** (GRCh38, TSS −1000/+500); variants from the **1000 Genomes Project 30× release** (3,202 genomes); conservation from **phyloP** (Zoonomia 470-way and 100-way); and the **BMAL1 (Arntl) matrix MA0603.2** from **JASPAR**. The same scripts reproduce the numbers quoted in the guide's checkpoints.

## Licences and credits

The workbench code and the guide are under the licences of the [bioinformatics workbench](https://github.com/comparativechrono/bioinformatics-workbench) (MIT for code; see `THIRD_PARTY.md` for the bundled programs and their licences). The bundled data are public: 1000 Genomes Project, Ensembl/GENCODE, UCSC/Zoonomia phyloP, JASPAR and Genomics England PanelApp — each under its own terms, cited in the guide's reading list. Built by the Chronomic Medicine Group, University of Cambridge.
