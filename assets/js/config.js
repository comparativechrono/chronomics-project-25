/* =====================================================================
   Settings for the Bioinformatics workbench.
   Edit and re-upload this file; nothing else needs to change.
   Never put an API key in this file: anyone could copy it from a public site.
   ===================================================================== */
window.MG_CONFIG = {
  courseTitle: 'Clock E-boxes in brain genes',
  courseSubtitle: 'A guided genomics research project · Chronomic Medicine Group',

  /* The prefix of what the page keeps in the browser (files, notebooks, settings). Change it only
     if two copies of the workbench on one web site must not share their saved work. */
  storePrefix: 'ebox25',
  hostname: 'workbench',

  /* Python in the browser (Pyodide). To host it with the site instead of the CDN, see
     "Python and R from your own server" in README.md and use 'assets/vendor/pyodide/'. */
  pyodideBase: 'https://cdn.jsdelivr.net/pyodide/v0.29.5/full/',
  /* Python packages loaded in the background once Python has started, so that the first
     import is quick. Any other Pyodide package loads when a cell imports it. */
  pythonPreload: ['numpy', 'pandas', 'matplotlib'],

  /* R in the browser (webR). R packages are downloaded when a notebook asks for them
     (install.packages("ggplot2")): CRAN packages from the webR repository, Bioconductor
     packages (the current release) from r-universe. */
  webrBase: 'https://webr.r-wasm.org/v0.6.0/',
  rRepos: ['https://repo.r-wasm.org', 'https://bioc-release.r-universe.dev', 'https://cran.r-universe.dev'],
  /* R packages that the repositories above have no browser build of, built for this site and
     looked in first: locfit, which DESeq2 and edgeR need (see the README.txt in that folder). */
  rLocalRepo: 'assets/vendor/r-packages',

  /* The WebAssembly bioinformatics programs of the terminal and of Galaxy. */
  biowasmBase: 'assets/vendor/biowasm',
  // the programs of the genomics-wasm project that run in a worker of their own: featureCounts, Salmon and PLINK 2 (on
  // several threads), Cutadapt, MACS3, LoFreq, SomaticSniper, Strelka2, Kraken2, Bracken, Prodigal, HMMER and bustools (see
  // tools-wasm.js: OWN_WORKER, somatic.js and metagenomics.js); and pyodide/314.0.7/, the Python that Cutadapt, MACS3 and
  // Bracken run in (with the packages that MACS3 needs)
  threadsBase: 'assets/vendor/genomics-wasm/',

  /* AI assistant. Each user enters their own API key in the assistant's settings; it is kept
     in the browser tab's memory (or, if the user ticks "Remember", in this tab's session
     storage, which is cleared when the tab closes) – never in this file or on the server.
     The service offered first ('gemini', 'anthropic' or 'openai') and the model of each.
     Models are retired from time to time: if the assistant reports that a model is not found,
     put a current one here (see ai.google.dev/gemini-api/docs/models). */
  liveProvider: 'gemini',
  geminiModel: 'gemini-3.8-flash',
  /* If that model is busy, over a limit or not found, the assistant asks these in turn. */
  geminiFallbackModels: ['gemini-3.7-flash', 'gemini-3.6-flash', 'gemini-3.5-flash', 'gemini-3.5-flash-lite', 'gemini-3.1-flash-lite'],
  aiWaitSeconds: 60,
  anthropicModel: 'claude-sonnet-5-5',
  /* The AI agent: after how many seconds a step that has not ended is stopped. */
  agentCommandSeconds: 150,
  /* The examples the AI tab offers (the workbench's own use example files that this project does not ship).
     They help a student find their way around; none of them does a step of the project for them. */
  aiChatStarters: [
    'What is an E-box, and why does CLOCK:BMAL1 bind CACGTG?',
    'How do I read the INFO field of /data/examples/variants/promoter_variants.vcf.gz with bcftools query?',
    'Why does the enrichment test need a CpG-aware expectation rather than 1 in 4^6?',
    'Which programs can I use in the terminal here?'
  ],
  aiAgentExample: 'Show me the first lines of each file in /data/examples and tell me what each file is',
  aiAgentTasks: [
    'Show me the first lines of each file in /data/examples and tell me what each file is',
    'Which chromosomes are the 21 genes of /data/examples/genes/brain_channelopathy.tsv on, and on which strand?',
    'How many promoter variants in /data/examples/variants/promoter_variants.vcf.gz does each of the five super-populations have with an allele frequency above 1%?'
  ],

  /* How much the workbench may keep in this browser before it warns (megabytes). */
  storageWarnMB: 1500,

  /* The Genome browser. This project does not ship the workbench's example genome (data/examples/reads/), so the
     browser starts on the human genome hosted by the IGV team (needs the internet), and the guide's genome-browser
     buttons go there too. The tracks below are shown with it: the project's promoters and variants (from this site)
     and phyloP 470-way conservation (the same bigWig the analysis read, from UCSC's mirror on Amazon S3). */
  igvGenome: 'hg38',
  igvHelpText: 'Type a gene (SCN1A) or a position (chr2:166,127,520-166,129,020) in the box to go there; drag sideways to move; + / − to zoom. Click a promoter, a variant or the conservation track for details.',
  igvTracks: {
    hg38: [
      { name: 'Project promoters (TSS −1000/+500)', type: 'annotation', format: 'bed', url: 'data/examples/genes/promoters.bed', indexed: false, displayMode: 'EXPANDED', color: '#2a7f62', height: 50 },
      { name: '1000 Genomes promoter SNVs', type: 'variant', format: 'vcf', url: 'data/examples/variants/promoter_variants.vcf.gz', indexURL: 'data/examples/variants/promoter_variants.vcf.gz.tbi', displayMode: 'EXPANDED', height: 60 },
      { name: 'phyloP 470-way (mammals)', type: 'wig', format: 'bigwig', url: 'https://genome-browser.s3.us-east-1.amazonaws.com/goldenPath/hg38/phyloP470way/hg38.phyloP470way.bw', color: '#3d5a99', altColor: '#c0504d', height: 60 }
    ]
  },

  /* Help pages of the workbench that this project's guide does not have: the help buttons of the tabs, menus and
     Databases cards that point at them open these pages of the guide instead ('page' or 'page/section'; '*' for any
     other). */
  helpAliases: {
    'h-files': 'h-setup/kept',
    'h-notebook': 'h-setup/warmups',
    'h-python': 'h-tools',
    'h-r': 'h-tools',
    'h-terminal': 'h-tools',
    'h-samtools': 'h-tools',
    'h-programs': 'h-tools',
    'h-databases': 'h-tools',
    'h-db-genes': 'h-tools',
    '*': 'h-tools'
  }
};
