# gene-oscillations

Analysis pipelines for detecting oscillatory / cyclic gene expression patterns in single-cell
RNA-seq data, using [Cyclum](https://github.com/KChen-lab/Cyclum) and
[scPrisma](https://github.com/nitzanlab/scPrisma).

## Repository layout

- `segmentation_clock/` — preprocessing and analysis scripts for the mouse E9.5 presomitic
  mesoderm segmentation clock dataset (GSE114186 / SPRING `Diaz2019/E95PM`).
- `myc/` — preprocessing and analysis scripts for the HCT116 cell line MYC
  oscillation dataset (GSM4286760).
- `nfkb/` — preprocessing and analysis scripts for the NF-kB oscillation dataset (GSE162992,
  WT / IkBaMM mouse BMDMs under TNF/LPS/PIC stimulation).
- `synthetic/` — generation (via [dyngen](https://github.com/dynverse/dyngen)) and evaluation
  scripts for the synthetic benchmark datasets, plus grid-search/tuning variants.
- `algorithms/` — one clearly-defined entry point per algorithm, each dataset's
  `run_cyclum.R`/`run_scPrisma.R`/`run_oscope.R` calls into these rather than
  re-inlining the invocation logic itself:
  - `run_cyclum.py`, `run_cyclum.R` — the actual Python Cyclum CLI, and an R
    `run_cyclum()` wrapper that shells out to it via `reticulate` (requires
    `use_condaenv("cyclum_env")` first).
  - `run_scPrisma.py`, `run_scPrisma.R` — same pattern for scPrisma
    (`use_condaenv("scPrisma_env")` first).
  - `run_oscope.R` — Oscope is an R/Bioconductor package, so `run_oscope()` runs
    in-process rather than shelling out to a separate script.
- `common/` — shared R helpers sourced by all of the above:
  - `hdfrw.R` — `.h5` read/write (`mat2hdf()`/`hdf2mat()`).
  - `utils.R` — cell QC (`qc_cell_mask()`), model-score loading (`get_cyclum_scores()`,
    `get_model_scores()`).
  - `gsea_functions.R` — Hallmark/GO/DoRothEA/TFEA.ChIP pathway analysis
    (`run_pathway_analyses()` and friends), used by every `analyse.R`.
  - `plotting_utils.R` — shared `theme_common()` plus score-comparison scatter plots
    (`plot_score_comparison()`, `plot_score_vs_expression()`).
  - `grid_search.R` — algorithm-agnostic hyperparameter-sweep drivers: `run_grid_search()`
    (one parameter grid) and `tune_per_combo()` (the same grid searched separately per
    (n_cells, n_genes) combination), used by `synthetic/gridsearch/cyclum_gs.R`/`scPrisma_gs.R`/`oscope_gs.R`.
  - `download_geo.R` — shared GEOquery bootstrap + skip-if-present download
    (`download_geo_supp_if_missing()`), used by every dataset's `download*.R`.
  - `run_step.sh` — shared `run_step()` bash helper (runs a pipeline step as its own
    `Rscript` process), sourced by every dataset's `run_pipeline.sh`.

## Environment setup

### Open this repo as an RStudio Project

Open `gene-oscillations.Rproj` in RStudio (not just the folder). This makes RStudio set the
working directory to the repo root automatically — every script here assumes that, via
`proj_root <- here::here(); setwd(proj_root)` — so relative paths resolve correctly regardless
of which script you run or click "Source" on. `*.Rproj` is gitignored, so this file has to exist
locally in every clone; without it, `here::here()` falls back to whatever directory RStudio
happened to start in (e.g. your home directory), and scripts fail with "file not found" errors
for paths like `data/...` or `algorithms/...`.

### R prerequisites

```r
install.packages(c("reticulate", "here", "data.table", "Matrix"))
install.packages("BiocManager")
BiocManager::install("rhdf5")
install.packages("hdf5r")
install.packages("testthat")
```

### Seurat (for nfkb preprocessing)

`nfkb/preprocess.R` uses Seurat for QC filtering and HTO demultiplexing - not needed for
myc/segmentation_clock's preprocessing, which read raw count matrices directly.

```r
install.packages("Seurat")
```

### GSEA / pathway-analysis packages (for analyse.R)

Every dataset's `analyse.R` (and the `common/gsea_functions.R`/`common/plotting_utils.R` it
sources) needs a further set of Bioconductor and CRAN packages, on top of the ones above - these
aren't required for preprocessing or running the algorithms themselves, only for the
Hallmark/GO/DoRothEA/TFEA.ChIP pathway analysis step:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "clusterProfiler", "enrichplot", "BiocParallel", "dorothea", "TFEA.ChIP",
  "GenomicRanges", "ExperimentHub",
  "org.Hs.eg.db",                      # human - myc, segmentation_clock (hIPSC)
  "org.Mm.eg.db",                      # mouse - nfkb, segmentation_clock (mme95/mESC)
  "TxDb.Hsapiens.UCSC.hg19.knownGene"  # myc only
))
install.packages(c(
  "dplyr", "readr", "msigdbr", "ggplot2", "ggpubr", "ggrepel", "patchwork"
))
```

`org.Hs.eg.db`/`org.Mm.eg.db` and `TxDb.Hsapiens.UCSC.hg19.knownGene` are large annotation
downloads, so this will take a while - only install the species-specific ones you actually need
if you want to skip the rest (e.g. `org.Hs.eg.db` alone covers `myc/analyse.R`).

### Python / conda environments

Cyclum and scPrisma need separate conda environments (`cyclum_env` and `scPrisma_env`) because
they pin different, sometimes conflicting, dependency versions. The R driver scripts select the
environment via `reticulate::use_condaenv("<env_name>", required = TRUE)`, so the **env names
must match exactly**: `cyclum_env` and `scPrisma_env`.

#### scPrisma_env

```bash
conda create -n scPrisma_env python=3.10 -y
conda activate scPrisma_env

pip install torch==2.10.0
pip install \
  numpy==2.2.6 \
  scipy==1.15.3 \
  pandas==2.3.3 \
  scikit-learn==1.7.2 \
  h5py==3.15.1 \
  matplotlib==3.10.8 \
  seaborn==0.13.2 \
  anndata==0.11.4 \
  scanpy==1.11.5 \
  umap-learn==0.5.11

# scPrisma isn't published on PyPI — install from source. setup.py on `master` currently
# declares version 0.0.5.
pip install "git+https://github.com/nitzanlab/scPrisma.git"
```

Verify the install:

```bash
python -c "
import scanpy, anndata, torch, scPrisma, scPrisma.algorithms_torch
print('scanpy', scanpy.__version__)
print('torch', torch.__version__)
print('OK')
"
```

The exact pinned versions above were captured from a working Windows install; on other
platforms (e.g. Apple Silicon macOS) pip will resolve matching platform-specific wheels, and
`torch` will automatically use the GPU backend available on that platform (CUDA on
Linux/Windows, MPS on Apple Silicon) — no separate install step needed for that.

#### cyclum_env

The reference install for this env (Windows) uses Python 3.7 and TensorFlow 2.4.0, because
Cyclum's own package metadata pins `python_requires='>=3.6, <3.8'`. **On Apple Silicon macOS,
that combination cannot run at all** — that generation of TensorFlow requires AVX CPU
instructions, and Rosetta 2 (Apple's x86_64-on-arm64 translator) does not support AVX. This
isn't a config issue; it's a hard ceiling, even running CPU-only.

The working alternative on Apple Silicon is to bypass the `<3.8` pin (it reflects what was
available in 2020, not an actual code constraint) and use a modern, natively-arm64 Python +
TensorFlow instead:

```bash
conda create -n cyclum_env python=3.10 -y
conda activate cyclum_env

# TF 2.16+ switches its default Keras to Keras 3, which breaks Cyclum's 2019-era tf.keras
# API usage. Pin to the 2.15 line, which still defaults to Keras 2 and matches Cyclum's code.
pip install tensorflow==2.15.0
pip install tensorflow-metal   # Apple GPU acceleration plugin
pip install pandas scikit-learn matplotlib

# Not published on PyPI — install from source, ignoring the outdated python_requires ceiling.
pip install --ignore-requires-python "git+https://github.com/KChen-lab/Cyclum.git"
```

**Known issue — must run CPU-only on Apple Silicon:** `tensorflow-metal`'s graph optimizer
currently crashes (`Mutation::Apply error ... fanout ... exist for missing node`) partway
through Cyclum's multi-phase training, once the model is rebuilt with a different number of
linear dimensions — this is a known TensorFlow-on-Metal issue with the modern Adam optimizer
graph, not a Cyclum bug (TF itself warns `tf.keras.optimizers.Adam` has problems on M1/M2
Macs). `algorithms/run_cyclum.py` already works around this by disabling the GPU on macOS
(`tf.config.set_visible_devices([], "GPU")`), so training runs correctly, just CPU-only, on
Apple Silicon. Other platforms (Linux/Windows with a real GPU) are unaffected by that check and
keep full GPU acceleration.

Verify the install:

```bash
python -c "
import cyclum.tuning, cyclum.models, tensorflow as tf
print('tensorflow', tf.__version__)
print('OK')
"
```

## Running tests

```r
source("tests/testthat.R")
```

or from the shell:

```bash
Rscript tests/testthat.R
```

Tests live in `tests/testthat/`, one file per module under test (`test-hdfrw.R`,
`test-utils.R`, `test-grid_search.R`, `test-plotting_utils.R`, `test-gsea_functions.R`,
`test-dyngen_utils.R`, `test-algorithms.R`, `test-run_oscope.R`), plus `test-syntax.R`, which
parses every `.R` file in the repo and catches syntax errors even in scripts that aren't
otherwise unit-testable (see below).

**What's covered**: the reusable function libraries in `common/`, `synthetic/utils/dyngen_utils.R`,
and `algorithms/run_cyclum.R`/`run_scPrisma.R`/`run_oscope.R` - argument-building, data
wrangling, file I/O, control flow, and error handling, all exercised with small in-memory
fixtures and mocked external calls (no conda envs, no network, no real training runs).

**What isn't covered, and why**: the per-dataset orchestration scripts (`preprocess.R`/
`preprocess_*.R`, `run_cyclum.R`/`run_scPrisma.R` in `myc/`/`nfkb`/`segmentation_clock`/`synthetic`,
every `analyse.R`) aren't unit-tested directly - their logic isn't factored into standalone
functions, and most of them `library()` heavy Bioconductor packages
(clusterProfiler/msigdbr/dorothea/TFEA.ChIP/Oscope/...) or drive real conda environments via
`reticulate`, neither of which is guaranteed to be available in every environment this suite
runs in. `test-syntax.R` still catches syntax regressions in all of them. Tests that need a
package not installed in the current environment (e.g. `org.Hs.eg.db`/`org.Mm.eg.db`) skip
cleanly with `skip_if_not_installed()` rather than failing, and will run for real in an
environment with the full package set installed (see "Python / conda environments" above and
the Bioconductor packages listed throughout this README).

Several test files rely on two non-obvious techniques to test code that would otherwise require
those unavailable packages - both are documented in `tests/testthat/helper-setup.R` and at the
point of use:
- Temporarily stubbing `library()` to a no-op so a file can be `source()`d (and its *logic*
  tested with mocked replacements for the specific functions actually exercised) without its
  own `library(clusterProfiler)`-style calls erroring.
- Explicitly targeting `.GlobalEnv` when stubbing/mocking a function, rather than just
  reassigning it inside a `test_that()` block - `source()` always evaluates the sourced code in
  `.GlobalEnv` regardless of where it's called from, and `testthat::test_file()` doesn't carry
  attached-package or global-binding state from a test file's own top-level code into its
  individual `test_that()` blocks, so a plain local reassignment often isn't visible to the code
  under test.

## Synthetic data generation (optional)

`synthetic/generate_datasets.R` produces the synthetic benchmark datasets (via
[dyngen](https://github.com/dynverse/dyngen)) used by `synthetic/evaluate.R` and the grid-search
scripts (`synthetic/gridsearch/cyclum_gs.R`, `synthetic/gridsearch/scPrisma_gs.R`,
`synthetic/gridsearch/oscope_gs.R`).
**This is optional** — the exact datasets used in this project will also be uploaded to Zenodo
(link TBD) for direct download, so you only need this section if you want to generate
new/different synthetic data yourself rather than using the pre-generated ones.

### Install dyngen and its dependencies

```r
install.packages(c("dyngen", "tidyverse", "dynwrap"))
```

`dynwrap` isn't a direct dependency of `dyngen` itself, but `synthetic/run_dyngen.R` calls
`dyngen::generate_dataset(model, format = "dyno")`, and that specific output format needs
`dynwrap` to wrap the result — without it you'll hit `there is no package called 'dynwrap'`
right at the last step of an otherwise-successful simulation.

### Install Oscope

`synthetic/run_oscope.R`/`gridsearch/oscope_gs.R` need Oscope, a Bioconductor package
(`algorithms/run_oscope.R` calls into it directly - Oscope runs in-process rather than via a
separate conda env, see "Repository layout" above):

```r
install.packages("BiocManager")
BiocManager::install("Oscope")
```

### Real network / real count reference data

`synthetic/run_dyngen.R` grounds each simulation in a real gene regulatory network topology and
a real scRNA-seq count distribution (rather than an arbitrary/unrealistic one), loaded from two
local files that aren't checked into the repo (`synthetic/data/` is gitignored):

- **Realnet**: `regulatorycircuits_26_neuron-associated_cells_cancer`
- **Realcount**: `zenodo_1443566_real_gold_psc-astrocyte-maturation-neuron_sloan`

These are the specific ones used for the datasets in this project — `synthetic/run_dyngen.R`
downloads them automatically (to `synthetic/data/realnet.rds` / `synthetic/data/realcount.rds`)
the first time it's sourced and reuses the local copy on subsequent runs, so nothing extra needs
to be done here. Browse `dyngen::realnets` / `dyngen::realcounts` for other available options if
you want to swap these out for a different simulation grounding.

### Generating data

```r
source("synthetic/generate_datasets.R")
```

For each of the 11 (n_cells, n_genes) combinations in the two cells/genes sweeps, this produces
**two independent 5-replicate sets**, each in its own `c<n_cells>g<n_genes>` subdirectory:

- **Tuning set** — `synthetic/data/dyngen/gridsearch/c<n_cells>g<n_genes>/`
- **Evaluation set** — `synthetic/data/dyngen/c<n_cells>g<n_genes>/`

The split exists so hyperparameters are never chosen and scored on the same data: the tuning set
is what `gridsearch/cyclum_gs.R`/`scPrisma_gs.R`/`oscope_gs.R` search over, and the evaluation set is what
`evaluate.R` reports performance on. Each combination gets its own subdirectory (rather than one
flat directory) because `algorithms/run_cyclum.py`/`run_scPrisma.py` apply one fixed
hyperparameter setting to every `.h5` file in whatever directory they're pointed at — different
combinations need different (tuned) hyperparameters, so they have to be physically separated.

Errors during generation (e.g. a simulation failing for a particular cells/genes/replicate
combination) are caught, logged to a `generation_errors.log` in the relevant combination's
subdirectory, and the script continues on to the next replicate rather than halting entirely.

### Hyperparameter tuning (optional but recommended)

```r
source("synthetic/gridsearch/cyclum_gs.R")
source("synthetic/gridsearch/scPrisma_gs.R")
source("synthetic/gridsearch/oscope_gs.R")
```

Each script searches its algorithm's hyperparameter grid **separately for every (n_cells,
n_genes) combination**, scoring each setting by mean AUC against dyngen's ground-truth cycling
labels averaged across that combination's 5 tuning replicates. The winning hyperparameters per
combination are saved to `synthetic/data/dyngen/gridsearch/best_hyperparams_<algorithm>.csv`.

If this has been run, `synthetic/run_cyclum.R`/`run_scPrisma.R`/`run_oscope.R` (below) pick up
each combination's own best hyperparameters automatically when run afterwards; otherwise they
fall back to fixed defaults. This step is substantially more compute than the rest of the
pipeline (the full hyperparameter grid × 11 combinations × 5 tuning replicates each), so it's
not included in `run_pipeline.sh` by default.

### Evaluating results

```r
source("synthetic/evaluate.R")  # AUROC/AUPRC/accuracy/F1/precision/recall vs. cells/genes sweeps
```

Scores the **evaluation set** (never the tuning set the grid search used) with whichever
hyperparameters `run_cyclum.R`/`run_scPrisma.R`/`run_oscope.R` were run with.

Figures are written to `synthetic/figures/`.

Or run generation, Cyclum, scPrisma, Oscope, and evaluation in one go with
`./synthetic/run_pipeline.sh` (hyperparameter tuning isn't included — see the script's own
comments for how to layer it in first).

## Running the segmentation clock pipeline

`preprocess_mme95.R` / `preprocess_mesc.R` / `preprocess_hipsc.R` each auto-download their
required raw data from GEO (accession `GSE114186`, via `GEOquery::getGEOSuppFiles()`) the first
time they're run, if it isn't already present under `segmentation_clock/data/GSE114186/` — no
separate download step needed. `preprocess_mme95.R` additionally needs
`segmentation_clock/data/mme95/cell_groupings.csv` and `original_cell_indices.txt`; these are a
manual lasso-selection export from the interactive
[SPRING viewer](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/Diaz2019/E95PM/full)
session, not a GEO supplementary file, so they can't be fetched automatically — the script fails
with a clear error if they're missing.

```r
# 1. Build per-cluster raw UMI count .h5 files from the raw GEO data (auto-downloaded if missing)
source("segmentation_clock/preprocess_mme95.R")

# 2. Run each algorithm (each activates its own conda env via reticulate)
source("segmentation_clock/run_cyclum.R")
source("segmentation_clock/run_scPrisma.R")

# 3. Compare scores against reference gene lists, plus per-cluster GSEA/DoRothEA/TFEA.ChIP
source("segmentation_clock/analyse.R")
```

Results are written to `segmentation_clock/results/cyclum/` and
`segmentation_clock/results/scPrisma/`.

`analyse.R` compares predicted cycling scores against three reference gene lists
(`segmentation_clock/data/{oscillating,housekeeping,cellcycle}_{mouse,human}.rds`), sourced from:
- **oscillating** — [Matsuda et al. 2020](https://www.nature.com/articles/s41586-020-2144-9)
  (Nature), the HES7-luciferase hIPSC segmentation clock paper this pipeline models
- **housekeeping** — [Eisenberg & Levanon 2013](https://www.cell.com/trends/genetics/fulltext/S0168-9525(13)00089-9)
  (Trends in Genetics)
- **cellcycle** — [Xue et al. 2020](https://www.nature.com/articles/s41586-019-1884-x) (Nature)

Or run every step above (plus evaluation) in one go with `./segmentation_clock/run_pipeline.sh`.

## Running the MYC oscillations (HCT116) pipeline

`preprocess.R` auto-downloads the raw GEO data (accession `GSM4286760`, via
`GEOquery::getGEOSuppFiles()`) the first time it's run, if it isn't already present under
`myc/data/GSM4286760/` — no separate download step needed.

```r
# 1. Build the raw UMI count .h5 file from the raw GEO data (auto-downloaded if missing)
source("myc/preprocess.R")

# 2. Run each algorithm (each activates its own conda env via reticulate)
source("myc/run_cyclum.R")
source("myc/run_scPrisma.R")

# 3. GSEA / DoRothEA / TFEA.ChIP analysis on each algorithm's results
source("myc/analyse.R")
```

Results are written to `myc/results/cyclum/` and
`myc/results/scPrisma/`; figures to `myc/figures/`.

Or run every step above in one go with `./myc/run_pipeline.sh`.

## Running the NF-kB oscillations pipeline

`preprocess.R` auto-downloads and unpacks the raw GEO data (accession `GSE162992`, via
`GEOquery::getGEOSuppFiles()`) the first time it's run, if it isn't already present under
`nfkb/data/GSE162992/` — no separate download step needed.

```r
# 1. Demultiplex by HTO, QC-filter, and split into per-condition raw UMI count .h5 files
#    (raw data auto-downloaded if missing)
source("nfkb/preprocess.R")

# 2. Run each algorithm (each activates its own conda env via reticulate)
source("nfkb/run_cyclum.R")
source("nfkb/run_scPrisma.R")

# 3. GSEA / DoRothEA / TFEA.ChIP analysis and cross-condition score comparison
source("nfkb/analyse.R")
```

Results are written to `nfkb/results/cyclum/` and `nfkb/results/scPrisma/`; figures to
`nfkb/figures/`. `nfkb/results/scPrisma0_1/` and `nfkb/results/scPrisma0_2/` are archived scPrisma
runs at alternative regularisation strengths (0.1 and 0.2), kept alongside the main run for
comparison.

Or run every step above in one go with `./nfkb/run_pipeline.sh`.
