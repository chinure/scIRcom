# scIRcom

A local adaptation of **scRepertoire: A toolkit for single-cell immune receptor
profiling**, originally authored by Nick Borcherding. The inherited functions
combine immune-receptor contigs, summarize clonotypes, and connect repertoire
measurements to Seurat or SingleCellExperiment objects. This repository maintains
selected import and compatibility changes; it does not introduce an independent
immune-repertoire algorithm.

Local maintainer: **Yao Zhi**, zhiyao@jlu.edu.cn.

## Origin, citation and license

- Upstream project: [ncborcherding/scRepertoire](https://github.com/ncborcherding/scRepertoire).
- Original paper: Borcherding N, Bormann NL, Kraus G. *scRepertoire: An R-based
  toolkit for single-cell immune receptor analysis*. F1000Research. 2020;9:47.
  [doi:10.12688/f1000research.22139.1](https://doi.org/10.12688/f1000research.22139.1).
- This is an older 1.7.x-derived development snapshot, not the current upstream or
  Bioconductor release. The exact original import commit was not recorded. The
  matching upstream development-era
  [DESCRIPTION, version 1.7.2](https://github.com/ncborcherding/scRepertoire/blob/43ae3a9cb903798f1f092dc5b15aeca7829f67c6/DESCRIPTION)
  declares **GPL-2**. This repository follows GPL-2; the previous unrelated
  Unlicense text has been replaced by the GPL version 2 text in [LICENSE](LICENSE).
  Upstream authorship is preserved in DESCRIPTION, with local maintenance
  separately attributed.
- STARTRAC and powerTCR are upstream dependencies/wrappers. Cite their original
  papers when using those methods; this repository does not claim their methods
  as local inventions. See the relevant function documentation.

## Installation

Use an R/Bioconductor combination compatible with the package dependencies.
DESCRIPTION lists the package requirements; dependency versions are not locked.
The obsolete `dev` branch previously listed here does not exist.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("powerTCR", "SingleCellExperiment", "SummarizedExperiment"))
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("chinure/scIRcom", ref = "main", build_vignettes = FALSE)
```

The R package name is `scIRcom`. Installing upstream `scRepertoire` from
Bioconductor installs a separate package and does not install this repository.
For current upstream tutorials, use the upstream project's documentation and
account for API changes between versions.

## Contig import

The maintained importer accepts sample directories or explicit file paths:

```r
library(scIRcom)
contigs <- loadContigs(c("sample1/outs", "sample2/outs"), format = "10X")
airr <- loadContigs("sample1/airr_rearrangement.tsv", format = "AIRR")
```

For 10X data, the canonical filename is `filtered_contig_annotations.csv`;
the older singular spelling is also accepted. AIRR directories use
`airr_rearrangement.tsv`. Input barcodes must still be reconciled with the
expression object before combining modalities.

## Checks and current limits

From a source checkout:

```bash
Rscript tests/smoke-load-contigs.R
R CMD build .
```

The regression check uses synthetic fixtures and exercises actual 10X/AIRR
import, productive-contig filtering, multiple samples, canonical/legacy filenames,
legacy WAT3R TRB/two-TRA mapping, TRUST4 six-field tuples, and missing-input
validation. It requires stringr, already listed in Imports. It does not reproduce a clinical paper.

The source files parse, the focused importer check passes, and a source package
build succeeds in the review environment. A complete `R CMD check` has **not**
passed: the required powerTCR dependency is absent in that environment. Bundled
example data and generated `man/` documentation from upstream are also absent.
WAT3R import accepts the legacy barcode_results.csv header with the first24
fields in upstream order, plus optional trailing annotation columns. A different
header/order is rejected explicitly. Its positional mapping was checked against
the [upstream scRepertoire parser](https://github.com/ncborcherding/scRepertoire/blob/a57be2f0491c022cda5a94301ae5736249a12bb6/R/loadContigs.R)
and [public WAT3R example header](https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv).
TRUST4 import accepts legacy six-field chain tuples and rejects unsupported
field counts explicitly; missing chains remain NA. A single-row fixture is
covered without matrix-dimension dropping. Neither
upstream wetlab/sequencing pipeline has been rerun on real study data.
Most inherited plotting and diversity functions have not received local
regression coverage. These limits should be resolved before presenting this as a
fully validated software release.

No source data from the maintainer's publications are bundled. Repertoire
software provenance does not establish which experiments or analyses an author
personally performed in any particular paper.
