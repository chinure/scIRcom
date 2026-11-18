# Run from repository root: Rscript tests/smoke-load-contigs.R
# Uses synthetic 10X/AIRR fixtures; no patient data or package installation needed.
source('R/loadContigs.R')
fixture_dir <- tempfile('contig-fixtures-')
dir.create(fixture_dir)
x <- data.frame(barcode = c('cell1', 'cell2', 'cell3'),
                chain = c('TRB', 'TRA', 'TRB'),
                productive = c(TRUE, TRUE, FALSE),
                cdr3 = c('CASSL', 'CAVR', 'CASST'), reads = c(20, 10, 30))
write.csv(x, file.path(fixture_dir, 'filtered_contig_annotations.csv'), row.names = FALSE)
loaded <- loadContigs(fixture_dir, format = '10X')
stopifnot(length(loaded) == 1L, nrow(loaded[[1]]) == 2L,
          identical(loaded[[1]]$barcode, c('cell2', 'cell1')),
          all(loaded[[1]]$productive))
by_file <- loadContigs(file.path(fixture_dir, 'filtered_contig_annotations.csv'))
stopifnot(identical(loaded, by_file))
# Older singular filenames remain readable.
legacy <- tempfile('legacy-contigs-'); dir.create(legacy)
file.copy(file.path(fixture_dir, 'filtered_contig_annotations.csv'),
          file.path(legacy, 'filtered_contig_annotation.csv'))
stopifnot(identical(loadContigs(legacy), loaded))
# Sample directories are read independently and retain supplied order.
stopifnot(length(loadContigs(c(legacy, fixture_dir))) == 2L)
airr <- data.frame(cell_id = 'airr-cell', locus = 'TRB', v_call = 'TRBV1',
                   d_call = 'TRBD1', j_call = 'TRBJ1', c_call = 'TRBC1',
                   junction = 'TGTGCCAGC', junction_aa = 'CAS')
write.table(airr, file.path(fixture_dir, 'airr_rearrangement.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)
y <- loadContigs(fixture_dir, format = 'AIRR')[[1]]
stopifnot(y$barcode == 'airr-cell', y$chain == 'TRB', y$cdr3 == 'CAS')
expect_error <- function(expr, pattern) {
  err <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(err, 'error'), grepl(pattern, conditionMessage(err), fixed = TRUE))
}
expect_error(loadContigs(character()), 'at least one')
expect_error(loadContigs(fixture_dir, format = 'unknown'), 'format')
expect_error(loadContigs(file.path(fixture_dir, 'missing')), 'does not exist')
# Legacy WAT3R wide output: TRB plus two possible TRA alleles.
wat3r <- as.data.frame(matrix(NA_character_, nrow = 1, ncol = 24))
names(wat3r) <- c("BC", "TCR_Recovery", "TRB_CDR3", "TRB_CDR3nuc",
  "TRB_CDR3_UMIcount", "TRB_CDR3_error", "TRB_nReads", "TRBV", "TRBD", "TRBJ",
  "TRA_CDR3", "TRA_CDR3nuc", "TRA_CDR3_UMIcount", "TRA_CDR3_error", "TRA_nReads",
  "TRAV", "TRAJ", "TRA.2_CDR3", "TRA.2_CDR3nuc", "TRA.2_CDR3_UMIcount",
  "TRA.2_CDR3_error", "TRA.2_nReads", "TRAV.2", "TRAJ.2")
wat3r[[1]] <- "wat3r-cell"
wat3r[[3]] <- "CASSB"; wat3r[[4]] <- "TGTGCCAGCAGC"
wat3r[[5]] <- 5L; wat3r[[7]] <- 30L
wat3r[[8]] <- "TRBV1"; wat3r[[9]] <- "TRBD1"; wat3r[[10]] <- "TRBJ1"
wat3r[[11]] <- "CAVRA"; wat3r[[12]] <- "TGTGCTGTG"
wat3r[[13]] <- 3L; wat3r[[15]] <- 20L
wat3r[[16]] <- "TRAV1"; wat3r[[17]] <- "TRAJ1"
wat3r[[18]] <- "CAVRB"; wat3r[[19]] <- "TGTGCTGTC"
wat3r[[20]] <- 2L; wat3r[[22]] <- 10L
wat3r[[23]] <- "TRAV2"; wat3r[[24]] <- "TRAJ2"
wat3r$InRNAseq <- TRUE # Optional annotation fields may follow the first24 core fields.
write.csv(wat3r, file.path(fixture_dir, "barcode_results.csv"), row.names = FALSE)
w <- loadContigs(fixture_dir, format = "WAT3R")[[1]]
stopifnot(identical(names(w), c("barcode", "chain", "v_gene", "d_gene", "j_gene",
  "c_gene", "cdr3_nt", "cdr3", "reads", "umis")), nrow(w) == 3L,
  identical(w$chain, c("TRA", "TRA", "TRB")),
  identical(w$v_gene, c("TRAV2", "TRAV1", "TRBV1")),
  identical(w$reads, c(10L, 20L, 30L)),
  identical(w$umis, c(2L, 3L, 5L)),
  w$cdr3[w$chain == "TRB"] == "CASSB", all(is.na(w$c_gene)),
  all(is.na(w$d_gene[w$chain == "TRA"])), w$d_gene[w$chain == "TRB"] == "TRBD1")
expect_error(parseWAT3R(list(wat3r[, 1:23])), "24-column")
misordered <- wat3r[, c(2, 1, 3:24)]
expect_error(parseWAT3R(list(misordered)), "header order")
# TRUST4 legacy tuples contain V,D,J,C,CDR3-nucleotide,CDR3-amino-acid.
# Source loading does not create a package namespace, so supply the declared import.
str_split <- stringr::str_split
trust <- data.frame(barcode = "trust-cell",
  chain1 = "TRBV1,TRBD1,TRBJ1,TRBC1,TGTGCCAGC,CASSB",
  chain2 = "TRAV1,*,TRAJ1,TRAC,TGTGCTGTG,CAVRA")
write.table(trust, file.path(fixture_dir, "barcode_report.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
t <- loadContigs(fixture_dir, format = "TRUST4")[[1]]
stopifnot(nrow(t) == 2L, identical(t$chain, c("TRA", "TRB")),
  identical(t$cdr3, c("CAVRA", "CASSB")), all(t$barcode == "trust-cell"),
  t$c_gene[t$chain == "TRB"] == "TRBC1")
invalid_trust <- trust; invalid_trust$chain1 <- "TRBV1,TRBD1"
expect_error(parseTRUST4(list(invalid_trust)), "six-field")
missing_chain <- trust; missing_chain$chain2 <- NA_character_
stopifnot(nrow(parseTRUST4(list(missing_chain))[[1]]) == 2L)
unlink(c(fixture_dir, legacy), recursive = TRUE)
cat('PASS: 10X/AIRR paths and filtering; legacy WAT3R TRB/two-TRA mapping; TRUST4 six-field tuples; validation.\n')
