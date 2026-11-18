#' Loading the contigs derived from single-cell sequencing
#'
#' This function generates a contig list and formats the data to allow for 
#' function with  \code{\link{combineTCR}} or \code{\link{combineTCR}}. If 
#' using data derived from filtered outputs of 10X Genomics, there is no 
#' need to use this function as the data is already compatible. The function
#' assumes if listing multiple directories, there are distinct outputs with
#' unmodified file names in them.
#' 
#' @examples
#' \dontrun{
#' dir <- c("./Sample1/outs/", "./Sample2/outs/", "./Sample3/outs/")
#' contig.list <- loadContigs(dir, format = "10X")
#' }
#' 
#' @param dir One or more sample directories, or paths to contig files
#' @param format The format of the single-cell contig, currently supporting: "10X", "AIRR", "TRUST4", and "WAT3R"
#' @export
#' @return List of contigs for further processing in scRepertoire
loadContigs <- function(dir, 
                        format = "10X") {
    
    format.list <- list("WAT3R" = "barcode_results.csv",
                        "10X" = "filtered_contig_annotations.csv",
                        "AIRR" = "airr_rearrangement.tsv",
                        "TRUST4" = "barcode_report.tsv")
    if (!is.character(format) || length(format) != 1L ||
        is.na(format) || !format %in% names(format.list)) {
        stop("Unsupported contig format; choose 10X, AIRR, TRUST4, or WAT3R.",
             call. = FALSE)
    }
    if (!is.character(dir) || length(dir) == 0L || anyNA(dir) || any(!nzchar(dir))) {
        stop("Supply at least one sample directory or contig file.", call. = FALSE)
    }
    contig.files <- vapply(dir, function(path) {
        if (!dir.exists(path)) return(path)
        candidates <- file.path(path, format.list[[format]])
        if (format == "10X") {
            candidates <- c(candidates, file.path(path, "filtered_contig_annotation.csv"))
        }
        existing <- candidates[file.exists(candidates)]
        if (!length(existing)) {
            stop(paste("Contig file does not exist in sample directory:", path),
                 call. = FALSE)
        }
        existing[[1L]]
    }, character(1), USE.NAMES = FALSE)
    if (any(!file.exists(contig.files))) {
        stop(paste("Contig path does not exist:",
                   paste(contig.files[!file.exists(contig.files)], collapse = ", ")),
             call. = FALSE)
    }

    if (format %in% c("10X", "WAT3R")) {
        df <- lapply(contig.files, read.csv)
        df <- if (format == "WAT3R") parseWAT3R(df) else parse10x(df)
    } else {
        df <- lapply(contig.files, read.delim)
        df <- if (format == "TRUST4") parseTRUST4(df) else parseAIRR(df)
    }
    return(df)
}

#Formats TRUST4 data
#' @importFrom stringr str_split
parseTRUST4 <- function(df) {
    split_chain <- function(value) {
        present <- !is.na(value) & nzchar(value)
        pieces <- str_split(value, ",")
        if (any(lengths(pieces)[present] != 6L)) {
            stop("TRUST4 import requires legacy six-field chain tuples (V,D,J,C,CDR3nt,CDR3aa).",
                 call. = FALSE)
        }
        parsed <- matrix(NA_character_, nrow = length(value), ncol = 6L)
        if (any(present)) parsed[present, ] <- do.call(rbind, pieces[present])
        parsed[parsed == "*"] <- NA
        colnames(parsed) <- c("v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
        parsed
    }
    for (i in seq_along(df)) {
        colnames(df[[i]])[1] <- "barcode"
        if (!all(c("chain1", "chain2") %in% names(df[[i]]))) {
            stop("TRUST4 input needs chain1 and chain2 columns.", call. = FALSE)
        }
        df[[i]][df[[i]] == "*"] <- NA
        chain2 <- data.frame(barcode = df[[i]]$barcode, split_chain(df[[i]]$chain1),
                             stringsAsFactors = FALSE)
        chain1 <- data.frame(barcode = df[[i]]$barcode, split_chain(df[[i]]$chain2),
                             stringsAsFactors = FALSE)
        data2 <- rbind(chain1, chain2)
        data2[data2 == ""] <- NA
        df[[i]] <- data2
    }
    df <- chain.parser(df)
    return(df)
}

# Formats the legacy WAT3R barcode_results.csv layout: one TRB and two TRA slots.
# The positional mapping is retained from scRepertoire 1.7.x (see README).
parseWAT3R <- function(df) {
    fields <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene",
                "cdr3_nt", "cdr3", "reads", "umis")
    for (i in seq_along(df)) {
        x <- df[[i]]
        if (!is.data.frame(x) || ncol(x) < 24L) {
            stop("WAT3R input requires the legacy 24-column barcode_results.csv layout.",
                 call. = FALSE)
        }
        expected <- c("BC", "TCR_Recovery", "TRB_CDR3", "TRB_CDR3nuc",
                      "TRB_CDR3_UMIcount", "TRB_CDR3_error", "TRB_nReads",
                      "TRBV", "TRBD", "TRBJ", "TRA_CDR3", "TRA_CDR3nuc",
                      "TRA_CDR3_UMIcount", "TRA_CDR3_error", "TRA_nReads",
                      "TRAV", "TRAJ", "TRA.2_CDR3", "TRA.2_CDR3nuc",
                      "TRA.2_CDR3_UMIcount", "TRA.2_CDR3_error", "TRA.2_nReads",
                      "TRAV.2", "TRAJ.2")
        if (!identical(names(x)[seq_len(24L)], expected)) {
            stop("WAT3R input requires the legacy 24-column header order; see README.",
                 call. = FALSE)
        }
        x[x == ""] <- NA
        make_chain <- function(chain, v, d, j, nt, aa, reads, umis) {
            data.frame(barcode = x[[1L]], chain = rep(chain, nrow(x)),
                       v_gene = x[[v]],
                       d_gene = if (is.null(d)) rep(NA_character_, nrow(x)) else x[[d]],
                       j_gene = x[[j]], c_gene = rep(NA_character_, nrow(x)),
                       cdr3_nt = x[[nt]], cdr3 = x[[aa]],
                       reads = x[[reads]], umis = x[[umis]],
                       stringsAsFactors = FALSE)
        }
        trb <- make_chain("TRB", 8L, 9L, 10L, 4L, 3L, 7L, 5L)
        tra1 <- make_chain("TRA", 16L, NULL, 17L, 12L, 11L, 15L, 13L)
        tra2 <- make_chain("TRA", 23L, NULL, 24L, 19L, 18L, 22L, 20L)
        parsed <- rbind(tra1, trb, tra2)
        stopifnot(identical(colnames(parsed), fields))
        df[[i]] <- parsed[order(parsed$reads, parsed$chain), , drop = FALSE]
    }
    return(df)
}

#Formats AIRR data
parseAIRR <- function(df) {
    for (i in seq_along(df)) {
        df[[i]] <- df[[i]][,c("cell_id", "locus", "v_call", "d_call", "j_call", "c_call", "junction", "junction_aa")]
        colnames(df[[i]]) <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
    }
    return(df)
}

#Loads 10x data
parse10x <- function(df) {
    for (i in seq_along(df)) {
        df[[i]] <- subset(df[[i]], chain != "Multi")
        df[[i]] <- subset(df[[i]], productive %in% c(TRUE, "TRUE", "True", "true"))
        if (nrow(df[[i]]) == 0) { stop(
            "There are 0 contigs after internal filtering -
            check the contig list to see if any issues exist 
            for productive chains", call. = FALSE) }
        df[[i]] <- subset(df[[i]], cdr3 != "None")
        df[[i]][df[[i]] == ""] <- NA
        df[[i]] <- df[[i]][with(df[[i]], order(reads, chain)),]
    }
    return(df)
}
#Grabs the chain info from v_gene
chain.parser <- function(df) {
  for (i in seq_along(df)) {
    df[[i]]$chain <- substr(df[[i]][,"v_gene"],1,3)
  }
  return(df)
}
    