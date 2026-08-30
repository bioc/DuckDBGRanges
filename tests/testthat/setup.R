# Pin DuckDB to a single thread so parallel float reductions accumulate in a
# fixed order. Otherwise a reduction can differ in its last ULP run-to-run and
# flake a tight-tolerance expectation. Test-harness only (real sessions use all
# cores); applied via configureOutOfCore() at connection setup, so it must be set
# before the first acquireDuckDBConn() call below.
options(DuckDBDataFrame.threads = 1L)

# GRanges dataset
granges_df <- data.frame(id = head(letters, 10L),
                         seqnames = rep.int(c("chr2", "chr2", "chr1", "chr3"), c(1L, 3L, 2L, 4L)),
                         start = 1:10, end = 10L, width = 10:1,
                         strand = strand(rep.int(c("-", "+", "*", "+", "-"), c(1L, 2L, 2L, 3L, 2L))),
                         score = 1:10,
                         GC = seq(1, 0, length = 10),
                         group = rep(c("gr1", "gr2", "gr3", "gr4"), 1:4))
granges_tf <- tempfile(fileext = ".parquet")
arrow::write_parquet(granges_df, granges_tf)


# GRangesList dataset
grlist_df <- data.frame(id = c("gr1", "gr2", "gr3", "gr4"),
                        seqnames = I(list("chr2", c("chr2", "chr2"),
                                          c("chr2", "chr1", "chr1"),
                                          c("chr3", "chr3", "chr3", "chr3"))),
                        start = I(list(1L, 2:3, 4:6, 7:10)),
                        end = I(list(10L, c(10L, 10L), c(10L, 10L, 10L), c(10L, 10L, 10L, 10L))),
                        width = I(list(10L, 9:8, 7:5, 4:1)),
                        strand = I(list("-", c("+", "+"), c("*", "*", "+"), c("+", "+", "-", "-"))),
                        score = I(list(1L, 2:3, 4:6, 7:10)),
                        GC = I(list(1, c(0.9, 0.8), c(0.7, 0.6, 0.5), c(0.4, 0.3, 0.2, 0.1))),
                        label = head(letters, 4L),
                        description = c("desc 1", "desc 2", "desc 3", "desc 4"))
grlist_tf <- tempfile(fileext = ".parquet")
arrow::write_parquet(grlist_df, grlist_tf)

grlist_grl <- GRangesList(gr1 = GRanges(seqnames = "chr2", IRanges(start = 1, end = 10), strand = "-"),
                          gr2 = GRanges(seqnames = c("chr2", "chr2"), IRanges(start = 2:3, end = 10), strand = c("+", "+")),
                          gr3 = GRanges(seqnames = c("chr2", "chr1", "chr1"), IRanges(start = 4:6, end = 10), strand = c("*", "*", "+")),
                          gr4 = GRanges(seqnames = c("chr3", "chr3", "chr3", "chr3"), IRanges(start = 7:10, end = 10), strand = c("+", "+", "-", "-")))
grlist_grl@unlistData@seqnames <- Rle(factor(as.character(grlist_grl@unlistData@seqnames)))
mcols(grlist_grl) <- DataFrame(score = I(list(1L, 2:3, 4:6, 7:10)),
                               GC = I(list(1, c(0.9, 0.8), c(0.7, 0.6, 0.5), c(0.4, 0.3, 0.2, 0.1))),
                               label = head(letters, 4L), description = c("desc 1", "desc 2", "desc 3", "desc 4"))


# Helper functions
checkDuckDBGRanges <- function(object, expected) {
    expect_s4_class(object, "DuckDBGRanges")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(length(object), length(expected))
    if (length(object) > 0L) {
        expect_identical(dbconn(object), acquireDuckDBConn())
        expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    }
    if (nkey(object@frame) > 0L) {
        expect_setequal(names(object), names(expected))
        object <- object[names(expected)]
        expect_identical(unname(as.vector(seqnames(object))), as.character(seqnames(expected)))
        expect_identical(unname(as.vector(start(object))), start(expected))
        expect_identical(unname(as.vector(end(object))), end(expected))
        expect_identical(unname(as.vector(width(object))), width(expected))
        expect_identical(unname(as.vector(strand(object))), as.character(strand(expected)))
        expect_setequal(seqlevels(object), seqlevels(expected))
        expect_setequal(seqlengths(object), seqlengths(expected))
        expect_setequal(isCircular(object), isCircular(expected))
        expect_setequal(genome(object), genome(expected))
        df <- as.data.frame(expected)
        # GenomicRanges >= 1.65.2 no longer puts names(expected) in row.names
        # of as.data.frame(expected); it uses default numeric row.names and
        # adds a "names" column instead.
        if (!is.null(names(expected))) {
            rownames(df) <- names(expected)
            df[["names"]] <- NULL
        }
        for (j in names(df)) {
            if (is.factor(df[[j]])) {
                df[[j]] <- as.character(df[[j]])
            }
        }
        expect_identical(as.data.frame(object)[names(expected), , drop=FALSE], df)
    }
}

checkDuckDBGRangesList <- function(object, expected) {
    expect_s4_class(object, "DuckDBGRangesList")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(length(object), length(expected))
    expect_identical(names(object), names(expected))
    expect_identical(elementNROWS(object), elementNROWS(expected))
    if (length(object) > 0L) {
        expect_identical(dbconn(object), acquireDuckDBConn())
        expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    }
    if (nkey(object@frame) > 0L) {
        expect_identical(as.list(seqnames(object)), lapply(as.list(seqnames(expected)), as.character))
        expect_identical(as.list(start(object)), as.list(start(expected)))
        expect_identical(as.list(end(object)), as.list(end(expected)))
        expect_identical(as.list(width(object)), as.list(width(expected)))
        expect_identical(as.list(strand(object)), lapply(as.list(strand(expected)), as.character))
        expect_setequal(seqlevels(object), seqlevels(expected))
        expect_setequal(seqlengths(object), seqlengths(expected))
        expect_setequal(isCircular(object), isCircular(expected))
        expect_identical(as(mcols(object), "DFrame"), mcols(expected))
    }
}
