# Genomic integration tests assisted by OpenAI Codex (October 2026).
hit_pairs <- function(hits) {
    sort(paste(S4Vectors::queryHits(hits), S4Vectors::subjectHits(hits)))
}

test_that("indexed genomic ranges retain original biological annotations", {
    subject <- GenomicRanges::GRanges(c("chr2", "chr1", "chr1"),
        IRanges::IRanges(c(5, 20, 1), c(10, 30, 10)),
        strand = c("-", "+", "*"), id = c("B", "C", "A"))
    names(subject) <- c("gene_B", "gene_C", "gene_A")
    index <- SuperIntervalMap(subject)
    expect_s4_class(index, "SuperIntervalMap")
    expect_identical(as_granges(index), subject)
    expect_identical(length(index), length(subject))
    expect_output(show(index), "3 ranges on 2 indexed sequences")
    query <- GenomicRanges::GRanges(c("chr1", "chr2"),
        IRanges::IRanges(c(5, 7), c(25, 8)), strand = c("+", "+"))
    for (ignore in c(FALSE, TRUE)) {
        actual <- IRanges::findOverlaps(query, index, ignore.strand = ignore)
        reference <- GenomicRanges::findOverlaps(query, subject,
            ignore.strand = ignore)
        expect_identical(hit_pairs(actual), hit_pairs(reference))
        expect_identical(IRanges::countOverlaps(query, index,
            ignore.strand = ignore), GenomicRanges::countOverlaps(query,
            subject, ignore.strand = ignore))
    }
    expect_error(SuperIntervalMap(1:3), "GRanges")
})

test_that("chromosome, strand, duplicates and unsorted input match GenomicRanges", {
    set.seed(100)
    make_ranges <- function(n) {
        starts <- sample(1:200, n, replace = TRUE)
        GenomicRanges::GRanges(sample(c("chr1", "chr2", "chrX"), n, TRUE),
            IRanges::IRanges(starts, width = sample(1:30, n, TRUE)),
            strand = sample(c("+", "-", "*"), n, TRUE))
    }
    subject <- make_ranges(200)
    subject <- c(subject, subject[1:10])
    query <- make_ranges(100)
    names(query) <- paste0("query_", seq_along(query))
    index <- SuperIntervalMap(subject)
    for (ignore in c(FALSE, TRUE)) {
        for (minoverlap in c(0L, 1L, 2L)) {
            actual <- IRanges::findOverlaps(query, index,
                ignore.strand = ignore, minoverlap = minoverlap)
            reference <- GenomicRanges::findOverlaps(query, subject,
                ignore.strand = ignore, minoverlap = minoverlap)
            expect_identical(hit_pairs(actual), hit_pairs(reference))
            expect_identical(IRanges::countOverlaps(query, index,
                ignore.strand = ignore, minoverlap = minoverlap),
                GenomicRanges::countOverlaps(query, subject,
                ignore.strand = ignore, minoverlap = minoverlap))
            for (select in c("first", "last")) {
                expect_identical(IRanges::findOverlaps(query, index,
                    select = select, ignore.strand = ignore,
                    minoverlap = minoverlap),
                    GenomicRanges::findOverlaps(query, subject,
                    select = select, ignore.strand = ignore,
                    minoverlap = minoverlap))
            }
        }
    }
})

test_that("zero widths and advanced overlap options use standard semantics", {
    subject <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(1, 10, 20), width = c(5, 0, 10)))
    query <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(5, 10, 15), width = c(0, 5, 10)))
    index <- SuperIntervalMap(subject)
    for (type in c("any", "start", "end", "within", "equal")) {
        actual <- IRanges::findOverlaps(query, index, type = type)
        reference <- GenomicRanges::findOverlaps(query, subject, type = type)
        expect_identical(hit_pairs(actual), hit_pairs(reference))
        expect_identical(IRanges::countOverlaps(query, index, type = type),
            GenomicRanges::countOverlaps(query, subject, type = type))
    }
    expect_identical(hit_pairs(IRanges::findOverlaps(query, index, maxgap = 3L)),
        hit_pairs(GenomicRanges::findOverlaps(query, subject, maxgap = 3L)))
    empty <- GenomicRanges::GRanges()
    expect_length(IRanges::findOverlaps(empty, index), 0)
    expect_identical(IRanges::countOverlaps(query, SuperIntervalMap(empty)),
        integer(length(query)))
    expect_error(IRanges::countOverlaps(query,
        SuperIntervalMap(subject[1]), ignore.strand = NA), "ignore.strand")
})


test_that("native genomic indexes survive collection and reject serialization", {
    subject <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(1, 1, 1), c(10, 10, 10)),
        strand = c("+", "-", "*"))
    index <- SuperIntervalMap(subject)
    rm(subject)
    invisible(gc())
    query <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(rep(5L, 3), rep(5L, 3)),
        strand = c("+", "-", "*"))
    expect_identical(unname(IRanges::countOverlaps(query, index)), c(2L, 2L, 3L))
    restored <- unserialize(serialize(index, NULL))
    expect_error(IRanges::countOverlaps(query, restored), "recreate maps")
})

test_that("queries on circular sequences match GenomicRanges", {
    si <- Seqinfo::Seqinfo(c("chr1", "chrM"), c(1000L, 100L), c(FALSE, TRUE))
    subject <- GenomicRanges::GRanges(c("chrM", "chrM", "chr1"),
        IRanges::IRanges(c(90, 10, 5), c(100, 20, 15)),
        strand = c("+", "-", "+"), seqinfo = si)
    # The first query wraps past the end of the circular chrM.
    query <- GenomicRanges::GRanges(c("chrM", "chrM", "chr1"),
        IRanges::IRanges(c(-5, 15, 1), c(5, 18, 10)),
        strand = c("+", "*", "+"), seqinfo = si)
    index <- SuperIntervalMap(subject)
    for (ignore in c(FALSE, TRUE)) {
        expect_identical(
            hit_pairs(IRanges::findOverlaps(query, index,
                                            ignore.strand = ignore)),
            hit_pairs(GenomicRanges::findOverlaps(query, subject,
                                                  ignore.strand = ignore)))
        expect_identical(
            IRanges::countOverlaps(query, index, ignore.strand = ignore),
            GenomicRanges::countOverlaps(query, subject,
                                         ignore.strand = ignore))
    }
})

test_that("wrapping queries and subjects on circular sequences match", {
    si <- Seqinfo::Seqinfo("chrM", 100L, TRUE)
    query <- GenomicRanges::GRanges("chrM",
        IRanges::IRanges(c(-20, 95, 50, -300, 1), c(3, 130, 60, 10, 100)),
        strand = c("+", "-", "*", "+", "*"), seqinfo = si)
    inside <- GenomicRanges::GRanges("chrM",
        IRanges::IRanges(c(1, 90, 40, 98), c(5, 100, 55, 99)),
        strand = c("+", "-", "+", "*"), seqinfo = si)
    # A subject range past the end of the circle uses the GenomicRanges path.
    wrapping <- c(inside, GenomicRanges::GRanges("chrM",
        IRanges::IRanges(95, 110), strand = "+", seqinfo = si))
    for (subject in list(inside, wrapping)) {
        index <- SuperIntervalMap(subject)
        for (ignore in c(FALSE, TRUE)) {
            expect_identical(
                hit_pairs(IRanges::findOverlaps(query, index,
                                                ignore.strand = ignore)),
                hit_pairs(GenomicRanges::findOverlaps(query, subject,
                                                      ignore.strand = ignore)))
            expect_identical(
                IRanges::countOverlaps(query, index, ignore.strand = ignore),
                GenomicRanges::countOverlaps(query, subject,
                                             ignore.strand = ignore))
        }
    }
})
