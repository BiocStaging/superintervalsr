# Assisted-by: OpenAI Codex (October 2026), genomic adapter and tests.
# Assisted-by: Anthropic Claude (November 2026), slot accessors and vectorized
# genomic dispatch.

#' An indexed collection of genomic ranges
#'
#' Index a GRanges once for repeated chromosome- and strand-aware overlap
#' queries.
#' The original ranges, sequence information, names and metadata are retained.
#'
#' @param x,object An indexed genomic map.
#' @slot ranges The original GRanges, in input order.
#' @slot index Internal per-sequence indexes. Do not edit slots;
#'   construct a new object when the subject ranges change.
#' @return The S4 class stores indexed GRanges. Its length method returns an
#'   integer feature count; show prints a summary and returns invisible NULL.
#' @importFrom methods show
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps countOverlaps
#' @exportClass SuperIntervalMap
#' @name SuperIntervalMap-class
#' @seealso [SuperIntervalMap()]
#' @examples
#' library(GenomicRanges)
#' genes <- GRanges(c("chr1", "chr2"), IRanges(c(1, 10), c(5, 20)))
#' index <- SuperIntervalMap(genes)
#' length(index)
#' index
methods::setClass("SuperIntervalMap",
    slots = c(ranges = "GRanges", index = "list")
)

# Internal slot accessors. Other code uses these rather than direct slot access.
.sim_ranges <- function(x) x@ranges
.sim_index <- function(x) x@index

#' Create an indexed genomic subject
#'
#' @param x A GRanges containing subject features, possibly on multiple
#' sequences
#'   and strands. Coordinates, names, seqinfo and metadata are preserved.
#' @return A SuperIntervalMap S4 object for use as the subject of
#'   `findOverlaps()` and `countOverlaps()`. `as_granges()` retrieves the
#' original
#'   GRanges, and `length()` returns its number of features.
#' @details Hits refer to the original subject row order. Strand matching
#' follows
#'   GenomicRanges: `"*"` matches either strand, and `ignore.strand = TRUE`
#'   ignores strand entirely. Ordinary positive-width overlaps use
#' SuperIntervals.
#'   Queries involving zero-width ranges, a non-default gap or overlap
#' threshold,
#'   or a type other than `"any"` delegate to GenomicRanges to preserve its
#'   semantics. The object stores one index per sequence, with strand filtering
#'   performed in C++. Queries on sequences marked circular in the seqinfo
#'   (such as the mitochondrial genome) also use GenomicRanges, which
#'   wraps coordinates around the sequence end. Query functions assume the
#'   indexes are
#'   current; do not alter the slots after construction.
#'
#'   Like IntervalMap, the indexes contain external pointers. Save
#'   `as_granges(x)` and reconstruct the index after loading serialized data.
#' @export
#' @examples
#' library(GenomicRanges)
#' genes <- GRanges(c("chr1", "chr2", "chr1"),
#'     IRanges(c(20, 1, 1), c(30, 10, 10)),
#'     strand = c("+", "-", "*"), gene_id = c("C", "B", "A"))
#' index <- SuperIntervalMap(genes)
#' query <- GRanges(c("chr1", "chr2"), IRanges(c(5, 5), c(25, 8)),
#'     strand = c("+", "+"))
#' findOverlaps(query, index)
#' countOverlaps(query, index, ignore.strand = TRUE)
#' as_granges(index)
SuperIntervalMap <- function(x) {
    if (!methods::is(x, "GRanges")) {
        stop("x must be a GRanges", call. = FALSE)
    }
    # Zero-width features are kept for the GenomicRanges fallback, but cannot
    # enter the closed-interval index.
    rows <- which(GenomicRanges::width(x) > 0L)
    by_sequence <- split(rows, as.character(GenomicRanges::seqnames(x))[rows])
    # GenomicRanges strand levels are +, -, * (codes 1, 2, 3).
    strands <- as.integer(GenomicRanges::strand(x))
    starts <- GenomicRanges::start(x)
    ends <- GenomicRanges::end(x)
    index <- lapply(by_sequence, function(ids) {
        create_genomic_index(starts[ids], ends[ids], ids, strands[ids])
    })
    attr(index, "zero.width") <- any(GenomicRanges::width(x) == 0L)
    # Per-sequence coordinate bounds, used to detect subject ranges that wrap
    # around a circular sequence.
    attr(index, "min.start") <- vapply(by_sequence,
        function(ids) min(starts[ids]), integer(1))
    attr(index, "max.end") <- vapply(by_sequence,
        function(ids) max(ends[ids]), integer(1))
    methods::new("SuperIntervalMap", ranges = x, index = index)
}

# Internal batch dispatch. No per-interval index-state checks are performed.
.genomic_query <- function(query, subject, ignore.strand, count.only = FALSE) {
    if (!is.logical(ignore.strand) || length(ignore.strand) != 1L ||
        is.na(ignore.strand)) {
        stop("ignore.strand must be TRUE or FALSE", call. = FALSE)
    }
    by_sequence <- split(seq_along(query),
                            as.character(GenomicRanges::seqnames(query)))
    strands <- as.integer(GenomicRanges::strand(query))
    starts <- GenomicRanges::start(query)
    ends <- GenomicRanges::end(query)
    index <- .sim_index(subject)
    ranges <- .sim_ranges(subject)
    by_sequence <- by_sequence[names(by_sequence) %in% names(index)]
    # Circular sequences (e.g. chrM) wrap coordinates around the sequence end.
    circle <- .circle_lengths(query, ranges, names(by_sequence))
    found <- lapply(names(by_sequence), function(seqname) {
        qi <- by_sequence[[seqname]]
        len <- circle[[seqname]]
        if (is.na(len)) {
            return(query_genomic_index(index[[seqname]], starts[qi],
                ends[qi], strands[qi], ignore.strand, count.only))
        }
        if (attr(index, "min.start")[[seqname]] < 1L ||
            attr(index, "max.end")[[seqname]] > len) {
            # Wrapping subject ranges are rare; let GenomicRanges handle them.
            return(.circular_query(query[qi], ranges, seqname,
                                    ignore.strand, count.only))
        }
        .wrapped_query(index[[seqname]], starts[qi], ends[qi], strands[qi],
                        len, ignore.strand, count.only)
    })
    query_rows <- as.integer(unlist(by_sequence, use.names = FALSE))
    if (count.only) {
        counts <- integer(length(query))
        names(counts) <- names(query)
        counts[query_rows] <- unlist(found, use.names = FALSE)
        return(counts)
    }
    found <- unlist(found, recursive = FALSE, use.names = FALSE)
    S4Vectors::Hits(rep(query_rows, lengths(found)),
                    as.integer(unlist(found, use.names = FALSE)),
                    length(query), length(.sim_ranges(subject)),
                    sort.by.query = TRUE)
}

# Circle length for each sequence name, or NA when the sequence is linear
# or its length is unknown (GenomicRanges then treats it as linear too).
.circle_lengths <- function(query, subject, seqnames) {
    circular <- Seqinfo::isCircular(query)[seqnames] %in% TRUE |
        Seqinfo::isCircular(subject)[seqnames] %in% TRUE
    len <- Seqinfo::seqlengths(query)[seqnames]
    len[is.na(len)] <- Seqinfo::seqlengths(subject)[seqnames][is.na(len)]
    len[!circular] <- NA_integer_
    len <- as.integer(len)
    names(len) <- seqnames
    len
}

# Queries on a circular sequence of length `len`, with every subject range
# inside [1, len]. Ranges overlap when they share a position modulo `len`.
# Queries inside [1, len] are ordinary linear queries. A query that runs off
# either end is mapped onto the circle as one or two linear pieces.
.wrapped_query <- function(map, starts, ends, strands, len, ignore.strand,
                            count.only) {
    w <- which(starts < 1L | ends > len)
    if (length(w) == 0L) {
        return(query_genomic_index(map, starts, ends, strands,
                                    ignore.strand, count.only))
    }
    width <- as.numeric(ends[w]) - starts[w] + 1
    whole <- width >= len
    first <- (starts[w] - 1L) %% len + 1L
    last <- first + width - 1
    split_piece <- !whole & last > len
    piece_query <- c(w, w[split_piece])
    piece_start <- c(ifelse(whole, 1L, first), rep(1L, sum(split_piece)))
    piece_end <- c(ifelse(whole, len, pmin(last, len)),
                    (last - len)[split_piece])
    found <- query_genomic_index(map, as.integer(piece_start),
        as.integer(piece_end), strands[piece_query], ignore.strand, FALSE)
    pairs <- unique(data.frame(
        query = rep(piece_query, lengths(found)),
        subject = as.integer(unlist(found, use.names = FALSE))))
    wrapped <- split(pairs$subject, factor(pairs$query, levels = w))
    linear <- query_genomic_index(map, starts[-w], ends[-w], strands[-w],
                                    ignore.strand, count.only)
    out <- if (count.only) integer(length(starts)) else
        vector("list", length(starts))
    out[-w] <- linear
    out[w] <- if (count.only) lengths(wrapped) else unname(wrapped)
    out
}

# Per-query results from GenomicRanges, in the same shape as
# query_genomic_index(): counts, or a list of subject rows per query.
# Only subject rows on `seqname` are passed.
.circular_query <- function(query, subject, seqname, ignore.strand,
                            count.only) {
    rows <- which(as.logical(GenomicRanges::seqnames(subject) == seqname))
    subject <- subject[rows]
    if (count.only) {
        return(GenomicRanges::countOverlaps(query, subject,
                                            ignore.strand = ignore.strand))
    }
    hits <- GenomicRanges::findOverlaps(query, subject,
                                        ignore.strand = ignore.strand)
    by_query <- factor(S4Vectors::queryHits(hits), levels = seq_along(query))
    unname(split(rows[S4Vectors::subjectHits(hits)], by_query))
}

.genomic_fast_path <- function(query, subject, maxgap, minoverlap, type) {
    identical(maxgap, -1L) && length(minoverlap) == 1L &&
        minoverlap %in% c(0L, 1L) && type == "any" &&
        all(GenomicRanges::width(query) > 0L) &&
        !isTRUE(attr(.sim_index(subject), "zero.width"))
}

#' Query an indexed genomic subject
#'
#' @param query A GRanges of query features.
#' @param subject A SuperIntervalMap constructed once from the subject
#' GRanges.
#' @param maxgap,minoverlap,type,select Overlap options, as in
#'   [GenomicRanges::findOverlaps()].
#' @param ignore.strand Whether to ignore strand (default FALSE).
#' @return `findOverlaps()` returns a Hits object with indices in the original
#'   input order, or a subject-index vector when `select` is not `"all"`.
#'   `countOverlaps()` returns one integer count per query feature.
#' @rdname genomic-overlaps
#' @aliases findOverlaps,GRanges,SuperIntervalMap-method
#' @export
#' @examples
#' library(GenomicRanges)
#' subject <- GRanges(c("chr1", "chr2"), IRanges(c(1, 10), c(5, 20)))
#' query <- GRanges(c("chr1", "chr2"), IRanges(c(3, 15), c(8, 18)))
#' index <- SuperIntervalMap(subject)
#' hits <- findOverlaps(query, index)
#' as_granges(index)[subjectHits(hits)]
#' countOverlaps(query, index)
methods::setMethod("findOverlaps",
    signature(query = "GRanges", subject = "SuperIntervalMap"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
                type = c("any", "start", "end", "within", "equal"),
                select = c("all", "first", "last", "arbitrary"),
                ignore.strand = FALSE) {
        type <- match.arg(type)
        select <- match.arg(select)
        if (!.genomic_fast_path(query, subject, maxgap, minoverlap, type)) {
            return(GenomicRanges::findOverlaps(query, .sim_ranges(subject),
                maxgap = maxgap, minoverlap = minoverlap, type = type,
                select = select, ignore.strand = ignore.strand))
        }
        hits <- .genomic_query(query, subject, ignore.strand)
        S4Vectors::selectHits(hits, select = select)
    }
)

#' @rdname genomic-overlaps
#' @aliases countOverlaps,GRanges,SuperIntervalMap-method
#' @export
methods::setMethod("countOverlaps",
    signature(query = "GRanges", subject = "SuperIntervalMap"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
                type = c("any", "start", "end", "within", "equal"),
                ignore.strand = FALSE) {
        type <- match.arg(type)
        if (!.genomic_fast_path(query, subject, maxgap, minoverlap, type)) {
            return(GenomicRanges::countOverlaps(query, .sim_ranges(subject),
                maxgap = maxgap, minoverlap = minoverlap, type = type,
                ignore.strand = ignore.strand))
        }
        .genomic_query(query, subject, ignore.strand, count.only = TRUE)
    }
)

#' @rdname SuperIntervalMap-class
#' @export
methods::setMethod("length", "SuperIntervalMap", function(x) {
    length(.sim_ranges(x))
})

#' @rdname SuperIntervalMap-class
#' @export
methods::setMethod("show", "SuperIntervalMap", function(object) {
    cat("SuperIntervalMap with", length(object), "ranges on",
        length(.sim_index(object)), "indexed sequences\n")
})
