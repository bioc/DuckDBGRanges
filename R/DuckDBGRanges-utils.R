#' Utilities for DuckDBGRanges objects
#'
#' @description
#' Various utility methods for DuckDBGRanges objects including overlap detection,
#' intra-range transformations, inter-range operations, and comparison methods.
#' These methods translate R operations into efficient SQL queries executed
#' directly in DuckDB.
#'
#' @details
#' The DuckDBGRanges methods translate genomic range operations into SQL queries,
#' enabling efficient manipulation of large disk-backed genomic annotations
#' without loading them into memory.
#'
#' @section Overlap Methods:
#' In the code snippets below, \code{query} and \code{subject} can be GRanges
#' or DuckDBGRanges objects, with at least one being a DuckDBGRanges:
#' \describe{
#'   \item{\code{findOverlaps(query, subject, maxgap=-1L, minoverlap=0L,
#'     type=c("any", "start", "end", "within", "equal"),
#'     select=c("all", "first", "last", "arbitrary"),
#'     ignore.strand=FALSE)}:}{
#'     Find overlapping intervals between query and subject via SQL JOIN.
#'     Returns a \linkS4class{Hits} object.
#'   }
#'   \item{\code{countOverlaps(query, subject, maxgap=-1L, minoverlap=0L,
#'     type=c("any", "start", "end", "within", "equal"),
#'     ignore.strand=FALSE)}:}{
#'     Count overlapping intervals for each query element.
#'     Returns an integer vector.
#'   }
#'   \item{\code{overlapsAny(query, subject, maxgap=-1L, minoverlap=0L,
#'     type=c("any", "start", "end", "within", "equal"), ...)}:}{
#'     Test if each query element overlaps any element in subject.
#'     Returns a logical vector.
#'   }
#'   \item{\code{subsetByOverlaps(x, ranges, maxgap=-1L, minoverlap=0L,
#'     type=c("any", "start", "end", "within", "equal"),
#'     invert=FALSE, ...)}:}{
#'     Subset x to elements that overlap ranges.
#'     When x is a DuckDBGRanges, returns a DuckDBGRanges object.
#'     When x is a GRanges, returns a GRanges object.
#'     Implemented using \code{overlapsAny} internally.
#'   }
#' }
#'
#' @section Intra-range Methods:
#' Intra-range methods transform individual ranges. The \code{x} parameter is
#' a DuckDBGRanges object, and these methods return a GRanges object:
#' \describe{
#'   \item{\code{shift(x, shift=0L, use.names=TRUE)}:}{
#'     Shift all ranges by the specified amount. Returns a GRanges object.
#'   }
#'   \item{\code{narrow(x, start=NA, end=NA, width=NA, use.names=TRUE)}:}{
#'     Narrow ranges by specifying new start/end/width relative to current ranges.
#'     Returns a GRanges object.
#'   }
#'   \item{\code{resize(x, width, fix="start", use.names=TRUE, ignore.strand=FALSE)}:}{
#'     Resize ranges to specified width, anchored by fix position.
#'     Returns a GRanges object.
#'   }
#'   \item{\code{flank(x, width, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE)}:}{
#'     Get flanking regions of specified width. Returns a GRanges object.
#'   }
#'   \item{\code{promoters(x, upstream=2000, downstream=200, use.names=TRUE)}:}{
#'
#'     Get promoter regions (upstream and downstream of TSS).
#'     Returns a GRanges object.
#'   }
#'   \item{\code{terminators(x, upstream=2000, downstream=200, use.names=TRUE)}:}{
#'     Get terminator regions (upstream and downstream of TES).
#'     Returns a GRanges object.
#'   }
#' }
#'
#' @section Inter-range Methods:
#' Inter-range methods operate across multiple ranges:
#' \describe{
#'   \item{\code{range(x, ..., with.revmap=FALSE, ignore.strand=FALSE, na.rm=FALSE)}:}{
#'     Returns the range (min start to max end) per seqname/strand combination.
#'     Computed via SQL MIN/MAX aggregation. Returns a DuckDBGRanges object.
#'   }
#'   \item{\code{reduce(x, drop.empty.ranges=FALSE, min.gapwidth=1L, with.revmap=FALSE,
#'     with.inframe.attrib=FALSE, ignore.strand=FALSE)}:}{
#'     Merge overlapping and adjacent ranges into single ranges per seqname/strand.
#'     Uses SQL window functions to identify merge groups. Returns a DuckDBGRanges object.
#'     Note: \code{with.revmap} and \code{with.inframe.attrib} are not supported.
#'   }
#'   \item{\code{gaps(x, start=1L, end=seqlengths(x), ignore.strand=FALSE)}:}{
#'     Find gaps (uncovered regions) between ranges per seqname/strand.
#'     Uses \code{reduce()} internally then computes gaps between consecutive ranges.
#'     Returns a DuckDBGRanges object.
#'     Note: When \code{ignore.strand=FALSE}, only returns gaps for seqname/strand
#'     combinations present in the input data, unlike GRanges which returns gaps
#'     for all possible combinations. Use \code{ignore.strand=TRUE} for exact
#'     matching behavior.
#'   }
#'   \item{\code{disjoin(x, with.revmap=FALSE, ignore.strand=FALSE)}:}{
#'     Break ranges into non-overlapping pieces at all breakpoints.
#'     Creates intervals between consecutive unique start/(end+1) positions.
#'     Returns a DuckDBGRanges object.
#'     Note: \code{with.revmap} is not supported.
#'   }
#'   \item{\code{coverage(x, shift=0L, width=NULL, weight=1L, method=c("auto", "sort", "hash"), ...)}:}{
#'     Compute per-base coverage depth per seqname (strand is ignored, as in
#'     \code{GenomicRanges::coverage}). Uses the same delta-event and
#'     window-function \code{cumsum()} SQL as \code{disjoin()}/\code{gaps()},
#'     collecting only the compact per-seqname breakpoint table before
#'     building a \code{SimpleRleList}. \code{shift} and \code{width} are
#'     recycled per seqlevel; a range shifted so that it falls (partly or
#'     entirely) before position 1 is clipped to position 1, matching
#'     \code{GenomicRanges::coverage} rather than erroring. \code{method} is
#'     accepted for API compatibility and ignored. Note: \code{weight} must be
#'     a single number; a per-range vector or an mcols column name is not
#'     supported.
#'   }
#' }
#'
#' @section Set Operations:
#' Set operations combine two DuckDBGRanges objects:
#' \describe{
#'   \item{\code{union(x, y, ignore.strand=FALSE)}:}{
#'     Union of ranges from x and y. Combines ranges then reduces overlapping ones.
#'     Returns a DuckDBGRanges object.
#'   }
#'   \item{\code{intersect(x, y, ignore.strand=FALSE)}:}{
#'     Intersection of ranges. Returns regions covered by both x and y.
#'     Uses SQL JOIN with overlap detection and GREATEST/LEAST for coordinates.
#'     Returns a DuckDBGRanges object.
#'   }
#'   \item{\code{setdiff(x, y, ignore.strand=FALSE)}:}{
#'     Set difference. Returns regions in x that are not covered by y.
#'     Computes fragments by subtracting overlapping y ranges from x.
#'     Returns a DuckDBGRanges object.
#'   }
#' }
#'
#' @section Nearest Neighbor Methods:
#' Methods for finding the nearest ranges in a subject:
#' \describe{
#'   \item{\code{precede(x, subject, select=c("first", "all"), ignore.strand=FALSE)}:}{
#'     For each range in x, find the index of the first range in subject that
#'     it precedes (i.e., the first subject range starting after x ends).
#'     Overlapping ranges are excluded. Returns an integer vector with NAs for
#'     ranges with no qualifying match, or a Hits object if select="all".
#'   }
#'   \item{\code{follow(x, subject, select=c("last", "all"), ignore.strand=FALSE)}:}{
#'     For each range in x, find the index of the last range in subject that
#'     x follows (i.e., the last subject range ending before x starts).
#'     Overlapping ranges are excluded. Returns an integer vector with NAs for
#'     ranges with no qualifying match, or a Hits object if select="all".
#'   }
#'   \item{\code{nearest(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE)}:}{
#'     For each range in x, find the index of the nearest range in subject.
#'     Distance is 0 for overlapping ranges. For ties, one match is selected
#'     arbitrarily (minimum index). Returns an integer vector or Hits object.
#'     When \code{subject} is omitted (\code{nearest(x)}), each range's nearest
#'     neighbour is found among the \emph{other} ranges of x -- a range is never
#'     its own nearest neighbour -- matching base GenomicRanges' \code{drop.self}.
#'   }
#'   \item{\code{distanceToNearest(x, subject, ignore.strand=FALSE, select=c("arbitrary", "all"))}:}{
#'     Find the nearest range and compute the distance. Returns a Hits object
#'     with a "distance" metadata column containing the distances. As with
#'     \code{nearest}, the \code{subject}-omitted form excludes self-hits.
#'   }
#' }
#'
#' @section Comparison Methods:
#' Methods for comparing DuckDBGRanges elements:
#' \describe{
#'   \item{\code{duplicated(x, incomparables=FALSE, fromLast=FALSE)}:}{
#'     Find duplicate ranges using SQL window functions. Returns a DuckDBColumn.
#'   }
#'   \item{\code{unique(x, incomparables=FALSE, fromLast=FALSE)}:}{
#'     Get unique ranges. Returns a DuckDBGRanges object.
#'   }
#'   \item{\code{match(x, table, nomatch=NA_integer_, incomparables=NULL)}:}{
#'     Find matches between DuckDBGRanges objects using SQL JOIN.
#'     Returns a DuckDBColumn.
#'   }
#' }
#'
#' @section Parallel Set Operations:
#' Parallel (element-wise) set operations between two DuckDBGRanges or
#' between DuckDBGRanges and GRanges objects. When one argument is a GRanges,
#' it is automatically converted to a DuckDBGRanges for optimized computation:
#' \describe{
#'   \item{\code{punion(x, y, fill.gap=FALSE, ignore.strand=FALSE)}:}{
#'     Parallel union of ranges. Returns a DuckDBGRanges object.
#'   }
#'   \item{\code{pintersect(x, y, drop.nohit.ranges=FALSE, ignore.strand=FALSE)}:}{
#'     Parallel intersection of ranges. Returns a DuckDBGRanges object.
#'   }
#'   \item{\code{psetdiff(x, y, ignore.strand=FALSE)}:}{
#'     Parallel set difference of ranges. Returns a DuckDBGRanges object.
#'   }
#' }
#'
#' @section Distance Methods:
#' Methods for computing distances between ranges:
#' \describe{
#'   \item{\code{distance(x, y, ignore.strand=FALSE)}:}{
#'     Compute pairwise distances between ranges.
#'     Returns an integer vector.
#'   }
#' }
#'
#' @section Ordering Methods:
#' Methods for ordering and sorting DuckDBGRanges objects:
#' \describe{
#'   \item{\code{sort(x, decreasing=FALSE, ignore.strand=FALSE)}:}{
#'     Sort ranges by seqnames, strand, start, and width using SQL ORDER BY.
#'     Returns a DuckDBGRanges object.
#'   }
#'   \item{\code{order(..., na.last=TRUE, decreasing=FALSE)}:}{
#'     Get the ordering permutation using SQL window functions.
#'     Returns an integer vector.
#'   }
#'   \item{\code{is.unsorted(x, na.rm=FALSE, strictly=FALSE, ignore.strand=FALSE)}:}{
#'     Check if ranges are unsorted with respect to genomic order.
#'     Returns a logical scalar.
#'   }
#'   \item{\code{rank(x, na.last=TRUE, ties.method=c("first", "min"), ignore.strand=FALSE)}:}{
#'     Get ranks of ranges using SQL window functions.
#'     Only \code{ties.method="first"} and \code{ties.method="min"} are supported
#'     for DuckDBGRanges (other methods like "average", "last", "random", "max"
#'     are not implemented).
#'     Returns an integer vector.
#'   }
#' }
#'
#' @section Range Restriction Methods:
#' Methods for restricting ranges to bounds:
#' \describe{
#'   \item{\code{trim(x, use.names=TRUE)}:}{
#'     Trim out-of-bound ranges on non-circular sequences to their seqlengths.
#'     Returns a DuckDBGRanges object.
#'   }
#'   \item{\code{restrict(x, start=NA, end=NA, keep.all.ranges=FALSE, use.names=TRUE)}:}{
#'     Restrict ranges to specified start/end bounds using SQL GREATEST/LEAST.
#'     Returns a DuckDBGRanges object.
#'   }
#' }
#'
#' @section Tiling Methods:
#' Methods for dividing ranges into sub-ranges:
#' \describe{
#'   \item{\code{tile(x, n, width)}:}{
#'     Divide each range into tiles. Specify either 'n' (number of tiles per
#'     range) or 'width' (tile width). Returns a GRangesList with one element
#'     per input range. Note: materializes data for processing.
#'   }
#'   \item{\code{slidingWindows(x, width, step=1L)}:}{
#'     Generate sliding windows of specified 'width' moving by 'step' positions.
#'     Returns a GRangesList. Note: materializes data for processing.
#'   }
#'   \item{\code{pgap(x, y, ignore.strand=FALSE)}:}{
#'     Compute pairwise gaps between ranges in x and y, entirely in SQL. For
#'     each pair, returns the gap region between the two ranges, or a
#'     zero-width range at the boundary if they overlap or are adjacent.
#'     Errors if a pair has incompatible seqnames, or (unless
#'     \code{ignore.strand}) incompatible strand. Returns a DuckDBGRanges
#'     object.
#'   }
#' }
#'
#' @param query,x A GRanges or DuckDBGRanges object.
#' @param subject,ranges A GRanges or DuckDBGRanges object.
#' @param y For parallel set operations, a GRanges or DuckDBGRanges object of
#'   the same length as \code{x}.
#' @param table For \code{match}, a DuckDBGRanges object to match against.
#' @param maxgap A non-negative integer. Intervals with a separation of
#'   \code{maxgap} or less are considered to be overlapping. Default is -1L
#'   (no gap allowed).
#' @param minoverlap A non-negative integer. The minimum number of base pairs
#'   that must be overlapping for two intervals to be considered overlapping.
#'   Default is 0L.
#' @param type The type of overlap. Currently only \code{"any"} is supported
#'   for DuckDBGRanges.
#' @param select When \code{"all"} (default), returns all overlapping pairs.
#'   When \code{"first"}, \code{"last"}, or \code{"arbitrary"}, returns a
#'   single match per query element.
#' @param ignore.strand If \code{TRUE}, strand information is ignored when
#'   computing operations. Default is \code{FALSE}.
#' @param invert For \code{subsetByOverlaps}, if \code{TRUE}, returns elements
#'   that do NOT overlap. Default is \code{FALSE}.
#' @param shift For \code{shift}, the amount to shift ranges (can be negative).
#' @param width For \code{resize} and \code{flank}, the target width.
#' @param fix For \code{resize}, where to anchor: "start", "end", or "center".
#' @param start,end For \code{narrow}, the new start/end positions relative to
#'   current ranges; a negative value counts back from the range end (\code{-1}
#'   is the last base), matching base \code{IRanges::narrow}. For \code{flank},
#'   if \code{TRUE} get upstream flanks.
#' @param both For \code{flank}, if \code{TRUE} get flanks on both sides.
#' @param upstream,downstream For \code{promoters}/\code{terminators}, distances.
#' @param use.names If \code{TRUE}, preserve names in output.
#' @param fill.gap For \code{punion}, if \code{TRUE} fill gaps between ranges.
#' @param drop.nohit.ranges For \code{pintersect}, if \code{TRUE} drop
#'   non-intersecting pairs.
#' @param with.revmap For \code{range}, if \code{TRUE} include reverse mapping.
#' @param na.rm,incomparables,fromLast,nomatch
#'   Standard R arguments for comparison functions.
#' @param decreasing For \code{sort} and \code{order}, if \code{TRUE} sort in
#'   decreasing order. Default is \code{FALSE}.
#' @param na.last For \code{order} and \code{rank}, handling of NA values.
#' @param strictly For \code{is.unsorted}, if \code{TRUE} check for strict
#'   ordering (no ties). Default is \code{FALSE}.
#' @param ties.method For \code{rank}, method for handling ties.
#'   Only \code{"first"} and \code{"min"} are supported for DuckDBGRanges.
#' @param keep.all.ranges For \code{restrict}, if \code{TRUE} keep ranges
#'   that become invalid after restriction. Default is \code{FALSE}.
#' @param method For \code{order}, sorting method (ignored, included for
#'   compatibility).
#' @param ... Additional arguments passed to methods.
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{DuckDBGRanges-class}} for the main class
#'   \item \code{\link[IRanges]{findOverlaps}} for the generic function
#'   \item \code{\link[GenomicRanges]{findOverlaps-methods}} for GRanges methods
#' }
#'
#' @examples
#' # Load required packages
#' library(GenomicRanges)
#'
#' # Create example DuckDBGRanges
#' df <- data.frame(
#'     seqnames = c("chr1", "chr1", "chr2", "chr2"),
#'     start = c(100, 200, 150, 300),
#'     end = c(150, 250, 200, 350),
#'     strand = c("+", "-", "+", "-"),
#'     score = 1:4
#' )
#' tf <- tempfile(fileext = ".parquet")
#' arrow::write_parquet(df, tf)
#' subject <- DuckDBGRanges(tf, seqnames = "seqnames", start = "start",
#'                          end = "end", strand = "strand", mcols = "score")
#'
#' # Create query GRanges
#' query <- GRanges(c("chr1:120-180", "chr2:100-400"))
#'
#' # Find overlaps
#' hits <- findOverlaps(query, subject)
#'
#' # Count overlaps
#' counts <- countOverlaps(query, subject)
#'
#' # Test for any overlaps
#' any_ov <- overlapsAny(query, subject)
#'
#' # Subset by overlaps (uses overlapsAny internally)
#' subset <- subsetByOverlaps(query, subject)
#'
#' @return
#' Overlap methods return a \link[S4Vectors]{Hits} object
#' (\code{findOverlaps()}), an integer vector (\code{countOverlaps()}), a
#' logical vector (\code{overlapsAny()}), or a DuckDBGRanges
#' (\code{subsetByOverlaps()}). Intra-range methods (\code{shift()},
#' \code{narrow()}, \code{resize()}, \code{flank()}, \code{promoters()},
#' \code{terminators()}) return a \link[GenomicRanges]{GRanges}. Inter-range,
#' set-operation, range-restriction, and tiling methods (such as \code{range()}
#' and \code{reduce()}) return a DuckDBGRanges (or DuckDBGRangesList where the
#' result is grouped). Nearest-neighbor, distance, comparison, and ordering
#' methods return the corresponding integer or logical vectors (or a
#' \link[S4Vectors]{Hits}).
#'
#' @aliases findOverlaps,GRanges,DuckDBGRanges-method
#' @aliases findOverlaps,DuckDBGRanges,DuckDBGRanges-method
#' @aliases findOverlaps,DuckDBGRanges,GRanges-method
#' @aliases countOverlaps,GRanges,DuckDBGRanges-method
#' @aliases countOverlaps,DuckDBGRanges,DuckDBGRanges-method
#' @aliases countOverlaps,DuckDBGRanges,GRanges-method
#' @aliases overlapsAny,GRanges,DuckDBGRanges-method
#' @aliases overlapsAny,DuckDBGRanges,DuckDBGRanges-method
#' @aliases overlapsAny,DuckDBGRanges,GRanges-method
#' @aliases subsetByOverlaps,DuckDBGRanges,DuckDBGRanges-method
#' @aliases subsetByOverlaps,DuckDBGRanges,GRanges-method
#' @aliases shift,DuckDBGRanges-method
#' @aliases narrow,DuckDBGRanges-method
#' @aliases resize,DuckDBGRanges-method
#' @aliases flank,DuckDBGRanges-method
#' @aliases promoters,DuckDBGRanges-method
#' @aliases terminators,DuckDBGRanges-method
#' @aliases range,DuckDBGRanges-method
#' @aliases reduce,DuckDBGRanges-method
#' @aliases gaps,DuckDBGRanges-method
#' @aliases disjoin,DuckDBGRanges-method
#' @aliases coverage,DuckDBGRanges-method
#' @aliases union,DuckDBGRanges,DuckDBGRanges-method
#' @aliases union,DuckDBGRanges,GRanges-method
#' @aliases union,GRanges,DuckDBGRanges-method
#' @aliases intersect,DuckDBGRanges,DuckDBGRanges-method
#' @aliases intersect,DuckDBGRanges,GRanges-method
#' @aliases intersect,GRanges,DuckDBGRanges-method
#' @aliases setdiff,DuckDBGRanges,DuckDBGRanges-method
#' @aliases setdiff,DuckDBGRanges,GRanges-method
#' @aliases setdiff,GRanges,DuckDBGRanges-method
#' @aliases duplicated,DuckDBGRanges-method
#' @aliases unique,DuckDBGRanges-method
#' @aliases match,DuckDBGRanges,DuckDBGRanges-method
#' @aliases punion,DuckDBGRanges,DuckDBGRanges-method
#' @aliases punion,DuckDBGRanges,GRanges-method
#' @aliases punion,GRanges,DuckDBGRanges-method
#' @aliases pintersect,DuckDBGRanges,DuckDBGRanges-method
#' @aliases pintersect,DuckDBGRanges,GRanges-method
#' @aliases pintersect,GRanges,DuckDBGRanges-method
#' @aliases psetdiff,DuckDBGRanges,DuckDBGRanges-method
#' @aliases psetdiff,DuckDBGRanges,GRanges-method
#' @aliases psetdiff,GRanges,DuckDBGRanges-method
#' @aliases distance,DuckDBGRanges,DuckDBGRanges-method
#' @aliases distance,DuckDBGRanges,GRanges-method
#' @aliases distance,GRanges,DuckDBGRanges-method
#' @aliases sort,DuckDBGRanges-method
#' @aliases order,DuckDBGRanges-method
#' @aliases is.unsorted,DuckDBGRanges-method
#' @aliases rank,DuckDBGRanges-method
#' @aliases trim,DuckDBGRanges-method
#' @aliases restrict,DuckDBGRanges-method
#' @aliases precede,DuckDBGRanges,DuckDBGRanges-method
#' @aliases precede,DuckDBGRanges,missing-method
#' @aliases precede,DuckDBGRanges,GRanges-method
#' @aliases precede,GRanges,DuckDBGRanges-method
#' @aliases follow,DuckDBGRanges,DuckDBGRanges-method
#' @aliases follow,DuckDBGRanges,missing-method
#' @aliases follow,DuckDBGRanges,GRanges-method
#' @aliases follow,GRanges,DuckDBGRanges-method
#' @aliases nearest,DuckDBGRanges,DuckDBGRanges-method
#' @aliases nearest,DuckDBGRanges,missing-method
#' @aliases nearest,DuckDBGRanges,GRanges-method
#' @aliases nearest,GRanges,DuckDBGRanges-method
#' @aliases distanceToNearest,DuckDBGRanges,DuckDBGRanges-method
#' @aliases distanceToNearest,DuckDBGRanges,missing-method
#' @aliases distanceToNearest,DuckDBGRanges,GRanges-method
#' @aliases distanceToNearest,GRanges,DuckDBGRanges-method
#' @aliases tile,DuckDBGRanges-method
#' @aliases slidingWindows,DuckDBGRanges-method
#' @aliases pgap,DuckDBGRanges,DuckDBGRanges-method
#' @aliases pgap,DuckDBGRanges,GRanges-method
#' @aliases pgap,GRanges,DuckDBGRanges-method
#'
#' @include DuckDBGRanges-class.R
#'
#' @keywords utilities methods
#'
#' @name DuckDBGRanges-utils
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions
###

# Create a temporary DuckDB table from a GRanges query
#' @importFrom dplyr copy_to tbl
#' @importFrom GenomicRanges seqnames
.granges_to_temp_table <- function(query, con, table_name = "query_ranges") {
    query_df <- data.frame(
        query_idx = seq_along(query),
        seqnames = as.character(seqnames(query)),
        query_start = as.integer(start(query)),
        query_end = as.integer(end(query)),
        query_strand = as.character(strand(query)),
        stringsAsFactors = FALSE
    )
    copy_to(con, query_df, name = table_name, temporary = TRUE, overwrite = TRUE)
    tbl(con, table_name)
}

# Build overlap filter conditions (excluding seqnames which is handled by join)
.build_overlap_conditions <-
function(q_start, s_start, q_end, s_end, q_strand, s_strand, maxgap = -1L,
         minoverlap = 0L, ignore.strand = FALSE)
{
    conditions <- list()

    # Overlap condition with maxgap adjustment
    # For maxgap = -1L (default), requires actual overlap
    # For maxgap >= 0, allows gaps up to maxgap bases
    gap_adj <- maxgap + 1L
    conditions$start_cond <- call("<=", call("-", q_start, gap_adj), s_end)
    conditions$end_cond <- call(">=", call("+", q_end, gap_adj), s_start)

    # Strand condition if not ignoring
    if (!ignore.strand) {
        # Overlap if: strand is * on either side, or strands match
        conditions$strand <- call("|",
            call("|",
                call("==", q_strand, "*"),
                call("==", s_strand, "*")),
            call("==", q_strand, s_strand))
    }

    # Minoverlap condition if needed
    if (minoverlap > 1L) {
        # overlap_width = min(q.end, s.end) - max(q.start, s.start) + 1
        overlap_width <- call("+",
            call("-",
                call("least", q_end, s_end),
                call("greatest", q_start, s_start)),
            1L)
        conditions$minoverlap <- call(">=", overlap_width, as.integer(minoverlap))
    }

    conditions
}

# Apply filter conditions to a dplyr connection
#' @importFrom dplyr filter
.apply_overlap_filters <- function(conn, conditions) {
    for (cond in conditions) {
        conn <- filter(conn, !!cond)
    }
    conn
}

# Build the (unmaterialized) overlap-hit join for two DuckDBGRanges.
#
# Split out of findOverlaps(DuckDBGRanges, DuckDBGRanges) so the query plan can
# be inspected with .explainQuery() without collecting it. The two sides are
# joined on seqnames (the equi-key) and the interval bounds applied as
# inequalities: this is an IEJoin / range-join, NOT an ASOF join (ASOF is
# nearest-match, which does not express interval overlap). Returns a lazy dplyr
# tbl of (query_idx, subject_idx).
#' @importFrom dplyr inner_join mutate select arrange
#' @importFrom dbplyr window_order
#' @importFrom DuckDBDataFrame tblconn
.overlap_join_tbl <-
function(query, subject, maxgap = -1L, minoverlap = 0L, ignore.strand = FALSE)
{
    query_select <- setNames(
        lapply(c("query_idx", "seqnames", "start", "end", "strand"), as.name),
        c("query_idx", "seqnames", "query_start", "query_end", "query_strand"))
    query_conn <- tblconn(query@frame)
    query_conn <- mutate(query_conn, query_idx = row_number())
    query_conn <- select(query_conn, !!!query_select)

    subject_select <- setNames(
        lapply(c("subject_idx", "seqnames", "start", "end", "strand"), as.name),
        c("subject_idx", "seqnames", "subject_start", "subject_end", "subject_strand"))
    subject_conn <- tblconn(subject@frame)
    subject_conn <- mutate(subject_conn, subject_idx = row_number())
    subject_conn <- select(subject_conn, !!!subject_select)

    conditions <- .build_overlap_conditions(
        q_start = as.name("query_start"),
        s_start = as.name("subject_start"),
        q_end = as.name("query_end"),
        s_end = as.name("subject_end"),
        q_strand = as.name("query_strand"),
        s_strand = as.name("subject_strand"),
        maxgap = maxgap, minoverlap = minoverlap,
        ignore.strand = ignore.strand)

    result_cols <- lapply(c("query_idx", "subject_idx"), as.name)
    result <- inner_join(query_conn, subject_conn, by = "seqnames")
    result <- .apply_overlap_filters(result, conditions)
    result <- select(result, !!!result_cols)
    arrange(result, !!!result_cols)
}

# Return DuckDB's query plan for a lazy dplyr tbl as a single string.
#
# A developer diagnostic used to confirm pushdown — most importantly that an
# interval-overlap join lowers to an IEJoin / piecewise-merge range join rather
# than a NESTED_LOOP_JOIN over the full cross product (which OOMs on skewed
# inputs). Note the shape is a range join, NOT an ASOF join (nearest-match).
#' @importFrom dbplyr remote_con sql_render
#' @importFrom DBI dbGetQuery
.explainQuery <- function(lazy_tbl) {
    con <- remote_con(lazy_tbl)
    sql <- as.character(sql_render(lazy_tbl))
    plan <- dbGetQuery(con, paste("EXPLAIN", sql))
    paste(unlist(plan, use.names = FALSE), collapse = "\n")
}

# Build a DuckDBGRanges from a dplyr connection with standard GRanges columns
# 
# This helper consolidates the common pattern of creating a new DuckDBGRanges
# from a dplyr connection that already has seqnames, start, end, width, strand.
# It adds the row_number keycol with proper window ordering and creates all
# necessary S4 objects.
#
# @param conn A dplyr connection with columns: seqnames, start, end, width, strand
# @param seqinfo The Seqinfo object to use for the new DuckDBGRanges
# @param order_by Optional list of column symbols for window_order before row_number.
#   Defaults to list(seqnames, strand, start). Pass NULL to skip ordering.
# @return A new DuckDBGRanges object
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importFrom DuckDBDataFrame set_row_number
#' @importFrom dbplyr window_order
#' @importFrom dplyr mutate
#' @importFrom S4Vectors new2
.build_DuckDBGRanges <- function(conn, seqinfo, order_by = NULL) {
    # Set default ordering if not specified
    if (is.null(order_by)) {
        conn <- window_order(conn, seqnames, strand, start)
    } else if (length(order_by) > 0) {
        conn <- window_order(conn, !!!order_by)
    }
    # If order_by is explicitly empty list(), skip ordering
    
    # Add row_number as keycol
    rownum_mutate <- list(row_number = call("row_number"))
    conn <- mutate(conn, !!!rownum_mutate)
    
    # Standard datacols for GRanges
    datacols <- expression(
        seqnames = seqnames,
        start = start,
        end = end,
        width = width,
        strand = strand
    )
    
    keycols <- list(row_number = set_row_number(conn))
    
    new_frame <- new2("DuckDBDataFrame",
                      conn = conn,
                      datacols = datacols,
                      keycols = keycols,
                      dimtbls = new.env(parent = emptyenv()),
                      check = FALSE)
    
    new2("DuckDBGRanges",
         frame = new_frame,
         seqinfo = seqinfo,
         elementMetadata = new("DFrame", nrows = nrow(new_frame)),
         check = FALSE)
}

# Recycle a start/end/shift/width-style per-seqlevel argument, matching a
# named vector up to seqlevels order or recycling an unnamed one. Shared by
# gaps() and coverage().
.recycle_per_seqlevel <- function(val, seqlevels) {
    if (!is.null(names(val)))
        val <- val[seqlevels]
    val <- rep_len(val, length(seqlevels))
    names(val) <- seqlevels
    val
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### findOverlaps methods
###

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors Hits selectHits
#' @importFrom dplyr arrange collect filter inner_join mutate row_number select
setMethod("findOverlaps", c("GRanges", "DuckDBGRanges"),
function(query, subject, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         select = c("all", "first", "last", "arbitrary"),
         ignore.strand = FALSE)
{
    type <- match.arg(type)
    select <- match.arg(select)

    if (type != "any") {
        stop("only type='any' is currently supported for DuckDBGRanges")
    }

    nq <- length(query)
    ns <- length(subject)

    if (nq == 0L || ns == 0L) {
        return(Hits(nLnode = nq, nRnode = ns))
    }

    # Create temp table for query
    con <- dbconn(subject)
    query_tbl <- .granges_to_temp_table(query, con)

    # Get subject connection with row indices
    subject_select <- setNames(
        lapply(c("subject_idx", "seqnames", "start", "end", "strand"), as.name),
        c("subject_idx", "seqnames", "subject_start", "subject_end", "subject_strand"))
    subject_conn <- tblconn(subject@frame)
    subject_conn <- mutate(subject_conn, subject_idx = row_number())
    subject_conn <- select(subject_conn, !!!subject_select)

    # Build overlap conditions (seqnames handled by inner_join)
    conditions <- .build_overlap_conditions(
        q_start = as.name("query_start"),
        s_start = as.name("subject_start"),
        q_end = as.name("query_end"),
        s_end = as.name("subject_end"),
        q_strand = as.name("query_strand"),
        s_strand = as.name("subject_strand"),
        maxgap = maxgap, minoverlap = minoverlap,
        ignore.strand = ignore.strand)

    # Join on seqnames and filter for overlaps
    result_cols <- lapply(c("query_idx", "subject_idx"), as.name)
    result <- inner_join(query_tbl, subject_conn, by = "seqnames")
    result <- .apply_overlap_filters(result, conditions)
    result <- select(result, !!!result_cols)
    result <- arrange(result, !!!result_cols)
    hits_df <- collect(result)

    # Build Hits object
    hits <- Hits(
        from = hits_df$query_idx,
        to = hits_df$subject_idx,
        nLnode = nq,
        nRnode = ns,
        sort.by.query = TRUE
    )

    if (select == "all") {
        hits
    } else {
        selectHits(hits, select = select)
    }
})

#' @export
#' @importFrom dplyr inner_join
#' @importFrom DuckDBDataFrame tblconn
setMethod("findOverlaps", c("DuckDBGRanges", "DuckDBGRanges"),
function(query, subject, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         select = c("all", "first", "last", "arbitrary"),
         ignore.strand = FALSE)
{
    type <- match.arg(type)
    select <- match.arg(select)

    if (type != "any") {
        stop("only type='any' is currently supported for DuckDBGRanges")
    }

    nq <- length(query)
    ns <- length(subject)

    if (nq == 0L || ns == 0L) {
        return(Hits(nLnode = nq, nRnode = ns))
    }

    # Build the overlap-hit range-join (seqnames equi-key + interval inequalities)
    # as a lazy tbl, then materialize. Shared with the .explainQuery plan guard.
    result <- .overlap_join_tbl(query, subject, maxgap = maxgap,
        minoverlap = minoverlap, ignore.strand = ignore.strand)
    hits_df <- collect(result)

    # Build Hits object
    hits <- Hits(
        from = hits_df$query_idx,
        to = hits_df$subject_idx,
        nLnode = nq,
        nRnode = ns,
        sort.by.query = TRUE
    )

    if (select == "all") {
        hits
    } else {
        selectHits(hits, select = select)
    }
})

#' @export
setMethod("findOverlaps", c("DuckDBGRanges", "GRanges"),
function(query, subject, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         select = c("all", "first", "last", "arbitrary"),
         ignore.strand = FALSE)
{
    type <- match.arg(type)
    select <- match.arg(select)

    # Swap query and subject, then transpose the Hits
    hits <- findOverlaps(subject, query, maxgap = maxgap,
        minoverlap = minoverlap, type = type, select = "all",
        ignore.strand = ignore.strand)

    hits <- t(hits)

    if (select == "all") {
        hits
    } else {
        selectHits(hits, select = select)
    }
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### countOverlaps methods
###

#' @export
#' @importFrom IRanges countOverlaps
#' @importFrom dplyr collect filter group_by inner_join mutate n row_number select summarize
#' @importFrom DuckDBDataFrame tblconn
setMethod("countOverlaps", c("GRanges", "DuckDBGRanges"),
function(query, subject, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         ignore.strand = FALSE)
{
    type <- match.arg(type)

    if (type != "any") {
        stop("only type='any' is currently supported for DuckDBGRanges")
    }

    nq <- length(query)
    if (nq == 0L) {
        return(integer(0L))
    }
    if (length(subject) == 0L) {
        return(rep(0L, nq))
    }

    # Create temp table for query
    con <- dbconn(subject)
    query_tbl <- .granges_to_temp_table(query, con)

    # Get subject connection
    subject_select <- setNames(
        lapply(c("seqnames", "start", "end", "strand"), as.name),
        c("seqnames", "subject_start", "subject_end", "subject_strand"))
    subject_conn <- tblconn(subject@frame)
    subject_conn <- select(subject_conn, !!!subject_select)

    # Build overlap conditions (seqnames handled by inner_join)
    conditions <- .build_overlap_conditions(
        q_start = as.name("query_start"),
        s_start = as.name("subject_start"),
        q_end = as.name("query_end"),
        s_end = as.name("subject_end"),
        q_strand = as.name("query_strand"),
        s_strand = as.name("subject_strand"),
        maxgap = maxgap, minoverlap = minoverlap,
        ignore.strand = ignore.strand)

    # Join on seqnames and count by query_idx
    groups <- list(query_idx = as.name("query_idx"))
    result <- inner_join(query_tbl, subject_conn, by = "seqnames")
    result <- .apply_overlap_filters(result, conditions)
    result <- group_by(result, !!!groups)
    result <- summarize(result, count = n(), .groups = "drop")
    counts <- collect(result)

    # Build result vector with zeros for non-overlapping queries
    ans <- rep(0L, nq)
    ans[counts$query_idx] <- as.integer(counts$count)
    ans
})

#' @export
#' @importFrom S4Vectors countQueryHits
setMethod("countOverlaps", c("DuckDBGRanges", "DuckDBGRanges"),
function(query, subject, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         ignore.strand = FALSE)
{
    hits <- findOverlaps(query, subject,
        maxgap = maxgap, minoverlap = minoverlap,
        type = type, select = "all",
        ignore.strand = ignore.strand)
    countQueryHits(hits)
})

#' @export
setMethod("countOverlaps", c("DuckDBGRanges", "GRanges"),
function(query, subject, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         ignore.strand = FALSE)
{
    hits <- findOverlaps(query, subject,
        maxgap = maxgap, minoverlap = minoverlap,
        type = type, select = "all",
        ignore.strand = ignore.strand)
    countQueryHits(hits)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### overlapsAny methods
###

#' @export
#' @importFrom IRanges overlapsAny
setMethod("overlapsAny", c("GRanges", "DuckDBGRanges"),
function(query, subject, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         ...)
{
    type <- match.arg(type)

    hits <- findOverlaps(query, subject,
        maxgap = maxgap, minoverlap = minoverlap,
        type = type, select = "first", ...)

    !is.na(hits)
})

# Helper to get logical vector of overlapping DuckDBGRanges elements
#' @importFrom bit64 as.integer64 is.integer64
#' @importFrom dplyr collect copy_to distinct inner_join pull select tbl
#' @importFrom DuckDBDataFrame has_row_number tblconn
.overlapsAny_DuckDBGRanges <-
function(query, subject, maxgap, minoverlap, type, ignore.strand = FALSE)
{
    type <- match.arg(type, c("any", "start", "end", "within", "equal"))

    if (type != "any") {
        stop("only type='any' is currently supported for DuckDBGRanges")
    }

    nq <- length(query)
    if (nq == 0L) {
        return(logical(0L))
    }

    if (length(subject) == 0L) {
        return(rep(FALSE, nq))
    }

    # Get keycol info from query
    frame <- query@frame
    keycol_name <- names(frame@keycols)
    has_row_num <- has_row_number(frame)

    # Get query connection with keycol and range columns
    query_select <- setNames(
        lapply(c(keycol_name, "seqnames", "start", "end", "strand"), as.name),
        c("query_keycol", "seqnames", "query_start", "query_end", "query_strand"))
    query_conn <- tblconn(frame)
    query_conn <- select(query_conn, !!!query_select)

    # Prepare subject connection/table
    if (is(subject, "DuckDBGRanges")) {
        subject_select <- setNames(
            lapply(c("seqnames", "start", "end", "strand"), as.name),
            c("seqnames", "subject_start", "subject_end", "subject_strand"))
        subject_conn <- tblconn(subject@frame)
        subject_conn <- select(subject_conn, !!!subject_select)
    } else {
        # GRanges: create temp table
        con <- dbconn(query)
        subject_df <- data.frame(
            seqnames = as.character(seqnames(subject)),
            subject_start = as.integer(start(subject)),
            subject_end = as.integer(end(subject)),
            subject_strand = as.character(strand(subject)),
            stringsAsFactors = FALSE
        )
        copy_to(con, subject_df, name = "subject_temp", temporary = TRUE, overwrite = TRUE)
        subject_conn <- tbl(con, "subject_temp")
    }

    # Build overlap conditions (seqnames handled by inner_join)
    conditions <- .build_overlap_conditions(
        q_start = as.name("query_start"),
        s_start = as.name("subject_start"),
        q_end = as.name("query_end"),
        s_end = as.name("subject_end"),
        q_strand = as.name("query_strand"),
        s_strand = as.name("subject_strand"),
        maxgap = maxgap, minoverlap = minoverlap,
        ignore.strand = ignore.strand)

    # Join on seqnames, filter to overlaps, get distinct keycol values
    query_keycol <- as.name("query_keycol")
    result <- inner_join(query_conn, subject_conn, by = "seqnames")
    result <- .apply_overlap_filters(result, conditions)
    result <- distinct(result, !!query_keycol)
    overlap_keys <- pull(collect(result), !!query_keycol)

    # Build logical result vector
    if (has_row_num) {
        # row_number mode: overlap_keys are row numbers (integer64)
        # Convert to integer if safe, otherwise keep as integer64
        if (is.integer64(overlap_keys)) {
            if (all(overlap_keys <= as.integer64(.Machine$integer.max))) {
                overlap_keys <- as.integer(overlap_keys)
            }
        }
        ans <- rep(FALSE, nq)
        ans[overlap_keys] <- TRUE
    } else {
        # Named keycols mode: overlap_keys are the actual key values
        all_keys <- frame@keycols[[1L]]
        ans <- all_keys %in% overlap_keys
    }

    ans
}

#' @export
setMethod("overlapsAny", c("DuckDBGRanges", "DuckDBGRanges"),
function(query, subject, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         ...)
{
    .overlapsAny_DuckDBGRanges(query, subject, maxgap, minoverlap, type, ...)
})

#' @export
setMethod("overlapsAny", c("DuckDBGRanges", "GRanges"),
function(query, subject, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         ...)
{
    .overlapsAny_DuckDBGRanges(query, subject, maxgap, minoverlap, type, ...)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subsetByOverlaps methods for DuckDBGRanges
###
### Note: These are needed because the default subsetByOverlaps uses x[logical]
### which doesn't work with row_number keycols. We use integer indexing instead.
###

#' @export
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors extractROWS
setMethod("subsetByOverlaps", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, ranges, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         invert = FALSE, ...)
{
    ov_any <- overlapsAny(x, ranges, maxgap = maxgap, minoverlap = minoverlap,
                          type = match.arg(type), ...)
    if (invert)
        ov_any <- !ov_any
    extractROWS(x, which(ov_any))
})

#' @export
setMethod("subsetByOverlaps", c("DuckDBGRanges", "GRanges"),
function(x, ranges, maxgap = -1L, minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         invert = FALSE, ...)
{
    ov_any <- overlapsAny(x, ranges, maxgap = maxgap, minoverlap = minoverlap,
                          type = match.arg(type), ...)
    if (invert)
        ov_any <- !ov_any
    extractROWS(x, which(ov_any))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Intra-range methods
###
### These methods modify the datacols expressions in the underlying
### DuckDBDataFrame to transform ranges via SQL operations.
###

# Helper to modify datacols and create a new DuckDBGRanges
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
.modify_DuckDBGRanges_datacols <- function(x, new_start = NULL, new_end = NULL,
                                            new_width = NULL)
{
    frame <- x@frame
    datacols <- frame@datacols

    # Update provided columns - only those explicitly provided
    if (!is.null(new_start))
        datacols[["start"]] <- new_start
    if (!is.null(new_end))
        datacols[["end"]] <- new_end
    if (!is.null(new_width))
        datacols[["width"]] <- new_width

    # Ensure consistency - only recompute if exactly one dimension is missing
    # after the explicit updates
    n_set <- sum(!is.null(new_start), !is.null(new_end), !is.null(new_width))

    if (n_set == 1L) {
        # Only one dimension changed - compute the other two from the original
        # This shouldn't happen in normal usage, but handle it
        if (!is.null(new_start)) {
            # Recompute width from new start and original end
            datacols[["width"]] <- call("-", call("+", 1L, datacols[["end"]]), datacols[["start"]])
        } else if (!is.null(new_end)) {
            # Recompute width from original start and new end
            datacols[["width"]] <- call("-", call("+", 1L, datacols[["end"]]), datacols[["start"]])
        } else {
            # new_width set - recompute end from original start
            datacols[["end"]] <- call("-", call("+", datacols[["start"]], datacols[["width"]]), 1L)
        }
    }
    # If 2 or 3 are set, the caller has already ensured consistency

    new_frame <- replaceSlots(frame, datacols = datacols, check = FALSE)

    mcols <- x@elementMetadata
    if (!is.null(mcols) && is(mcols, "DuckDBDataFrame") && ncol(mcols) > 0L) {
        mcols <- replaceSlots(new_frame, datacols = mcols@datacols, check = FALSE)
    } else {
        mcols <- NULL
    }

    replaceSlots(x, frame = new_frame, elementMetadata = mcols, check = FALSE)
}

#' @export
#' @importFrom IRanges shift
setMethod("shift", "DuckDBGRanges",
function(x, shift = 0L, use.names = TRUE)
{
    if (length(x) == 0L)
        return(x)

    frame <- x@frame
    datacols <- frame@datacols

    # new_start = old_start + shift
    # new_end = old_end + shift
    # width stays the same
    new_start <- call("+", datacols[["start"]], as.integer(shift))
    new_end <- call("+", datacols[["end"]], as.integer(shift))

    .modify_DuckDBGRanges_datacols(x, new_start = new_start, new_end = new_end)
})

#' @export
#' @importFrom IRanges narrow
setMethod("narrow", "DuckDBGRanges",
function(x, start = NA, end = NA, width = NA, use.names = TRUE)
{
    if (length(x) == 0L)
        return(x)

    frame <- x@frame
    datacols <- frame@datacols
    old_start <- datacols[["start"]]
    old_end <- datacols[["end"]]
    old_width <- datacols[["width"]]

    # narrow adjusts positions relative to each range. Base IRanges
    # (solveUserSEW) treats a NEGATIVE start OR end as a position counting back
    # from the range end (-1 is the last base); a positive one counts from the
    # start.
    s <- if (is.na(start)) NA_integer_ else as.integer(start)
    e <- if (is.na(end)) NA_integer_ else as.integer(end)
    w <- if (is.na(width)) NA_integer_ else as.integer(width)

    if (is.na(s) && is.na(e) && is.na(w))
        return(x)

    # New start/end as `column +/- scalar` (NULL = unchanged).
    new_start <- NULL
    new_end <- NULL
    if (!is.na(s))
        new_start <- if (s < 0L) call("+", old_end, s + 1L) else
            call("+", old_start, s - 1L)
    if (!is.na(e))
        new_end <- if (e < 0L) call("+", old_end, e + 1L) else
            call("-", call("+", old_start, e), 1L)

    # A given width fixes the missing side; anchor at the already-resolved start
    # (or end). Adding a scalar to `new_start`, or subtracting a scalar from
    # `new_end`, keeps the expression `column +/- scalar` (safe).
    if (!is.na(w)) {
        if (!is.na(s)) {
            new_end <- call("+", new_start, w - 1L)   # new_start + (w-1)
        } else if (!is.na(e)) {
            new_start <- call("-", new_end, w - 1L)   # new_end - (w-1)
        }
    }

    # width column. `old_width` is itself the expression `end - start + 1`, so it
    # may only be ADDED to (never subtracted, and never subtract any compound).
    # Each case below is a scalar, `old_width + scalar`, or `(old_start -
    # old_end) + scalar` — where old_start/old_end are plain coordinate columns,
    # matching the subtractions the rest of this method already uses.
    if (!is.na(w)) {
        new_width <- w
    } else if (!is.na(s) && !is.na(e)) {
        if (s >= 0L && e >= 0L) {
            new_width <- e - s + 1L
        } else if (s >= 0L && e < 0L) {
            new_width <- call("+", old_width, e - s + 2L)
        } else if (s < 0L && e >= 0L) {
            # width = new_end - new_start + 1
            #       = (old_start + e - 1) - (old_end + s + 1) + 1
            #       = (old_start - old_end) + (e - s - 1)
            new_width <- call("+", call("-", old_start, old_end), e - s - 1L)
        } else {
            new_width <- e - s + 1L
        }
    } else if (!is.na(s)) {
        new_width <- if (s < 0L) -s else call("-", old_width, s - 1L)
    } else {  # end only
        new_width <- if (e < 0L) call("+", old_width, e + 1L) else e
    }

    .modify_DuckDBGRanges_datacols(x, new_start = new_start, new_end = new_end,
                                    new_width = new_width)
})

#' @export
#' @importFrom IRanges resize
setMethod("resize", "DuckDBGRanges",
function(x, width, fix = "start", use.names = TRUE, ignore.strand = FALSE)
{
    if (length(x) == 0L)
        return(x)

    width <- as.integer(width)
    frame <- x@frame
    datacols <- frame@datacols
    old_start <- datacols[["start"]]
    old_end <- datacols[["end"]]
    strand_col <- datacols[["strand"]]

    # For strand-aware resize, we need SQL CASE WHEN
    # fix="start" anchors at start for + strand, end for - strand
    # fix="end" anchors at end for + strand, start for - strand
    # fix="center" anchors at center

    if (fix == "start") {
        if (ignore.strand) {
            new_start <- old_start
            new_end <- call("-", call("+", old_start, width), 1L)
        } else {
            # CASE WHEN strand = '-' THEN end - width + 1 ELSE start END
            new_start <- call("if_else",
                call("==", strand_col, "-"),
                call("+", call("-", old_end, width), 1L),
                old_start)
            # CASE WHEN strand = '-' THEN end ELSE start + width - 1 END
            new_end <- call("if_else",
                call("==", strand_col, "-"),
                old_end,
                call("-", call("+", old_start, width), 1L))
        }
    } else if (fix == "end") {
        if (ignore.strand) {
            new_end <- old_end
            new_start <- call("+", call("-", old_end, width), 1L)
        } else {
            # Opposite of start
            new_start <- call("if_else",
                call("==", strand_col, "-"),
                old_start,
                call("+", call("-", old_end, width), 1L))
            new_end <- call("if_else",
                call("==", strand_col, "-"),
                call("-", call("+", old_start, width), 1L),
                old_end)
        }
    } else if (fix == "center") {
        # new_start = start + (width(x) - width) %/% 2
        # new_end   = new_start + width - 1
        old_width <- call("+", call("-", old_end, old_start), 1L)
        delta <- call("as.integer", call("%/%", call("-", old_width, width), 2L))
        new_start <- call("+", old_start, delta)
        new_end <- call("-", call("+", new_start, width), 1L)
    } else {
        stop("'fix' must be 'start', 'end', or 'center'")
    }

    .modify_DuckDBGRanges_datacols(x, new_start = new_start, new_end = new_end,
                                    new_width = width)
})

#' @export
#' @importFrom IRanges flank
setMethod("flank", "DuckDBGRanges",
function(x, width, start = TRUE, both = FALSE, use.names = TRUE,
         ignore.strand = FALSE)
{
    if (length(x) == 0L)
        return(x)

    width <- as.integer(width)
    frame <- x@frame
    datacols <- frame@datacols
    old_start <- datacols[["start"]]
    old_end <- datacols[["end"]]
    strand_col <- datacols[["strand"]]

    if (both) {
        # Flanking region on both sides of specified end
        # For both=TRUE, center at the specified position and extend width on each side
        # Total width is 2*width
        if (start) {
            if (ignore.strand) {
                # Center at old_start, extend width on each side
                new_start <- call("-", old_start, width)
                new_end <- call("-", call("+", old_start, width), 1L)
            } else {
                # + strand: center at old_start; - strand: center at old_end + 1
                # For - strand: start = (end+1) - width, end = (end+1) + width - 1
                new_start <- call("if_else",
                    call("==", strand_col, "-"),
                    call("+", call("-", old_end, width), 1L),
                    call("-", old_start, width))
                new_end <- call("if_else",
                    call("==", strand_col, "-"),
                    call("+", old_end, width),
                    call("-", call("+", old_start, width), 1L))
            }
        } else {
            if (ignore.strand) {
                # Center at old_end + 1, extend width on each side
                new_start <- call("+", call("-", old_end, width), 1L)
                new_end <- call("+", old_end, width)
            } else {
                # + strand: center at old_end + 1; - strand: center at old_start
                new_start <- call("if_else",
                    call("==", strand_col, "-"),
                    call("-", old_start, width),
                    call("+", call("-", old_end, width), 1L))
                new_end <- call("if_else",
                    call("==", strand_col, "-"),
                    call("-", call("+", old_start, width), 1L),
                    call("+", old_end, width))
            }
        }
        new_width <- call("*", width, 2L)
    } else {
        if (start) {
            if (ignore.strand) {
                new_start <- call("-", old_start, width)
                new_end <- call("-", old_start, 1L)
            } else {
                # + strand: region before start; - strand: region after end
                new_start <- call("if_else",
                    call("==", strand_col, "-"),
                    call("+", old_end, 1L),
                    call("-", old_start, width))
                new_end <- call("if_else",
                    call("==", strand_col, "-"),
                    call("+", old_end, width),
                    call("-", old_start, 1L))
            }
        } else {
            if (ignore.strand) {
                new_start <- call("+", old_end, 1L)
                new_end <- call("+", old_end, width)
            } else {
                # + strand: region after end; - strand: region before start
                new_start <- call("if_else",
                    call("==", strand_col, "-"),
                    call("-", old_start, width),
                    call("+", old_end, 1L))
                new_end <- call("if_else",
                    call("==", strand_col, "-"),
                    call("-", old_start, 1L),
                    call("+", old_end, width))
            }
        }
        new_width <- width
    }

    .modify_DuckDBGRanges_datacols(x, new_start = new_start, new_end = new_end,
                                    new_width = new_width)
})

#' @export
#' @importFrom IRanges promoters
setMethod("promoters", "DuckDBGRanges",
function(x, upstream = 2000, downstream = 200, use.names = TRUE)
{
    if (length(x) == 0L)
        return(x)

    upstream <- as.integer(upstream)
    downstream <- as.integer(downstream)
    frame <- x@frame
    datacols <- frame@datacols
    old_start <- datacols[["start"]]
    old_end <- datacols[["end"]]
    strand_col <- datacols[["strand"]]

    # Promoter: upstream of TSS (transcription start site)
    # + strand: TSS = start; promoter = [start - upstream, start + downstream - 1]
    # - strand: TSS = end; promoter = [end - downstream + 1, end + upstream]
    new_start <- call("if_else",
        call("==", strand_col, "-"),
        call("+", call("-", old_end, downstream), 1L),
        call("-", old_start, upstream))
    new_end <- call("if_else",
        call("==", strand_col, "-"),
        call("+", old_end, upstream),
        call("-", call("+", old_start, downstream), 1L))

    new_width <- upstream + downstream

    .modify_DuckDBGRanges_datacols(x, new_start = new_start, new_end = new_end,
                                    new_width = new_width)
})

#' @export
#' @importFrom GenomicRanges terminators
setMethod("terminators", "DuckDBGRanges",
function(x, upstream = 2000, downstream = 200, use.names = TRUE)
{
    if (length(x) == 0L)
        return(x)

    upstream <- as.integer(upstream)
    downstream <- as.integer(downstream)
    frame <- x@frame
    datacols <- frame@datacols
    old_start <- datacols[["start"]]
    old_end <- datacols[["end"]]
    strand_col <- datacols[["strand"]]

    # Terminator: upstream of TES (transcription end site)
    # + strand: TES = end; terminator = [end - upstream, end + downstream - 1]
    # - strand: TES = start; terminator = [start - downstream + 1, start + upstream]
    new_start <- call("if_else",
        call("==", strand_col, "-"),
        call("+", call("-", old_start, downstream), 1L),
        call("-", old_end, upstream))
    new_end <- call("if_else",
        call("==", strand_col, "-"),
        call("+", old_start, upstream),
        call("-", call("+", old_end, downstream), 1L))

    new_width <- upstream + downstream

    .modify_DuckDBGRanges_datacols(x, new_start = new_start, new_end = new_end,
                                    new_width = new_width)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Inter-range methods
###

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importFrom DuckDBDataFrame set_row_number tblconn
#' @importFrom dplyr distinct group_by mutate row_number summarize
#' @importFrom Seqinfo seqinfo
#' @importFrom S4Vectors new2
setMethod("range", "DuckDBGRanges",
function(x, ..., with.revmap = FALSE, ignore.strand = FALSE, na.rm = FALSE)
{
    if (!identical(na.rm, FALSE))
        warning("'na.rm' argument is ignored")

    if (with.revmap)
        warning("'with.revmap' is not supported for DuckDBGRanges, ignoring")

    args <- list(x, ...)
    if (length(args) > 1L) {
        stop("multiple arguments not yet supported for DuckDBGRanges range()")
    }

    if (length(x) == 0L)
        return(x)

    # Create a new DuckDBGRanges with aggregated ranges via SQL
    # We need to create a new connection with GROUP BY + MIN/MAX
    frame <- x@frame
    conn <- tblconn(frame)

    if (ignore.strand) {
        # Group by seqnames only, set strand to "*"
        groups <- list(seqnames = as.name("seqnames"))
        aggr <- list(
            start = call("min", as.name("start"), na.rm = TRUE),
            end = call("max", as.name("end"), na.rm = TRUE)
        )
        new_conn <- conn |>
            group_by(!!!groups) |>
            summarize(!!!aggr, .groups = "drop") |>
            mutate(strand = "*",
                   width = !!call("+", call("-", as.name("end"), as.name("start")), 1L))
    } else {
        # Group by seqnames and strand
        groups <- list(seqnames = as.name("seqnames"), strand = as.name("strand"))
        aggr <- list(
            start = call("min", as.name("start"), na.rm = TRUE),
            end = call("max", as.name("end"), na.rm = TRUE)
        )
        new_conn <- conn |>
            group_by(!!!groups) |>
            summarize(!!!aggr, .groups = "drop") |>
            mutate(width = !!call("+", call("-", as.name("end"), as.name("start")), 1L))
    }

    # Add row_number as keycol
    new_conn <- mutate(new_conn, row_number = row_number())

    # Create new datacols expression
    datacols <- expression(
        seqnames = seqnames,
        start = start,
        end = end,
        width = width,
        strand = strand
    )

    # Create keycols
    keycols <- list(row_number = set_row_number(new_conn))

    # Create new DuckDBDataFrame
    new_frame <- new2("DuckDBDataFrame",
                      conn = new_conn,
                      datacols = datacols,
                      keycols = keycols,
                      dimtbls = new.env(parent = emptyenv()),
                      check = FALSE)

    # Create new DuckDBGRanges with proper empty DataFrame for elementMetadata
    new2("DuckDBGRanges",
         frame = new_frame,
         seqinfo = seqinfo(x),
         elementMetadata = new("DFrame", nrows = nrow(new_frame)),
         check = FALSE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparison methods
###

#' @importClassesFrom DuckDBDataFrame DuckDBColumn
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom GenomicRanges duplicated
#' @importFrom dbplyr window_order
#' @importFrom dplyr row_number mutate select group_by ungroup arrange desc inner_join copy_to
#' @importFrom S4Vectors new2
.duplicated.DuckDBGRanges <- function(x, incomparables = FALSE, fromLast = FALSE,
                                       nmax = NA, method = c("auto", "quick", "hash"))
{
    if (!identical(incomparables, FALSE))
        stop("\"duplicated\" method for DuckDBGRanges objects ",
             "only accepts 'incomparables=FALSE'")

    frame <- x@frame
    n <- length(x)

    if (n == 0L) {
        # Return empty DuckDBColumn
        empty_mutate <- list(duplicated = FALSE)
        empty_select <- lapply(c(names(frame@keycols), "duplicated"), as.name)
        empty_conn <- mutate(tblconn(frame), !!!empty_mutate)
        empty_conn <- select(empty_conn, !!!empty_select)
        datacols <- as.expression(setNames(list(as.name("duplicated")), "duplicated"))
        table <- new2("DuckDBTable", conn = empty_conn, datacols = datacols,
                      keycols = frame@keycols, check = FALSE)
        return(new2("DuckDBColumn", table = table, check = FALSE))
    }

    # Use SQL WINDOW function to identify duplicates in-database
    conn <- tblconn(frame)
    keycol_name <- names(frame@keycols)
    all_keys <- frame@keycols[[1L]]

    # Create position mapping: keycol value -> position in keycols
    # Position determines row order (position 1 = row 1, etc.)
    key_positions <- seq_along(all_keys)
    names(key_positions) <- all_keys

    # Add row position to the connection for proper ordering
    position_df <- setNames(
        data.frame(all_keys, key_positions, stringsAsFactors = FALSE),
        c(keycol_name, "row_pos"))

    # Upload position lookup to DuckDB
    position_tbl <- copy_to(dbconn(frame), position_df,
                            name = paste0("pos_lookup_", sample.int(1e6, 1)),
                            temporary = TRUE, overwrite = TRUE)

    # Join to get positions, then compute duplicates
    conn_with_pos <- inner_join(conn, position_tbl, by = keycol_name)

    # Column symbols for group_by
    group_cols <- lapply(c("seqnames", "start", "end", "strand"), as.name)

    # Build expressions using call() and setNames() to avoid R CMD check NOTEs
    dup_rank_expr <- setNames(list(call("row_number")), "dup_rank")
    dup_cond_expr <- setNames(list(call(">", as.name("dup_rank"), 1L)), "duplicated")
    select_cols <- lapply(c(keycol_name, "duplicated"), as.name)

    # Use group_by + window_order + row_number() for window function in dbplyr
    # window_order() sets ordering for window functions without causing
    # "ORDER BY ignored in subqueries" warnings
    # Rows with dup_rank > 1 are duplicates
    if (fromLast) {
        # For fromLast, we want the last occurrence to be first (dup_rank = 1)
        result_conn <- group_by(conn_with_pos, !!!group_cols)
        result_conn <- window_order(result_conn, desc(!!as.name("row_pos")))
        result_conn <- mutate(result_conn, !!!dup_rank_expr)
        result_conn <- ungroup(result_conn)
        result_conn <- mutate(result_conn, !!!dup_cond_expr)
        result_conn <- select(result_conn, !!!select_cols)
    } else {
        result_conn <- group_by(conn_with_pos, !!!group_cols)
        result_conn <- window_order(result_conn, !!as.name("row_pos"))
        result_conn <- mutate(result_conn, !!!dup_rank_expr)
        result_conn <- ungroup(result_conn)
        result_conn <- mutate(result_conn, !!!dup_cond_expr)
        result_conn <- select(result_conn, !!!select_cols)
    }

    # Create DuckDBColumn with compatible keycols
    datacols <- as.expression(setNames(list(as.name("duplicated")), "duplicated"))
    table <- new2("DuckDBTable", conn = result_conn, datacols = datacols,
                  keycols = frame@keycols, check = FALSE)
    new2("DuckDBColumn", table = table, check = FALSE)
}

#' @export
setMethod("duplicated", "DuckDBGRanges", .duplicated.DuckDBGRanges)

#' @export
setMethod("unique", "DuckDBGRanges",
function(x, incomparables = FALSE, fromLast = FALSE, ...)
{
    if (!identical(incomparables, FALSE))
        stop("\"unique\" method for DuckDBGRanges objects ",
             "only accepts 'incomparables=FALSE'")

    if (length(x) == 0L)
        return(x)

    # Find duplicates and remove them
    # Convert DuckDBColumn to vector for subsetting
    dups <- as.vector(duplicated(x, fromLast = fromLast))
    x[!dups]
})

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBColumn
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom GenomicRanges match
#' @importFrom S4Vectors new2
#' @importFrom dplyr left_join coalesce summarize all_of rename
setMethod("match", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, table, nomatch = NA_integer_, incomparables = NULL,
         method = c("auto", "quick", "hash"), ignore.strand = FALSE)
{
    if (!is.null(incomparables))
        stop("\"match\" method for DuckDBGRanges objects ",
             "only accepts 'incomparables=NULL'")

    frame_x <- x@frame
    frame_table <- table@frame
    n_x <- length(x)
    n_table <- length(table)
    keycol_name_x <- names(frame_x@keycols)

    # Handle empty cases - return DuckDBColumn with nomatch values
    if (n_x == 0L) {
        empty_mutate <- list(match_result = NA_integer_)
        empty_select <- lapply(c(keycol_name_x, "match_result"), as.name)
        empty_conn <- mutate(tblconn(frame_x), !!!empty_mutate)
        empty_conn <- select(empty_conn, !!!empty_select)
        datacols <- as.expression(setNames(list(as.name("match_result")), "match_result"))
        table_obj <- new2("DuckDBTable", conn = empty_conn, datacols = datacols,
                          keycols = frame_x@keycols, check = FALSE)
        return(new2("DuckDBColumn", table = table_obj, check = FALSE))
    }

    if (n_table == 0L) {
        # All rows get nomatch value
        nomatch_mutate <- list(match_result = nomatch)
        nomatch_select <- lapply(c(keycol_name_x, "match_result"), as.name)
        result_conn <- mutate(tblconn(frame_x), !!!nomatch_mutate)
        result_conn <- select(result_conn, !!!nomatch_select)
        datacols <- as.expression(setNames(list(as.name("match_result")), "match_result"))
        table_obj <- new2("DuckDBTable", conn = result_conn, datacols = datacols,
                          keycols = frame_x@keycols, check = FALSE)
        return(new2("DuckDBColumn", table = table_obj, check = FALSE))
    }

    # Get frame info
    keycol_name_table <- names(frame_table@keycols)
    all_keys_table <- frame_table@keycols[[1L]]

    # Create position mapping for table: keycol value -> position
    table_positions <- seq_along(all_keys_table)
    names(table_positions) <- all_keys_table
    position_df <- setNames(
        data.frame(all_keys_table, table_positions, stringsAsFactors = FALSE),
        c(keycol_name_table, "table_pos"))

    # Upload position lookup to DuckDB
    position_tbl <- copy_to(dbconn(frame_table), position_df,
                            name = paste0("match_pos_", sample.int(1e6, 1)),
                            temporary = TRUE, overwrite = TRUE)

    # Get connections
    conn_x <- tblconn(frame_x)
    conn_table <- tblconn(frame_table, select = TRUE, filter = TRUE)

    # Add positions to table
    conn_table_with_pos <- inner_join(conn_table, position_tbl, by = keycol_name_table)

    # Build join conditions
    if (ignore.strand) {
        join_cols <- c("seqnames", "start", "end")
    } else {
        join_cols <- c("seqnames", "start", "end", "strand")
    }

    # Rename table columns for join to avoid conflicts using setNames pattern
    table_select_cols <- lapply(c(join_cols, "table_pos"), as.name)
    table_renamed <- select(conn_table_with_pos, all_of(c(join_cols, "table_pos")))
    rename_list <- setNames(
        lapply(c("seqnames", "start", "end"), as.name),
        c("t_seqnames", "t_start", "t_end"))
    table_renamed <- rename(table_renamed, !!!rename_list)
    if (!ignore.strand) {
        strand_rename <- setNames(list(as.name("strand")), "t_strand")
        table_renamed <- rename(table_renamed, !!!strand_rename)
    }

    # Prepare x for join
    x_select_cols <- lapply(c(keycol_name_x, "seqnames", "start", "end", "strand"), as.name)
    x_prepared <- select(conn_x, !!!x_select_cols)

    # Perform left join to find matches
    # For each x row, we want the minimum table_pos (first match)
    if (ignore.strand) {
        joined <- left_join(x_prepared, table_renamed,
            by = c("seqnames" = "t_seqnames", "start" = "t_start", "end" = "t_end"))
    } else {
        joined <- left_join(x_prepared, table_renamed,
            by = c("seqnames" = "t_seqnames", "start" = "t_start",
                   "end" = "t_end", "strand" = "t_strand"))
    }

    # Get first match (minimum table_pos) for each x row
    # Use COALESCE to replace NULL with nomatch value
    group_cols <- list(as.name(keycol_name_x))
    summ_expr <- list(match_pos = call("min", as.name("table_pos"), na.rm = FALSE))
    coalesce_expr <- list(match_result = call("coalesce",
        call("as.integer", as.name("match_pos")), nomatch))
    final_select <- lapply(c(keycol_name_x, "match_result"), as.name)

    result_conn <- group_by(joined, !!!group_cols)
    result_conn <- summarize(result_conn, !!!summ_expr, .groups = "drop")
    result_conn <- mutate(result_conn, !!!coalesce_expr)
    result_conn <- select(result_conn, !!!final_select)

    # Create DuckDBColumn with compatible keycols
    datacols <- as.expression(setNames(list(as.name("match_result")), "match_result"))
    table_obj <- new2("DuckDBTable", conn = result_conn, datacols = datacols,
                      keycols = frame_x@keycols, check = FALSE)
    new2("DuckDBColumn", table = table_obj, check = FALSE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Parallel set operations
###
### These require both objects to have the same length and operate element-wise.
### We join the two tables on row position and compute new coordinates.
###

# Helper to convert a GRanges to a DuckDBGRanges using the same connection
# as an existing DuckDBGRanges object. This enables DuckDB-optimized operations
# when one argument is GRanges and the other is DuckDBGRanges.
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importFrom GenomicRanges seqnames start end strand seqinfo
#' @importFrom dplyr copy_to
.granges_to_duckdb <- function(gr, ddb_gr, preserve_names = FALSE)
{
    # Get the DuckDB connection from the existing DuckDBGRanges
    conn <- dbconn(ddb_gr@frame)

    # Convert GRanges to data.frame
    n <- length(gr)

    # Determine keycol name and values
    if (preserve_names && !is.null(names(gr))) {
        keycol_name <- "name"
        keycol_vals <- names(gr)
        gr_df <- data.frame(
            name = keycol_vals,
            seqnames = as.character(seqnames(gr)),
            start = as.integer(start(gr)),
            end = as.integer(end(gr)),
            width = as.integer(end(gr) - start(gr) + 1L),
            strand = as.character(strand(gr)),
            stringsAsFactors = FALSE
        )
    } else {
        keycol_name <- "row_number"
        keycol_vals <- seq_len(n)
        gr_df <- data.frame(
            row_number = keycol_vals,
            seqnames = as.character(seqnames(gr)),
            start = as.integer(start(gr)),
            end = as.integer(end(gr)),
            width = as.integer(end(gr) - start(gr) + 1L),
            strand = as.character(strand(gr)),
            stringsAsFactors = FALSE
        )
    }

    # Upload to DuckDB as temporary table
    table_name <- paste0("gr_temp_", sample.int(1e6, 1))
    ddb_tbl <- copy_to(conn, gr_df, name = table_name, temporary = TRUE, overwrite = TRUE)

    # Create datacols expression
    datacols <- expression(
        seqnames = seqnames,
        start = start,
        end = end,
        width = width,
        strand = strand
    )

    # Create keycols
    keycols <- setNames(list(keycol_vals), keycol_name)

    # Create the DuckDBDataFrame
    new_frame <- new2("DuckDBDataFrame",
                      conn = ddb_tbl,
                      datacols = datacols,
                      keycols = keycols,
                      dimtbls = new.env(parent = emptyenv()),
                      check = FALSE)

    # Create and return DuckDBGRanges
    new2("DuckDBGRanges",
         frame = new_frame,
         seqinfo = seqinfo(gr),
         elementMetadata = new("DFrame", nrows = n),
         check = FALSE)
}

# Helper to create a DuckDBGRanges from joined x and y connections
#' @importFrom DuckDBDataFrame tblconn
.parallel_set_op <- function(x, y, ignore.strand, new_start_expr, new_end_expr)
{
    if (length(x) != length(y))
        stop("'x' and 'y' must have the same length")

    if (length(x) == 0L)
        return(x)

    # Get connections with row indices
    row_idx_mutate <- setNames(list(call("row_number")), ".row_idx")
    x_conn <- tblconn(x@frame)
    x_conn <- mutate(x_conn, !!!row_idx_mutate)

    y_conn <- tblconn(y@frame)
    y_conn <- mutate(y_conn, !!!row_idx_mutate)
    y_select_list <- setNames(
        lapply(c(".row_idx", "seqnames", "start", "end", "strand"), as.name),
        c(".row_idx", "y_seqnames", "y_start", "y_end", "y_strand"))
    y_conn <- select(y_conn, !!!y_select_list)

    # Join on row index
    joined <- inner_join(x_conn, y_conn, by = ".row_idx", copy = TRUE)

    # Apply the new start/end expressions
    # new_width = new_end - new_start + 1
    new_width_expr <- call("+", call("-", as.name("new_end"), as.name("new_start")), 1L)
    coord_mutate <- list(new_start = new_start_expr, new_end = new_end_expr,
                         new_width = new_width_expr)
    joined <- mutate(joined, !!!coord_mutate)

    # Select and rename to standard columns
    select_rename_list <- setNames(
        lapply(c("seqnames", "strand", "new_start", "new_end", "new_width"), as.name),
        c("seqnames", "strand", "start", "end", "width"))
    joined <- select(joined, !!!select_rename_list)

    .build_DuckDBGRanges(joined, seqinfo(x))
}

#' @export
#' @importFrom IRanges punion
#' @importFrom dplyr inner_join mutate row_number select
#' @importFrom S4Vectors new2
setMethod("punion", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, y, fill.gap = FALSE, ignore.strand = FALSE)
{
    if (!fill.gap) {
        # Without fill.gap, ranges must overlap or be adjacent
        # For simplicity, we still compute the union but don't validate overlap
        # (the standard method also doesn't strictly enforce this for all cases)
    }

    # punion: new_start = min(x_start, y_start), new_end = max(x_end, y_end)
    new_start_expr <- call("least", as.name("start"), as.name("y_start"))
    new_end_expr <- call("greatest", as.name("end"), as.name("y_end"))

    .parallel_set_op(x, y, ignore.strand, new_start_expr, new_end_expr)
})

#' @export
setMethod("punion", c("DuckDBGRanges", "GRanges"),
function(x, y, fill.gap = FALSE, ignore.strand = FALSE)
{
    # Convert GRanges to DuckDBGRanges and use optimized method
    y_ddb <- .granges_to_duckdb(y, x)
    punion(x, y_ddb, fill.gap = fill.gap, ignore.strand = ignore.strand)
})

#' @export
setMethod("punion", c("GRanges", "DuckDBGRanges"),
function(x, y, fill.gap = FALSE, ignore.strand = FALSE)
{
    # Convert GRanges to DuckDBGRanges and use optimized method
    x_ddb <- .granges_to_duckdb(x, y)
    punion(x_ddb, y, fill.gap = fill.gap, ignore.strand = ignore.strand)
})

#' @export
#' @importFrom IRanges pintersect
#' @importFrom dplyr inner_join mutate row_number select
setMethod("pintersect", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, y, drop.nohit.ranges = FALSE, ignore.strand = FALSE,
         strict.strand = FALSE)
{
    # pintersect: new_start = max(x_start, y_start), new_end = min(x_end, y_end)
    new_start_expr <- call("greatest", as.name("start"), as.name("y_start"))
    new_end_expr <- call("least", as.name("end"), as.name("y_end"))

    .parallel_set_op(x, y, ignore.strand, new_start_expr, new_end_expr)
})

#' @export
setMethod("pintersect", c("DuckDBGRanges", "GRanges"),
function(x, y, drop.nohit.ranges = FALSE, ignore.strand = FALSE,
         strict.strand = FALSE)
{
    # Convert GRanges to DuckDBGRanges and use optimized method
    y_ddb <- .granges_to_duckdb(y, x)
    pintersect(x, y_ddb, drop.nohit.ranges = drop.nohit.ranges,
               ignore.strand = ignore.strand, strict.strand = strict.strand)
})

#' @export
setMethod("pintersect", c("GRanges", "DuckDBGRanges"),
function(x, y, drop.nohit.ranges = FALSE, ignore.strand = FALSE,
         strict.strand = FALSE)
{
    # Convert GRanges to DuckDBGRanges and use optimized method
    x_ddb <- .granges_to_duckdb(x, y)
    pintersect(x_ddb, y, drop.nohit.ranges = drop.nohit.ranges,
               ignore.strand = ignore.strand, strict.strand = strict.strand)
})

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom IRanges psetdiff
#' @importFrom dplyr case_when inner_join mutate row_number select
setMethod("psetdiff", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE)
{
    # psetdiff is complex: x minus y
    # Result depends on relationship between x and y
    # For ranges where y overlaps x:
    # - If y starts at or before x start: result is y_end+1 to x_end
    # - If y ends at or after x end: result is x_start to y_start-1
    # This implementation assumes y is fully contained overlap scenario
    # or edge-aligned (which is the common use case)

    if (length(x) != length(y))
        stop("'x' and 'y' must have the same length")

    if (length(x) == 0L)
        return(x)

    row_idx_mutate <- setNames(list(call("row_number")), ".row_idx")
    x_conn <- tblconn(x@frame)
    x_conn <- mutate(x_conn, !!!row_idx_mutate)

    y_conn <- tblconn(y@frame)
    y_conn <- mutate(y_conn, !!!row_idx_mutate)
    y_select_list <- setNames(
        lapply(c(".row_idx", "seqnames", "start", "end", "strand"), as.name),
        c(".row_idx", "y_seqnames", "y_start", "y_end", "y_strand"))
    y_conn <- select(y_conn, !!!y_select_list)

    joined <- inner_join(x_conn, y_conn, by = ".row_idx", copy = TRUE)

    # Use CASE WHEN logic for psetdiff
    # If y_start <= start: new_start = y_end + 1, new_end = end
    # Else: new_start = start, new_end = y_start - 1
    # Build case_when expression using call()
    cond1 <- call("<=", as.name("y_start"), as.name("start"))
    case_start_expr <- call("case_when",
        call("~", cond1, call("+", as.name("y_end"), 1L)),
        call("~", TRUE, as.name("start")))
    case_end_expr <- call("case_when",
        call("~", cond1, as.name("end")),
        call("~", TRUE, call("-", as.name("y_start"), 1L)))
    new_width_expr <- call("+", call("-", as.name("new_end"), as.name("new_start")), 1L)

    coord_mutate <- list(new_start = case_start_expr, new_end = case_end_expr,
                         new_width = new_width_expr)
    joined <- mutate(joined, !!!coord_mutate)

    select_rename_list <- setNames(
        lapply(c("seqnames", "strand", "new_start", "new_end", "new_width"), as.name),
        c("seqnames", "strand", "start", "end", "width"))
    joined <- select(joined, !!!select_rename_list)

    .build_DuckDBGRanges(joined, seqinfo(x))
})

#' @export
setMethod("psetdiff", c("DuckDBGRanges", "GRanges"),
function(x, y, ignore.strand = FALSE)
{
    # Convert GRanges to DuckDBGRanges and use optimized method
    y_ddb <- .granges_to_duckdb(y, x)
    psetdiff(x, y_ddb, ignore.strand = ignore.strand)
})

#' @export
setMethod("psetdiff", c("GRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE)
{
    # Convert GRanges to DuckDBGRanges and use optimized method
    x_ddb <- .granges_to_duckdb(x, y)
    psetdiff(x_ddb, y, ignore.strand = ignore.strand)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Distance methods
###
### distance() returns an integer vector, so we compute via SQL and collect.
###

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom IRanges distance
#' @importFrom dplyr arrange collect if_else inner_join mutate row_number select
setMethod("distance", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE, ...)
{
    if (length(x) != length(y))
        stop("'x' and 'y' must have the same length")

    n <- length(x)
    if (n == 0L)
        return(integer(0L))

    # SQL-based distance computation. Carry seqnames + strand through the join so
    # incompatible pairs can be set to NA (base GenomicRanges::distance returns NA
    # for a pair on different seqnames, or — unless ignore.strand — on '+' vs '-';
    # '*' matches any strand).
    row_idx_mutate <- setNames(list(call("row_number")), ".row_idx")
    x_conn <- tblconn(x@frame)
    x_conn <- mutate(x_conn, !!!row_idx_mutate)
    x_select_list <- setNames(
        lapply(c(".row_idx", "seqnames", "strand", "start", "end"), as.name),
        c(".row_idx", "x_seqnames", "x_strand", "x_start", "x_end"))
    x_conn <- select(x_conn, !!!x_select_list)

    y_conn <- tblconn(y@frame)
    y_conn <- mutate(y_conn, !!!row_idx_mutate)
    y_select_list <- setNames(
        lapply(c(".row_idx", "seqnames", "strand", "start", "end"), as.name),
        c(".row_idx", "y_seqnames", "y_strand", "y_start", "y_end"))
    y_conn <- select(y_conn, !!!y_select_list)

    joined <- inner_join(x_conn, y_conn, by = ".row_idx", copy = TRUE)

    # Distance is 0 if overlapping, otherwise gap between ranges
    # distance = max(0, max(x_start, y_start) - min(x_end, y_end) - 1)
    # Use call() for greatest/least which are SQL functions translated by dbplyr
    greatest_call <- call("greatest", as.name("x_start"), as.name("y_start"))
    least_call <- call("least", as.name("x_end"), as.name("y_end"))
    gap_expr <- call("greatest", 0L,
                     call("-", call("-", greatest_call, least_call), 1L))

    # Valid pair: same seqname, and (unless ignore.strand) compatible strand.
    same_seq <- call("==", as.name("x_seqnames"), as.name("y_seqnames"))
    if (ignore.strand) {
        valid_expr <- same_seq
    } else {
        strand_ok <- call("|",
                          call("|",
                               call("==", as.name("x_strand"), as.name("y_strand")),
                               call("==", as.name("x_strand"), "*")),
                          call("==", as.name("y_strand"), "*"))
        valid_expr <- call("&", same_seq, strand_ok)
    }
    dist_expr <- call("if_else", valid_expr, gap_expr, NA_integer_)
    dist_mutate <- list(dist = dist_expr)
    joined <- mutate(joined, !!!dist_mutate)

    arrange_cols <- list(as.name(".row_idx"))
    select_cols <- list(as.name("dist"))
    result <- arrange(joined, !!!arrange_cols)
    result <- select(result, !!!select_cols)
    result <- collect(result)

    as.integer(result$dist)
})

#' @export
setMethod("distance", c("DuckDBGRanges", "GRanges"),
function(x, y, ignore.strand = FALSE, ...)
{
    # Convert GRanges to DuckDBGRanges and use optimized method
    y_ddb <- .granges_to_duckdb(y, x)
    distance(x, y_ddb, ignore.strand = ignore.strand, ...)
})

#' @export
setMethod("distance", c("GRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE, ...)
{
    # Convert GRanges to DuckDBGRanges and use optimized method
    x_ddb <- .granges_to_duckdb(x, y)
    distance(x_ddb, y, ignore.strand = ignore.strand, ...)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ordering methods
###
### sort(), order(), is.unsorted(), and rank() use SQL ORDER BY and window
### functions for efficient in-database computation.
###

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dplyr arrange collect mutate row_number select
#' @importFrom S4Vectors new2
setMethod("order", "DuckDBGRanges",
function(..., na.last = TRUE, decreasing = FALSE,
         method = c("auto", "shell", "radix"))
{
    if (!isTRUEorFALSE(decreasing))
        stop("'decreasing' must be TRUE or FALSE")

    args <- list(...)
    if (length(args) != 1L)
        stop("only single argument ordering is currently supported for DuckDBGRanges")

    x <- args[[1L]]
    n <- length(x)
    if (n == 0L)
        return(integer(0L))

    frame <- x@frame

    # Get data in natural order and assign row positions
    conn <- tblconn(frame)
    pos_mutate <- list(orig_pos = call("row_number"))
    conn_with_pos <- mutate(conn, !!!pos_mutate)

    # Order by seqnames, strand, start, end (natural genomic order)
    # Note: ordering by end is equivalent to ordering by width when start is fixed
    if (decreasing) {
        arrange_cols <- list(
            call("desc", as.name("seqnames")),
            call("desc", as.name("strand")),
            call("desc", as.name("start")),
            call("desc", as.name("end")))
    } else {
        arrange_cols <- list(
            as.name("seqnames"), as.name("strand"),
            as.name("start"), as.name("end"))
    }
    result_conn <- arrange(conn_with_pos, !!!arrange_cols)

    select_cols <- list(as.name("orig_pos"))
    result_conn <- select(result_conn, !!!select_cols)
    result <- collect(result_conn)

    as.integer(result$orig_pos)
})

#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dbplyr window_order
#' @importFrom dplyr arrange mutate row_number select
#' @importFrom S4Vectors new2
.sort.DuckDBGRanges <- function(x, decreasing = FALSE, ignore.strand = FALSE, by)
{
    if (!isTRUEorFALSE(decreasing))
        stop("'decreasing' must be TRUE or FALSE")
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    if (!missing(by))
        stop("'by' argument is not yet supported for DuckDBGRanges")

    n <- length(x)
    if (n == 0L)
        return(x)

    frame <- x@frame
    conn <- tblconn(frame)

    # Use window_order() to set ordering for row_number() window function
    # This avoids "ORDER BY ignored in subqueries" warnings
    # Note: ordering by end is equivalent to ordering by width when start is fixed
    if (ignore.strand) {
        if (decreasing) {
            sorted_conn <- window_order(conn, desc(seqnames), desc(start), desc(end))
        } else {
            sorted_conn <- window_order(conn, seqnames, start, end)
        }
    } else {
        if (decreasing) {
            sorted_conn <- window_order(conn, desc(seqnames), desc(strand), desc(start), desc(end))
        } else {
            sorted_conn <- window_order(conn, seqnames, strand, start, end)
        }
    }

    # Window ordering already set above, use empty list to skip helper's default
    .build_DuckDBGRanges(sorted_conn, seqinfo(x), order_by = list())
}

#' @export
setMethod("sort", "DuckDBGRanges", .sort.DuckDBGRanges)

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dplyr collect filter mutate lag summarize
setMethod("is.unsorted", "DuckDBGRanges",
function(x, na.rm = FALSE, strictly = FALSE, ignore.strand = FALSE, ...)
{
    if (!isTRUEorFALSE(strictly))
        stop("'strictly' must be TRUE or FALSE")
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    n <- length(x)
    if (n <= 1L)
        return(FALSE)

    frame <- x@frame

    # Get data in its natural order and compare with LAG window function
    conn <- tblconn(frame)

    # Add LAG values for comparison
    lag_mutate <- list(
        prev_seqnames = call("lag", as.name("seqnames")),
        prev_strand = call("lag", as.name("strand")),
        prev_start = call("lag", as.name("start")),
        prev_end = call("lag", as.name("end")))
    conn_with_lag <- mutate(conn, !!!lag_mutate)

    # Filter out first row (no previous value)
    filter_expr <- list(call("!", call("is.na", as.name("prev_seqnames"))))
    conn_filtered <- filter(conn_with_lag, !!!filter_expr)

    # Build the "is_unsorted" condition in SQL
    # Unsorted if current < prev (seqnames < prev_seqnames, or equal and start < prev_start, etc.)
    # String comparison in SQL: seqnames < prev_seqnames
    seq_lt <- call("<", as.name("seqnames"), as.name("prev_seqnames"))
    seq_eq <- call("==", as.name("seqnames"), as.name("prev_seqnames"))
    strand_lt <- call("<", as.name("strand"), as.name("prev_strand"))
    strand_eq <- call("==", as.name("strand"), as.name("prev_strand"))
    start_lt <- call("<", as.name("start"), as.name("prev_start"))
    start_eq <- call("==", as.name("start"), as.name("prev_start"))
    end_lt <- call("<", as.name("end"), as.name("prev_end"))
    end_eq <- call("==", as.name("end"), as.name("prev_end"))

    if (ignore.strand) {
        if (strictly) {
            # Strictly unsorted: seq < prev OR (seq == prev AND start < prev_start)
            # OR (seq == prev AND start == prev_start AND end <= prev_end)
            unsorted_expr <- call("|", seq_lt,
                call("|", call("&", seq_eq, start_lt),
                    call("&", call("&", seq_eq, start_eq),
                        call("|", end_lt, end_eq))))
        } else {
            # Non-strictly unsorted: seq < prev OR (seq == prev AND start < prev_start)
            # OR (seq == prev AND start == prev_start AND end < prev_end)
            unsorted_expr <- call("|", seq_lt,
                call("|", call("&", seq_eq, start_lt),
                    call("&", call("&", seq_eq, start_eq), end_lt)))
        }
    } else {
        if (strictly) {
            unsorted_expr <- call("|", seq_lt,
                call("|", call("&", seq_eq, strand_lt),
                    call("|", call("&", call("&", seq_eq, strand_eq), start_lt),
                        call("&", call("&", call("&", seq_eq, strand_eq), start_eq),
                            call("|", end_lt, end_eq)))))
        } else {
            unsorted_expr <- call("|", seq_lt,
                call("|", call("&", seq_eq, strand_lt),
                    call("|", call("&", call("&", seq_eq, strand_eq), start_lt),
                        call("&", call("&", call("&", seq_eq, strand_eq), start_eq), end_lt))))
        }
    }

    # Add is_unsorted column
    unsorted_mutate <- list(is_unsorted = unsorted_expr)
    conn_with_flag <- mutate(conn_filtered, !!!unsorted_mutate)

    # Aggregate: check if ANY row is unsorted using MAX (TRUE = 1, FALSE = 0)
    # In SQL: MAX(CAST(is_unsorted AS INT)) > 0
    agg_expr <- list(any_unsorted = call("max", call("as.integer", as.name("is_unsorted"))))
    result <- collect(summarize(conn_with_flag, !!!agg_expr))

    # Return TRUE if any row is unsorted
    isTRUE(result$any_unsorted > 0)
})

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dbplyr window_order
#' @importFrom dplyr arrange collect dense_rank group_by mutate row_number select ungroup
setMethod("rank", "DuckDBGRanges",
function(x, na.last = TRUE, ties.method = c("first", "min"), ignore.strand = FALSE, ...)
{
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    ties.method <- match.arg(ties.method)

    n <- length(x)
    if (n == 0L)
        return(integer(0L))

    frame <- x@frame

    # Get data in natural order and assign row positions
    conn <- tblconn(frame)
    pos_mutate <- list(orig_pos = call("row_number"))
    conn_with_pos <- mutate(conn, !!!pos_mutate)

    # Use SQL window functions for efficient in-database ranking
    # Use window_order() to set ordering for rank window functions
    # Note: use end instead of width (ordering by end is equivalent when start is fixed)
    if (ignore.strand) {
        ranked_conn <- window_order(conn_with_pos, seqnames, start, end)
    } else {
        ranked_conn <- window_order(conn_with_pos, seqnames, strand, start, end)
    }

    if (ties.method == "first") {
        # Use row_number() for unique ranks
        rank_mutate <- list(rnk = call("row_number"))
    } else {
        # ties.method == "min": use dense_rank() for tied ranks
        rank_mutate <- list(rnk = call("dense_rank"))
    }
    ranked_conn <- mutate(ranked_conn, !!!rank_mutate)

    # Reorder by original position and extract ranks
    reorder_cols <- list(as.name("orig_pos"))
    ranked_conn <- arrange(ranked_conn, !!!reorder_cols)
    select_cols <- list(as.name("rnk"))
    result <- collect(select(ranked_conn, !!!select_cols))
    as.integer(result$rnk)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Range restriction methods
###
### trim() and restrict() use SQL arithmetic for in-database computation.
###

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dplyr copy_to filter left_join mutate select tbl
#' @importFrom IRanges trim
#' @importFrom S4Vectors new2
#' @importFrom Seqinfo seqlengths isCircular
setMethod("trim", "DuckDBGRanges",
function(x, use.names = TRUE, ...)
{
    n <- length(x)
    if (n == 0L)
        return(x)

    # Get seqlengths and circularity
    seqlens <- seqlengths(x)
    is_circ <- isCircular(x)

    # If all seqlengths are NA or all are circular, nothing to trim
    if (all(is.na(seqlens)) || all(is_circ %in% TRUE))
        return(x)

    # Identify seqnames that need trimming (non-NA seqlength, non-circular)
    needs_trim <- !is.na(seqlens) & !(is_circ %in% TRUE)
    if (!any(needs_trim))
        return(x)

    # Create seqlengths data frame for upload
    seqlen_df <- data.frame(
        seqnames = names(seqlens)[needs_trim],
        seqlength = as.integer(seqlens[needs_trim]),
        stringsAsFactors = FALSE
    )

    frame <- x@frame
    conn <- tblconn(frame)
    db_conn <- dbconn(frame)

    # Upload seqlengths as temporary table
    table_name <- paste0("seqlen_trim_", sample.int(1e6, 1))
    seqlen_tbl <- copy_to(db_conn, seqlen_df, name = table_name,
                          temporary = TRUE, overwrite = TRUE)

    # Left join to get seqlength for each range
    conn_with_len <- left_join(conn, seqlen_tbl, by = "seqnames")

    # Apply trimming:
    # - start: max(start, 1) - but only if seqlength is not NA
    # - end: min(end, seqlength) - but only if seqlength is not NA
    # Use CASE WHEN to handle NULL seqlength (keep original)
    new_start_expr <- call("case_when",
        call("~", call("is.na", as.name("seqlength")), as.name("start")),
        call("~", TRUE, call("greatest", as.name("start"), 1L)))
    new_end_expr <- call("case_when",
        call("~", call("is.na", as.name("seqlength")), as.name("end")),
        call("~", TRUE, call("least", as.name("end"), as.name("seqlength"))))

    coord_mutate <- list(new_start = new_start_expr, new_end = new_end_expr)
    conn_trimmed <- mutate(conn_with_len, !!!coord_mutate)

    # Calculate new width
    width_mutate <- list(
        new_width = call("+", call("-", as.name("new_end"), as.name("new_start")), 1L))
    conn_trimmed <- mutate(conn_trimmed, !!!width_mutate)

    # Rename columns and drop seqlength
    select_list <- setNames(
        lapply(c("seqnames", "new_start", "new_end", "new_width", "strand"), as.name),
        c("seqnames", "start", "end", "width", "strand"))
    conn_trimmed <- select(conn_trimmed, !!!select_list)

    .build_DuckDBGRanges(conn_trimmed, seqinfo(x))
})

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importFrom DuckDBDataFrame set_row_number tblconn
#' @importFrom dplyr copy_to filter left_join mutate select
#' @importFrom IRanges restrict
#' @importFrom S4Vectors new2
setMethod("restrict", "DuckDBGRanges",
function(x, start = NA, end = NA, keep.all.ranges = FALSE, use.names = TRUE)
{
    n <- length(x)
    if (n == 0L)
        return(x)

    frame <- x@frame
    conn <- tblconn(frame)

    # Handle named start/end (per seqname) using table upload + join
    if (!is.null(names(start)) || !is.null(names(end))) {
        db_conn <- dbconn(frame)

        # Build constraint data frame
        all_seqnames <- unique(c(names(start), names(end)))
        constraint_df <- data.frame(
            seqnames = all_seqnames,
            restrict_start = as.integer(start[all_seqnames]),
            restrict_end = as.integer(end[all_seqnames]),
            stringsAsFactors = FALSE
        )
        # Replace NA with sentinel values that won't affect GREATEST/LEAST
        # For start: NA means no lower bound, use -Inf equivalent (very small number)
        # For end: NA means no upper bound, use Inf equivalent (very large number)
        constraint_df$restrict_start[is.na(constraint_df$restrict_start)] <- -.Machine$integer.max
        constraint_df$restrict_end[is.na(constraint_df$restrict_end)] <- .Machine$integer.max

        # Upload constraints as temporary table
        table_name <- paste0("restrict_", sample.int(1e6, 1))
        constraint_tbl <- copy_to(db_conn, constraint_df, name = table_name,
                                  temporary = TRUE, overwrite = TRUE)

        # Left join to get per-seqname constraints
        conn <- left_join(conn, constraint_tbl, by = "seqnames")

        # For seqnames not in constraint table, use original values
        # COALESCE handles NULL from left join
        new_start_expr <- call("case_when",
            call("~", call("is.na", as.name("restrict_start")), as.name("start")),
            call("~", TRUE, call("greatest", as.name("start"), as.name("restrict_start"))))
        new_end_expr <- call("case_when",
            call("~", call("is.na", as.name("restrict_end")), as.name("end")),
            call("~", TRUE, call("least", as.name("end"), as.name("restrict_end"))))
    } else {
        # Simple case: scalar start/end
        old_start <- as.name("start")
        old_end <- as.name("end")

        if (!is.na(start)) {
            new_start_expr <- call("greatest", old_start, as.integer(start))
        } else {
            new_start_expr <- old_start
        }

        if (!is.na(end)) {
            new_end_expr <- call("least", old_end, as.integer(end))
        } else {
            new_end_expr <- old_end
        }
    }

    # Apply mutations
    coord_mutate <- list(
        new_start = new_start_expr,
        new_end = new_end_expr)
    conn <- mutate(conn, !!!coord_mutate)

    # Calculate new width
    width_mutate <- list(
        new_width = call("+", call("-", as.name("new_end"), as.name("new_start")), 1L))
    conn <- mutate(conn, !!!width_mutate)

    # Filter out ranges that became invalid (new_end < new_start) unless keep.all.ranges
    if (!keep.all.ranges) {
        filter_expr <- list(call(">=", as.name("new_end"), as.name("new_start")))
        conn <- filter(conn, !!!filter_expr)
    }

    # Rename columns
    select_list <- setNames(
        lapply(c("seqnames", "new_start", "new_end", "new_width", "strand"), as.name),
        c("seqnames", "start", "end", "width", "strand"))
    conn <- select(conn, !!!select_list)

    # Add row_number
    rownum_mutate <- list(row_number = call("row_number"))
    conn <- mutate(conn, !!!rownum_mutate)

    datacols <- expression(
        seqnames = seqnames,
        start = start,
        end = end,
        width = width,
        strand = strand
    )

    keycols <- list(row_number = set_row_number(conn))

    new_frame <- new2("DuckDBDataFrame",
                      conn = conn,
                      datacols = datacols,
                      keycols = keycols,
                      dimtbls = new.env(parent = emptyenv()),
                      check = FALSE)

    new2("DuckDBGRanges",
         frame = new_frame,
         seqinfo = seqinfo(x),
         elementMetadata = new("DFrame", nrows = nrow(new_frame)),
         check = FALSE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reduce()
###
### Merge overlapping and adjacent ranges into single ranges per seqname/strand.
### Uses SQL window functions to identify merge groups via cumulative max of end.
###

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom IRanges reduce
#' @importFrom dbplyr window_order
#' @importFrom dplyr mutate filter select arrange group_by summarize ungroup
#' @importFrom S4Vectors new2
setMethod("reduce", "DuckDBGRanges",
function(x, drop.empty.ranges = FALSE, min.gapwidth = 1L,
         with.revmap = FALSE, with.inframe.attrib = FALSE,
         ignore.strand = FALSE)
{
    if (!isTRUEorFALSE(drop.empty.ranges))
        stop("'drop.empty.ranges' must be TRUE or FALSE")
    if (!isSingleNumber(min.gapwidth))
        stop("'min.gapwidth' must be a single integer")
    min.gapwidth <- as.integer(min.gapwidth)
    if (min.gapwidth < 0L)
        stop("'min.gapwidth' must be non-negative")
    if (!isTRUEorFALSE(with.revmap))
        stop("'with.revmap' must be TRUE or FALSE")
    if (with.revmap)
        warning("'with.revmap' is not supported for DuckDBGRanges, ignoring")
    if (!identical(with.inframe.attrib, FALSE))
        stop("'with.inframe.attrib' is not supported for DuckDBGRanges")
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    if (length(x) == 0L) {
        # Return empty DuckDBGRanges
        return(x)
    }

    frame <- x@frame
    conn <- tblconn(frame)

    # Drop empty ranges if requested (width <= 0)
    if (drop.empty.ranges) {
        filter_expr <- list(call(">", as.name("end"), call("-", as.name("start"), 1L)))
        conn <- filter(conn, !!!filter_expr)
    }

    # Optimized reduce algorithm using consolidated window function passes:
    #
    # The algorithm uses cumulative max to correctly handle nested ranges:
    #   Range A: 100-500, Range B: 150-200, Range C: 400-600
    #   With LAG(end): C sees prev_end=200 (from B), thinks 400 > 200 = new group. WRONG!
    #   With cummax: C sees max_end_so_far=500 (from A), 400 <= 500 = same group. CORRECT!
    #
    # Optimizations applied:
    # 1. Single code path for ignore.strand=TRUE/FALSE (only group_cols differ)
    # 2. Consolidated window passes: cummax → lag+case_when → cumsum (3 passes instead of 4)
    # 3. Minimized group_by/ungroup cycles by keeping grouping where possible

    # Determine grouping columns based on ignore.strand
    if (ignore.strand) {
        strand_mutate <- list(strand = "*")
        conn <- mutate(conn, !!!strand_mutate)
        group_cols <- list(as.name("seqnames"))
    } else {
        group_cols <- list(as.name("seqnames"), as.name("strand"))
    }

    # Pass 1: Compute cumulative max of end (running max including current row)
    conn <- group_by(conn, !!!group_cols)
    conn <- window_order(conn, start, end)
    cummax_expr <- list(cummax_end = call("cummax", as.name("end")))
    conn <- mutate(conn, !!!cummax_expr)

    # Pass 2: LAG the cummax + compute new_group flag in single mutate
    # (LAG is a window function, case_when is not - they can share a mutate)
    # Must be separate from cummax because DuckDB doesn't allow nested window functions
    lag_and_group_expr <- list(
        max_end_so_far = call("lag", as.name("cummax_end")),
        new_group = call("case_when",
            call("~", call("is.na", call("lag", as.name("cummax_end"))), 1L),
            call("~", call(">", as.name("start"),
                           call("+", call("lag", as.name("cummax_end")), min.gapwidth)), 1L),
            call("~", TRUE, 0L))
    )
    conn <- mutate(conn, !!!lag_and_group_expr)

    # Pass 3: Cumulative sum to get group IDs
    cumsum_expr <- list(grp_id = call("cumsum", as.name("new_group")))
    conn <- mutate(conn, !!!cumsum_expr)
    conn <- ungroup(conn)

    # Aggregate by seqnames + strand + grp_id
    agg_group_cols <- list(as.name("seqnames"), as.name("strand"), as.name("grp_id"))
    agg_expr <- list(
        start = call("min", as.name("start"), na.rm = TRUE),
        end = call("max", as.name("end"), na.rm = TRUE)
    )
    conn <- group_by(conn, !!!agg_group_cols)
    conn <- summarize(conn, !!!agg_expr, .groups = "drop")

    # Compute width
    width_mutate <- list(
        width = call("+", call("-", as.name("end"), as.name("start")), 1L)
    )
    conn <- mutate(conn, !!!width_mutate)

    # Select final columns and build result
    select_cols <- list(
        seqnames = as.name("seqnames"),
        start = as.name("start"),
        end = as.name("end"),
        width = as.name("width"),
        strand = as.name("strand")
    )
    conn <- select(conn, !!!select_cols)

    .build_DuckDBGRanges(conn, seqinfo(x))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### gaps()
###
### Find gaps (uncovered regions) between ranges per seqname/strand.
### Uses coverage tracking: gaps are intervals where coverage == 0.
###

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom IRanges gaps
#' @importFrom dbplyr window_order
#' @importFrom dplyr arrange copy_to filter group_by inner_join mutate select summarize ungroup union_all
#' @importFrom S4Vectors new2
#' @importFrom Seqinfo seqlengths seqlevels
setMethod("gaps", "DuckDBGRanges",
function(x, start = 1L, end = seqlengths(x), ignore.strand = FALSE)
{
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    x_seqlevels <- seqlevels(x)
    x_seqlengths <- seqlengths(x)

    # Normalize start and end vectors
    start <- as.integer(.recycle_per_seqlevel(start, x_seqlevels))
    end <- as.integer(.recycle_per_seqlevel(end, x_seqlevels))
    names(start) <- x_seqlevels
    names(end) <- x_seqlevels

    # Build bounds data frame
    if (ignore.strand) {
        bounds_df <- data.frame(
            seqnames = x_seqlevels,
            strand = "*",
            bound_start = start,
            bound_end = end,
            stringsAsFactors = FALSE
        )
    } else {
        # Replicate bounds for each strand (+, -, *)
        bounds_list <- lapply(c("+", "-", "*"), function(str) {
            data.frame(
                seqnames = x_seqlevels,
                strand = str,
                bound_start = start,
                bound_end = end,
                stringsAsFactors = FALSE
            )
        })
        bounds_df <- do.call(rbind, bounds_list)
    }
    bounds_df <- bounds_df[!is.na(bounds_df$bound_start) & !is.na(bounds_df$bound_end) &
                           bounds_df$bound_end >= bounds_df$bound_start, ]

    if (nrow(bounds_df) == 0L) {
        # No valid bounds - return empty
        return(x[integer(0)])
    }

    # Handle empty input case
    if (length(x) == 0L) {
        # No ranges - gaps are the full bounds for each seqname/strand
        gap_data <- bounds_df
        gap_data$start <- gap_data$bound_start
        gap_data$end <- gap_data$bound_end
        gap_data$width <- gap_data$end - gap_data$start + 1L
        gap_data <- gap_data[, c("seqnames", "start", "end", "width", "strand")]

        gap_tf <- tempfile(fileext = ".parquet")
        arrow::write_parquet(gap_data, gap_tf)
        return(DuckDBGRanges(gap_tf, seqnames = "seqnames", start = "start",
                            end = "end", strand = "strand", seqinfo = seqinfo(x)))
    }

    # Optimized gaps algorithm using coverage tracking:
    #
    # 1. Create boundary events: bound_start = 0 (marker), bound_end+1 = 0 (marker)
    # 2. Create range events: range_start = +1, range_end+1 = -1
    # 3. Compute cumulative coverage at each position
    # 4. Create intervals between consecutive positions
    # 5. Keep intervals where coverage == 0 AND within bounds
    #
    # This avoids calling reduce() and eliminates the need for
    # three separate UNION ALL operations.

    frame <- x@frame
    conn <- tblconn(frame)
    db_conn <- dbconn(frame)

    # Upload bounds table
    bounds_tbl <- copy_to(db_conn, bounds_df,
                          name = paste0("bounds_", sample.int(1e6, 1)),
                          temporary = TRUE, overwrite = TRUE)

    # Handle strand for grouping
    if (ignore.strand) {
        strand_mutate <- list(strand = "*")
        conn <- mutate(conn, !!!strand_mutate)
    }

    # Create range events: start = +1, end+1 = -1
    start_events <- select(conn, !!!list(
        seqnames = as.name("seqnames"),
        strand = as.name("strand"),
        bp = as.name("start")
    ))
    start_events <- mutate(start_events, !!!list(delta = 1L))

    end_events <- mutate(conn, !!!list(bp = call("+", as.name("end"), 1L)))
    end_events <- select(end_events, !!!list(
        seqnames = as.name("seqnames"),
        strand = as.name("strand"),
        bp = as.name("bp")
    ))
    end_events <- mutate(end_events, !!!list(delta = -1L))

    # Create boundary marker events (delta = 0, just to mark positions)
    bound_start_events <- select(bounds_tbl, !!!list(
        seqnames = as.name("seqnames"),
        strand = as.name("strand"),
        bp = as.name("bound_start")
    ))
    bound_start_events <- mutate(bound_start_events, !!!list(delta = 0L))

    bound_end_events <- mutate(bounds_tbl, !!!list(bp = call("+", as.name("bound_end"), 1L)))
    bound_end_events <- select(bound_end_events, !!!list(
        seqnames = as.name("seqnames"),
        strand = as.name("strand"),
        bp = as.name("bp")
    ))
    bound_end_events <- mutate(bound_end_events, !!!list(delta = 0L))

    # Combine all events
    all_events <- union_all(start_events, end_events)
    all_events <- union_all(all_events, bound_start_events)
    all_events <- union_all(all_events, bound_end_events)

    # Aggregate deltas at each position
    agg_group_cols <- list(as.name("seqnames"), as.name("strand"), as.name("bp"))
    all_events <- group_by(all_events, !!!agg_group_cols)
    all_events <- summarize(all_events, !!!list(delta = call("sum", as.name("delta"))), .groups = "drop")

    # Compute cumulative coverage using window function
    all_events <- group_by(all_events, !!!list(as.name("seqnames"), as.name("strand")))
    all_events <- window_order(all_events, !!as.name("bp"))
    all_events <- mutate(all_events, !!!list(
        coverage = call("cumsum", as.name("delta")),
        next_bp = call("lead", as.name("bp"))
    ))
    all_events <- ungroup(all_events)

    # Create intervals: [bp, next_bp - 1] where coverage == 0
    # Filter out: last position in group (next_bp is NULL), non-zero coverage
    result <- filter(all_events, !!!list(
        call("!", call("is.na", as.name("next_bp"))),
        call("==", as.name("coverage"), 0L)
    ))

    result <- mutate(result, !!!list(
        interval_start = as.name("bp"),
        interval_end = call("-", as.name("next_bp"), 1L)
    ))

    # Join with bounds to filter to valid gap regions
    result <- inner_join(result, bounds_tbl, by = c("seqnames", "strand"))

    # Clip gaps to bounds and filter valid ones
    result <- mutate(result, !!!list(
        start = call("greatest", as.name("interval_start"), as.name("bound_start")),
        end = call("least", as.name("interval_end"), as.name("bound_end"))
    ))

    # Keep only valid gaps (start <= end and within bounds)
    result <- filter(result, !!!list(
        call("<=", as.name("start"), as.name("end")),
        call(">=", as.name("start"), as.name("bound_start")),
        call("<=", as.name("end"), as.name("bound_end"))
    ))

    # Select final columns
    result <- select(result, !!!list(
        seqnames = as.name("seqnames"),
        start = as.name("start"),
        end = as.name("end"),
        strand = as.name("strand")
    ))

    # Compute width
    width_mutate <- list(
        width = call("+", call("-", as.name("end"), as.name("start")), 1L)
    )
    result <- mutate(result, !!!width_mutate)

    .build_DuckDBGRanges(result, seqinfo(x))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### disjoin()
###
### Break ranges into non-overlapping pieces at all breakpoints.
### Creates new ranges between consecutive unique start/end positions.
###

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom IRanges disjoin
#' @importFrom dbplyr window_order
#' @importFrom dplyr arrange distinct filter group_by mutate select summarize ungroup union_all
#' @importFrom S4Vectors new2
setMethod("disjoin", "DuckDBGRanges",
function(x, with.revmap = FALSE, ignore.strand = FALSE)
{
    if (!isTRUEorFALSE(with.revmap))
        stop("'with.revmap' must be TRUE or FALSE")
    if (with.revmap)
        warning("'with.revmap' is not supported for DuckDBGRanges, ignoring")
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    if (length(x) == 0L) {
        return(x)
    }

    frame <- x@frame
    conn <- tblconn(frame)
    db_conn <- dbconn(frame)

    # Optimized disjoin algorithm using coverage tracking (no O(n²) join):
    #
    # 1. Create events: each range start is +1, each end+1 is -1
    # 2. At each unique position, compute net change in coverage
    # 3. Compute cumulative sum to get coverage at each position
    # 4. An interval is valid when coverage > 0 at its start
    #
    # This avoids the expensive join between intervals and original ranges.

    # Determine grouping columns based on ignore.strand
    if (ignore.strand) {
        strand_mutate <- list(strand = "*")
        conn <- mutate(conn, !!!strand_mutate)
        group_cols <- list(as.name("seqnames"))
    } else {
        group_cols <- list(as.name("seqnames"), as.name("strand"))
    }

    # Create start events (+1 coverage)
    start_events <- select(conn, !!!list(
        seqnames = as.name("seqnames"),
        strand = as.name("strand"),
        bp = as.name("start")
    ))
    start_events <- mutate(start_events, !!!list(delta = 1L))

    # Create end events (-1 coverage at end + 1)
    end_events <- mutate(conn, !!!list(bp = call("+", as.name("end"), 1L)))
    end_events <- select(end_events, !!!list(
        seqnames = as.name("seqnames"),
        strand = as.name("strand"),
        bp = as.name("bp")
    ))
    end_events <- mutate(end_events, !!!list(delta = -1L))

    # Combine all events
    all_events <- union_all(start_events, end_events)

    # Aggregate deltas at each position (multiple ranges may start/end at same pos)
    # Group by seqnames, strand, bp (strand always included for final select)
    agg_group_cols <- list(as.name("seqnames"), as.name("strand"), as.name("bp"))
    all_events <- group_by(all_events, !!!agg_group_cols)
    all_events <- summarize(all_events, !!!list(delta = call("sum", as.name("delta"))), .groups = "drop")

    # Compute cumulative coverage using window function
    all_events <- group_by(all_events, !!!group_cols)
    all_events <- window_order(all_events, !!as.name("bp"))
    all_events <- mutate(all_events, !!!list(
        coverage = call("cumsum", as.name("delta")),
        next_bp = call("lead", as.name("bp"))
    ))
    all_events <- ungroup(all_events)

    # Create intervals: [bp, next_bp - 1] where coverage > 0
    # Filter out: last position in group (next_bp is NULL) and zero coverage
    result <- filter(all_events, !!!list(
        call("!", call("is.na", as.name("next_bp"))),
        call(">", as.name("coverage"), 0L)
    ))

    result <- mutate(result, !!!list(
        start = as.name("bp"),
        end = call("-", as.name("next_bp"), 1L)
    ))

    result <- select(result, !!!list(
        seqnames = as.name("seqnames"),
        start = as.name("start"),
        end = as.name("end"),
        strand = as.name("strand")
    ))

    # Compute width
    width_mutate <- list(
        width = call("+", call("-", as.name("end"), as.name("start")), 1L)
    )
    result <- mutate(result, !!!width_mutate)

    .build_DuckDBGRanges(result, seqinfo(x))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### coverage()
###
### Per-base coverage depth per seqname. Strand is ignored, matching
### GenomicRanges::coverage() (a position's read depth does not split by
### strand). Uses the same delta-event + window-function cumsum() sweep-line
### pattern as disjoin()/gaps(), but keeps the running sum itself instead of
### thresholding it, and turns the compact per-breakpoint result (bounded by
### the number of ranges, never by genome length) into a run-length encoding.
###

# Reconstruct a SimpleRleList from coverage()'s compact per-seqname breakpoint
# table (columns: seqnames, bp, depth; one row per position where the summed
# delta is non-zero-or-not, sorted by bp within seqname). By construction the
# last breakpoint for each seqname always has depth 0 (every start event has a
# matching end event), so it is used only to close the final run, never
# emitted itself; a leading zero run is added when the first breakpoint is
# after position 1. `width` is a named-by-seqlevel vector (NA entries mean
# "natural extent": stop exactly at the last covered position). `int_values`
# controls whether runs are stored as integer or numeric, matching
# GenomicRanges::coverage()'s own contract (integer unless a non-integer
# 'weight' is used) -- DuckDB's SUM()/cumsum() always return a double
# regardless of the input type, so this is not just cosmetic.
#' @importFrom IRanges RleList
#' @importFrom S4Vectors Rle window
.breakpoints_to_rlelist <- function(breakpoints, seqlevels, width, int_values) {
    zero <- if (int_values) 0L else 0
    rles <- lapply(seqlevels, function(sn) {
        w <- width[[sn]]
        rows <- breakpoints[breakpoints$seqnames == sn, , drop = FALSE]
        if (nrow(rows) == 0L)
            return(Rle(zero, if (is.na(w)) 0L else w))

        bp <- rows$bp
        depth <- rows$depth
        if (int_values) depth <- as.integer(depth)
        n <- length(bp)

        # bp is guaranteed >= 1 by the clipping done in .coverage_events_tbl().
        values <- depth[-n]
        lengths <- diff(bp)
        if (bp[1L] > 1L) {
            values <- c(zero, values)
            lengths <- c(bp[1L] - 1L, lengths)
        }
        rle <- Rle(values = values, lengths = lengths)

        if (!is.na(w)) {
            if (length(rle) < w) {
                rle <- c(rle, Rle(zero, w - length(rle)))
            } else if (length(rle) > w) {
                rle <- window(rle, start = 1L, end = w)
            }
        }
        rle
    })
    names(rles) <- seqlevels
    RleList(rles, compress = FALSE)
}

# Build the lazy per-(seqname, breakpoint) running-depth tbl for coverage(),
# without collecting it -- exposed separately so tests can inspect its query
# plan the same way .overlap_join_tbl() does for findOverlaps(). Ranges are
# clipped to positions >= 1 after 'shift' is applied (matching
# GenomicRanges::coverage(), which silently drops the portion of a shifted
# range that falls before position 1, rather than erroring): a range that
# lands entirely before 1 is dropped, and one that straddles 1 has its start
# event clipped up to 1.
.coverage_events_tbl <- function(x, x_seqlevels, shift, weight) {
    frame <- x@frame
    conn <- tblconn(frame)
    db_conn <- dbconn(frame)

    shift_df <- data.frame(seqnames = x_seqlevels, .shift = shift,
                           stringsAsFactors = FALSE)
    shift_tbl <- copy_to(db_conn, shift_df,
                         name = paste0("coverage_shift_", sample.int(1e6, 1)),
                         temporary = TRUE, overwrite = TRUE)

    conn <- inner_join(conn, shift_tbl, by = "seqnames")

    # Drop ranges that land entirely before position 1 after shifting.
    conn <- filter(conn, !!!list(
        call(">=", call("+", as.name("end"), as.name(".shift")), 1L)
    ))

    # Start events: +weight at max(start + shift, 1); end events: -weight at
    # (end + shift + 1).
    start_events <- mutate(conn, !!!list(
        bp = call("greatest", call("+", as.name("start"), as.name(".shift")), 1L)
    ))
    start_events <- select(start_events, !!!list(
        seqnames = as.name("seqnames"),
        bp = as.name("bp")
    ))
    start_events <- mutate(start_events, !!!list(delta = weight))

    end_events <- mutate(conn, !!!list(
        bp = call("+", call("+", as.name("end"), as.name(".shift")), 1L)
    ))
    end_events <- select(end_events, !!!list(
        seqnames = as.name("seqnames"),
        bp = as.name("bp")
    ))
    end_events <- mutate(end_events, !!!list(delta = -weight))

    all_events <- union_all(start_events, end_events)

    agg_group_cols <- list(as.name("seqnames"), as.name("bp"))
    all_events <- group_by(all_events, !!!agg_group_cols)
    all_events <- summarize(all_events,
        !!!list(delta = call("sum", as.name("delta"), na.rm = TRUE)),
        .groups = "drop")

    all_events <- group_by(all_events, !!!list(as.name("seqnames")))
    all_events <- window_order(all_events, !!as.name("bp"))
    all_events <- mutate(all_events,
        !!!list(depth = call("cumsum", as.name("delta"))))
    all_events <- ungroup(all_events)

    arrange(all_events, !!!list(as.name("seqnames"), as.name("bp")))
}

#' @export
#' @importFrom dbplyr window_order
#' @importFrom dplyr arrange collect copy_to filter group_by inner_join mutate
#' @importFrom dplyr select summarize ungroup union_all
#' @importFrom DuckDBDataFrame dbconn tblconn
#' @importFrom IRanges coverage
#' @importFrom S4Vectors Rle
#' @importFrom Seqinfo seqlengths seqlevels
setMethod("coverage", "DuckDBGRanges",
function(x, shift = 0L, width = NULL, weight = 1L,
         method = c("auto", "sort", "hash"), ...)
{
    if (!is.numeric(weight) || length(weight) != 1L || is.na(weight))
        stop("'weight' must be a single number for DuckDBGRanges ",
             "(a per-range or mcols-column 'weight' is not supported)")
    int_values <- is.integer(weight)

    x_seqlevels <- seqlevels(x)
    x_seqlengths <- seqlengths(x)

    shift <- as.integer(.recycle_per_seqlevel(shift, x_seqlevels))
    width <- if (is.null(width)) x_seqlengths
             else as.integer(.recycle_per_seqlevel(width, x_seqlevels))
    names(width) <- x_seqlevels

    if (length(x) == 0L) {
        zero <- if (int_values) 0L else 0
        rles <- lapply(width, function(w) Rle(zero, if (is.na(w)) 0L else w))
        names(rles) <- x_seqlevels
        return(RleList(rles, compress = FALSE))
    }

    all_events <- .coverage_events_tbl(x, x_seqlevels, shift, weight)

    breakpoints <- collect(select(all_events, !!!list(
        seqnames = as.name("seqnames"),
        bp = as.name("bp"),
        depth = as.name("depth")
    )))

    .breakpoints_to_rlelist(breakpoints, x_seqlevels, width, int_values)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Set operations: union(), intersect(), setdiff()
###
### These methods combine two DuckDBGRanges objects using SQL operations.
###

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importFrom DuckDBDataFrame set_row_number tblconn
#' @importFrom dplyr union_all
#' @importFrom S4Vectors new2
setMethod("union", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE)
{
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    # Union is simply: reduce(c(x, y), drop.empty.ranges=TRUE)
    # We need to combine the two connections and then reduce

    frame_x <- x@frame
    frame_y <- y@frame
    conn_x <- tblconn(frame_x)
    conn_y <- tblconn(frame_y)

    if (ignore.strand) {
        strand_mutate <- list(strand = "*")
        conn_x <- mutate(conn_x, !!!strand_mutate)
        conn_y <- mutate(conn_y, !!!strand_mutate)
    }

    # Select only the core columns from both
    core_cols <- lapply(c("seqnames", "start", "end", "strand"), as.name)
    conn_x <- select(conn_x, !!!core_cols)
    conn_y <- select(conn_y, !!!core_cols)

    # Combine using union_all
    combined <- union_all(conn_x, conn_y)

    # Compute width
    width_mutate <- list(
        width = call("+", call("-", as.name("end"), as.name("start")), 1L)
    )
    combined <- mutate(combined, !!!width_mutate)

    # Add row_number as keycol
    rownum_mutate <- list(row_number = call("row_number"))
    combined <- mutate(combined, !!!rownum_mutate)

    datacols <- expression(
        seqnames = seqnames,
        start = start,
        end = end,
        width = width,
        strand = strand
    )

    keycols <- list(row_number = set_row_number(combined))

    new_frame <- new2("DuckDBDataFrame",
                      conn = combined,
                      datacols = datacols,
                      keycols = keycols,
                      dimtbls = new.env(parent = emptyenv()),
                      check = FALSE)

    # Merge seqinfo
    merged_seqinfo <- merge(seqinfo(x), seqinfo(y))

    temp_gr <- new2("DuckDBGRanges",
                    frame = new_frame,
                    seqinfo = merged_seqinfo,
                    elementMetadata = new("DFrame", nrows = nrow(new_frame)),
                    check = FALSE)

    # Reduce to merge overlapping ranges
    reduce(temp_gr, drop.empty.ranges = TRUE, ignore.strand = ignore.strand)
})

#' @export
setMethod("union", c("DuckDBGRanges", "GRanges"),
function(x, y, ignore.strand = FALSE)
{
    y_ddb <- .granges_to_duckdb(y, x)
    union(x, y_ddb, ignore.strand = ignore.strand)
})

#' @export
setMethod("union", c("GRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE)
{
    x_ddb <- .granges_to_duckdb(x, y)
    union(x_ddb, y, ignore.strand = ignore.strand)
})

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dplyr arrange distinct filter group_by inner_join mutate select summarize
#' @importFrom S4Vectors new2
setMethod("intersect", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE)
{
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    # Intersect approach:
    # 1. Find all overlapping pairs between reduced(x) and reduced(y)
    # 2. For each overlapping pair, compute intersection: [MAX(x.start, y.start), MIN(x.end, y.end)]
    # 3. Remove duplicates with distinct()
    #
    # Note: No final reduce() is needed because:
    # - x_reduced ranges are non-overlapping with each other
    # - y_reduced ranges are non-overlapping with each other
    # - Intersections from different x ranges can't overlap (bounded by non-overlapping x ranges)
    # - Intersections from different y ranges can't overlap (bounded by non-overlapping y ranges)
    # - The only potential duplicates are exact coordinate matches, handled by distinct()

    # First reduce both inputs
    x_reduced <- reduce(x, ignore.strand = ignore.strand)
    y_reduced <- reduce(y, ignore.strand = ignore.strand)

    if (length(x_reduced) == 0L || length(y_reduced) == 0L) {
        return(x_reduced[integer(0)])  # Return empty
    }

    frame_x <- x_reduced@frame
    frame_y <- y_reduced@frame
    conn_x <- tblconn(frame_x)
    conn_y <- tblconn(frame_y)

    # Rename columns for join
    x_rename <- list(
        x_start = as.name("start"),
        x_end = as.name("end")
    )
    conn_x <- mutate(conn_x, !!!x_rename)
    conn_x <- select(conn_x, !!!lapply(c("seqnames", "strand", "x_start", "x_end"), as.name))

    y_rename <- list(
        y_start = as.name("start"),
        y_end = as.name("end")
    )
    conn_y <- mutate(conn_y, !!!y_rename)
    conn_y <- select(conn_y, !!!lapply(c("seqnames", "strand", "y_start", "y_end"), as.name))

    # Join on seqnames (and strand if not ignore.strand)
    if (ignore.strand) {
        joined <- inner_join(conn_x, conn_y, by = "seqnames")
        # Use strand from x (which is "*" if ignore.strand was applied in reduce)
        strand_fix <- list(strand = as.name("strand.x"))
        joined <- mutate(joined, !!!strand_fix)
    } else {
        joined <- inner_join(conn_x, conn_y, by = c("seqnames", "strand"))
    }

    # Filter for actual overlaps: x_start <= y_end AND x_end >= y_start
    overlap_filter <- list(
        call("<=", as.name("x_start"), as.name("y_end")),
        call(">=", as.name("x_end"), as.name("y_start"))
    )
    joined <- filter(joined, !!!overlap_filter)

    # Compute intersection coordinates: [MAX(x_start, y_start), MIN(x_end, y_end)]
    intersect_mutate <- list(
        start = call("greatest", as.name("x_start"), as.name("y_start")),
        end = call("least", as.name("x_end"), as.name("y_end"))
    )
    joined <- mutate(joined, !!!intersect_mutate)

    # Compute width
    width_mutate <- list(
        width = call("+", call("-", as.name("end"), as.name("start")), 1L)
    )
    joined <- mutate(joined, !!!width_mutate)

    # Select final columns and remove duplicates
    final_select <- lapply(c("seqnames", "start", "end", "width", "strand"), as.name)
    result <- select(joined, !!!final_select)
    result <- distinct(result)

    merged_seqinfo <- merge(seqinfo(x), seqinfo(y))
    .build_DuckDBGRanges(result, merged_seqinfo)
})

#' @export
setMethod("intersect", c("DuckDBGRanges", "GRanges"),
function(x, y, ignore.strand = FALSE)
{
    y_ddb <- .granges_to_duckdb(y, x)
    intersect(x, y_ddb, ignore.strand = ignore.strand)
})

#' @export
setMethod("intersect", c("GRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE)
{
    x_ddb <- .granges_to_duckdb(x, y)
    intersect(x_ddb, y, ignore.strand = ignore.strand)
})

#' @export
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dbplyr window_order
#' @importFrom dplyr arrange distinct filter group_by lag left_join mutate select summarize ungroup union_all
#' @importFrom S4Vectors new2
setMethod("setdiff", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE)
{
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    # setdiff algorithm using event-based approach:
    # For each x range, we need to subtract ALL overlapping y ranges.
    # The correct approach is:
    # 1. For each x range, collect all overlapping y ranges
    # 2. Sort y ranges by start within each x
    # 3. Compute gaps: [x_start, first_y_start-1], [y_end+1, next_y_start-1], [last_y_end+1, x_end]
    # 4. Clip these to x bounds

    x_reduced <- reduce(x, ignore.strand = ignore.strand)
    y_reduced <- reduce(y, ignore.strand = ignore.strand)

    if (length(x_reduced) == 0L) {
        return(x_reduced)
    }

    if (length(y_reduced) == 0L) {
        return(x_reduced)
    }

    frame_x <- x_reduced@frame
    frame_y <- y_reduced@frame
    conn_x <- tblconn(frame_x)
    conn_y <- tblconn(frame_y)

    # Add unique IDs to x ranges
    x_id_mutate <- list(x_id = call("row_number"))
    conn_x <- mutate(conn_x, !!!x_id_mutate)

    # Rename x columns
    x_cols <- list(
        x_start = as.name("start"),
        x_end = as.name("end"),
        x_strand = as.name("strand")
    )
    conn_x <- mutate(conn_x, !!!x_cols)
    conn_x <- select(conn_x, !!!lapply(c("seqnames", "x_strand", "x_id", "x_start", "x_end"), as.name))

    # Rename y columns
    y_cols <- list(
        y_start = as.name("start"),
        y_end = as.name("end")
    )
    conn_y <- mutate(conn_y, !!!y_cols)
    if (ignore.strand) {
        conn_y <- select(conn_y, !!!lapply(c("seqnames", "y_start", "y_end"), as.name))
    } else {
        y_strand <- list(y_strand = as.name("strand"))
        conn_y <- mutate(conn_y, !!!y_strand)
        conn_y <- select(conn_y, !!!lapply(c("seqnames", "y_strand", "y_start", "y_end"), as.name))
    }

    # Left join x with y
    if (ignore.strand) {
        joined <- left_join(conn_x, conn_y, by = "seqnames")
    } else {
        joined <- left_join(conn_x, conn_y, by = c("seqnames", x_strand = "y_strand"))
    }

    # Identify overlaps: y overlaps x if y_start <= x_end AND y_end >= x_start
    overlap_cond <- call("case_when",
        call("~", call("is.na", as.name("y_start")), FALSE),
        call("~", call("&",
                       call("<=", as.name("y_start"), as.name("x_end")),
                       call(">=", as.name("y_end"), as.name("x_start"))), TRUE),
        call("~", TRUE, FALSE))
    joined <- mutate(joined, has_overlap = !!overlap_cond)

    # Part 1: X ranges with NO overlapping y at all - keep them as-is
    # An x range has no overlap if ALL its joined rows have has_overlap = FALSE
    # We identify these by grouping by x_id and checking if max(has_overlap) = 0
    x_has_any_overlap <- group_by(joined, !!!lapply(c("seqnames", "x_strand", "x_id", "x_start", "x_end"), as.name))
    # Use sum() on the boolean (TRUE=1, FALSE=0) to count overlaps
    any_overlap_flag <- list(
        overlap_count = call("sum", call("if_else", as.name("has_overlap"), 1L, 0L))
    )
    x_has_any_overlap <- summarize(x_has_any_overlap, !!!any_overlap_flag)
    x_has_any_overlap <- ungroup(x_has_any_overlap)

    # X ranges with no overlap at all (overlap_count = 0)
    no_overlap <- filter(x_has_any_overlap, !!!list(call("==", as.name("overlap_count"), 0L)))
    no_overlap <- mutate(no_overlap,
                         start = !!as.name("x_start"),
                         end = !!as.name("x_end"),
                         strand = !!as.name("x_strand"))
    no_overlap <- select(no_overlap, !!!lapply(c("seqnames", "start", "end", "strand"), as.name))

    # Part 2: X ranges WITH at least one overlapping y - compute gaps
    # Filter to keep only rows where y actually overlaps x
    has_overlap <- filter(joined, !!!list(as.name("has_overlap")))

    # Clip y ranges to x bounds for proper gap computation
    clip_mutate <- list(
        y_start_clipped = call("greatest", as.name("y_start"), as.name("x_start")),
        y_end_clipped = call("least", as.name("y_end"), as.name("x_end"))
    )
    has_overlap <- mutate(has_overlap, !!!clip_mutate)

    # For each x_id, we need to find gaps between the clipped y ranges
    # Group by x_id and seqnames, order by y_start_clipped
    has_overlap <- group_by(has_overlap, !!!lapply(c("seqnames", "x_strand", "x_id", "x_start", "x_end"), as.name))
    has_overlap <- window_order(has_overlap, !!as.name("y_start_clipped"))

    # Compute the previous y_end using lag, and running max of y_end
    # We need cummax to handle nested y ranges
    prev_mutate <- list(
        cummax_y_end = call("cummax", as.name("y_end_clipped")),
        y_row = call("row_number")
    )
    has_overlap <- mutate(has_overlap, !!!prev_mutate)

    prev_end_mutate <- list(
        prev_y_end = call("lag", as.name("cummax_y_end"))
    )
    has_overlap <- mutate(has_overlap, !!!prev_end_mutate)
    has_overlap <- ungroup(has_overlap)

    # Gap before first y: [x_start, min(y_start_clipped) - 1]
    # We identify first y by y_row == 1
    first_y <- filter(has_overlap, !!!list(call("==", as.name("y_row"), 1L)))
    gap_before <- mutate(first_y,
                         gap_start = !!as.name("x_start"),
                         gap_end = !!call("-", as.name("y_start_clipped"), 1L))
    gap_before <- filter(gap_before, !!!list(call("<=", as.name("gap_start"), as.name("gap_end"))))
    gap_before <- mutate(gap_before,
                         start = !!as.name("gap_start"),
                         end = !!as.name("gap_end"),
                         strand = !!as.name("x_strand"))
    gap_before <- select(gap_before, !!!lapply(c("seqnames", "start", "end", "strand"), as.name))

    # Gaps between consecutive y ranges: [prev_y_end + 1, y_start_clipped - 1]
    # These occur when y_row > 1 and there's a gap
    gaps_between <- filter(has_overlap, !!!list(call(">", as.name("y_row"), 1L)))
    gaps_between <- mutate(gaps_between,
                           gap_start = !!call("+", as.name("prev_y_end"), 1L),
                           gap_end = !!call("-", as.name("y_start_clipped"), 1L))
    gaps_between <- filter(gaps_between, !!!list(call("<=", as.name("gap_start"), as.name("gap_end"))))
    gaps_between <- mutate(gaps_between,
                           start = !!as.name("gap_start"),
                           end = !!as.name("gap_end"),
                           strand = !!as.name("x_strand"))
    gaps_between <- select(gaps_between, !!!lapply(c("seqnames", "start", "end", "strand"), as.name))

    # Gap after last y: [max(y_end_clipped) + 1, x_end]
    # We need to find the last y for each x (max cummax_y_end per x_id)
    last_y_info <- group_by(has_overlap, !!!lapply(c("seqnames", "x_strand", "x_id", "x_start", "x_end"), as.name))
    max_y_end <- list(max_y_end = call("max", as.name("cummax_y_end")))
    last_y_info <- summarize(last_y_info, !!!max_y_end)
    last_y_info <- ungroup(last_y_info)

    gap_after <- mutate(last_y_info,
                        gap_start = !!call("+", as.name("max_y_end"), 1L),
                        gap_end = !!as.name("x_end"))
    gap_after <- filter(gap_after, !!!list(call("<=", as.name("gap_start"), as.name("gap_end"))))
    gap_after <- mutate(gap_after,
                        start = !!as.name("gap_start"),
                        end = !!as.name("gap_end"),
                        strand = !!as.name("x_strand"))
    gap_after <- select(gap_after, !!!lapply(c("seqnames", "start", "end", "strand"), as.name))

    # Combine all parts
    result <- union_all(no_overlap, gap_before)
    result <- union_all(result, gaps_between)
    result <- union_all(result, gap_after)
    result <- distinct(result)

    # Compute width
    width_mutate <- list(
        width = call("+", call("-", as.name("end"), as.name("start")), 1L)
    )
    result <- mutate(result, !!!width_mutate)

    merged_seqinfo <- merge(seqinfo(x), seqinfo(y))
    .build_DuckDBGRanges(result, merged_seqinfo)
})

#' @export
setMethod("setdiff", c("DuckDBGRanges", "GRanges"),
function(x, y, ignore.strand = FALSE)
{
    y_ddb <- .granges_to_duckdb(y, x)
    setdiff(x, y_ddb, ignore.strand = ignore.strand)
})

#' @export
setMethod("setdiff", c("GRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE)
{
    x_ddb <- .granges_to_duckdb(x, y)
    setdiff(x_ddb, y, ignore.strand = ignore.strand)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Nearest neighbor methods
###
### precede(), follow(), nearest(), and distanceToNearest() find the
### nearest range in a subject for each query range. These methods
### return integer vectors or Hits objects, requiring data collection.
###

# Helper function to add keycol-based row indices to a connection
# This ensures row indices match the iteration order of DuckDBGRanges
#' @importFrom DuckDBDataFrame has_row_number
#' @importFrom dbplyr window_order
#' @importFrom dplyr mutate
#' @importFrom stats setNames
.add_keycol_indices <- function(conn, frame, db_conn, idx_name, prefix = "") {
    keycol_name <- names(frame@keycols)[1L]
    keycol_vals <- frame@keycols[[1L]]

    if (has_row_number(frame)) {
        rank_expr <- call("row_number", as.name(keycol_name))
        return(mutate(conn, !!!setNames(list(rank_expr), idx_name)))
    }

    # Create mapping dataframe
    mapping_df <- data.frame(
        keycol_val = keycol_vals,
        idx = seq_along(keycol_vals),
        stringsAsFactors = FALSE
    )
    names(mapping_df) <- c(keycol_name, idx_name)
    
    # Copy mapping to DuckDB
    temp_name <- paste0(prefix, "_keycol_map_", sample(1e6, 1))
    mapping_tbl <- copy_to(db_conn, mapping_df, temp_name, overwrite = TRUE)
    
    # Join to add the index
    joined <- inner_join(conn, mapping_tbl, by = keycol_name)
    joined
}

# Helper function for nearest neighbor methods (precede, follow, nearest, distanceToNearest)
# Consolidates common setup: connection handling, subject renaming, joining
#
# @param x DuckDBGRanges query object
# @param subject DuckDBGRanges subject object
# @param ignore.strand Whether to ignore strand
# @return List with: joined (connection), n_x (length of x), n_subj (length of subject)
#' @importFrom DuckDBDataFrame tblconn
.setup_nearest_neighbor_join <- function(x, subject, ignore.strand) {
    # Get connections and database connection
    x_frame <- x@frame
    x_conn <- tblconn(x_frame)
    db_conn <- dbconn(x_frame)
    
    subj_frame <- subject@frame
    subj_conn <- tblconn(subj_frame)
    
    # Add keycol-based indices for proper row ordering
    x_conn <- .add_keycol_indices(x_conn, x_frame, db_conn, "x_idx", "x")
    subj_conn <- .add_keycol_indices(subj_conn, subj_frame, db_conn, "subj_idx", "subj")
    
    # Rename subject columns to avoid conflicts
    subj_rename <- list(
        subj_seqnames = as.name("seqnames"),
        subj_start = as.name("start"),
        subj_end = as.name("end"),
        subj_strand = as.name("strand")
    )
    subj_conn <- mutate(subj_conn, !!!subj_rename)
    subj_select <- lapply(c("subj_idx", "subj_seqnames", "subj_start", "subj_end", "subj_strand"), as.name)
    subj_conn <- select(subj_conn, !!!subj_select)
    
    # Join on seqnames only, then (unless ignore.strand) apply a strand
    # COMPATIBILITY filter rather than an equi-join on strand. Base
    # GenomicRanges treats '*' as matching any strand, so '+' pairs with {+,*},
    # '-' with {-,*}, and '*' with all -- an equi-join on strand would silently
    # drop every '*' pair. Joining on seqnames alone also keeps 'subj_strand' as
    # a column, which precede()/follow() need to choose the transcription
    # direction per row.
    joined <- left_join(x_conn, subj_conn,
                        by = c("seqnames" = "subj_seqnames"),
                        suffix = c("", "_y"))
    if (!ignore.strand) {
        joined <- filter(joined, !!!list(
            call("|",
                 call("|",
                      call("==", as.name("strand"), "*"),
                      call("==", as.name("subj_strand"), "*")),
                 call("==", as.name("strand"), as.name("subj_strand")))))
    }

    list(joined = joined, n_x = length(x), n_subj = length(subject))
}

# Helper to return empty results for nearest neighbor methods
#' @importFrom S4Vectors Hits
.empty_nearest_result <- function(n_x, n_subj, select, select_all_value) {
    if (select == select_all_value) {
        return(Hits(nLnode = n_x, nRnode = n_subj))
    }
    if (n_x == 0L) {
        return(integer(0))
    }
    return(rep(NA_integer_, n_x))
}

# Helper to build result for nearest neighbor methods with single selection
.build_nearest_single_result <- function(result_df, n_x) {
    out <- rep(NA_integer_, n_x)
    if (nrow(result_df) > 0) {
        out[as.integer(result_df$x_idx)] <- as.integer(result_df$subj_idx)
    }
    out
}

# Helper to build Hits result for nearest neighbor methods with "all" selection
#' @importFrom S4Vectors Hits
.build_nearest_hits_result <- function(result_df, n_x, n_subj) {
    if (nrow(result_df) == 0) {
        return(Hits(nLnode = n_x, nRnode = n_subj))
    }
    Hits(
        from = as.integer(result_df$x_idx),
        to = as.integer(result_df$subj_idx),
        nLnode = n_x,
        nRnode = n_subj
    )
}

# Core nearest-neighbour engine shared by nearest(x, subject) and nearest(x).
# `drop.self = TRUE` excludes a range from being its own neighbour, matching base
# GenomicRanges' nearest(x, missing) (which calls .nearest(x, x, drop.self=TRUE)).
# Self must be removed BEFORE the minimum-distance step so a range with a self
# overlap still resolves to its next-nearest neighbour rather than to itself.
#' @importFrom dplyr mutate filter select arrange group_by summarize collect inner_join distinct ungroup
.nearest_ddb <- function(x, subject, select, ignore.strand, drop.self = FALSE) {
    n_x <- length(x)
    n_subj <- length(subject)
    if (n_x == 0L || n_subj == 0L) {
        return(.empty_nearest_result(n_x, n_subj, select, "all"))
    }

    setup <- .setup_nearest_neighbor_join(x, subject, ignore.strand)
    joined <- setup$joined

    if (drop.self) {
        joined <- filter(joined,
                         !!!list(call("!=",
                                      as.name("x_idx"),
                                      as.name("subj_idx"))))
    }

    # Compute distance: max(subj_start - x_end - 1, x_start - subj_end - 1, 0);
    # overlapping ranges have distance 0.
    dist_expr <- list(
        dist = call("greatest",
                    call("-", call("-", as.name("subj_start"), as.name("end")), 1L),
                    call("-", call("-", as.name("start"), as.name("subj_end")), 1L),
                    0L)
    )
    joined <- mutate(joined, !!!dist_expr)
    joined <- filter(joined,
                     !!!list(call("!", call("is.na", as.name("subj_idx")))))
    joined <- group_by(joined, !!!list(as.name("x_idx")))
    min_dists <- summarize(joined,
                           !!!list(min_dist = call("min", as.name("dist"), na.rm = TRUE)))
    min_dists <- ungroup(min_dists)

    joined <- ungroup(joined)
    joined <- inner_join(joined, min_dists, by = "x_idx", suffix = c("", "_min"))
    joined <- filter(joined, !!!list(call("==", as.name("dist"), as.name("min_dist"))))

    if (select == "arbitrary") {
        joined <- group_by(joined, !!!list(as.name("x_idx")))
        result <- summarize(joined,
                            !!!list(subj_idx = call("min", as.name("subj_idx"), na.rm = TRUE)))
        result <- ungroup(result)
        return(.build_nearest_single_result(collect(result), n_x))
    } else {
        result <- select(joined, !!!lapply(c("x_idx", "subj_idx"), as.name))
        result <- distinct(result)
        result <- arrange(result, !!!lapply(c("x_idx", "subj_idx"), as.name))
        return(.build_nearest_hits_result(collect(result), n_x, n_subj))
    }
}

# Core distanceToNearest engine shared by the (x, subject) and (x) forms; see
# .nearest_ddb for the drop.self rationale (base uses drop.self=TRUE for x alone).
#' @importFrom dplyr mutate filter select arrange group_by summarize collect inner_join distinct ungroup
#' @importFrom S4Vectors DataFrame Hits mcols<-
.distanceToNearest_ddb <- function(x, subject, ignore.strand, select, drop.self = FALSE) {
    n_x <- length(x)
    n_subj <- length(subject)
    if (n_x == 0L || n_subj == 0L) {
        hits <- Hits(nLnode = n_x, nRnode = n_subj)
        mcols(hits) <- DataFrame(distance = integer(0))
        return(hits)
    }

    setup <- .setup_nearest_neighbor_join(x, subject, ignore.strand)
    joined <- setup$joined

    if (drop.self) {
        joined <- filter(joined,
                         !!!list(call("!=", as.name("x_idx"), as.name("subj_idx"))))
    }

    dist_expr <- list(
        dist = call("greatest",
                    call("-", call("-", as.name("subj_start"), as.name("end")), 1L),
                    call("-", call("-", as.name("start"), as.name("subj_end")), 1L),
                    0L)
    )
    joined <- mutate(joined, !!!dist_expr)

    # Filter out rows with NULL subj_idx (no match on seqnames/strand)
    joined <- filter(joined,
                     !!!list(call("!", call("is.na", as.name("subj_idx")))))

    joined <- group_by(joined, !!!list(as.name("x_idx")))
    min_dists <- summarize(joined,
                           !!!list(min_dist = call("min", as.name("dist"), na.rm = TRUE)))
    min_dists <- ungroup(min_dists)

    joined <- ungroup(joined)
    joined <- inner_join(joined, min_dists, by = "x_idx", suffix = c("", "_min"))
    joined <- filter(joined,
                     !!!list(call("==", as.name("dist"), as.name("min_dist"))))

    if (select == "arbitrary") {
        joined <- group_by(joined, !!!list(as.name("x_idx")))
        result <- summarize(joined, !!!list(
            subj_idx = call("min", as.name("subj_idx"), na.rm = TRUE),
            distance = call("min", as.name("dist"), na.rm = TRUE)
        ))
        result <- ungroup(result)
    } else {
        result <- select(joined,
                         !!!lapply(c("x_idx", "subj_idx", "dist"), as.name))
        result <- mutate(result, distance = as.name("dist"))
        result <- select(result,
                         !!!lapply(c("x_idx", "subj_idx", "distance"), as.name))
        result <- distinct(result)
    }

    result <- arrange(result, !!!lapply(c("x_idx", "subj_idx"), as.name))
    result_df <- collect(result)

    if (nrow(result_df) == 0) {
        hits <- Hits(nLnode = n_x, nRnode = n_subj)
        mcols(hits) <- DataFrame(distance = integer(0))
        return(hits)
    }

    hits <- Hits(from = as.integer(result_df$x_idx),
                 to = as.integer(result_df$subj_idx),
                 nLnode = n_x,
                 nRnode = n_subj)
    mcols(hits) <- DataFrame(distance = as.integer(result_df$distance))
    hits
}

#' @export
#' @importFrom IRanges precede
#' @importFrom dplyr mutate filter select arrange group_by summarize collect left_join inner_join copy_to distinct ungroup if_else
#' @importFrom S4Vectors isTRUEorFALSE
setMethod("precede", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, subject, select = c("first", "all"), ignore.strand = FALSE)
{
    select <- match.arg(select)
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    
    n_x <- length(x)
    n_subj <- length(subject)
    
    # Handle empty cases
    if (n_x == 0L || n_subj == 0L) {
        return(.empty_nearest_result(n_x, n_subj, select, "all"))
    }
    
    # Setup join
    setup <- .setup_nearest_neighbor_join(x, subject, ignore.strand)
    joined <- setup$joined

    # precede is strand-DIRECTIONAL. Base picks the convention per row by the
    # QUERY strand (and, for a '*' query, by the SUBJECT strand): for '+'/'*' the
    # preceded subject lies to the right (subj_start > x_end); for '-' it lies to
    # the left (subj_end < x_start). ignore.strand collapses everything to '+'.
    # `pf_dist` is the (positive) gap in the transcription direction, so the
    # nearest is min(pf_dist) whichever side applies.
    if (ignore.strand) {
        valid_expr <- call(">", as.name("subj_start"), as.name("end"))
        dist_expr  <- call("-", as.name("subj_start"), as.name("end"))
    } else {
        use_minus <- call("|",
                          call("==", as.name("strand"), "-"),
                          call("&",
                               call("==", as.name("strand"), "*"),
                               call("==", as.name("subj_strand"), "-")))
        valid_expr <- call("if_else", use_minus,
                           call("<", as.name("subj_end"), as.name("start")),
                           call(">", as.name("subj_start"), as.name("end")))
        dist_expr  <- call("if_else", use_minus,
                           call("-", as.name("start"), as.name("subj_end")),
                           call("-", as.name("subj_start"), as.name("end")))
    }
    joined <- mutate(joined, !!!list(pf_valid = valid_expr, pf_dist = dist_expr))
    joined <- filter(joined, !!!list(as.name("pf_valid")))

    # Group by x_idx
    joined <- group_by(joined, !!!list(as.name("x_idx")))

    if (select == "first") {
        # Rank by transcription-direction distance (nearest first)
        joined <- mutate(joined, !!!list(rk = call("dense_rank", as.name("pf_dist"))))
        joined <- ungroup(joined)
        joined <- filter(joined, !!!list(call("==", as.name("rk"), 1L)))

        # For ties, take minimum subj_idx (base select = "first")
        joined <- group_by(joined, !!!list(as.name("x_idx")))
        result <- summarize(joined, !!!list(subj_idx = call("min", as.name("subj_idx"), na.rm = TRUE)))
        result <- ungroup(result)

        return(.build_nearest_single_result(collect(result), n_x))
    } else {
        # select == "all": base returns every subject at the NEAREST distance
        # (the ties), not all directional subjects — so filter to the minimum
        # pf_dist per query first (mirrors .nearest_ddb/.distanceToNearest_ddb).
        joined <- mutate(joined, !!!list(rk = call("dense_rank", as.name("pf_dist"))))
        joined <- ungroup(joined)
        joined <- filter(joined, !!!list(call("==", as.name("rk"), 1L)))
        result <- select(joined, !!!lapply(c("x_idx", "subj_idx"), as.name))
        result <- distinct(result)
        result <- arrange(result, !!!lapply(c("x_idx", "subj_idx"), as.name))

        return(.build_nearest_hits_result(collect(result), n_x, n_subj))
    }
})

#' @export
setMethod("precede", c("DuckDBGRanges", "missing"),
function(x, subject, select = c("first", "all"), ignore.strand = FALSE)
{
    precede(x, x, select = select, ignore.strand = ignore.strand)
})

#' @export
setMethod("precede", c("DuckDBGRanges", "GRanges"),
function(x, subject, select = c("first", "all"), ignore.strand = FALSE)
{
    subject_ddb <- .granges_to_duckdb(subject, x)
    precede(x, subject_ddb, select = select, ignore.strand = ignore.strand)
})

#' @export
setMethod("precede", c("GRanges", "DuckDBGRanges"),
function(x, subject, select = c("first", "all"), ignore.strand = FALSE)
{
    x_ddb <- .granges_to_duckdb(x, subject)
    precede(x_ddb, subject, select = select, ignore.strand = ignore.strand)
})

#' @export
#' @importFrom IRanges follow
#' @importFrom dplyr mutate filter select arrange group_by summarize collect left_join inner_join copy_to distinct ungroup desc if_else
#' @importFrom S4Vectors isTRUEorFALSE
setMethod("follow", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, subject, select = c("last", "all"), ignore.strand = FALSE)
{
    select <- match.arg(select)
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    
    n_x <- length(x)
    n_subj <- length(subject)
    
    # Handle empty cases
    if (n_x == 0L || n_subj == 0L) {
        return(.empty_nearest_result(n_x, n_subj, select, "all"))
    }
    
    # Setup join
    setup <- .setup_nearest_neighbor_join(x, subject, ignore.strand)
    joined <- setup$joined

    # follow is strand-DIRECTIONAL (the mirror of precede). Base picks the
    # convention per row by the QUERY strand (and, for a '*' query, by the SUBJECT
    # strand): for '+'/'*' the followed subject lies to the left (subj_end <
    # x_start); for '-' it lies to the right (subj_start > x_end). ignore.strand
    # collapses everything to '+'. `pf_dist` is the positive gap in the
    # transcription direction, so the nearest is min(pf_dist).
    if (ignore.strand) {
        valid_expr <- call("<", as.name("subj_end"), as.name("start"))
        dist_expr  <- call("-", as.name("start"), as.name("subj_end"))
    } else {
        use_minus <- call("|",
                          call("==", as.name("strand"), "-"),
                          call("&",
                               call("==", as.name("strand"), "*"),
                               call("==", as.name("subj_strand"), "-")))
        valid_expr <- call("if_else", use_minus,
                           call(">", as.name("subj_start"), as.name("end")),
                           call("<", as.name("subj_end"), as.name("start")))
        dist_expr  <- call("if_else", use_minus,
                           call("-", as.name("subj_start"), as.name("end")),
                           call("-", as.name("start"), as.name("subj_end")))
    }
    joined <- mutate(joined, !!!list(pf_valid = valid_expr, pf_dist = dist_expr))
    joined <- filter(joined, !!!list(as.name("pf_valid")))

    # Group by x_idx
    joined <- group_by(joined, !!!list(as.name("x_idx")))

    if (select == "last") {
        # Rank by transcription-direction distance (nearest first)
        joined <- mutate(joined, !!!list(rk = call("dense_rank", as.name("pf_dist"))))
        joined <- ungroup(joined)
        joined <- filter(joined, !!!list(call("==", as.name("rk"), 1L)))

        # For ties, take maximum subj_idx (base select = "last")
        joined <- group_by(joined, !!!list(as.name("x_idx")))
        result <- summarize(joined, !!!list(subj_idx = call("max", as.name("subj_idx"), na.rm = TRUE)))
        result <- ungroup(result)

        return(.build_nearest_single_result(collect(result), n_x))
    } else {
        # select == "all": base returns every subject at the NEAREST distance
        # (the ties), not all directional subjects — so filter to the minimum
        # pf_dist per query first (mirrors .nearest_ddb/.distanceToNearest_ddb).
        joined <- mutate(joined, !!!list(rk = call("dense_rank", as.name("pf_dist"))))
        joined <- ungroup(joined)
        joined <- filter(joined, !!!list(call("==", as.name("rk"), 1L)))
        result <- select(joined, !!!lapply(c("x_idx", "subj_idx"), as.name))
        result <- distinct(result)
        result <- arrange(result, !!!lapply(c("x_idx", "subj_idx"), as.name))

        return(.build_nearest_hits_result(collect(result), n_x, n_subj))
    }
})

#' @export
setMethod("follow", c("DuckDBGRanges", "missing"),
function(x, subject, select = c("last", "all"), ignore.strand = FALSE)
{
    follow(x, x, select = select, ignore.strand = ignore.strand)
})

#' @export
setMethod("follow", c("DuckDBGRanges", "GRanges"),
function(x, subject, select = c("last", "all"), ignore.strand = FALSE)
{
    subject_ddb <- .granges_to_duckdb(subject, x)
    follow(x, subject_ddb, select = select, ignore.strand = ignore.strand)
})

#' @export
setMethod("follow", c("GRanges", "DuckDBGRanges"),
function(x, subject, select = c("last", "all"), ignore.strand = FALSE)
{
    x_ddb <- .granges_to_duckdb(x, subject)
    follow(x_ddb, subject, select = select, ignore.strand = ignore.strand)
})

#' @export
#' @importFrom IRanges nearest
#' @importFrom dplyr mutate filter select arrange group_by summarize collect left_join inner_join copy_to distinct ungroup
#' @importFrom S4Vectors isTRUEorFALSE
setMethod("nearest", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)
{
    select <- match.arg(select)
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    # An explicit subject keeps self-matches (matches base GenomicRanges).
    .nearest_ddb(x, subject, select, ignore.strand, drop.self = FALSE)
})

#' @export
setMethod("nearest", c("DuckDBGRanges", "missing"),
function(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)
{
    select <- match.arg(select)
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    # No subject: a range is never its own nearest neighbour (base uses
    # drop.self=TRUE for nearest(x, missing)).
    .nearest_ddb(x, x, select, ignore.strand, drop.self = TRUE)
})

#' @export
setMethod("nearest", c("DuckDBGRanges", "GRanges"),
function(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)
{
    subject_ddb <- .granges_to_duckdb(subject, x)
    nearest(x, subject_ddb, select = select, ignore.strand = ignore.strand)
})

#' @export
setMethod("nearest", c("GRanges", "DuckDBGRanges"),
function(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)
{
    x_ddb <- .granges_to_duckdb(x, subject)
    nearest(x_ddb, subject, select = select, ignore.strand = ignore.strand)
})

#' @export
#' @importFrom IRanges distanceToNearest
#' @importFrom dplyr mutate filter select arrange group_by summarize collect left_join inner_join copy_to distinct ungroup
#' @importFrom S4Vectors DataFrame Hits isTRUEorFALSE mcols<-
setMethod("distanceToNearest", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, subject, ignore.strand = FALSE, select = c("arbitrary", "all"))
{
    select <- match.arg(select)
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    # An explicit subject keeps self-matches (matches base GenomicRanges).
    .distanceToNearest_ddb(x, subject, ignore.strand, select, drop.self = FALSE)
})

#' @export
setMethod("distanceToNearest", c("DuckDBGRanges", "missing"),
function(x, subject, ignore.strand = FALSE, select = c("arbitrary", "all"))
{
    select <- match.arg(select)
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    # No subject: exclude self-hits (base uses drop.self=TRUE for x alone).
    .distanceToNearest_ddb(x, x, ignore.strand, select, drop.self = TRUE)
})

#' @export
setMethod("distanceToNearest", c("DuckDBGRanges", "GRanges"),
function(x, subject, ignore.strand = FALSE, select = c("arbitrary", "all"))
{
    subject_ddb <- .granges_to_duckdb(subject, x)
    distanceToNearest(x, subject_ddb, ignore.strand = ignore.strand, select = select)
})

#' @export
setMethod("distanceToNearest", c("GRanges", "DuckDBGRanges"),
function(x, subject, ignore.strand = FALSE, select = c("arbitrary", "all"))
{
    x_ddb <- .granges_to_duckdb(x, subject)
    distanceToNearest(x_ddb, subject, ignore.strand = ignore.strand, select = select)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tiling methods
###
### tile() and slidingWindows() divide ranges into sub-ranges.
### These return GRangesList objects since they create multiple output
### ranges per input range.
###

#' @export
#' @importFrom IRanges tile
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom S4Vectors isTRUEorFALSE isSingleNumber
setMethod("tile", "DuckDBGRanges",
function(x, n, width)
{
    if (length(x) == 0L)
        return(GRangesList())
    
    # Validate arguments
    has_n <- !missing(n) && !is.null(n)
    has_width <- !missing(width) && !is.null(width)
    
    if (!has_n && !has_width)
        stop("either 'n' or 'width' must be specified")
    if (has_n && has_width)
        stop("only one of 'n' or 'width' can be specified")
    
    if (has_n) {
        if (!isSingleNumber(n) || n < 1)
            stop("'n' must be a single positive integer")
        n <- as.integer(n)
    }
    if (has_width) {
        if (!isSingleNumber(width) || width < 1)
            stop("'width' must be a single positive integer")
        width <- as.integer(width)
    }
    
    # For DuckDBGRanges, materialize and use GRanges tile implementation
    # This is a pragmatic approach since tile creates a 1-to-many relationship
    # that would require complex SQL generation
    gr <- as(x, "GRanges")
    if (has_n) {
        tile(gr, n = n)
    } else {
        tile(gr, width = width)
    }
})

#' @export
#' @importFrom IRanges slidingWindows
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom S4Vectors isSingleNumber
setMethod("slidingWindows", "DuckDBGRanges",
function(x, width, step = 1L)
{
    if (length(x) == 0L)
        return(GRangesList())
    
    # Validate arguments
    if (!isSingleNumber(width) || width < 1)
        stop("'width' must be a single positive integer")
    if (!isSingleNumber(step) || step < 1)
        stop("'step' must be a single positive integer")
    
    width <- as.integer(width)
    step <- as.integer(step)
    
    # For DuckDBGRanges, materialize and use GRanges implementation
    # This creates a 1-to-many relationship best handled by GRanges
    gr <- as(x, "GRanges")
    slidingWindows(gr, width = width, step = step)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pgap method
###
### pgap() computes pairwise gaps between ranges.
###

#' @export
#' @importFrom DuckDBDataFrame dbconn has_row_number tblconn
#' @importFrom IRanges pgap
#' @importFrom S4Vectors isTRUEorFALSE
#' @importFrom dplyr inner_join mutate pull select summarize
setMethod("pgap", c("DuckDBGRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE, ...)
{
    if (length(x) != length(y))
        stop("'x' and 'y' must have the same length")

    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    if (length(x) == 0L)
        return(x)

    if (has_row_number(x@frame) || has_row_number(y@frame))
        stop("pgap() requires 'x' and 'y' to have an explicit keycol; ",
             "row-number-keyed DuckDBGRanges are not supported")

    x_frame <- x@frame
    x_conn <- tblconn(x_frame)
    db_conn <- dbconn(x_frame)

    y_frame <- y@frame
    y_conn <- tblconn(y_frame)

    x_conn <- .add_keycol_indices(x_conn, x_frame, db_conn, ".row_idx", "pgap_x")
    y_conn <- .add_keycol_indices(y_conn, y_frame, db_conn, ".row_idx", "pgap_y")
    y_select_list <- setNames(
        lapply(c(".row_idx", "seqnames", "start", "end", "strand"), as.name),
        c(".row_idx", "y_seqnames", "y_start", "y_end", "y_strand"))
    y_conn <- select(y_conn, !!!y_select_list)

    joined <- inner_join(x_conn, y_conn, by = ".row_idx", copy = TRUE)

    # pgap()'s gap formula (matching IRanges:::pgap,IntegerRanges,IntegerRanges):
    # new_start = min(end(x), end(y)) + 1; the naive gap end, max(start(x),
    # start(y)) - 1, is clamped up to (new_start - 1) so an overlapping or
    # adjacent pair collapses to a zero-width range at the boundary instead of
    # a negative-width one.
    new_start_expr <- call("+", call("least", as.name("end"), as.name("y_end")), 1L)
    naive_end_expr <- call("-", call("greatest", as.name("start"), as.name("y_start")), 1L)
    new_end_expr <- call("greatest", naive_end_expr,
                         call("least", as.name("end"), as.name("y_end")))

    # pgap() requires x[i]/y[i] to have compatible seqnames (equal) and, unless
    # ignore.strand, compatible strand ('*' matches any strand). The strand
    # sub-expression is wrapped in an explicit call("(", .) because dbplyr's
    # SQL translation does not parenthesize a '|' operand of '&' on its own,
    # which would otherwise silently change AND/OR precedence in the rendered
    # SQL (verified against dbplyr::sql_render()).
    compat_expr <- call("==", as.name("seqnames"), as.name("y_seqnames"))
    if (!ignore.strand) {
        strand_ok <- call("(", call("|",
            call("|", call("==", as.name("strand"), "*"),
                     call("==", as.name("y_strand"), "*")),
            call("==", as.name("strand"), as.name("y_strand"))))
        compat_expr <- call("&", compat_expr, strand_ok)
    }

    coord_mutate <- list(new_start = new_start_expr, new_end = new_end_expr,
                         compatible = compat_expr)
    joined <- mutate(joined, !!!coord_mutate)

    n_bad_expr <- call("sum", call("as.integer", call("!", as.name("compatible"))), na.rm = TRUE)
    n_bad <- pull(summarize(joined, !!!setNames(list(n_bad_expr), "n_bad")))
    if (isTRUE(n_bad > 0L))
        stop("'x' and 'y' elements must have compatible 'seqnames' and 'strand' values")

    new_width_expr <- call("+", call("-", as.name("new_end"), as.name("new_start")), 1L)
    joined <- mutate(joined, new_width = !!new_width_expr)

    # Retain .row_idx so the result can be ordered to match the original
    # x[i]/y[i] pairing; .build_DuckDBGRanges()'s datacols expression below
    # doesn't name it, so it never becomes a visible column of the result.
    select_rename_list <- setNames(
        lapply(c(".row_idx", "seqnames", "strand", "new_start", "new_end", "new_width"), as.name),
        c(".row_idx", "seqnames", "strand", "start", "end", "width"))
    joined <- select(joined, !!!select_rename_list)

    # pgap() is positional (result[i] is the gap of x[i]/y[i]): order by the
    # pairing index rather than .build_DuckDBGRanges()'s default coordinate
    # sort, which would silently permute the result relative to x/y.
    .build_DuckDBGRanges(joined, seqinfo(x), order_by = list(as.name(".row_idx")))
})

#' @export
setMethod("pgap", c("DuckDBGRanges", "GRanges"),
function(x, y, ignore.strand = FALSE, ...)
{
    y_ddb <- .granges_to_duckdb(y, x)
    pgap(x, y_ddb, ignore.strand = ignore.strand, ...)
})

#' @export
setMethod("pgap", c("GRanges", "DuckDBGRanges"),
function(x, y, ignore.strand = FALSE, ...)
{
    x_ddb <- .granges_to_duckdb(x, y)
    pgap(x_ddb, y, ignore.strand = ignore.strand, ...)
})
