# DuckDBGRanges 0.99.6

## Bug fixes

- `nearest()`, `precede()`, `follow()`, and `distanceToNearest()` returned
  all-`NA` (or an empty `Hits`) for any `DuckDBGRanges` that does not carry an
  explicit keycol, which is the default for an object built from a plain file.
  Three separate defects were stacked here; the first two masked the third.
  - `.add_keycol_indices()` treated the `keycols` slot as a set of literal key
    values and inner-joined the table against them. A row-number-keyed frame
    does not store key values there: it stores the `c(NA, -n)` sentinel
    returned by `set_row_number()`. The join therefore matched nothing and
    silently emptied the result rather than erroring. The row position is now
    derived directly instead. Note that the frame's own `row_number` column
    cannot be used as the index as-is, since a subset keeps the parent's
    numbering (`g[3:4]` carries `3`, `4` for what are now positions `1` and
    `2`); a row-number-keyed object always iterates in ascending `row_number`
    order, so ranking over that column yields the required `1..n` position.
  - `.build_nearest_single_result()` used the collected index directly as a
    subscript. Now that the index is computed in SQL it arrives as
    `integer64`, and subsetting with an `integer64` vector silently returns
    `NA` for every position instead of indexing. (The `Hits`-building path was
    unaffected only because it already coerced with `as.integer()`.)
  - `nearest()` lacked the `is.na()` guard on the joined subject index that
    `distanceToNearest()` already had. The setup left-joins on seqnames, so a
    query range on a seqname with no subject produces a row whose subject
    bounds are `NULL`; because DuckDB's `greatest()` *skips* nulls rather than
    propagating them, `greatest(NULL, NULL, 0)` is `0`, so that row scored a
    perfect distance, won its own minimum-distance filter, and yielded a hit
    to a `NULL` subject. The directional `precede()`/`follow()` paths were
    never affected, as their validity comparisons against a `NULL` bound
    evaluate to `NULL` and drop the row.

  Every pre-existing test in this family supplied an explicit keycol, which is
  why none of this was covered. Reported in part as item R-G2 of the
  six-architect cross-engine parity review; the other two defects were found
  while reproducing it. All four generics are now checked against
  `GenomicRanges` for both keying modes.

- `psetdiff()` handled only one of the four ways `y[i]` can overlap `x[i]`,
  as its own header comment conceded ("assumes y is fully contained overlap
  scenario or edge-aligned"). A pair that does not overlap at all came back
  **wider** than `x`, and a pair where `y` covers `x` came back with a
  **negative width**, which is not merely a wrong number but unmaterializable
  (`each range must have a non-negative width`). A `y` sitting strictly
  inside `x` cannot be expressed as a single range, so base refuses it;
  this returned nonsense instead. The four shapes now match base:
  `(y_end + 1, x_end)` for a left-edge overlap, `(x_start, y_start - 1)` for
  a right-edge one, a zero-width range at `(x_start, x_start - 1)` when `y`
  covers `x`, an unchanged `x` when the two do not overlap, and an error when
  `y` falls strictly inside.

  Three further defects in the same method:
  - `seqnames` and `strand` compatibility were ignored entirely, so ranges
    were subtracted from one another across different chromosomes. Base
    leaves `x[i]` untouched whenever the pair cannot be compared.
  - `ignore.strand` was accepted and then never used. It is now honoured, and
    validated.
  - The result was passed through `.build_DuckDBGRanges()`'s default
    coordinate sort, which silently permuted it: `result[i]` was not
    `x[i]` minus `y[i]`. As with `pgap()`, the result is now ordered by the
    pairing index.

  Pairing also moves off a bare `row_number()` over DuckDB's undefined scan
  order onto `.add_keycol_indices()`, matching the `pgap()` rewrite. That is
  separately tracked as review item R-G4 for `punion()`/`pintersect()`, but
  correct pairing is a precondition for this method being correct at all.

  Reported as item R-G1 of the six-architect cross-engine parity review. The
  pre-existing test covered only left-edge overlaps, which is the single
  shape the old code got right; all four shapes are now checked against
  `GenomicRanges` for both keying modes and both `ignore.strand` settings.

## Internal changes

- `case_when` is now declared in the `@importFrom` tags of `trim()`,
  `restrict()`, `reduce()`, and `setdiff()`, which all use it. They had been
  relying on the import declared by `psetdiff()`, which no longer needs it.

# DuckDBGRanges 0.99.5

## New features

- Added `coverage()` for `DuckDBGRanges`, computed lazily via the same
  delta-event + window-function `cumsum()` sweep-line pattern already used by
  `disjoin()`/`gaps()`: only the compact per-seqname breakpoint table is
  collected (bounded by the number of ranges, never by genome length), and a
  real `SimpleRleList` is built from it, matching
  `GenomicRanges::coverage()`'s value type (integer unless a non-integer
  `weight` is used). `shift` and `width` (per-seqlevel, truncating or
  zero-padding) are supported, including clipping a range that lands (partly
  or entirely) before position 1 after `shift`, exactly as
  `GenomicRanges::coverage()` does, rather than erroring. `weight` is
  scalar-only for now (a per-range vector or mcols-column-name weight is not
  supported, and errors clearly rather than risking a silent row-order
  misalignment).

## Bug fixes

- Previously, calling `coverage()` on a `DuckDBGRanges` object crashed with a
  C stack overflow: with no `coverage()` method of its own, it inherited
  `coverage,GenomicRanges-method`, which calls
  `split(ranges(x), seqnames(x))`, and `ranges,DuckDBGRanges-method` returns
  a `DuckDBDataFrame` rather than real `IRanges`. The new method above takes
  priority via S4 dispatch, so this no longer happens.

# DuckDBGRanges 0.99.4

## Bug fixes

- `pgap(x="DuckDBGRanges", y="DuckDBGRanges")` is rewritten as a genuine, lazy
  SQL implementation (no `collect()`). Previously the method built a SQL
  pipeline for the non-overlapping-gap case, discarded it unused, and always
  fell back to materializing both inputs and delegating to
  `IRanges::pgap()`; every call, including through the `DuckDBGRanges,GRanges`
  and `GRanges,DuckDBGRanges` methods, paid for a wasted SQL round trip and
  never stayed lazy despite otherwise participating in the package's SQL
  pushdown design. Also fixes a latent return-type inconsistency where a
  zero-length call returned a `DuckDBGRanges` while every other call returned
  a `GRanges`.
  - The gap formula (`start = min(end(x), end(y)) + 1`,
    `end = max(max(start(x), start(y)) - 1, min(end(x), end(y)))`) matches
    `IRanges:::pgap,IntegerRanges,IntegerRanges`, including its convention of
    a zero-width range at the boundary for overlapping or adjacent pairs.
  - `x[i]`/`y[i]` are paired by each object's own recorded key values via
    `.add_keycol_indices()` (the same mechanism `.setup_nearest_neighbor_join()`
    uses for `nearest()`/`precede()`/`follow()`), not a plain `row_number()`;
    a plain `row_number()` has no guaranteed correspondence to a table's
    logical row order across separate query executions and can silently
    mispair rows. Row-number-keyed (no explicit `keycol`) `DuckDBGRanges`
    inputs now error clearly rather than risk a silently wrong pairing.
  - `x[i]`/`y[i]` compatibility (matching `seqnames`, and unless
    `ignore.strand`, compatible `strand`) is validated in SQL, erroring like
    `IRanges::pgap()` on a mismatch; the new `ignore.strand` argument is
    threaded through the `DuckDBGRanges,GRanges`/`GRanges,DuckDBGRanges`
    methods as well.
  - The result is ordered to match the original `x`/`y` pairing rather than
    the lazy-builder helper's default coordinate sort, which would otherwise
    silently permute the result.
  - The result is now a lazy `DuckDBGRanges`, consistent with `range()` and
    `reduce()`, rather than a materialized `GRanges`.

## Testing

- Updated the `pgap()` tests for the `DuckDBGRanges` return type and added
  coverage for the seqnames/strand compatibility errors (including
  `ignore.strand`) and the row-number-keyed guard.

# DuckDBGRanges 0.99.3

## Follow-up review changes

- Adopted the undotted names for the shared BiocDuckDB internal generics, per the
  Bioconductor guideline that `.`-prefixed names should not be exported:
  `.keycols` -> `keycols`, `.has_row_number` -> `has_row_number`,
  `.set_row_number` -> `set_row_number`, and
  `.makePrettyCharacterMatrixForDisplay` -> `makePrettyCharacterMatrixForDisplay`.
  Updated the `importFrom`/`exportMethods` entries, the `DuckDBGRanges` and
  `DuckDBGRangesList` method definitions, and all call sites. Requires the
  renamed version of DuckDBDataFrame.
- Replaced `sapply(..., simplify = FALSE)` with `setNames(lapply(...), ...)` in
  the `mcols` normalization (preserving the element names), and `1:h` with
  `seq_len(h)` in the display method, resolving the corresponding `R CMD
  BiocCheck` NOTEs.
- Used `system2()` instead of `system()` in `inst/scripts/run_vignette_benchmarks.R`.

## Documentation

- Benchmarking vignette: wrapped the live comparison in `system.time()` so the
  runtimes are visible, and rewrote the summary section as prose.
- Replaced em dashes with commas (or parentheses) in the vignettes.

# DuckDBGRanges 0.99.2

## Documentation

- Added `URL` and `BugReports` fields to DESCRIPTION.
- Added a package-level man page (`?DuckDBGRanges`).
- Removed the redundant `library(BiocStyle)` call from the vignettes (the
  `BiocStyle::html_document` output and `::`-qualified helpers apply the style
  without loading the package). Applied suite-wide for consistency with the
  Bioconductor review of DuckDBDataFrame.

# DuckDBGRanges 0.99.1

## Documentation

- Added `\value` sections (roxygen `@return`) to the `DuckDBGRanges` and
  `DuckDBGRangesList` class and utility man pages, documenting the values
  returned by the constructors, accessors, replacement, coercion, and
  range-operation methods. Resolves the `R CMD BiocCheck` "missing \value"
  WARNING.

# DuckDBGRanges 0.9.8

## Testing and diagnostics

- Added a query-plan **regression guard** for interval overlaps (no change to
  overlap results or their current speed). A new internal `.explainQuery()`
  returns DuckDB's plan for a lazy query, and a test uses it to assert the
  overlap join stays an IEJoin / range join rather than degrading to a
  `NESTED_LOOP_JOIN` over the full cross product (which OOMs on skewed inputs).
  Interval overlap is a range join, not an ASOF (nearest-match) join. The
  overlap join was factored into an internal `.overlap_join_tbl()` helper so its
  plan can be inspected without collecting.

# DuckDBGRanges 0.9.7

## Changes

- Relicensed under the MIT License.

# DuckDBGRanges 0.9.6

## Bug fixes

- `precede()`/`follow()` with `select = "all"` now return only the subject(s) at
  the **nearest** distance (the ties), matching base GenomicRanges — previously
  they returned every directional subject. The `select = "first"`/`"last"` paths
  were already correct; the `"all"` branch was missing the min-distance filter.
- `distance()` now returns `NA` for a paired range on a **different seqname** or
  (unless `ignore.strand = TRUE`) an **incompatible strand** (`+` vs `-`; `*`
  matches any strand), matching base `GenomicRanges::distance`. Previously it
  ignored seqnames and strand and returned a bare coordinate gap. Added an
  `ignore.strand` argument.

# DuckDBGRanges 0.9.5

## Bug fixes

- `narrow()` now treats a **negative `start`** as a position counting back from
  the range end (`-1` is the last base), matching base `IRanges::narrow`
  (`solveUserSEW`) and the existing negative-`end` behaviour. Previously a
  negative `start` was applied as `start + (start - 1)`, walking left of the
  range start and producing an invalid interval. The start/end/width expressions
  are rebuilt so every sign-and-`width` combination is correct; because the
  `datacols` expressions flatten to SQL without parentheses, each is kept as a
  `column +/- scalar` (never a subtraction of a compound). Oracle-tested against
  base `narrow` across positive/negative `start`/`end` and their `width`
  combinations.
- `nearest(x)` and `distanceToNearest(x)` (no `subject`) now exclude self-hits,
  matching base GenomicRanges' `drop.self = TRUE` for the missing-subject form:
  a range is never its own nearest neighbour, so its next-nearest is returned
  instead of itself (distance 0). Previously both delegated to the two-argument
  form against `x` itself, so every range matched itself. The explicit
  `nearest(x, x)` / `distanceToNearest(x, x)` forms still keep self-hits, as in
  base. (`precede`/`follow` were already correct — their strict inequalities
  never admit a self-match.) Added a self-query test comparing to the base oracle.
- `nearest`/`precede`/`follow`/`distanceToNearest` now treat strand `"*"` as
  compatible with any strand, matching base GenomicRanges: `"+"` pairs with
  `{+,*}`, `"-"` with `{-,*}`, and `"*"` with all. The neighbour join used a
  strict equi-join on strand, which silently dropped every `"*"` pair (a `"*"`
  query against `"+"`/`"-"` subjects returned `NA`). The join is now on
  `seqnames` only, followed by a strand-compatibility filter.
- `precede()` and `follow()` are now strand-directional. Base defines them in the
  transcription direction — for a `"-"` strand query the roles of upstream and
  downstream are reversed, and for a `"*"` query the direction is chosen per
  subject strand. The previous implementation used a fixed genomic-coordinate
  direction regardless of strand, so results were inverted on the `"-"` strand.
  The convention is now selected per row (`subj_start > end` vs `subj_end <
  start`) via `use_minus = strand == "-" | (strand == "*" & subj_strand == "-")`,
  ranking by the transcription-direction gap; `ignore.strand = TRUE` collapses to
  the `"+"` convention as before. `follow()` keeps its base `select = "last"`
  (largest index) tie-break, `precede()` its `select = "first"` (smallest index).
  Added oracle tests against base GenomicRanges across `+`/`-`/`*` queries with
  mixed-strand subjects, plus tie-break tests.

# DuckDBGRanges 0.9.4

## Bug fixes

- `resize(x, width, fix = "center")` now matches base `IRanges::resize()`:
  `new_start = start + (width(x) - width) %/% 2`, anchored at the (strand-independent)
  center with an exact integer width. The previous implementation used the midpoint
  `(start + end) / 2` with true division --- yielding fractional coordinates for
  even-width ranges --- and `center + width %/% 2 - 1`, which returned `width - 1`
  for odd target widths. Added `fix = "center"` cases to the `resize` tests.

# DuckDBGRanges 0.9.3

## Documentation

- Restructured the vignettes into a user-first set, replacing the single
  internals-heavy *Architecture of the DuckDBGRanges Package*:
  - *Introduction to DuckDBGRanges* --- motivation, construction, and the common
    operations (accessors, subsetting, range operations, the filter-then-materialize
    workflow, `DuckDBGRangesList`).
  - *Benchmarking DuckDBGRanges* --- a best-effort comparison against in-memory
    `GRanges` on scATAC-seq (1M peaks) and variant (10M variants) scenarios,
    rendered from precomputed results so the vignette builds quickly.
  - *Design and extension of DuckDBGRanges* --- the five coordinate columns, class
    structure, SQL translation, and the `LIST[]` representation of grouped ranges,
    for developers.
- Added `inst/scripts/` with the offline benchmark generator
  (`run_vignette_benchmarks.R`) and the vignette table helpers
  (`make_timings_table.R`), following the `HDF5Array` performance-vignette
  precompute pattern.
- Rewrote the README.
