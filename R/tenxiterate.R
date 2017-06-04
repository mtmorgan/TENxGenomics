#' @rdname tenxiterate
#'
#' @title Visit and 'reduce()' all values in a 10xGenomics hdf5 file.
#'
#' This function visits all records in the hdf5 file, passing a
#' 'chunk' containing a list of value, row index, and column index to
#' a user-provided \code{FUN()}. The chunks represent one or more
#' complete cells (columns); the number of cells is controled by
#' \code{yieldSize}. \code{FUN()} manipulates the arguments (typically
#' a data reduction) in an arbitary way, returning the result.
#'
#' All records in the hdf5 file are visited, but only those correspond
#' to rows or columns of \code{X} are reported. The implementation
#' uses [BiocParallel::bpiterate()].
#'
#' @param X A \code{TENxGenomics-class} instance. Rows and columns of
#'     \code{X} cannot be duplicated.
#'
#' @param FUN A function of at least one argument \code{VALUE}, used
#'     to summarize successive chunks of the hdf5 file. \code{VALUE}
#'     is a list with elements \code{ridx} (row index), \code{cidx}
#'     (column index) and \code{value} (count value). See
#'     [BiocParallel::bpiterate()].
#' 
#' @param ... Additional arguments, passed to \code{FUN()}.
#'
#' @param yieldSize numeric(1) Number of columns (cells) to process
#'     during each iteration.
#'
#' @importFrom BiocParallel bpiterate
#'
#' @export
tenxiterate <-
    function(X, FUN, ..., yieldSize = 10000)
{
    stopifnot(
        is(X, "TENxGenomics"),
        !anyDuplicated(.rowidx(X)), !anyDuplicated(.colidx(X)),
        is.function(FUN), length(formals(FUN)) > 0,
        is.numeric(yieldSize), length(yieldSize) == 1L, !is.na(yieldSize),
            yieldSize > 0
    )

    h5f <- H5Fopen(.h5path(X))
    h5indptr <- H5Dopen(h5f, .indptr(X))

    datalen <- H5Sget_simple_extent_dims(H5Dget_space(h5indptr))$size

    H5Dclose(h5indptr)
    H5Fclose(h5f)

    ITER <- local({
        X <- X
        start <- 1
        datalen <- datalen
        yieldSize <- yieldSize

        function() {
            if (start >= datalen)
                return(NULL)
            block <- min(start + yieldSize, datalen) - start
            result <- list(X=X, start=start, block=block)
            start <<- start + block
            result
        }
    })

    .USER_MAP <- FUN
    bpiterate(
        ITER, FUN = .MAP, .YIELD_MAP = .YIELD_MAP, .USER_MAP = .USER_MAP, ...
    )
}

## forward yield coordinates (*not* data) to worker
.YIELD_MAP <- function(X, start, block) {
    h5f <- H5Fopen(.h5path(X))
    on.exit(H5Fclose(h5f))

    offsets <- h5read(
        h5f, .indptr(X), start = start, block = block + 1, count = 1L,
        bit64conversion = "double"
    ) + 1
    start0 <- min(offsets)
    block0 <- diff(range(offsets))
    cidx <- rep(start - 1 + seq_len(block), diff(offsets))

    ridx <- h5read(
        h5f, .indices(X), start = start0, block = block0, count = 1L
    ) + 1

    value <- h5read(
        h5f, .dataname(X), start = start0, block = block0, count = 1L
    )

    keep <- (ridx %in% .rowidx(X)) & (cidx %in% .colidx(X))

    list(
        ridx = as.vector(ridx[keep]),
        cidx = cidx[keep],
        value  = as.vector(value[keep])
    )
}

.MAP <-
    function(YIELD, .YIELD_MAP, .USER_MAP, ...)
{
    x <- do.call(.YIELD_MAP, YIELD)     # workers read data, then...
    .USER_MAP(x, ...)                   # do work
}
