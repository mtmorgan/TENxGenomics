#' @rdname tenxreduce
#'
#' @title Visit and 'reduce()' all values in the hdf5 file.
#'
#' @description This function visits all records in the hdf5 file,
#'     passing a list of value, row index, and column index to a
#'     user-provided \code{FUN}. \code{FUN} reduces the arguments in
#'     an arbitary way, returning the result. All records in the hdf5
#'     file are visited, but only those correspond to rows or columns
#'     of \code{X} are reported.
#'
#' @param X A \code{TENxGenomics-class} instance. Rows and columns of
#'     \code{X} cannot be duplicated.
#'
#' @param FUN A function of two arguments. The first argument is a
#'     list with elements \code{ridx} (row index), \code{cidx} (column
#'     index) and \code{value} (value). The second argument is NULL
#'     (on first iteration) or the return value of the previous call
#'     to \code{FUN} (all subsequent iterations).
#'
#' @param ... Additional arguments, passed to \code{FUN}.
#'
#' @param verbose logical(1) when TRUE, print a progress bar.
#'
#' @param BLOCK numeric(1) Number of values to process during each
#'     iteration.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
tenxreduce <-
    function(X, FUN, ..., verbose = interactive(), BLOCK = 100000000)
{
    stopifnot(
        is(X, "TENxGenomics"),
        !anyDuplicated(.rowidx(X)), !anyDuplicated(.colidx(X)),
        is.logical(verbose), length(verbose) == 1L, !is.na(verbose),
        is.numeric(BLOCK), length(BLOCK) == 1L, !is.na(BLOCK), BLOCK > 0
    )

    h5f <- H5Fopen(.h5path(X))
    on.exit(H5Fclose(h5f))

    h5data <- H5Dopen(h5f, "/mm10/data")
    datalen <- H5Sget_simple_extent_dims(H5Dget_space(h5data))$size
    H5Dclose(h5data)

    indptr <- h5read(h5f, "/mm10/indptr", bit64conversion = "double") + 1

    ## initialize
    result <- NULL
    start <- 1
    if (verbose) {
        updt <- local({
            progress <- txtProgressBar(min = 0, max = datalen, style=3)
            function(...) setTxtProgressBar(progress, ...)
        })
    } else {
        updt <- function(...) {}
    }

    ## iterate
    repeat {
        block <- min(start + BLOCK, datalen) - start + 1
        if (block <= 0)
            break

        value <- h5read(
            h5f, "/mm10/data", start =start, block = block, count = 1L
        )
        ridx <- h5read(
            h5f, "/mm10/indices", start =start, block = block, count = 1L
        ) + 1
        cidx <- findInterval(seq(start, start + block - 1), indptr)
        keep <- (ridx %in% .rowidx(X)) & (cidx %in% .colidx(X))

        result <- FUN(
            list(
                ridx=as.vector(ridx[keep]),
                cidx = cidx[keep],
                value  = as.vector(value[keep])
            ),
            result
        )

        updt(start + block)
        start <- start + block
    }
    result
}
