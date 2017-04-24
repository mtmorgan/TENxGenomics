#' @rdname tenxreduce
#'
#' @title Visit and 'reduce()' all values in a 10xGenomics hdf5 file.
#'
#' @description This function visits all records in the hdf5 file,
#'     passing a list of value, row index, and column index to a
#'     user-provided \code{f()}. \code{f()} reduces the arguments in
#'     an arbitary way, returning the result. All records in the hdf5
#'     file are visited, but only those correspond to rows or columns
#'     of \code{x} are reported. This is modeled after the base
#'     function \code{Reduce()}.
#'
#' @param f A function of two arguments. The first argument \code{x}
#'     is a list with elements \code{ridx} (row index), \code{cidx}
#'     (column index) and \code{value} (value). The second argument
#'     \code{y} is \code{init} (on first iteration) or the return
#'     value of the previous call to \code{f()} (all subsequent
#'     iterations).
#'
#' @param x A \code{TENxGenomics-class} instance. Rows and columns of
#'     \code{x} cannot be duplicated.
#'
#' @param init Initial value, used as \code{y} on first call to
#'     \code{f()}.
#'
#' @param size numeric(1) Number of columns to process during each
#'     iteration.
#'
#' @param verbose logical(1) when TRUE, print a progress bar.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
tenxreduce <-
    function(f, x, init = NULL, size = 10000, verbose = interactive())
{
    stopifnot(
        is(x, "TENxGenomics"),
        !anyDuplicated(.rowidx(x)), !anyDuplicated(.colidx(x)),
        is.logical(verbose), length(verbose) == 1L, !is.na(verbose),
        is.numeric(size), length(size) == 1L, !is.na(size), size > 0
    )

    h5f <- H5Fopen(.h5path(x))
    on.exit(H5Fclose(h5f))

    indptr <- h5read(h5f, .indptr(x), bit64conversion = "double") + 1
    datalen <- length(indptr)

    ## initialize
    result <- init
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
        if (start >= datalen)
            break
        block <- min(start + size, datalen) - start

        offsets <- h5read(
            h5f, .indptr(x), start = start, block = block + 1, count = 1L
        ) + 1L
        start0 <- min(offsets)
        block0 <- diff(range(offsets))
        cidx <- rep(start - 1L + seq_len(block), diff(offsets))

        ridx <- h5read(
            h5f, .indices(x), start = start0, block = block0, count = 1L
        ) + 1

        value <- h5read(
            h5f, .dataname(x), start = start0, block = block0, count = 1L
        )

        keep <- (ridx %in% .rowidx(x)) & (cidx %in% .colidx(x))

        result <- f(
            list(
                ridx = as.vector(ridx[keep]),
                cidx = cidx[keep],
                value  = as.vector(value[keep])
            ),
            result
        )

        updt(start + block)
        start <- start + block
    }
    if (verbose) message()              # newline after progress bar
    result
}
