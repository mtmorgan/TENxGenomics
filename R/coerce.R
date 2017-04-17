##
## coerce
##

# Helper for common logic underlying .as.matrix() and .as.sparseMatrix()
.values_and_indices <-
    function(x, sparse = FALSE, ..., withDimnames=TRUE)
{
    stopifnot(
        is.logical(withDimnames), length(withDimnames) == 1L,
        !is.na(withDimnames)
    )

    h5f <- H5Fopen(.h5path(x))
    on.exit(H5Fclose(h5f))

    ## maximum index; needed when selecting last individual
    h5indptr <- H5Dopen(h5f, "/mm10/indptr")
    on.exit(H5Dclose(h5indptr), add=TRUE)
    indlen <- H5Sget_simple_extent_dims(H5Dget_space(h5indptr))$size

    ## get all rows for selected columns
    cidx <- .colidx(x)
    startidx <- h5read(
        h5f, "/mm10/indptr", list(cidx), bit64conversion = "double"
    )

    endidx <- h5read(                   # indptr contains the last index, too
        h5f, "/mm10/indptr", list(cidx + 1L),
        bit64conversion = "double"
    ) - 1L

    idx <- Map(seq, startidx, endidx, MoreArgs=list(by = 1))
    lens <- lengths(idx)
    idx <- unlist(idx) + 1L
    ridx <- h5read(
        h5f, "/mm10/indices", list(idx), bit64conversion = "double"
    ) + 1L

    ## get values for rows of interest
    keep <- ridx %in% .rowidx(x)
    idx <- idx[keep]
    values <-as.vector(h5read(h5f, "/mm10/data", index=list(idx)))
    ridx <- match(ridx, .rowidx(x))[keep]
    cidx <- rep(seq_along(cidx), lens)[keep]

    if (sparse) {
        ## It is more (memory) efficient to supply i and p, rather
        ## than i and j, to Matrix::sparseMatrix():
        ##
        ## i = ridx: length(i) == length(values)
        ## j = cidx: length(j) == length(values)
        ## p = c(startidx, endidx[length(endidx)] + 1L):
        ##     length(p) == (ncol(x) + 1)
        p <- c(0, cumsum(tabulate(cidx)))
        list(values = values, i = ridx, p = p)
    } else {
        list(values = values, ridx = ridx, cidx = cidx)
    }
}

#' @rdname TENxGenomics-class
#'
#' @param withDimnames logical(1) Include dimnames on returned matrix?
#'
#' @return \code{as.matrix(tenx)} and \code{as(tenx, "matrix")} return
#'     a matrix with dim and dimnames equal to \code{tenx}, and values
#'     the read counts overlapping corresponding genes and
#'     samples. Use \code{as.matrix(withDimnames=FALSE)} to suppress
#'     dimnames on the returned matrix. NOTE: consider the size of the
#'     matrix, \code{prod(as.numeric(dim(tenx)))} before invoking this
#'     function.
#'
#' @method as.matrix TENxGenomics
#'
#' @export
as.matrix.TENxGenomics <-
    function(x, ..., withDimnames=TRUE)
{
    values_and_indices <- .values_and_indices(
        x=x, sparse=FALSE, ..., withDimnames=withDimnames
    )
    values <- values_and_indices[["values"]]
    ridx <- values_and_indices[["ridx"]]
    cidx <- values_and_indices[["cidx"]]

    ## formulate result as matrix
    m <- matrix(
        0L, nrow(x), ncol(x),
        dimnames = if (withDimnames) dimnames(x) else list(NULL, NULL)
    )
    m[cbind(ridx, cidx)] <- values
    m
}

#' @rdname TENxGenomics-class
#'
#' @name coerce,TENxGenomics,matrix-method
#'
#' @exportMethod coerce
setAs("TENxGenomics", "matrix", function(from) as.matrix.TENxGenomics(from))

#' @rdname TENxGenomics-class
#'
#' @return \code{as.dgCMatrix(tenx)} and \code{as(tenx, "dgCMatrix")}
#'     return a sparse matrix (from the Matrix package) with dim and
#'     dimnames equal to \code{tenx}, and values the read counts
#'     overlapping corresponding genes and samples. Use
#'     \code{as.matrix(withDimnames=FALSE)} to suppress dimnames on
#'     the returned matrix.
#'
#' @export
as.dgCMatrix <-
    function(x, ..., withDimnames=TRUE)
{
    # TODO: Support withDimnames
    stopifnot(withDimnames)
    values_and_indices <- .values_and_indices(
        x, sparse = TRUE, ..., withDimnames = withDimnames
    )

    sparseMatrix(
        i = values_and_indices[["i"]],
        p = values_and_indices[["p"]],
        x = values_and_indices[["values"]],
        dims = dim(x), dimnames = dimnames(x),
        giveCsparse = TRUE
    )
}

## NOTE: This uses a dgCMatrix, a compressed, sparse, column-oriented
##     numeric (double) matrix. What we really want to use to store
##     these data is a igCMatrix, a compressed, sparse,
##     column-oriented *integer* matrix.  However, the igCMatrix
##     class, while defined in the Matrix package, is not actually
##     implemented.' @rdname TENxGenomics-class
##
#' @rdname TENxGenomics-class
#'
#' @name coerce,TENxGenomics,dgCMatrix-method
#'
#' @exportMethod coerce
setAs("TENxGenomics", "dgCMatrix", function(from) as.dgCMatrix(from))
