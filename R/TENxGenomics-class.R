#' @import methods
#' @import rhdf5
#' @import Matrix

.TENxGenomics <- setClass(
    "TENxGenomics",
    slots = c(
        h5path = "character",
        colidx = "integer",
        rowidx = "integer",
        ## internal
        dataname = "character",
        rowname = "character",
        colname = "character"
    ),
    prototype = prototype(
        dataname = "/mm10/data",
        rowname = "/mm10/genes",
        colname = "/mm10/barcodes"
    )
)

.h5path <- function(x) x@h5path

.colidx <- function(x) x@colidx

.rowidx <- function(x) x@rowidx

.dataname <- function(x) x@dataname

.rowname <- function(x) x@rowname

.colname <- function(x) x@colname

.h5_dimidx <-
    function(h5f, name)
{
    h5d <- H5Dopen(h5f, name)
    on.exit(H5Dclose(h5d))
    n <- H5Sget_simple_extent_dims(H5Dget_space(h5d))$size
    seq_len(n)
}

#' @rdname TENxGenomics-class
#'
#' @title Create and manipulate a reference to 10xGenomics data.
#'
#' @description The TENxGenomics class provides a simple interface to
#'     subset and input 10xGenomis'
#'     1M_neurons_filtered_gene_bc_matrices_h5.h5 file. Subsetting is
#'     a light-weight operation; input (typically of the subset
#'     matrix) is as a dense matrix and hence consumes memory.
#'
#' @param h5path character(1) file path to the
#'     1M_neurons_filtered_gene_bc_matrices_h5.h5 file.
#'
#' @return \code{TENxGenomics()} returns a \code{TENxGenomics-class}
#'     instance.
#'
#' @export
TENxGenomics <-
    function(h5path)
{
    stopifnot(
        is.character(h5path), length(h5path) == 1L, !is.na(h5path),
        nzchar(h5path), file.exists(h5path)
    )

    h5f <- H5Fopen(h5path)
    on.exit(H5Fclose(h5f))

    tmpl <- .TENxGenomics()
    rowidx <- .h5_dimidx(h5f, .rowname(tmpl))
    colidx <- .h5_dimidx(h5f, .colname(tmpl))

    .TENxGenomics(tmpl, h5path = h5path, rowidx = rowidx, colidx = colidx)
}

##
## dim, dimnames
##

#' @rdname TENxGenomics-class
#'
#' @aliases dim,TENxGenomics-method
#'
#' @exportMethod dim
setMethod("dim", "TENxGenomics",
    function(x)
{
    c(length(.rowidx(x)), length(.colidx(x)))
})

#' @rdname TENxGenomics-class
#'
#' @aliases dimnames,TENxGenomics-method
#'
#' @exportMethod dimnames
setMethod("dimnames", "TENxGenomics",
    function(x)
{
    dim <- dim(x)

    h5f <- H5Fopen(.h5path(x))
    on.exit(H5Fclose(h5f))

    list(
        as.character(h5read(h5f, .rowname(x), index=list(.rowidx(x)))),
        as.character(h5read(h5f, .colname(x), index=list(.colidx(x))))
    )
})

##
## Subset
##

.subset_from_character <-
    function(x, name, idx, i)
{
    h5f <- H5Fopen(.h5path(x))
    on.exit(H5Fclose(h5f))

    names <- h5read(h5f, name, index = list(idx))
    match(i, names)
}

.subset_from_logical <-
    function(x, idx, i)
{
    if (length(i) > length(idx))
        stop("logical subscript too long")
    which(rep(i, length=length(idx)))
}

.subset_as_integer <-
    function(x, name, idx, i)
{
    if (missing(i))
        return(idx)
    switch(
        class(i),
        integer = i,
        numeric = as.integer(i),
        character = .subset_from_character(x, name, idx, i),
        logical = .subset_from_logical(x, idx, i),
        default = stop("unsupport subset class ", class(i))
    )
}

#' @rdname TENxGenomics-class
#'
#' @aliases [,TENxGenomics,ANY,ANY-method
#'
#' @param x A \code{TENxGenomics-class} instance.
#'
#' @param i integer(), numeric(), character(), or logical() index into
#'     rows of \code{x}.
#'
#' @param j integer(), numeric(), character(), or logical() index into
#'     columns of \code{x}.
#'
#' @param drop logical(1) TRUE only.
#'
#' @param ... Additional arguments, ignored.
#'
#' @exportMethod [
setMethod("[", c("TENxGenomics", "ANY", "ANY"),
    function(x, i, j, ..., drop = TRUE)
{
    stopifnot(identical(drop, TRUE))

    i <- .subset_as_integer(x, .rowname(x), .rowidx(x), i)
    j <- .subset_as_integer(x, .colname(x), .colidx(x), j)
    stopifnot(all(i %in% .rowidx(x)), all(j %in% .colidx(x)))

    initialize(x, rowidx = .rowidx(x)[i], colidx = .colidx(x)[j])
})

##
## coerce
##

# Helper for common logic underlying .as.matrix() and .as.sparseMatrix()
.getValuesAndIndex <-
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
    values <- h5read(h5f, "/mm10/data", index=list(idx))
    ridx <- match(ridx, .rowidx(x))[keep]
    if (sparse) {
      # It is more (memory) efficient to supply i and p, rather than i and j,
      # to Matrix::sparseMatrix():
      # i = ridx: length(i) == length(values)
      # j = cidx: length(j) == length(values)
      # p = c(startidx, endidx[length(endidx)] + 1L): length(p) == (ncol(x) + 1)
      i <- ridx
      p <- c(startidx, endidx[length(endidx)] + 1L)
      return(list(values = values, i = i, p = p))
    }
    cidx <- rep(seq_along(cidx), lens)[keep]

    list(values = values, ridx = ridx, cidx = cidx)
}


.as.matrix <-
    function(x, ..., withDimnames=TRUE)
{
      val_and_idx <- .getValuesAndIndex(x=x, sparse=FALSE, ...,
                                        withDimnames=withDimnames)
      values <- val_and_idx[["values"]]
      ridx <- val_and_idx[["ridx"]]
      cidx <- val_and_idx[["cidx"]]

      ## formulate result as matrix
      m <- matrix(
        0L, nrow(x), ncol(x),
        dimnames = if (withDimnames) dimnames(x) else list(NULL, NULL)
      )
      m[cbind(ridx, cidx)] <- values
}

.as.CsparseMatrix <-
  function(x, ..., withDimnames=TRUE)
{
    # TODO: Support withDimnames
    stopifnot(withDimnames)
    val_and_idx <- .getValuesAndIndex(x, sparse = TRUE, ...,
                                      withDimnames = withDimnames)
    values <- val_and_idx[["values"]]
    i <- val_and_idx[["i"]]
    p <- val_and_idx[["p"]]
    sparseMatrix(i = i,
                 p = p,
                 x = values,
                 dims = dim(x),
                 dimnames = dimnames(x),
                 giveCsparse = TRUE)
}

#' @rdname TENxGenomics-class
#'
#' @aliases as.matrix.TENxGenomics
#'
#' @param withDimnames logical(1) Include dimnames on returned matrix?
#'
#' @return \code{as.matrix(tenx)} and \code{as(tenx, "matrix")} return
#'     a matrix with dim and dimnames equal to \code{tenx}, and values
#'     the read counts overlapping corresponding genes and
#'     samples. Use \code{as.matrix(withDimnames=FALSE)} to suppress
#'     dimnames on the returned matrix. NOTE: consider the size of the
#'     matrix, \code{prod(as.numeric(dim(tenx)))} before invoking this
#'     function.  #' @export
#'
#' @export
as.matrix.TENxGenomics <- .as.matrix

#' @rdname TENxGenomics-class
#'
#' @name coerce,TENxGenomics,matrix-method
#'
#' @exportMethod coerce
setAs("TENxGenomics", "matrix", function(from) .as.matrix(from))

#' @rdname TENxGenomics-class
#'
#' @name coerce,TENxGenomics,CsparseMatrix-method
#'
#' @exportMethod coerce
setAs("TENxGenomics", "CsparseMatrix", function(from) .as.CsparseMatrix(from))

##
## show
##

#' @rdname TENxGenomics-class
#'
#' @aliases show,TENxGenomics-method
#'
#' @param object A \code{TENxGenomics-class} instance.
#'
#' @exportMethod show
setMethod("show", "TENxGenomics",
    function(object)
{
    cat(
        "class:", class(object),
        "\nh5path:", .h5path(object),
        "\ndim():", nrow(object), "x", ncol(object),
        "\n"
    )
})
