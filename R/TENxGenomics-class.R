#' @import methods
#' @import rhdf5

.TENxGenomics <- setClass(
    "TENxGenomics",
    slots = c(
        h5path = "character",
        colidx = "integer",
        rowidx = "integer",
        ## internal
        dataname = "character",
        rowname = "character",
        genenames = "character",
        colname = "character",
        indices = "character",
        indptr = "character"
    ),
    prototype = prototype(
        dataname = "/mm10/data",
        rowname = "/mm10/genes",
        genenames = "/mm10/gene_names",
        colname = "/mm10/barcodes",
        indices = "/mm10/indices",
        indptr = "/mm10/indptr"
    )
)

.h5path <- function(x) x@h5path

.colidx <- function(x) x@colidx

.rowidx <- function(x) x@rowidx

.dataname <- function(x) x@dataname

.rowname <- function(x) x@rowname

.genenames <- function(x) x@genenames

.colname <- function(x) x@colname

.indptr <- function(x) x@indptr

.indices <- function(x) x@indices

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
#'     subset and input 10xGenomics'
#'     1M_neurons_filtered_gene_bc_matrices_h5.h5 file. Subsetting is
#'     a light-weight operation; input (typically of the subset
#'     matrix) is as a dense matrix and hence consumes memory.
#'
#' @param h5path character(1) file path to a 1M_neurons_*.h5 file.
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
#' @aliases [,TENxGenomics,ANY,ANY-method [,TENxGenomics-method
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
    i <- .subset_as_integer(x, .rowname(x), .rowidx(x), i)
    j <- .subset_as_integer(x, .colname(x), .colidx(x), j)
    stopifnot(all(i %in% .rowidx(x)), all(j %in% .colidx(x)))

    initialize(x, rowidx = .rowidx(x)[i], colidx = .colidx(x)[j])
})

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
