#' @import methods
#' @import rhdf5

.TENxGenomics <- setClass(
    "TENxGenomics",
    slots = c(
        h5path = "character",
        colidx = "integer",
        rowidx = "integer"
    )
)

.h5path <- function(x) x@h5path

.colidx <- function(x) x@colidx

.rowidx <- function(x) x@rowidx

.h5_dimidx <-
    function(h5f, name)
{
    h5d <- H5Dopen(h5f, name)
    on.exit(H5Dclose(h5d))
    n <- H5Sget_simple_extent_dims(H5Dget_space(h5d))$size
    seq_len(n)
}

.h5_colidx <- function(h5f) .h5_dimidx(h5f, "/mm10/barcodes")

.h5_rowidx <- function(h5f) .h5_dimidx(h5f, "/mm10/genes")

.h5_dimnames <-
    function(h5f, rowidx, colidx)
{
    h5d <- H5Dopen(h5f, "/mm10/barcodes")
    on.exit(H5Dclose(h5d))
browser()
    list(
        h5read(h5f, "/mm10/genes", colidx),
        h5read(h5f, "/mm10/barcodes", colidx)
    )
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
    colidx <- .h5_colidx(h5f)
    rowidx <- .h5_rowidx(h5f)

    .TENxGenomics(h5path = h5path, colidx = colidx, rowidx = rowidx)
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
        h5read(h5f, "/mm10/genes", index=list(.rowidx(x))),
        h5read(h5f, "/mm10/barcodes", index=list(.colidx(x)))
    )
})

##
## Subset
##

.subset_i_from_character <-
    function(x, i)
{
}

.subset_i_as_integer <-
    function(x, idx)
{
    if (missing(idx))
        return(.rowidx(x))
    switch(
        class(idx),
        integer = idx,
        numeric = as.integer(idx),
        character = .subset_i_from_character(x, idx),
        logical = .subset_i_from_logical(x, idx)
    )
}

.subset_j_as_integer <-
    function(x, idx)
{
    if (missing(idx))
        return(.colidx(x))

    switch(
        class(idx),
        integer = idx,
        numeric = as.integer(idx),
        character = .subset_i_from_character(x, idx),
        logical = .subset_i_from_logical(x, idx)
    )
}

#' @rdname TENxGenomics-class
#'
#' @aliases [,TENxGenomics,ANY,ANY-method
#'
#' @param x A \code{TENxGenomics-class} instance.
#'
#' @param i integer() or numeric() index into rows of \code{x}.
#'
#' @param j integer() or numeric() index into columns of \code{x}.
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

    i <- .subset_i_as_integer(x, i)
    j <- .subset_j_as_integer(x, j)
    stopifnot(all(i %in% .rowidx(x)), all(j %in% .colidx(x)))

    initialize(x, rowidx = .rowidx(x)[i], colidx = .colidx(x)[j])
})

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
        "\ndim:", nrow(object), "x", ncol(object),
        "\n"
    )
})
