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

.h5_colidx <-
    function(h5f)
{
    .h5_dimidx(h5f, "/mm10/barcodes")
}

.h5_rowidx <-
    function(h5f)
{
    .h5_dimidx(h5f, "/mm10/genes")
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
    colidx <- .h5_colidx(h5f)
    rowidx <- .h5_rowidx(h5f)
    H5Fclose(h5f)

    .TENxGenomics(h5path = h5path, colidx = colidx, rowidx = rowidx)
}

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
        "\ndim:", length(.rowidx(object)), "x", length(.colidx(object)),
        "\n"
    )
})
