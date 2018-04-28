### =========================================================================
### TENxMatrix objects
### -------------------------------------------------------------------------

#' @rdname TENxMatrix-class
#'
#' @title TENxMatrix objects
#'
#' @name TENxMatrix-class
#'
#' @aliases TENxMatrixSeed-class
#'
#' @aliases DelayedArray,TENxMatrixSeed-method
#'
#' @description A container for representing 10xGenomics data.
#'     TENxMatrix extends \link{DelayedArray} so all the operations
#'     available on \link{DelayedArray} objects work on TENxMatrix objects.
#'
#' @import methods BiocGenerics S4Vectors IRanges DelayedArray HDF5Array
#'
#' @importFrom tools file_path_as_absolute
#'
#' @exportClass TENxMatrixSeed
setClass("TENxMatrixSeed",
    representation(
        file="character",     # Absolute path to the HDF5 file so the object
                              # doesn't break when the user changes the working
                              # directory (e.g. with setwd()).
        group="character",    # Name of the group in the HDF5 file containing
                              # 10xGenomics data.
        dim="integer",
        dimnames="list",
        col_ranges="data.frame"  # Can't use an IRanges object for this at the
                                 # moment because they don't support Linteger
                                 # start/end values yet.
    )
)

### All the 10xGenomics components are monodimensional.
.get_TENx_component <- function(file, group, name, idx=NULL)
{
    name <- paste0(group, "/", name)
    if (!is.null(idx))
        idx <- list(idx)
    as.vector(h5read(file, name, index=idx))
}

.get_barcodes <- function(file, group, idx=NULL)
    .get_TENx_component(file, group, "barcodes", idx=idx)

.get_data <- function(file, group, idx=NULL)
    .get_TENx_component(file, group, "data", idx=idx)

.get_gene_names <- function(file, group, idx=NULL)
    .get_TENx_component(file, group, "gene_names", idx=idx)

.get_genes <- function(file, group, idx=NULL)
    .get_TENx_component(file, group, "genes", idx=idx)

.get_indices <- function(file, group, idx=NULL)
    .get_TENx_component(file, group, "indices", idx=idx)

.get_indptr <- function(file, group, idx=NULL)
{
    name <- paste0(group, "/indptr")
    if (!is.null(idx))
        idx <- list(idx)
    as.vector(h5read(file, name, index=idx, bit64conversion="double"))
}

.get_shape <- function(file, group, idx=NULL)
    .get_TENx_component(file, group, "shape", idx=idx)

TENxMatrixSeed <- function(h5path, group="mm10")
{
    if (!isSingleString(h5path))
        stop(wmsg("'h5path' must be a single string specifying the path to ",
                  "the HDF5 file where the 10xGenomics data is located"))
    h5path <- file_path_as_absolute(h5path)
    if (!isSingleString(group))
        stop(wmsg("'group' must be a single string specifying the name ",
                  "of the group in the HDF5 file containing the ",
                  "10xGenomics data"))
    if (group == "")
        stop(wmsg("'group' cannot be the empty string"))

    ## dim
    dim <- .get_shape(h5path, group)
    stopifnot(length(dim) == 2L)

    ## dimnames
    rownames <- .get_genes(h5path, group)
    stopifnot(length(rownames) == dim[[1L]])
    colnames <- .get_barcodes(h5path, group)
    stopifnot(length(colnames) == dim[[2L]])
    dimnames <- list(rownames, colnames)

    ## col_ranges
    ## "/data" and "/indices" are monodimensional arrays of length >= 2^31 so
    ## we need to call h5dim() with 'check=FALSE'.
    data_len <- HDF5Array:::h5dim(h5path, paste0(group, "/data"),
                                  check=FALSE)
    stopifnot(length(data_len) == 1L)
    indices_len <- HDF5Array:::h5dim(h5path, paste0(group, "/indices"),
                                  check=FALSE)
    stopifnot(identical(data_len, indices_len))
    indptr <- .get_indptr(h5path, group)
    stopifnot(length(indptr) == dim[[2L]] + 1L,
              indptr[[1L]] == 0L,
              indptr[[length(indptr)]] == data_len)
    col_ranges <- data.frame(start=indptr[-length(indptr)] + 1,
                             width=as.integer(diff(indptr)))

    new2("TENxMatrixSeed", file=h5path,
                           group=group,
                           dim=dim,
                           dimnames=dimnames,
                           col_ranges=col_ranges)
}

#' @rdname TENxMatrix-class
#'
#' @aliases dim,TENxMatrixSeed-method
#'
#' @param x A \code{TENxMatrix-class} or \code{TENxMatrixSeed-class}
#'     object.
#'
#' @exportMethod dim
setMethod("dim", "TENxMatrixSeed", function(x) x@dim)

#' @rdname TENxMatrix-class
#'
#' @aliases dimnames,TENxMatrixSeed-method
#'
#' @exportMethod dimnames
setMethod("dimnames", "TENxMatrixSeed", function(x) x@dimnames)

### S4Vectors:::fancy_mseq() does not accept 'offset' of type double yet so
### we implement a version that does.
.fancy_mseq <- function(lengths, offset=0)
{
    lengths_len <- length(lengths)
    if (lengths_len == 0L)
        return(numeric(0))
    offsets <- offset - cumsum(c(0L, lengths[-lengths_len]))
    seq_len(sum(lengths)) + rep.int(offsets, lengths)
}

.extract_array_from_TENxMatrixSeed <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    ans <- array(0L, dim=ans_dim)

    j <- index[[2L]]
    col_ranges <- x@col_ranges
    if (is.null(j)) {
        j <- seq_len(ncol(ans))
    } else {
        col_ranges <- S4Vectors:::extract_data_frame_rows(col_ranges, j)
    }
    start2 <- col_ranges[ , "start"]
    width2 <- col_ranges[ , "width"]
    idx2 <- .fancy_mseq(width2, offset=start2 - 1)
    i2 <- .get_indices(x@file, x@group, idx=idx2) + 1L
    j2 <- rep.int(seq_along(j), width2)

    i <- index[[1L]]
    if (!is.null(i)) {
        m <- findMatches(i2, i)
        idx2 <- idx2[from(m)]
        i2 <- to(m)
        j2 <- j2[from(m)]
    }

    ans[cbind(i2, j2)] <- .get_data(x@file, x@group, idx=idx2)
    ans
}

setMethod("extract_array", "TENxMatrixSeed",
    .extract_array_from_TENxMatrixSeed
)

#' @rdname TENxMatrix-class
#'
#' @exportClass TENxMatrix
setClass("TENxMatrix", contains="DelayedMatrix")

.validate_TENxMatrix <- function(x)
{
    if (!is(x@seed, "TENxMatrixSeed"))
        return(wmsg("'x@seed' must be a TENxMatrixSeed object"))
    if (!DelayedArray:::is_pristine(x))
        return(wmsg("'x' carries delayed operations"))
    TRUE
}

setValidity2("TENxMatrix", .validate_TENxMatrix)

setMethod("DelayedArray", "TENxMatrixSeed",
    function(seed) DelayedArray:::new_DelayedArray(seed, Class="TENxMatrix")
)

### Works directly on a TENxMatrixSeed object, in which case it must be called
### with a single argument.
#' @rdname TENxMatrix-class
#'
#' @param h5path character(1) path to the 10xGenomics hdf5 file, or a
#'     \code{TENxMatrixSeed} instance.
#'
#' @param group character(1) group (e.g., \dQuote{mm10}) containing
#'     scRNA-seq data; only allowed when \code{h5path} is a
#'     character(1).
#'
#' @return \code{TENxMatrix()} returns a TENxMatrix instance;
#'     \code{dim()}, \code{dimnames()} and all other
#'     \code{DelayedArray} methods operate on this instance.
#'
#' @export
TENxMatrix <- function(h5path, group="mm10")
{
    if (is(h5path, "TENxMatrixSeed")) {
        if (!missing(group))
            stop(wmsg("TENxMatrix() must be called with a single argument ",
                      "when passed a TENxMatrixSeed object"))
        seed <- h5path
    } else {
        seed <- TENxMatrixSeed(h5path, group)
    }
    DelayedArray(seed)
}
