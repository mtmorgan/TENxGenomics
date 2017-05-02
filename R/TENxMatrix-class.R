### =========================================================================
### TENxMatrix objects
### -------------------------------------------------------------------------

#' @rdname TENxMatrix-class
#'
#' @title TENxMatrix objects
#'
#' @description A container for representing 10xGenomics data.
#'     TENxMatrix extends \link{DelayedArray} so all the operations
#'     available on \link{DelayedArray} objects work on TENxMatrix objects.
#'
#' @exportClass TENxMatrixSeed TENxMatrix
#' @export TENxMatrix
#' @exportMethod dim dimnames coerce
#' @import methods BiocGenerics S4Vectors IRanges DelayedArray HDF5Array
#' @importFrom tools file_path_as_absolute
#' @aliases dim,TENxMatrixSeed-method
#' @aliases dimnames,TENxMatrixSeed-method


setClass("TENxMatrixSeed",
    representation(
        file="character",     # Absolute path to the HDF5 file so the object
                              # doesn't break when the user changes the working
                              # directory (e.g. with setwd()).
        group="character",    # Name of the group in the HDF5 file containing
                              # 10xGenomics data.
        dim="integer",
        dimnames="list",
        col_ranges="IRanges"
    )
)

### All the 10xGenomics components are monodimensional.
.get_TENx_component <- function(file, group, name, idx=NULL)
{
    name <- paste0(group, "/", name)
    if (!is.null(idx))
        idx <- list(idx)
    as.vector(HDF5Array:::h5read2(file, name, index=idx))
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
    .get_TENx_component(file, group, "indptr", idx=idx)

.get_shape <- function(file, group, idx=NULL)
    .get_TENx_component(file, group, "shape", idx=idx)

TENxMatrixSeed <- function(file, group="mm10")
{
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the HDF5 file where the 10xGenomics data is located"))
    file <- file_path_as_absolute(file)
    if (!isSingleString(group))
        stop(wmsg("'group' must be a single string specifying the name ",
                  "of the group in the HDF5 file containing the ",
                  "10xGenomics data"))
    if (group == "")
        stop(wmsg("'group' cannot be the empty string"))

    ## dim
    dim <- .get_shape(file, group)
    stopifnot(length(dim) == 2L)

    ## dimnames
    rownames <- .get_genes(file, group)
    stopifnot(length(rownames) == dim[[1L]])
    colnames <- .get_barcodes(file, group)
    stopifnot(length(colnames) == dim[[2L]])
    dimnames <- list(rownames, colnames)

    ## col_ranges
    data_len <- HDF5Array:::h5dim(file, paste0(group, "/data"))
    stopifnot(length(data_len) == 1L)
    indices_len <- HDF5Array:::h5dim(file, paste0(group, "/indices"))
    stopifnot(identical(data_len, indices_len))
    indptr <- .get_indptr(file, group)
    stopifnot(length(indptr) == dim[[2L]] + 1L,
              indptr[[1L]] == 0L,
              indptr[[length(indptr)]] == data_len)
    col_ranges <- as(PartitioningByEnd(indptr[-1L]), "IRanges")

    new2("TENxMatrixSeed", file=file,
                           group=group,
                           dim=dim,
                           dimnames=dimnames,
                           col_ranges=col_ranges)
}

setMethod("dim", "TENxMatrixSeed", function(x) x@dim)

setMethod("dimnames", "TENxMatrixSeed", function(x) x@dimnames)

.subset_TENxMatrixSeed_as_array <- function(seed, index)
{
    ans_dim <- DelayedArray:::get_subscripts_lengths(index, dim(seed))
    ans <- array(0L, dim=ans_dim)

    i <- index[[2L]]
    if (is.name(i))
        i <- seq_len(ncol(ans))
    col_ranges2 <- seed@col_ranges[i]
    idx2 <- as.integer(col_ranges2)
    i2 <- .get_indices(seed@file, seed@group, idx=idx2) + 1L
    j2 <- rep.int(seq_along(i), width(col_ranges2))

    j <- index[[1L]]
    if (!is.name(j)) {
        m <- findMatches(i2, j)
        idx2 <- idx2[from(m)]
        i2 <- to(m)
        j2 <- j2[from(m)]
    }

    ans[cbind(i2, j2)] <- .get_data(seed@file, seed@group, idx=idx2)
    ans
}

setMethod("subset_seed_as_array", "TENxMatrixSeed",
    .subset_TENxMatrixSeed_as_array
)

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
TENxMatrix <- function(file, group="mm10")
{
    if (is(file, "TENxMatrixSeed")) {
        if (!missing(group))
            stop(wmsg("TENxMatrix() must be called with a single argument ",
                      "when passed a TENxMatrixSeed object"))
        seed <- file
    } else {
        seed <- TENxMatrixSeed(file, group)
    }
    DelayedArray(seed)
}

