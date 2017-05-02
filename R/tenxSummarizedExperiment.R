#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment

.tenxGenomics <-
    function(h5path, i, j)
{
    if (missing(i) && missing(j)) {
        TENxGenomics(h5path)
    } else if (missing(i)) {
        TENxGenomics(h5path)[, j]
    } else if (missing(j)) {
        TENxGenomics(h5path)[i, ]
    } else {
        TENxGenomics(h5path)[i, j]
    }
}

.rowData <-
    function(x, rowData)
{
    if (!missing(rowData))
        return(DataFrame(rowData))

    DataFrame(
        Ensembl = rownames(x),
        Symbol = h5read(.h5path(x), .genenames(x), list(.rowidx(x)))
    )
}

.colData_1M_neurons <-
    function(x)
{
    barcode <- colnames(x)
    sequence <- sub("-.*", "", barcode)
    lib <- as.integer(sub(".*-", "", barcode))
    mouse <- ifelse(lib <= 69, "A", "B")
    DataFrame(
        Barcode = barcode,
        Sequence = sequence,
        Library = lib,
        Mouse = factor(mouse)
    )
}

.colData <-
    function(x, colData)
{
    if (!missing(colData))
        DataFrame(colData)
    else {
        tryCatch({
            .colData_1M_neurons(x)
         }, error = function(e) {
             warning(
                 "failed to create colData:",
                 "\n  ", conditionMessage(e)
             )
             DataFrame(i=seq_len(ncol(x)))[, FALSE]
        })
    }
}

.as.SummarizedExperiment <-
    function(tenx, rowData, colData, how)
{
    rowData <- .rowData(tenx, rowData)
    colData <- .colData(tenx, colData)

    SummarizedExperiment(
        assays = list(how(tenx)), rowData = rowData, colData = colData
    )
}

#' @rdname tenxSummarizedExperiment
#'
#' @title Create a 'SummarizedExperiment' from an 10xGenomics hdf5 file.
#'
#' @description The SummarizedExperiment \code{assay()} contains the
#'     \code{TENxGenomics} object corresponding to the underlying hdf5
#'     file. It also contains row and column annotations provided by
#'     the user or inferred from the hdf5 file. Inferred data requires
#'     a simple match between the file name prefix and
#'     \dQuote{1M_neurons_}; if the file name prefix does not match,
#'     row and column annotations are not created.
#'
#' @param h5path character(1) file path to a 1M_neurons_*.h5 file.
#'
#' @param x A \code{TENxGenomics-class} instance.
#'
#' @param i Optional integer(), character(), or logical() index used to subset
#'     rows (genes) of the \code{TENxGenomics} object.
#'
#' @param j Optional integer(), character(), or logical() index used to subset
#'     columns (samples) of the \code{TENxGenomics} object.
#'
#' @param rowData Optional \code{DataFrame()} with as many rows as
#'     there are genes in the 10xGenomics file or object. If missing,
#'     an object is created with 'gene' and 'genename' fields from the
#'     hdf5 file.
#'
#' @param colData Optional \code{DataFrame()} with as many rows as
#'     there are samples in the 10xGenomics file or object. If
#'     missing, and object is constructed from the barcodes of the
#'     hdf5 file. The sequence and library portions are separated, and
#'     the mouse is inferred (libraries <= 69 are from mouse "A",
#'     others are from mouse "B").
#'
#' @return \code{tenxSummarizedExperiment()} and
#'     \code{as.tenxSummarizedExperiment()} return a
#'     \code{SummarizedExperiment} instance where the assay() data are
#'     represented as a \code{TENxGenomics} object. Down-stream
#'     analysis will typically extract this object from (a subset of)
#'     the SummarizedExperiment, and coerce it to a, e.g,. matrix,
#'     \code{as.matrix(assay(se[, 1:100]))}.
#'
#' @export
tenxSummarizedExperiment <-
    function(h5path, i, j, rowData, colData)
{
    tenx <- .tenxGenomics(h5path, i, j)
    as.tenxSummarizedExperiment(tenx, rowData, colData)
}

#' @rdname tenxSummarizedExperiment
#'
#' @export
as.tenxSummarizedExperiment <-
    function(x, rowData, colData)
{
    stopifnot(is(x, "TENxGenomics"))
    .as.SummarizedExperiment(x, rowData, colData, identity)
}

#' @rdname tenxSummarizedExperiment
#'
#' @return \code{matrixSummarizedExperiment()} and
#'     \code{as.matrixSummarizedExperiment()} return a
#'     \code{SummarizedExperiment} instance where the assay() data are
#'     represented as a \code{Matrix::dgCMatrix} object. There are
#'     practical limits to the size of this object (e.g., 20k
#'     samples); the code is most efficient when consecutive samples
#'     are selected.
#'
#' @export
matrixSummarizedExperiment <-
    function(h5path, i, j, rowData, colData)
{
    tenx <- .tenxGenomics(h5path, i, j)
    as.matrixSummarizedExperiment(tenx, rowData, colData)
}

#' @rdname tenxSummarizedExperiment
#'
#' @export
as.matrixSummarizedExperiment <-
    function(x, rowData, colData)
{
    stopifnot(is(x, "TENxGenomics"))
    .as.SummarizedExperiment(x, rowData, colData, as.matrix)
}

#' @rdname tenxSummarizedExperiment
#'
#' @return \code{dgCMatrixSummarizedExperiment()} and
#'     \code{as.dgCMatrixSummarizedExperiment()} return a
#'     \code{SummarizedExperiment} instance where the assay() data are
#'     represented as a \code{Matrix::dgCMatrix} object. There are
#'     practical limits to the size of this object; the code is most
#'     efficient when consecutive samples are selected.
#' @export
dgCMatrixSummarizedExperiment <-
    function(h5path, i, j, rowData, colData)
{
    tenx <- .tenxGenomics(h5path, i, j)
    as.dgCMatrixSummarizedExperiment(tenx, rowData, colData)
}

#' @rdname tenxSummarizedExperiment
#'
#' @export
as.dgCMatrixSummarizedExperiment <-
    function(x, rowData, colData)
{
    .as.SummarizedExperiment(x, rowData, colData, as.dgCMatrix)
}
