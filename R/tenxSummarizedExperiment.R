.rowData <-
    function(x)
{
    S4Vectors::DataFrame(
        Ensembl = rownames(x),
        Symbol = h5read(.h5path(x), .genenames(x))
    )
}

.colData_1M_neurons <-
    function(x)
{
    barcode <- colnames(x)
    sequence <- sub("-.*", "", barcode)
    lib <- as.integer(sub(".*-", "", barcode))
    mouse <- ifelse(lib <= 69, "A", "B")
    S4Vectors::DataFrame(
        Barcode = barcode,
        Sequence = sequence,
        Library = lib,
        Mouse = factor(mouse)
    )
}

.colData <-
    function(x)
{
    if (startsWith(basename(.h5path(x)), "1M_neurons_")) {
        .colData_1M_neurons(x)
    } else {
        S4Vectors::DataFrame(i=seq_len(ncol(x)))[, FALSE]
    }
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
#' @param rowData A \code{DataFrame()} with as many rows as there are
#'     genes in the 10xGenomics file. If missing, an object is created
#'     with 'gene' and 'genename' fields from the hdf5 file.
#'
#' @param colData A \code{DataFrame()} with as many rows as there are
#'     samples in the 10xGenomics file. If missing, and object is
#'     constructed from the barcodes of the hdf5 file. The sequence
#'     and library portions are separated, and the mouse is inferred
#'     (libraries <= 69 are from mouse "A", others are from mouse
#'     "B").
#'
#' @return A \code{SummarizedExperiment} object. For effective use of
#'     the return object, load the SummarizedExperiment library into
#'     the R session.
#' 
#' @export
tenxSummarizedExperiment <-
    function(h5path, rowData, colData)
{
    requireNamespace("S4Vectors", quietly=TRUE)
    requireNamespace("SummarizedExperiment", quietly=TRUE)
    if (!"SummarizedExperiment" %in% sub("package:", "", search()))
        message("use 'library(\"SummarizedExperiment\")' to access this object")

    tenx <- TENxGenomics(h5path)

    if (missing(rowData))
        rowData <- .rowData(tenx)
    if (missing(colData))
        colData <- .colData(tenx)

    SummarizedExperiment::SummarizedExperiment(
        assays = list(tenx), rowData = rowData, colData = colData
    )
}
