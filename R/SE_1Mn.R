#' @import SummarizedExperiment
#' @title creates dense SummarizedExperiment with the 1Mneuron data
#' @param samps a vector of sample indices to use (might be able to use barcodes)
#' @param h5fname path to hdf5 for the 10x genomics data as distributed
#' @author VJ Carey 
#' @export
SE_1Mn = function(samps=1:100,
      h5fname = "./1M_neurons_filtered_gene_bc_matrices_h5.h5"
) {
  stopifnot(file.exists(h5fname))
  report = h5ls(h5fname)
  nrec = as.numeric(report$dim)[3] # length of data component
  bc = h5read(h5fname, "mm10/barcodes")
  gn = h5read(h5fname, "mm10/gene_names")
  ensg = h5read(h5fname, "mm10/genes")
  indptr = h5read(h5fname, "mm10/indptr",
    bit64conversion='double')
  nsamp = length(indptr)
  starts = as.numeric(indptr)+1
  ends = c(as.numeric(indptr)[-1], nrec)
  ngenes = length(ensg)
  assay = matrix(0, nr=ngenes, nc=length(samps))
  rownames(assay) = as.character(ensg)
  colnames(assay) = as.character(bc[samps])
  for (cursamp in samps) {  # could be parallel
    dat = as.numeric(h5read(h5fname, "mm10/data", index=list(starts[cursamp]:ends[cursamp])))
    curinds = 1+as.numeric(h5read(h5fname, "mm10/indices", index=list(starts[cursamp]:ends[cursamp])))
    assay[ensg[curinds],cursamp] = dat
    }
  SummarizedExperiment(assay)
}
    
