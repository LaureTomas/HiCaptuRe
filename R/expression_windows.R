#' Extract counts for windows flanking a set of regions
#'
#' This function extract the counts for the genes in windows flanking a set of given regions
#'
#' @param regions a GenomicRanges object
#' @param windows a sequence of windows starting from 0
#' @param dds a DESeqDataSet or DESeqTransform object
#' @param promoters a GenomicRanges with the promoter regions and the gene associated
#'
#' @return dataframe with 7 columns: gene_id, sample, count, cond, windows, total_genes, all_zero
#' 1. 'gene_id' gene_id from the dds object
#' 2. 'sample' name of the sample in the dds object
#' 3. 'count' count value for the specific gene
#' 4. 'cond' condition (name of sample without _1/_2)
#' 5. 'windows' number of the window
#' 6. 'total_genes' number of total genes in the specific windows
#' 7. 'all_zero' number of genes that have 0 counts in all samples
#'
#' @importFrom IRanges mergeByOverlaps
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors queryHits
#' @importFrom GenomicRanges disjoin findOverlaps reduce
#' @importFrom regioneR extendRegions
#' @importFrom reshape2 melt
#'
#' @export
expression_windows <- function(regions,windows,dds,promoters)
{
  reg_prom <- IRanges::mergeByOverlaps(regions, promoters)
  genes <- unique(reg_prom$gene_id)
  net_counts <- as.data.frame(SummarizedExperiment::assay(dds[genes]))
  all_zero <- apply(net_counts, 1, function(x)all(x==0))
  net_counts$gene_id <- rownames(net_counts)

  # first <- reshape2::melt(net_counts[!all_zero])
  first <- suppressMessages(reshape2::melt(net_counts))
  first$cond <- gsub("_[12]","",first$variable)
  first$windows <- windows[1]
  first$total_genes <- length(genes)
  first$all_zero <- sum(all_zero)
  colnames(first)[2:3] <- c("sample","count")

  windows_counts <- data.frame()
  progress_bar = txtProgressBar(min=2, max=length(windows), style = 3, char="=")

  for (i in 2:length(windows))
  {
    prev <- regioneR::extendRegions(regions, extend.start = windows[i-1],extend.end = windows[i-1])

    a <- regioneR::extendRegions(regions, extend.start = windows[i],extend.end = windows[i])
    a2 <- GenomicRanges::disjoin(a,ignore.strand=T)

    b <- c(prev,a2)
    b <- GenomicRanges::disjoin(b,ignore.strand=T)

    c <- GenomicRanges::findOverlaps(b,prev)

    flank_regions <- b[!1:length(b) %in% unique(queryHits(c))]
    flank_regions <- GenomicRanges::reduce(flank_regions)

    flank_prom <- IRanges::mergeByOverlaps(flank_regions, promoters)
    flank_genes <- unique(flank_prom$gene_id)

    net_counts <- as.data.frame(assay(dds[flank_genes]))
    all_zero <- apply(net_counts, 1, function(x)all(x==0))
    net_counts$gene_id <- rownames(net_counts)

    # net_melt <- reshape2::melt(net_counts[!all_zero,])
    net_melt <- reshape2::melt(net_counts)

    net_melt$cond <- gsub("_[12]","",net_melt$variable)
    net_melt$windows <- windows[i]
    net_melt$total_genes <- length(flank_genes)
    net_melt$all_zero <- sum(all_zero)
    colnames(net_melt)[2:3] <- c("sample","count")

    windows_counts <- rbind(windows_counts,net_melt)
    setTxtProgressBar(progress_bar, value = i)
  }

  return(rbind(first,windows_counts))
  close(progress_bar)
}
