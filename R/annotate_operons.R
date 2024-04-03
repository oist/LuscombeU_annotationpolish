#' Calculate Operons
#'
#' Calculates operons based on genomic coordinates of genes. Operons are identified by grouping genes on the same strand that are within a specified distance from each other.
#'
#' @param genes A `GRanges` object containing genomic coordinates of genes.
#' @param window The maximum distance between consecutive genes in base pairs to consider them as part of the same operon. Default is 500 bp.
#'
#' @return A `GRanges` object representing the genomic coordinates of identified operons.
#' @export
#' @import GenomicRanges
#' @import BiocGenerics
#'
#' @examples
calcOperons <- function(genes, window = 500) {
  # Calculating operons on each strand separately so that nested genes do not
  # interrupt collinearity
  strands <- split(genes, BiocGenerics::strand(genes))
  ops.plus  <- calcOperons.onestrand(strands[['+']], window = window)
  ops.minus <- calcOperons.onestrand(strands[['-']], window = window)
  sort(c(ops.plus, ops.minus), ignore.strand = TRUE)
}

#' Calculate Operons on a Single Strand
#'
#' Helper function for `calcOperons` to calculate operons on a single genomic strand. This is used internally by `calcOperons` to handle plus and minus strands separately.
#'
#' @param genes A `GRanges` object containing genomic coordinates of genes on a single strand.
#' @param window The maximum distance between consecutive genes in base pairs to consider them as part of the same operon. Default is 500 bp.
#'
#' @return A `GRanges` object representing the genomic coordinates of identified operons on the given strand.
#' @export
#' @import GenomicRanges
#' @import BiocGenerics
#'
#' @examples
calcOperons.onestrand <- function(genes, window = 500) {
  if (identical(genes, GenomicRanges::GRanges())) return(GRanges())
  if(S4Vectors::nrun(BiocGenerics::strand(genes)) != 1)
    stop("Do not use this function on GRanges objects with more than one strand")
  # Divide the window by two because we will expand both sides of the range.
  halfwin <- round(window / 2)
  # Sort the genes ignoring strand.
  g <- sort(genes, ignore.strand = TRUE)
  # Expand by `halfwin` nucleotides
  # We suppress warnings because coordinates at the edges of seqfeatures can
  # become transiently negative.
  g <- (g + halfwin) |> suppressWarnings()
  # Merge into operons
  operons   <- GenomicRanges::reduce(g, min.gapwidth = 0L)
  # Remove the flanking sequences
  operons <- operons - halfwin
  # Count number of genes in operons
  operons$n <- GenomicRanges::countOverlaps(operons, genes)
  # Annotate operons with gene names.
  if(!is.null(genes$gene_id)){
    ov <- GenomicRanges::findOverlaps(operons, genes)
    operons$gene_id <- IRanges::CharacterList(split(genes$gene_id[S4Vectors::subjectHits(ov)], S4Vectors::queryHits(ov))) |> unname()
  }
  # Remove non-merged genes and return
  operons[operons$n > 1]
}


#' Split Operons and Assign Genes Based on CAGE Peaks
#'
#' Splits operons and assigns genes based on CAGE peaks, specifically excluding operons overlapping with intergenic CAGE peaks.
#'
#' @param operons A `GRanges` object representing operons.
#' @param genes A `GRanges` object representing genes.
#' @param cage A `GRanges` object representing CAGE peaks with thick metadata.
#'
#' @return A modified `GRanges` object of genes, with operon assignments split using CAGE peaks.
#' @export
#' @import tidyverse
#' @import plyranges
#' @import S4Vectors
#'
#' @examples
split_operons_assign_genes <- function(operons, genes, cage){
  # resize gene so that CAGE peaks that are assigned to 5' UTR is included
  genes_resized <- resize(genes, width(genes) - 1, fix="end")

  GenomicRanges::ranges(cage) <- cage$thick
  cage <- cage |> plyranges::filter(tss_type == "nosl")

  # remove CAGE peaks overlapping with genes (retain only intergenic CAGE peaks, including 5' UTRs)
  cage_intergenic <- cage[-(S4Vectors::subjectHits(GenomicRanges::findOverlaps(cage, genes_resized)))]

  # new operon boundaries split by CAGE non sl peaks
  operon_boundaries <- GenomicRanges::setdiff(operons, cage_intergenic)

  # overlap between genes and new operons
  new_overlaps <- findOverlaps(genes, operon_boundaries)

  operon_index <- subjectHits(new_overlaps)

  # remove non-operons if it is unique
  operonic_genes_index <- unique(operon_index[operon_index %in% names(table(operon_index)[table(operon_index) > 1])])

  # filter by this
  new_overlaps <- new_overlaps[subjectHits(new_overlaps) %in% operonic_genes_index]

  genes <- genes[queryHits(new_overlaps)] |> mutate(operon_id = subjectHits(new_overlaps))

  # assign new ID
  genes$operon_id <- paste("operon", genes$operon_id |> as.factor() |> as.numeric(), sep = "_")

  return(genes)

}


#' Summarize Operonic Genes with Leader and Trailing Assignment
#'
#' Summarizes genes within operons, identifying leader genes (first gene on + strand or last gene on - strand) and trailing genes.
#'
#' @param operonic_genes A tibble or data frame containing operonic genes with columns for operon_id, gene_id, and possibly others.
#'
#' @return A summarized `GRanges` with operon boundaries, leader, and trailing gene IDs.
#' @export
#' @import tidyverse
#'
#' @examples
summarize_operons <- function(operonic_genes){
  # label the genes if they are leader or not leader in an operon
  operonic_genes_labeled <- operonic_genes |> as_tibble()|>
    mutate(strand = as.character(strand)) |>
    group_by(operon_id, strand) |>
    mutate(
      status = case_when(
        strand == '+' & start == min(start) ~ 'leader',
        strand == '-' & end == max(end) ~ 'leader',
        TRUE ~ 'trailing'
      )
    ) |>
    ungroup()

  # summarise to get operon boundaries
  operon_summarized <- operonic_genes_labeled |>
    group_by(operon_id) |>
    summarise(
      seqnames = first(seqnames), # Assuming all entries in an operon have the same seqnames
      start = min(start),
      end = max(end),
      strand = first(strand), # Assuming all entries in an operon share the same strand
      source = first(source), # Assuming all entries in an operon come from the same source
      leader = list(gene_id[status == "leader"]), # Aggregating leader gene_id
      trailing = list(gene_id[status == "trailing"]), # Aggregating trailing gene_id
      .groups = 'drop'
    ) |> mutate(type = "operon") |> makeGRangesFromDataFrame(keep.extra.columns = T)

  return(list(summary = operon_summarized, gene_details = operonic_genes_labeled))
}


#' Create and Summarize Operons Based on Genes and CAGE Peaks
#'
#' This function integrates the process of calculating operons from genomic coordinates of genes,
#' splitting these operons based on CAGE peak information, and summarizing the results. It first
#' calculates operons based on proximity of genes on the same strand. Then, it refines these operons
#' by excluding regions overlapping with specified CAGE peaks, effectively splitting operons based
#' on potential transcription start sites. Finally, it summarizes the operons, identifying leader
#' and trailing genes within each operon based on their genomic positions.
#'
#' @param genes A `GRanges` object containing genomic coordinates of genes to be used for initial
#' operon calculation.
#' @param window An integer specifying the maximum distance in base pairs between genes to consider
#' them as part of the same operon. Defaults to 500.
#' @param cage A `GRanges` object containing CAGE peaks used to refine operon boundaries by
#' identifying potential transcription start sites.
#'
#' @return A summarized data frame containing the boundaries of each operon, the `seqnames`,
#' `strand`, and `source` of the operon, and lists of gene IDs categorized as either leader or
#' trailing within their respective operons.
#' @export
#'
#' @examples
#' # Assume 'genes' and 'cage' are predefined GRanges objects
#' annotation_gr <- system.file("extdata", "genes_5utr.gff3", package="annotationpolish") |> rtracklayer::import()
#' cage_peaks <- read_cage("https://raw.githubusercontent.com/oist/LuscombeU-CAGE_libraries/main/2021-12-17_Okinawa_Oik/CrossAlignments/consensus_clusters_no_OKItoOKI.bed", "https://raw.githubusercontent.com/oist/LuscombeU-CAGE_libraries/main/2021-12-17_Okinawa_Oik/CrossAlignments/consensus_clusters_sl_OKItoOKI.bed" )
#' make_operons(annotation_gr |> filter(type == "gene"), window = 500, cage_peaks |> filter(tss_type == "nosl"))
make_operons <- function(genes, window = 500, cage){
  operons_unsplit <- calcOperons(genes, window)
  operonic_genes <- split_operons_assign_genes(operons_unsplit, genes, cage)
  operon_summary <- summarize_operons(operonic_genes)
  return(operon_summary)
}
