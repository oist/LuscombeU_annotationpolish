#' Title
#'
#' This function processes genomic annotations by splitting genes based on their overlap with provided CAGE peaks and filtering based on gene size. It is designed to integrate new gene models from a 'liftoff' annotation with an original annotation, prioritizing genes that intersect with CAGE peaks and are above a minimum size threshold.
#'
#' @param original_annotation A GenomicRanges object containing the original gene annotations.
#' @param liftoff A GenomicRanges object with 'liftoff' gene models to be integrated into the original annotation.
#' @param cage_peaks A GenomicRanges object containing CAGE peaks used for identifying genes to be split.
#' @param min_gene_size Integer, minimum gene size (in base pairs) for a gene to be retained. Default is 1000 bp.
#' @param remove_no_ORF Logical, indicating whether to remove liftoff gene models without a valid ORF. Default is FALSE.
#'
#' @return A GenomicRanges object containing the modified gene annotations, with some genes split or replaced based on the overlap with CAGE peaks and the size threshold.
#' @export
#' @import tidyverse
#'
#' @examples
#' cage_peaks <- read_cage("https://raw.githubusercontent.com/oist/LuscombeU-CAGE_libraries/main/2021-12-17_Okinawa_Oik/CrossAlignments/consensus_clusters_no_OKItoOKI.bed", "https://raw.githubusercontent.com/oist/LuscombeU-CAGE_libraries/main/2021-12-17_Okinawa_Oik/CrossAlignments/consensus_clusters_sl_OKItoOKI.bed" )
#' original_annotation <- system.file("extdata", "split_example_target.gff3", package="annotationpolish") |> rtracklayer::import()
#' liftoff <- system.file("extdata", "split_example_liftoff.gff3", package="annotationpolish") |> rtracklayer::import() |> fix_liftoff_oikobase()
#'
#' split_transcript(original_annotation, liftoff, cage_peaks)
split_transcript <- function(original_annotation, liftoff, cage_peaks, min_gene_size = 1000, remove_no_ORF = F){

  # Remove genes smaller than the min_gene_size threshold from the original annotation
  old_annotation <- original_annotation |> filter(gene_id %in% (original_annotation |> filter(type == "gene") |> filter(width >= 1000) %>% .$gene_id))

  # remove liftoff gene models if filtering is TRUE
  # if remove_no_ORF = T
  if (remove_no_ORF == T) {
    liftoff_orf_genes <- (liftoff |> filter(valid_ORF == "True"))$gene_id |> unique()
    liftoff <- liftoff |> filter(gene_id %in% liftoff_orf_genes)
  }

  # Remove liftoff gene models that do not intersect with the original annotation
  liftoff_annot_overlaps <- plyranges::join_overlap_intersect(liftoff |> filter(type == "gene"),
                                                              old_annotation |> filter(type == "gene"))


  liftoff <- liftoff |>
    filter(gene_id %in% (liftoff_annot_overlaps %>%
                           .$gene_id.x |> unique()))


  old_annotation <- old_annotation |>
    filter(gene_id %in% (liftoff_annot_overlaps %>%
                           .$gene_id.y |> unique()))

  # Identify genes to be split based on overlap
  gene_to_split_index <- findOverlaps(liftoff |> filter(type == "gene"), old_annotation |> filter(type == "gene")) |>
    subjectHits()


  annotation_to_split <- old_annotation |> filter(
    gene_id %in% (old_annotation |> filter(type == "gene"))[gene_to_split_index[duplicated(gene_to_split_index)]]$gene_id)

  # Identify gaps within the target genes for potential splitting
  gaps_within_target <- subsetByOverlaps(gaps(liftoff |> filter(type == "gene")),
                                         annotation_to_split |> filter(type == "gene"), type = "within")


  # Adjust CAGE peaks to their thick range and identify overlaps within gene gaps
  GenomicRanges::ranges(cage_peaks) <- cage_peaks$thick

  cage_peaks_within_gaps <- join_overlap_inner_directed(cage_peaks, gaps_within_target)

  # Find original genes overlapping with CAGE peaks for replacement
  genes_to_replace <- unique((annotation_to_split |> filter(type == "gene"))[subjectHits(findOverlaps(cage_peaks_within_gaps, annotation_to_split |> filter(type == "gene")))]$gene_id)

  message(paste(length(genes_to_replace), "Genes to replace", paste(genes_to_replace, collapse = ", ")))

  drop_mcols <- function(gr){
    mcols(gr) <- NULL
    return(gr)
  }

  genes_replacement <- join_overlap_inner_directed(drop_mcols(original_annotation |> filter(type == "gene", gene_id %in% genes_to_replace)),
                                                   liftoff |> filter(type == "gene"))$gene_id

  message(paste(length(genes_replacement), "New genes", paste(genes_replacement, collapse = ", ")))


  # Compile the final GenomicRanges object with modified annotations
  final_gr <- c(liftoff |> filter(gene_id %in% genes_replacement),
                original_annotation |> filter(!gene_id %in% genes_to_replace)) |> sort()


  final_gr <- final_gr |> select(-ID) |> add_unique_id()


  # unlist Parent and remove_transcript_id if gene
  mcols(final_gr)$Parent <- sapply(mcols(final_gr)$Parent, function(x) {
    if (is.null(x) || length(x) == 0) {
      return(NA)
    } else {
      return(x[[1]])
    }
  })


  return(final_gr)

}















