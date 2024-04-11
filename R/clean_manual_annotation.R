#' Clean Manual Annotations
#'
#' This function imports manual annotation data from a specified GFF3 file. It processes
#' the annotations to select specific fields, including 'source', 'type', 'score', 'phase',
#' 'Name', 'gene_name', 'gene_id', 'transcript_id', 'Parent', and optionally 'tss_type',
#' if present. The 'ID' field is explicitly removed to clean the annotation data for further
#' analysis. Terminal exon and gene/transcript boundaries are adjusted based on the terminal features.
#'
#' @param file A character string specifying the path to the GFF3 file containing
#' manual annotation data.
#'
#' @return An object of class `GRanges` representing the cleaned manual annotation data.
#' This object includes selected fields relevant for further analysis, with 'tss_type'
#' included only if it was present in the original data.
#'
#' @export
#'
#' @examples
#' # Load and clean manual annotation data
#' cleaned_data <- clean_manual_anno(system.file("extdata", "manual_data.gff3", package="annotationpolish"))
#'
clean_manual_anno <- function(file) {
  # Import the annotation file
  anno <- rtracklayer::import.gff3(file)

  # Determine if 'tss_type' is present in the data
  all_fields <- names(GenomicRanges::mcols(anno))
  if ("tss_type" %in% all_fields) {
    # Select fields including 'tss_type' if present
    anno <- anno |>
      plyranges::select(source, type, score, phase, Name = gene_name, gene_name, gene_id, transcript_id, -ID, Parent, tss_type)
  } else {
    # Select fields without 'tss_type'
    anno <- anno |>
      plyranges::select(source, type, score, phase, Name = gene_name, gene_name, gene_id, transcript_id, -ID, Parent)
  }

  # adjust transcript and gene boundaries
  anno <- adjust_transcript_gene_boundaries(anno)

  # adjust exons
  anno <- adjust_terminal_exons(anno)

  return(anno)
}





#' Adjust Transcript and Gene Boundaries in GRanges Object
#'
#' This function takes a GRanges object and adjusts the start and end positions
#' for transcripts and genes based on their associated features. For transcripts,
#' the start and end positions are adjusted to encompass all its features. Similarly,
#' gene boundaries are adjusted to encompass all associated transcripts and features.
#' The function is designed to work with genomic annotation data where accurate
#' boundary information is crucial for downstream analysis.
#'
#' @param gr A GRanges object containing genomic annotations, including transcripts
#' and genes. The object should have metadata columns for 'transcript_id' and 'gene_id'.
#'
#' @return A modified GRanges object with adjusted start and end positions for
#' transcripts and genes.
#'
#' @export
#' @import GenomicRanges
#' @import dplyr
#'
#' @examples
#' # Assuming `gr` is a GRanges object with transcript and gene annotations
#' adjusted_gr <- adjust_transcript_gene_boundaries(gr)
adjust_transcript_gene_boundaries <- function(gr){
  gr <- gr |> unlist_parent() |> as_tibble() |>
    # adjust transcript boundaries
    group_by(transcript_id) |> mutate(start = case_when(type == "transcript" ~ min(start), T ~ start),
                                      end = case_when(type == "transcript" ~ max(end), T ~ end)) |>
    # adjust gene boundaries
    group_by(gene_id) |> mutate(start = case_when(type == "gene" ~ min(start), T ~ start),
                             end = case_when(type == "gene" ~ max(end), T ~ end)) |>
    makeGRangesFromDataFrame(keep.extra.columns = T)

  return(gr)

}



#' Adjust Terminal Exons in GRanges Object
#'
#' This function processes a GRanges object to identify and adjust terminal exons
#' for each transcript. It first identifies the terminal (first and last) exons based
#' on their positions within the transcript. Then, it adjusts the start and end positions
#' of these terminal exons to ensure accurate representation of exon boundaries. This
#' adjustment is critical for analyses that rely on precise exon-intron structure,
#' such as differential exon usage or transcript assembly.
#'
#' @param gr A GRanges object containing genomic annotations with exon features.
#' The object should have a metadata column for 'transcript_id' and must distinguish
#' exon types in a 'type' metadata column.
#'
#' @return A modified GRanges object with adjusted positions for terminal exons
#' within each transcript.
#' @import GenomicRanges
#' @import dplyr
#'
#' @examples
adjust_terminal_exons <- function(gr){
  everything_else <- gr |> unlist_parent() |> as_tibble() |>
    group_by(transcript_id) |> filter(type != "exon")

  # mark terminal exons
  exons <- gr |> unlist_parent() |> as_tibble() |>
    filter(type == "exon") |>
    group_by(transcript_id) |>
    mutate(is_terminal_start = case_when(start == min(start) ~ TRUE, T ~ F),
           is_terminal_end = case_when(end == max(end) ~ TRUE, T ~ F)) |>
    ungroup()

  models <- bind_rows(exons, everything_else)

  # adjust the exons
  models <- models |> group_by(transcript_id) |>
    mutate(start = case_when(is_terminal_start == T ~ min(start), T ~ start),
           end = case_when(is_terminal_end == T ~ max(end), T ~ end))

  models_gr <- models |> makeGRangesFromDataFrame(keep.extra.columns = T)

  return(models_gr)
}

