#' Clean Manual Annotations
#'
#' This function imports manual annotation data from a specified GFF3 file. It processes
#' the annotations to select specific fields, including 'source', 'type', 'score', 'phase',
#' 'Name', 'gene_name', 'gene_id', 'transcript_id', 'Parent', and optionally 'tss_type',
#' if present. The 'ID' field is explicitly removed to clean the annotation data for further
#' analysis.
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

  return(anno)
}

