#' Extract and Refine Transcripts with Valid ORFs
#'
#' This function processes genomic range data to select and refine transcripts that contain valid open reading frames (ORFs). It filters the genomic ranges to include only those transcripts marked with a valid ORF and ensures that the gene boundaries are adjusted accordingly.
#'
#' @param gr A GenomicRanges object, usually from liftoff containing the genomic data with fields for `valid_ORF`, `transcript_id`, `type`, and `gene_id`. This should include annotations that distinguish different genomic features such as genes and transcripts.
#'
#' @return A GenomicRanges object that has been filtered and adjusted to include only transcripts with valid ORFs, along with their corresponding genes. The boundaries of genes are adjusted to encompass all associated transcripts.
#'
#' @export
#'
#' @examples
#' # Assume `gr` is a GenomicRanges object pre-loaded with necessary metadata
#' # Example usage:
#' refined_gr <- get_valid_orf_transcript(gr)
#'
get_valid_orf_transcript <- function(gr){

  valid_orf_tx <- (gr |> plyranges::filter(valid_ORF == "True"))$transcript_id |> unique()

  gr_tx <- gr |> plyranges::filter(type != "gene", transcript_id %in% valid_orf_tx)

  gr_gene <- gr |> plyranges::filter(type == "gene", gene_id %in% gr_tx$gene_id)


  #combine transcript and other features with gene and adjust gene boundaries
  gr <- c(gr_tx, gr_gene) |> annotationpolish::adjust_transcript_gene_boundaries()

  gr
}
