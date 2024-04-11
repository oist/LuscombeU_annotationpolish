#' Identifies Overlapping Genomic Ranges
#'
#' This function finds and returns genomic ranges within a `GRanges` object that overlap with other ranges,
#' excluding self-overlaps. It's designed to work with genomic range objects defined in the GenomicRanges package,
#' leveraging the `findOverlaps` method for the identification of overlaps.
#'
#' @param gr A `GRanges` object representing genomic ranges, where each range specifies a region on a genome.
#'
#' @return A `GRanges` object containing only those ranges that overlap with others, excluding any self-overlaps.
#'
#' @export
#'
#' @import GenomicRanges
#'
#' @examples
#' # Assuming `gr` is a GRanges object with genomic ranges:
#' overlapping_granges <- overlaps_within(gr)
#' # `overlapping_granges` will contain genomic ranges from `gr` that overlap with others.
overlaps_within <- function(gr){
  # Find overlaps, excluding self-overlaps
  overlaps <- findOverlaps(gr, drop.self = T)

  # Extract indices of the GRanges object that overlap with others
  overlapping_indices <- unique(c(queryHits(overlaps), subjectHits(overlaps)))

  # Subset the GRanges object to only those ranges that overlap with others
  overlapping_granges <- gr[overlapping_indices]

  return(overlapping_granges)
}



#' Checks for Various Types of Overlapping Genomic Ranges
#'
#' This function examines a `GRanges` object for overlaps between different genomic features,
#' such as exons, CDS (Coding Sequences), and UTRs (Untranslated Regions). It identifies overlaps
#' within exons, between CDS and 5' UTRs, and between CDS and 3' UTRs, returning a detailed list of these overlaps.
#' The function is intended to work on annotated genomic ranges that include information on feature types.
#'
#' @param annot_gr A `GRanges` object with annotations for genomic features, including but not limited to,
#' types like "CDS", "five_prime_UTR", and "three_prime_UTR". Each range is expected to include metadata specifying
#' the feature type.
#'
#' @return A list containing `GRanges` objects for each type of overlap identified: exons overlapping with other exons,
#' CDS overlapping with 5' UTRs, and CDS overlapping with 3' UTRs. If no overlaps are found for a specific category,
#' a message is printed, and the category is excluded from the output list.
#'
#' @export
#' @import GenomicRanges
#' @import dplyr
#'
#' @examples
#' # Assuming `annot_gr` is a GRanges object with annotations:
#' overlaps_list <- check_overlaps(annot_gr)
#' # `overlaps_list` will contain lists of overlapping genomic ranges for specified features.
check_overlaps <- function(annot_gr){


  CDS <- annot_gr |> filter(type == "CDS")

  fiveutr_cds <- c(annot_gr |> filter(type == "five_prime_UTR"), CDS)

  threeutr_cds <- c(annot_gr |> filter(type == "three_prime_UTR"), CDS)

  # find overlapping ranges
  exons_overlaps <- exons |>
    split(exons$transcript_id) |>
    lapply(overlaps_within) |> GRangesList() |> unlist()

  fiveutr_overlaps <- fiveutr_cds |>
    split(fiveutr_cds$transcript_id) |>
    lapply(overlaps_within) |> GRangesList() |> unlist()

  threeutr_overlaps <- threeutr_cds |>
    split(threeutr_cds$transcript_id) |>
    lapply(overlaps_within) |> GRangesList() |> unlist()

  # create a list for output if there is anything wrong
  output_list <- c()

  if (length(exons_overlaps) == 0 ) {
    message("No exon overlaps")
  }else{
    output_list$exons <- exons_overlaps
    message(paste("check exons of ", paste(unique(exons_overlaps$transcript_id), collapse = ", "), "!!!"))
  }

  if (length(fiveutr_overlaps) == 0 ) {
    message("No 5' UTR overlapping CDS")
  }else{
    output_list$fiveutr_overlaps <- fiveutr_overlaps
    message(paste("check 5' UTR of ", paste(unique(fiveutr_overlaps$transcript_id), collapse = ", "), "!!!"))
  }

  if (length(threeutr_overlaps) == 0 ) {
    message("No 3' UTR overlapping CDS")
  }else{
    output_list$threeutr_overlaps <- threeutr_overlaps
    message(paste("check 3' UTR of ", paste(unique(threeutr_overlaps$transcript_id), collapse = ", "), "!!!"))
  }


  return(output_list)

}


