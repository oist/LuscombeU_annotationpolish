#' Prioritize and Combine Transcript Annotations
#'
#' This function takes two sets of transcript annotations and prioritizes one over the other.
#' When transcripts from the two sets overlap, those from the priority set are retained,
#' and the overlapping transcripts from the secondary set are removed. The result is a combined
#' set of annotations that includes all non-overlapping transcripts from both sets, with precedence
#' given to transcripts from the priority set in case of overlaps.
#'
#' @param anno_priority A `GRanges` object containing the priority set of transcript annotations.
#' These annotations are given precedence, and any overlapping transcripts in the secondary set
#' (`anno_second`) are removed.
#' @param anno_second A `GRanges` object containing the secondary set of transcript annotations.
#' Overlaps with the priority set result in these transcripts being excluded from the final combined set.
#'
#' @return A `GRanges` object representing the combined set of annotations from both input sets,
#' after resolving overlaps by retaining priority set transcripts and excluding overlapping transcripts
#' from the secondary set. The resulting annotations are sorted by genomic coordinates.
#'
#' @export
#'
#' @import GenomicRanges
#' @import plyranges
#'
#' @examples
#' library(GenomicRanges)
#' library(plyranges)
#' # Assume `anno_priority` and `anno_second` are `GRanges` objects loaded or created previously
#' combined_annotations <- pick_transcripts(anno_priority, anno_second)
#' print(combined_annotations)
pick_transcripts <- function(anno_priority, anno_second){
  # Ensure annotations are only of type "transcript"
  anno_priority_transcript <- anno_priority |> filter(type == "transcript")
  anno_second_transcript <- anno_second |> filter(type == "transcript")

  # Find overlaps between priority and secondary transcript annotations
  overlaps <- findOverlaps(anno_priority_transcript, anno_second_transcript)

  # Get indices of overlapping transcripts in the secondary set for removal
  indices_to_remove <- subjectHits(overlaps) |> unique()

  # Remove overlapping transcripts from the secondary set
  anno_second_transcript <- anno_second_transcript[-indices_to_remove]

  # Filter the secondary annotations to include only non-overlapping transcripts
  anno_second_filtered <- anno_second |> filter(transcript_id %in% anno_second_transcript$transcript_id)

  # Combine priority annotations with filtered secondary annotations and sort
  combined_annotations <- c(anno_priority, anno_second_filtered) |> sort()

  return(combined_annotations)
}
