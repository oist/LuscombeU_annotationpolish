#' Check for Unique IDs in a GRanges Object
#'
#' This function checks if the 'ID' values in a GRanges object are unique. It is useful
#' for ensuring data integrity after assigning unique identifiers to genomic features.
#'
#' @param gr A GRanges object that has been processed to include an 'ID' column with
#' identifiers for genomic features.
#'
#' @return Returns `TRUE` if all IDs in the 'ID' column are unique; otherwise, it returns
#' `FALSE` and emits a warning indicating that some IDs are not unique.
#'
#' @export
#'
#' @examples
#' # Assuming 'gr' is a GRanges object with an 'ID' column
#' result <- check_unique_ids(gr)
#' if (result) {
#'   message("IDs validation passed.")
#' } else {
#'   warning("Duplicate IDs found.")
#' }
check_unique_ids <- function(gr) {
  ids <- gr$ID
  unique_ids <- length(unique(ids))
  total_ids <- length(ids)

  if (unique_ids == total_ids) {
    message("All IDs are unique.")
    return(TRUE)
  } else {
    warning("Some IDs are not unique.")
    return(FALSE)
  }
}


#' Add Unique IDs to GRanges Entries
#'
#' This function adds a unique ID to each entry in a GRanges object based on its feature
#' type. For 'gene' and 'transcript' types, it uses 'gene_id' and 'transcript_id', respectively.
#' For other feature types (e.g., 'exon', 'CDS'), it creates a unique ID by concatenating the
#' feature type, 'transcript_id', and a unique number for each feature associated with the same
#' transcript. It also performs a check to ensure all IDs are unique.
#'
#' @param gr A GRanges object containing genomic annotations, which must include metadata columns
#' for 'type', 'gene_id', and 'transcript_id'. The 'type' column should categorize entries as
#' 'gene', 'transcript', or other feature types such as 'exon', 'CDS', etc.
#'
#' @return A GRanges object identical to the input but with an added 'ID' metadata column containing
#' unique identifiers for each entry. It also prints a message indicating whether all IDs are unique
#' after processing.
#'
#' @export
#'
#' @examples
#' # Assuming you have a GRanges object `gr` with the necessary columns
#' gr <- add_unique_id(gr)
#' # `gr` now includes an 'ID' column with unique identifiers for each entry.
add_unique_id <- function(gr) {
  gr$ID <- NA_character_
  gr$ID[gr$type == "gene"] <- as.character(gr$gene_id[gr$type == "gene"])
  gr$ID[gr$type == "transcript"] <- as.character(gr$transcript_id[gr$type == "transcript"])

  other_types <- (gr$type |> levels())[!(gr$type |> levels()) %in% c("gene", "transcript")]
  for (type in other_types) {
    indices <- which(gr$type == type)
    if (length(indices) > 0) {
      for (idx in indices) {
        transcript_id <- gr$transcript_id[idx]
        same_transcript_features <- which(gr$transcript_id == transcript_id & gr$type == type)
        number <- match(idx, same_transcript_features)
        gr$ID[idx] <- paste(type, transcript_id, number, sep = "_")
      }
    }
  }

  # Check if IDs are unique after assignment
  check_unique_ids(gr)

  return(gr)
}
