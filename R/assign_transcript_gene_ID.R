#' Assign IDs Based on Parent and ID Fields from GFF Files
#'
#' This function processes a GFF file or `GRanges` object and assigns transcript and gene IDs based on the "Parent" and "ID" fields. It ensures that IDs are properly assigned to the appropriate annotations.
#'
#' @param gff A `GRanges` object or a file path (character) to a GFF3 file.
#'
#' @return A GRanges object with the assigned `transcript_id` and `gene_id`.
#' @export
#'
#' @import tidyverse
#' @import plyranges
#'
#' @examples
#' # Example 1: Using a GRanges object
#' gr <- rtracklayer::import.gff3("example.gff3")
#' result <- assign_IDs_from_ParentandID(gr)
#'
#' # Example 2: Using a GFF3 file directly
#' result <- assign_IDs_from_ParentandID("example.gff3")
assign_IDs_from_ParentandID <- function(gff){

  # Check the input type and import GFF data accordingly
  if (class(gff) == "GRanges") {
    gr <- gff
  } else if (class(gff) == "character") {
    gr <- rtracklayer::import.gff3(gff)
  } else {
    stop("Input must be either a GRanges object or a GFF3 file path")
  }

  # Process the GRanges object, extracting relevant fields
  gr_df <- gr |>
    plyranges::select(source, type, score, phase, ID, Name, Parent) |>
    unlist_parent() |>
    as_tibble()

  # Assign transcript IDs
  gr_df_with_txid <- gr_df |>
    mutate(transcript_id = case_when(
      (!type %in% c("gene", "transcript")) ~ Parent,
      type == "transcript" ~ ID,
      TRUE ~ NA_character_
    ))

  # Assign gene IDs based on transcript IDs
  gr_df_with_geneid_1 <- gr_df_with_txid |>
    mutate(gene_id = case_when(
      type == "gene" ~ ID,
      type == "transcript" ~ Parent,
      TRUE ~ NA_character_
    ))

  # Remove gene features and assign gene IDs to remaining features
  gr_df_with_geneid <- gr_df_with_geneid_1 |>
    filter(type != "gene") |>
    group_by(transcript_id) |>
    mutate(gene_id = na.omit(unique(gene_id)))

  gr_df_genes <- gr_df_with_geneid_1 |>
    filter(type == "gene")

  # convert to gr
  final_gr <- bind_rows(gr_df_with_geneid, gr_df_genes) |>
    makeGRangesFromDataFrame(keep.extra.columns = T) |> sort() |>
    add_unique_id() |> rename_Name_based_on_gene_id()

  return(final_gr)
}



#' Rename "Name" Field Based on Gene ID
#'
#' This function updates the `Name` column of a `GRanges` object or data frame to the `Name` corresponding to the "gene" type within each `gene_id` group. If the `Name` field is missing (NA), it can optionally be filled with the `gene_id` value.
#'
#' @param gr A `GRanges` object or a data frame containing the genomic annotation data.
#' @param convert_NA_to_gene_id Logical value indicating whether to convert missing (NA) `Name` values to `gene_id`. Defaults to `TRUE`.
#'
#' @return If `convert_NA_to_gene_id` is `TRUE`, returns a tibble with the updated `Name` values. Otherwise, returns a `GRanges` object with extra columns.
#' @export
#'
#' @examples
#' # Example with GRanges object
#' gr <- rtracklayer::import.gff3("example.gff3")
#' result <- rename_Name_based_on_gene_id(gr)
#'
#' # Example with data frame input
#' df <- data.frame(
#'   gene_id = c("g1", "g1", "g2"),
#'   type = c("gene", "exon", "gene"),
#'   Name = c("gene1", "exon1", "gene2")
#' )
#' result <- rename_Name_based_on_gene_id(df, convert_NA_to_gene_id = FALSE)
rename_Name_based_on_gene_id <- function(gr, convert_NA_to_gene_id = TRUE) {
  # Convert the GRanges or input to a tibble and group by gene_id
  gr_df <- gr |> as_tibble() |>
    group_by(gene_id) |>
    mutate(Name = first(Name[type == "gene"])) |>
    ungroup()

  # Conditionally handle missing (NA) values in the Name column
  if (convert_NA_to_gene_id) {
    gr_df <- gr_df |>
      mutate(Name = case_when(
        is.na(Name) ~ gene_id,
        TRUE ~ Name
      ))
  }

  # Return either as a tibble or as a GRanges object
  return(makeGRangesFromDataFrame(gr_df, keep.extra.columns = TRUE))
}

