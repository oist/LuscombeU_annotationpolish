#' Create new isoforms from overlapping gene models
#'
#' This function identifies overlaps between gene models annotated by different sources
#' (e.g., AUGUSTUS and Liftoff) and creates new isoforms by merging overlapping genes.
#' The new isoforms are assigned unique identifiers. The function returns a GenomicRanges
#' object containing the updated gene and isoform annotations.
#'
#' @param gr GenomicRanges object containing gene annotations with metadata columns.
#'           Expected to have columns like `source`, `type`, `gene_id`, and `Parent`.
#'
#' @return GenomicRanges object with updated annotations including new isoforms.
#' @export
#' @import plyranges
#' @import tidyverse
#'
#' @examples
#' gr <- makeGRangesFromDataFrame(your_dataframe)
#' result <- make_new_isoforms(gr)
make_new_isoforms <- function(gr){
  # Identify overlaps between AUGUSTUS and Liftoff gene annotations
  overlap_1 <- join_overlap_intersect_within_directed(gr |> filter(source == "AUGUSTUS", type == "gene"),
                                                      gr |> filter(source == "Liftoff", type == "gene")) |>
    as_tibble() |> select(ID.x, ID.y) |> mutate(isoform_ID = paste(ID.x, ID.y, sep = "_"))



  overlap_2 <- join_overlap_intersect_within_directed(gr |> filter(source == "Liftoff", type == "gene"),
                                                      gr |> filter(source == "AUGUSTUS", type == "gene")) |>
    as_tibble() |> select(ID.x, ID.y) |> mutate(isoform_ID = paste(ID.x, ID.y, sep = "_"))



  # Return the original data if no overlaps are found
  all_overlaps <- bind_rows(overlap_1, overlap_2)

  if (length(all_overlaps) == 0) {
    return(gr)
  }

  # Convert GenomicRanges to a data frame for further manipulation
  gr_df <- gr |> as_tibble()


  # Transform and subset data to focus on genes involved in overlaps
  df1_long <- all_overlaps %>%
    pivot_longer(cols = c(ID.x, ID.y), names_to = "ID_type", values_to = "gene_id")


  # Separate data into those with new isoforms and those without
  with_isoforms <- gr_df |> filter(gene_id %in% df1_long$gene_id)
  remaining <- gr_df |> filter(!gene_id %in% df1_long$gene_id)
  remaining <- remaining |> makeGRangesFromDataFrame(keep.extra.columns = T) |>
    unlist_parent() |> as_tibble()

  # Integrate new isoform IDs and adjust the 'Parent' column for transcripts
  df2_mutated <- with_isoforms |> mutate(Parent = unlist(Parent)) %>%
    left_join(df1_long, by = "gene_id") %>%
    mutate(Parent = if_else(type == "transcript" & !is.na(isoform_ID), isoform_ID, Parent),
           gene_id = if_else(type == "gene" & !is.na(isoform_ID), isoform_ID, gene_id)) %>%
    select(-ID_type, -isoform_ID)  # Remove the extra columns from the join

  # Create new gene features from transcripts and ensure correct IDs
  df2_tx <- df2_mutated |> filter(type == "transcript")
  df2_non_tx <- df2_mutated |> filter(type != "transcript")


  # change the gene_id for df2_tx and df2_non_tx
  df2_tx <- df2_tx |> left_join(df1_long |> select(1, gene_id)) |> mutate(gene_id = isoform_ID) |> select(-isoform_ID)
  df2_non_tx <- df2_non_tx |> left_join(df1_long |> select(1, gene_id)) |> mutate(gene_id = isoform_ID) |> select(-isoform_ID)


  # Consolidate gene records from transcripts to define new gene boundaries
  df2_gene <- df2_tx |> group_by(Parent) |> mutate(start = min(start), end = max(end)) |>
    slice_head(n = 1) |> ungroup() |>
    mutate(source = "annotationpolish", type = "gene", transcript_id = NA, gene_id = Parent, ID = Parent) |>
    select(seqnames, start, end, width, strand, source, type, gene_id, ID)

  # Combine all modifications into a final data frame
  result_df <- bind_rows(df2_gene, df2_tx, df2_non_tx, remaining)

  # Convert the result back to a GenomicRanges object for output
  result_gr <- result_df |> makeGRangesFromDataFrame(keep.extra.columns = T) |>
    filter(!is.na(gene_id)) |>
    adjust_transcript_gene_boundaries()

  return(result_gr |> sort())

}




