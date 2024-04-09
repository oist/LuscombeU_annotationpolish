#' Title: Fix Liftoff data from Oikobase gff
#'
#' Description: This function processes a GenomicRanges object (or compatible tibble) to standardize and fix IDs for genes and transcripts based on specified naming conventions. It performs several operations including removing descriptions from names, changing transcript types, assigning gene and transcript IDs, converting transcript IDs to gene IDs (and vice versa), selecting important columns, and dealing with parent-child relationships in genomic data. The output is a cleaned and standardized GenomicRanges object ready for further analysis.
#'
#' @param gr A GenomicRanges object or a tibble that represents genomic features with columns like type, Name, ID, and potentially others related to genomic annotations. This input is expected to be the result of a Liftoff process or similar annotation transfer tool.
#'
#' @return A GenomicRanges object with standardized and fixed IDs, and cleaned genomic feature annotations. It includes operations like renaming, ID conversion, and the establishment of gene-transcript relationships. The function also ensures the retention of important columns and adds a unique ID for each feature.
#'
#' @export
#'
#' @examples
#' liftoff <- system.file("extdata", "split_example_liftoff.gff3", package="annotationpolish") |> rtracklayer::import()
#' gr_cleaned <- fix_liftoff_oikobase(liftoff)
#' print(gr_cleaned)
fix_liftoff_oikobase <- function(gr){
  # work with tibble because it is much easier
  gr <- gr |> as_tibble() |>
    # remove description
    mutate(ID = case_when(type == "gene" ~ Name |> stringr::str_remove(":.*") |> stringr::str_remove(" .*"), T ~ ID)) |>
    mutate(Name = Name |> stringr::str_remove(":.*") |> stringr::str_remove(" .*")) |>
    # change transcript type
    mutate(type = case_when(type == "mRNA" ~ "transcript", T ~ type)) |>
    # assign gene_id and transcirpt_id
    mutate(gene_id = Name, transcript_id = Name) |> select(-Name) |>
    # use string manipulation to convert transcript id to gene id and vv
    mutate(gene_id = case_when(type == "gene" ~ stringr::str_replace(gene_id, "GSOIDT", "GSOIDG"))) |>
    mutate(transcript_id = case_when(type == "transcript" ~ stringr::str_replace(transcript_id, "GSOIDG", "GSOIDT"))) |>
    # select important columns
    select(1:10, coverage, Parent, valid_ORF, gene_id, transcript_id) |>
    # assign ID
    mutate(ID = case_when(type == "gene" ~ ID |> stringr::str_replace("GSOIDT", "GSOIDG"), T ~ stringr::str_replace(ID, "GSOIDG", "GSOIDT")))

  # unlist Parent
  gr$Parent <- sapply(gr$Parent, function(x) {
    if (is.null(x) || length(x) == 0) {
      return(NA)
    } else {
      return(x[[1]])
    }
  })


  gr <- gr |> mutate(Parent = case_when(type == "gene" ~ NA,
                                 type == "transcript" ~ stringr::str_replace(ID, "GSOIDT", "GSOIDG"),
                                 T ~ stringr::str_replace(Parent, "GSOIDG", "GSOIDT"))) |>
    mutate(gene_id = case_when(type != "gene" ~ stringr::str_replace(Parent, "GSOIDT", "GSOIDG"), T ~ gene_id)) |>
    mutate(transcript_id = case_when(type != "gene" ~ stringr::str_replace(gene_id, "GSOIDG", "GSOIDT"), T ~ gene_id)) |>
    distinct() |>
    makeGRangesFromDataFrame(keep.extra.columns = T)

  # assign new ID
  gr <- gr |> select(-ID) |> add_unique_id()

  return(gr)
}

