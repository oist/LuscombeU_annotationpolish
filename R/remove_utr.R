#' Remove UTRs and Adjust Transcript and Gene Boundaries
#'
#' This function takes a `GRanges` object containing genomic annotations and removes UTR annotations from it.
#' For transcripts associated with UTRs, it adjusts their start and end positions based on the boundaries of the CDS.
#' The function aims to retain only the coding sequences (CDS) and other relevant annotations while discarding UTRs
#' and adjusting the genomic ranges of affected transcripts and genes accordingly.
#'
#' @param gr A `GRanges` object containing genomic annotations with at least the columns 'type' and 'transcript_id'.
#'          The 'type' column should indicate the annotation type (e.g., "utr", "transcript", "CDS"). The function
#'          identifies UTRs by detecting 'utr' in the 'type' field, removes these annotations, and then adjusts
#'          the start and end positions of transcripts and genes based on the boundaries of the remaining annotations.
#' @return A `GRanges` object similar to the input but with UTR annotations removed and transcript as well as
#'         gene boundaries adjusted based on the remaining annotations. If no UTRs are found, it returns the original `GRanges` object.
#' @export
#' @import GenomicRanges
#' @import dplyr
#' @import plyranges
#' @examples
#' # Assuming `gr` is your GRanges object with genomic annotations:
#' gr_updated <- remove_utr(gr)
#'
#' # This will return a GRanges object with UTRs removed and adjusted annotations.
remove_utr <- function(gr){
  transcript_ids_utr <- gr %>%
    filter(type %>% stringr::str_detect(regex('utr', ignore_case = TRUE))) %>%
    .$transcript_id

  if (length(transcript_ids_utr) != 0) {
    gene_list_utr <- gr %>%
      filter(transcript_id %in% transcript_ids_utr) %>%
      as_tibble() %>%
      dplyr::group_split(transcript_id)

    noutr_genes <- gene_list_utr %>%
      lapply(function(x) {
        no_utr <- x %>% filter(!type %>% stringr::str_detect(regex('utr', ignore_case = TRUE)))

        if (nrow(no_utr) > 0) {
          start_transcript <- no_utr %>% filter(type == "transcript") %>% pull(start) %>% min(na.rm = TRUE)
          end_transcript <- no_utr %>% filter(type == "transcript") %>% pull(end) %>% max(na.rm = TRUE)

          start_cds <- no_utr %>% filter(type == "CDS") %>% pull(start) %>% min(na.rm = TRUE)
          end_cds <- no_utr %>% filter(type == "CDS") %>% pull(end) %>% max(na.rm = TRUE)

          # mutate if the start and end are at transcript
          no_utr <- no_utr %>%
            mutate(start = ifelse(type == "transcript", start_cds, start),
                   end = ifelse(type == "transcript", end_cds, end))
        }

        return(no_utr)
      }) %>%
      bind_rows() %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    # Return combined GRanges object, excluding original UTRs
    return(c(gr %>% filter(!transcript_id %in% transcript_ids_utr), noutr_genes) %>% sort())
  } else {
    message("No UTR to remove, returning the original GRanges object.")
    return(gr)
  }
}
