#' Add Five Prime UTR Annotations Based on CAGE Peaks
#'
#' This function enhances genomic annotations by identifying the closest upstream CAGE peak within a specified distance
#' to the transcription start site (TSS) of each transcript. It adjusts transcript and exon starts and adds five prime UTR
#' annotations accordingly. The function operates on a `GRanges` object containing genomic annotations and another `GRanges`
#' object representing CAGE peaks.
#'
#' @param annotation_gr A `GRanges` object containing genomic annotations with columns for `type`, `transcript_id`, and
#' strand information. It should include annotations for transcripts, genes, and exons.
#' @param cage_peaks A `GRanges` object containing CAGE peaks, with the `thick` column representing the peak regions.
#' @param dist The maximum distance upstream of the transcription start site (TSS) within which to search for CAGE peaks.
#' Default is 500 base pairs.
#'
#' @return A `GRanges` object containing the original genomic annotations plus the newly added five prime UTR annotations
#' and adjusted transcript and exon start positions based on the nearest upstream CAGE peaks within the specified distance.
#'
#' @export
#' @import GenomicRanges
#' @import dplyr
#' @import plyranges
#'
#' @examples
#' # Assuming `annotation_gr` is your GRanges object with transcript, gene, and exon annotations,
#' # and `cage_peaks` is another GRanges object with CAGE peak data:
#' annotation_gr <- system.file("extdata", "OKI.I69_v2_minimized.gff3", package="annotationpolish") |> rtracklayer::import()
#' cage_peaks <- read_cage("https://raw.githubusercontent.com/oist/LuscombeU-CAGE_libraries/main/2021-12-17_Okinawa_Oik/CrossAlignments/consensus_clusters_no_OKItoOKI.bed", "https://raw.githubusercontent.com/oist/LuscombeU-CAGE_libraries/main/2021-12-17_Okinawa_Oik/CrossAlignments/consensus_clusters_sl_OKItoOKI.bed" )
#' updated_annotations <- add_five_prime_utr(annotation_gr, cage_peaks, dist = 500)
#'
#' # This will return a GRanges object with adjusted annotations and added five prime UTRs.
add_five_prime_utr <- function(annotation_gr, cage_peaks, dist = 500){

  # dist should be minus 1
  dist <- dist - 1

  # add unique ID for manipulation
  annotation_gr <- add_unique_id(annotation_gr)

  atx <- annotation_gr |> plyranges::filter(type == "transcript")

  # get the current TSS
  atx_tss <- GenomicRanges::promoters(atx, upstream = 0, downstream = 1)

  # get the max peak
  GenomicRanges::ranges(cage_peaks) <- cage_peaks$thick

  # When there are both SL and non-SL CAGE peaks at the same coordinate for one
  # transcript, pick the non-SL peak
  cage_peaks <- cage_peaks[-subjectHits(findOverlaps(cage_peaks |> filter(tss_type == "nosl"),
                                                     cage_peaks |> filter(tss_type == "sl")))]

  # get the closest one upstream within 500bp
  nearest_upstream_peak <- plyranges::join_nearest_downstream(cage_peaks, atx_tss, distance = T) |> plyranges::filter(distance <= dist)

  # start or end does not matter because it is a 1bp range
  nearest_upstream_peak_df <- nearest_upstream_peak |> dplyr::as_tibble() |> dplyr::select(tss = start, tss_type, transcript_id, distance)

  nearest_upstream_peak_df <- nearest_upstream_peak_df |> group_by(transcript_id) |> filter(distance == min(distance))

  # adjust first exon, transcript, and add 5' UTR
  transcripts <- annotation_gr |> tibble::as_tibble() |> dplyr::filter(type == "transcript") |> dplyr::right_join(nearest_upstream_peak_df)

  # create 5' UTR by adjusting transcript models
  five_prime_utr <- transcripts |> mutate(temp_start = start, temp_end = end) |>
    dplyr::mutate(end = dplyr::case_when(strand == "+" ~ temp_start - 1, T ~ tss),
                  start = dplyr::case_when(strand == "+" ~ tss, T ~ temp_end + 1)) |>
    dplyr::mutate(type = "five_prime_UTR", source = "CAGE", score = NA) |>
    dplyr::select(-tss, -temp_start, -temp_end)

  # adjust transcripts
  transcripts <- transcripts |> dplyr::mutate(start = dplyr::case_when(strand == "+" ~ tss, T ~ start),
                                              end = dplyr::case_when(strand == "-" ~ tss, T ~ end)) |>
    dplyr::select(-tss, -tss_type)

  # adjust genes
  genes <- right_join(annotation_gr |> as_tibble() |> select(-start, -end) |> filter(type == "gene"),
                      transcripts |> group_by(gene_id) |> mutate(start = min(start), end = max(end)) |>
                        select(start, end) |> distinct())

  # adjust exons
  exon_sense <- annotation_gr |> plyranges::filter(type == "exon", strand == "+") |> dplyr::as_tibble() |>
    dplyr::group_by(transcript_id) |> arrange(start) |> slice_head(n = 1)

  exon_antisense <- annotation_gr |> plyranges::filter(type == "exon", strand == "-") |> dplyr::as_tibble() |>
    dplyr::group_by(transcript_id) |> arrange(desc(end)) |> slice_head(n = 1)


  # get the new TSS for both antisense and sense strands
  transcripts_sense <- transcripts |> filter(strand == "+") |> select(tss = start, transcript_id)
  transcripts_antisense <- transcripts |> filter(strand == "-") |> select(tss = end, transcript_id)

  # modify the start and ends to reflect the new TSS
  exon_sense <- exon_sense |> right_join(transcripts_sense) |> mutate(start = tss) |> select(-tss)
  exon_antisense <- exon_antisense |> right_join(transcripts_antisense) |> mutate(end = tss) |> select(-tss)

  exons <- bind_rows(exon_sense, exon_antisense)

  # return those removed before (non first exons)
  other_features <- annotation_gr |> filter(transcript_id %in% transcripts$transcript_id, !type %in% c("gene", "transcript"), !ID %in% exons$ID) |> as_tibble()

  all_features <- bind_rows(exons, other_features) |> ungroup()

  # combine everything together again
  granges_with_5utr <- bind_rows(transcripts, genes, five_prime_utr, all_features) |>
    select(seqnames, start, end, width, strand, source, type, score, phase, gene_id, transcript_id, Parent, tss_type) |>
    makeGRangesFromDataFrame(keep.extra.columns = T)

  # combine eveerything and sort
  final_gr <- c(granges_with_5utr, annotation_gr |> plyranges::filter(!transcript_id %in% granges_with_5utr$transcript_id)) |>
    plyranges::select(-ID) |>
    sort()

  # unlist Parent and remove_transcript_id if gene
  mcols(final_gr)$Parent <- sapply(mcols(final_gr)$Parent, function(x) {
    if (is.null(x) || length(x) == 0) {
      return(NA)
    } else {
      return(x[[1]])
    }
  })

  final_gr$transcript_id[final_gr$type == "gene"] <- NA

  # fix stupid bug where Parent of five prime UTR is not transcript ID
  final_gr$Parent[final_gr$type == "five_prime_UTR"] <- final_gr$transcript_id[final_gr$type == "five_prime_UTR"]

  return(final_gr)
}









