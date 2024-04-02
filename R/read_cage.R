#' Read CAGE (Cap Analysis Gene Expression) Data
#'
#' This function imports CAGE data from provided BED files. It supports importing data
#' with or without SL (Spliced Leader) sequences. For files with SL sequences, the function
#' can combine and annotate the two types of sequences (SL and non-SL) for comprehensive analysis.
#'
#' @param no_sl A character string specifying the path or URL to the BED file containing
#' CAGE data without SL sequences.
#' @param sl An optional character string specifying the path or URL to the BED file containing
#' CAGE data with SL sequences. If provided, data from this file will be combined with the
#' non-SL data, and each sequence will be annotated with its type.
#'
#' @return An object of class `GRanges` representing the imported CAGE data. If both SL and
#' non-SL data are provided, the return value will be a combined `GRanges` object with an
#' additional column `tss_type` indicating whether a sequence is 'sl' (with Spliced Leader)
#' or 'nosl' (without Spliced Leader).
#'
#' @export
#'
#' @examples
#' no_sl_url <- "https://raw.githubusercontent.com/oist/LuscombeU-CAGE_libraries/main/2021-12-17_Okinawa_Oik/CrossAlignments/consensus_clusters_no_OKItoOKI.bed"
#' sl_url <- "https://raw.githubusercontent.com/oist/LuscombeU-CAGE_libraries/main/2021-12-17_Okinawa_Oik/CrossAlignments/consensus_clusters_sl_OKItoOKI.bed"
#'
#' # Reading non-SL CAGE data only
#' read_cage(no_sl_url)
#'
#' # Reading and combining non-SL and SL CAGE data
#' read_cage(no_sl_url, sl_url)
#'
read_cage <- function(no_sl, sl = NULL) {
  message(paste("reading", no_sl))

  no_sl_gr <- rtracklayer::import.bed(no_sl) |>
    GenomicRanges::GRanges() |> plyranges::select(thick)

  # if there is SL file, also import and combine
  if (!is.null(sl)) {
    message(paste("SL file present, reading", no_sl))
    sl_gr <- rtracklayer::import.bed(sl) |> GenomicRanges::GRanges() |>
      plyranges::select(thick) |> plyranges::mutate(tss_type = "sl")
    no_sl_gr <- no_sl_gr |> plyranges::mutate(tss_type = "nosl")

    # return the combined granges
    return(c(sl_gr, no_sl_gr))
  } else {
    return(no_sl_gr)
  }
}
