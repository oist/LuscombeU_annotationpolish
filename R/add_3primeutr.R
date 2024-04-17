#' Scan for Motif Occurrences in a Genome
#'
#' This function searches for all occurrences of a specified DNA motif and its reverse complement across all chromosomes of a provided genome. It utilizes the `BSgenome` package to access genomic sequences and the `Biostrings` package for pattern matching. The function is designed to be strand-aware, identifying motif occurrences on both the forward and reverse strands.
#'
#' @param genome A `BSgenome` object representing the genome to be scanned. This object should be loaded and passed to the function. The genome object must be from the `BSgenome` package, which provides comprehensive genomic sequences for various organisms.
#' @param motif A character string representing the DNA motif to be searched for. The default motif is "AATAAA", a common polyadenylation signal in many genomes. The function will also search for the reverse complement of this motif.
#'
#' @return A `GRanges` object containing the genomic locations of all occurrences of the specified motif and its reverse complement. Each location is annotated with the chromosome (seqnames), start and end positions (ranges), and the strand ("+" for the forward strand and "-" for the reverse strand) on which the motif occurs.
#'
#' @export
#' @import Biostrings
#' @import GenomicRanges
#'
#' @examples
#' # Load the BSgenome object for a specific organism (Example: Okinawa O. dioica)
#' genome <- BSgenome::getBSgenome("BSgenome.Oidioi.OIST.OKI2018.I69")
#'
#' # Scan for the "AATAAA" motif and its reverse complement in the human genome
#' motif_occurrences <- scan_motif(genome, motif = "AATAAA")
#' print(motif_occurrences)
scan_motif <- function(genome, motif = "AATAAA", ...){
  # Define the motif and its reverse complement
  motif_revcomp <- reverseComplement(DNAString(motif))

  # Function to search for the motif in both directions on a given chromosome
  search_motif_both_strands <- function(chromosome) {
    # Search for motif in the forward direction
    forward_matches <- matchPattern(motif, genome[[chromosome]], ...)
    # Search for motif in the reverse direction (reverse complement)
    reverse_matches <- matchPattern(motif_revcomp, genome[[chromosome]], ...)

    # Check if there are no matches in both directions
    if (length(forward_matches) == 0 && length(reverse_matches) == 0) {
      return(NULL)  # Return NULL to indicate no matches found
    }

    # Create GRanges objects for matches, checking for non-zero lengths
    forward_gr <- if (length(forward_matches) > 0) {
      GRanges(seqnames = Rle(chromosome), ranges = ranges(forward_matches), strand = "+")
    } else NULL

    reverse_gr <- if (length(reverse_matches) > 0) {
      GRanges(seqnames = Rle(chromosome), ranges = ranges(reverse_matches), strand = "-")
    } else NULL

    # Combine the results, excluding NULL values
    combined_gr <- c(forward_gr, reverse_gr)
    return(combined_gr)
  }


  # Apply the search function across chromosomes (example uses a subset for speed)
  # Consider parallelizing or running on a high-performance computing environment for full genomes
  all_matches <- lapply(seqnames(genome), search_motif_both_strands)

  # remove NULL
  all_matches <- all_matches[!sapply(all_matches, is.null)]

  # unlist some that are somehow structured
  for (i in seq_along(all_matches)) {
    all_matches[i] <- all_matches[i] |> unlist()
  }

  all_matches <- all_matches |>
    GRangesList() |> unlist()

  return(all_matches)
}



#' Calculate Regions Around Feature Ends
#'
#' Generates new genomic regions centered around the ends of features in a given `GRanges` object.
#' The function adjusts these regions based on specified upstream and downstream distances,
#' taking into account the orientation of each feature (strand). The resulting regions aim to
#' capture areas immediately adjacent to the feature ends, potentially useful for analyzing
#' regulatory elements or extending annotations.
#'
#' @param gr A `GRanges` object containing genomic features for which to calculate end-centered regions.
#' @param upstream The distance upstream of the feature end (towards the 5' end for positive strands,
#' towards the 3' end for negative strands) to include in the new region. Default is 500 bases.
#' @param downstream The distance downstream of the feature end (towards the 3' end for positive strands,
#' towards the 5' end for negative strands) to include in the new region. Default is 500 bases.
#'
#' @return A `GRanges` object containing the new regions around the ends of the original features,
#' including all metadata (elementMetadata) from the input `GRanges` object. For positive-strand features,
#' the regions extend upstream and downstream from the feature's end; for negative-strand features,
#' the regions extend upstream and downstream from the feature's start.
#'
#' @export
#' @import IRanges
#'
#' @examples
#' library(GenomicRanges)
#' # Create a GRanges object with example features
#' gr <- GRanges(seqnames = c("chr1", "chr1"),
#'               ranges = IRanges(start = c(100, 500), end = c(200, 600)),
#'               strand = c("+", "-"),
#'               mcols = DataFrame(gene_id = c("gene1", "gene2")))
#' # Calculate regions around feature ends
#' regions_around_end_gr <- get_regions_around_end(gr, upstream = 500, downstream = 500)
#' print(regions_around_end_gr)
get_regions_around_end <- function(gr, upstream = 500, downstream = 500) {
  # Ensure 'gr' is a GRanges object
  stopifnot(is(gr, "GRanges"))

  # Initialize a list to collect results
  regions_list <- list()

  # Process positive strand
  pos_strands <- gr[strand(gr) == "+"]
  if (length(pos_strands) > 0) {
    pos_regions <- GRanges(seqnames = seqnames(pos_strands),
                           ranges = IRanges(start = end(pos_strands) - upstream + 1,
                                            end = end(pos_strands) + downstream),
                           strand = strand(pos_strands))
    mcols(pos_regions) <- mcols(pos_strands)
    regions_list[[1]] <- pos_regions
  }

  # Process negative strand
  neg_strands <- gr[strand(gr) == "-"]
  if (length(neg_strands) > 0) {
    neg_regions <- GRanges(seqnames = seqnames(neg_strands),
                           ranges = IRanges(start = start(neg_strands) - downstream,
                                            end = start(neg_strands) + upstream - 1),
                           strand = strand(neg_strands))
    mcols(neg_regions) <- mcols(neg_strands)
    regions_list[[2]] <- neg_regions
  }

  # Combine results from both strands into one GRanges object
  regions_around_end <- do.call(c, regions_list)

  return(regions_around_end)
}



#' Extend Transcripts to Nearest Downstream UTR Motif
#'
#' This function adjusts the 3' ends of transcripts in an annotation set to extend to the nearest downstream UTR motif, based on a specified maximum distance. It is designed to refine gene models by leveraging external data on UTR motifs, potentially improving annotations for downstream analyses.
#'
#' @param annotation_gr A `GRanges` object containing genomic annotations, specifically transcripts that need to be adjusted. The annotations should include at least the basic features like sequence names, start and end positions, and strand information.
#' @param utr_motif_gr A `GRanges` object representing UTR motifs that serve as potential extension points for the 3' ends of transcripts. These motifs are used to identify the nearest downstream motif for each transcript.
#' @param dist An integer specifying the maximum distance (in base pairs) within which to search for a downstream UTR motif relative to the end of each transcript. The search distance is adjusted internally by subtracting one to ensure precise boundary matching. Default value is 500 bp.
#'
#' @return A modified `GRanges` object containing the original annotations with adjusted 3' ends of transcripts. Each transcript's 3' end is extended to the nearest downstream UTR motif within the specified maximum distance, when possible. The function also integrates new 3' UTR annotations into the dataset and performs necessary adjustments to related genomic features such as exons and genes.
#'
#' @export
#'
#' @examples
add_three_prime_utr <- function(annotation_gr, utr_motif_gr, dist = 500){
  # dist should be minus 1
  dist <- dist - 1

  # add unique ID for manipulation
  annotation_gr <- add_unique_id(annotation_gr)

  atx <- annotation_gr |> plyranges::filter(type == "transcript")

  # get the current gene ends
  atx_end <- get_regions_around_end(atx, upstream = 1, downstream = 0)

  # get the 3' ends of canonical utr motif
  utr_end <- get_regions_around_end(utr_motif_gr, upstream = 1, downstream = 0)

  # get the closest one downstream within 500bp
  nearest_downstream_utr <- plyranges::join_nearest_upstream(utr_end, atx_end, distance = T) |> plyranges::filter(distance <= dist)

  nearest_downstream_utr_df <- nearest_downstream_utr |> dplyr::as_tibble() |> dplyr::select(txend = start, transcript_id, distance)

  nearest_downstream_utr_df <- nearest_downstream_utr_df |> group_by(transcript_id) |> filter(distance == min(distance))

  # adjust first exon, transcript, and add 3' UTR
  transcripts <- annotation_gr |> tibble::as_tibble() |> dplyr::filter(type == "transcript") |> dplyr::right_join(nearest_downstream_utr_df)

  # create 3' UTR by adjusting transcript models
  three_prime_utr <- transcripts |> mutate(temp_start = start, temp_end = end) |>
    dplyr::mutate(start = dplyr::case_when(strand == "+" ~ temp_end + 1, T ~ txend),
                  end = dplyr::case_when(strand == "+" ~ txend, T ~ temp_start - 1)) |>
    dplyr::mutate(type = "three_prime_UTR", source = "genome", score = NA) |>
    dplyr::select(-txend, -temp_start, -temp_end)

  # adjust transcripts
  transcripts <- transcripts |> dplyr::mutate(start = dplyr::case_when(strand == "-" ~ txend, T ~ start),
                                              end = dplyr::case_when(strand == "+" ~ txend, T ~ end)) |>
    dplyr::select(-txend)

  # adjust genes based on min max of each gene
  genes <- right_join(annotation_gr |> as_tibble() |> select(-start, -end) |> filter(type == "gene"),
                      transcripts |> group_by(gene_id) |> mutate(start = min(start), end = max(end)) |>
                        select(start, end) |> distinct()) |> mutate(transcript_id = NA) |> distinct()

  # adjust exons should be OK
  exon_sense <- annotation_gr |> plyranges::filter(type == "exon", strand == "+") |> dplyr::as_tibble() |>
    dplyr::group_by(transcript_id) |> arrange(desc(end)) |> slice_head(n = 1)

  exon_antisense <- annotation_gr |> plyranges::filter(type == "exon", strand == "-") |> dplyr::as_tibble() |>
    dplyr::group_by(transcript_id) |> arrange(start) |> slice_head(n = 1)


  # get the new TSS for both antisense and sense strands
  transcripts_sense <- transcripts |> filter(strand == "+") |> select(txend = end, transcript_id)
  transcripts_antisense <- transcripts |> filter(strand == "-") |> select(txend = start, transcript_id)

  # modify the start and ends to reflect the new TSS
  exon_sense <- exon_sense |> right_join(transcripts_sense) |> mutate(end = txend) |> select(-txend)
  exon_antisense <- exon_antisense |> right_join(transcripts_antisense) |> mutate(start = txend) |> select(-txend)

  exons <- bind_rows(exon_sense, exon_antisense)

  # return those removed before (non first exons)
  other_features <- annotation_gr |> filter(transcript_id %in% transcripts$transcript_id, !type %in% c("gene", "transcript"), !ID %in% exons$ID) |> as_tibble()

  all_features <- bind_rows(exons, other_features) |> ungroup()

  # combine everything together again
  granges_with_3utr <- bind_rows(transcripts, genes, three_prime_utr, all_features) |>
    select(seqnames, start, end, width, strand, source, type, score, phase, gene_id, transcript_id, Parent, tss_type) |>
    makeGRangesFromDataFrame(keep.extra.columns = T)

  # get those genes without detectable 3' UTR
  # when there are isoforms with 3' UTR, remove it from the noutrgenes
  noutrgenes <- annotation_gr |> plyranges::filter(!transcript_id %in% granges_with_3utr$transcript_id)
  noutrgenes <- annotation_gr |> filter(gene_id %in% unique(noutrgenes$gene_id)) |>
    filter(!transcript_id %in% granges_with_3utr$transcript_id)

  # combine eveerything and sort
  final_gr <- c(granges_with_3utr,
                noutrgenes) |>
    plyranges::select(-ID) |>
    sort()

  # remove all gene entries and add again
  final_gr <- c(final_gr |> filter(type != "gene"), annotation_gr |> filter(type == "gene"))

  # fix genes again by adjusting the start and ends
  final_gr <- bind_rows(
    final_gr |> as_tibble() |> group_by(gene_id) |> mutate(start = min(start), end = max(end)) |> filter(type == "gene"),
    final_gr |> as_tibble() |> filter(type != "gene")
  ) |> makeGRangesFromDataFrame(keep.extra.columns = T) |> sort()

  # unlist Parent and remove_transcript_id if gene
  mcols(final_gr)$Parent <- sapply(mcols(final_gr)$Parent, function(x) {
    if (is.null(x) || length(x) == 0) {
      return(NA)
    } else {
      return(x[[1]])
    }
  })

  final_gr$transcript_id[final_gr$type == "gene"] <- NA

  # fix stupid bug where Parent of three prime UTR is not transcript ID
  final_gr$Parent[final_gr$type == "three_prime_UTR"] <- final_gr$transcript_id[final_gr$type == "three_prime_UTR"]

  return(final_gr)
}





#' Extend 3' UTRs to Nearest Downstream Motif
#'
#' This function extends the 3' untranslated regions (UTRs) of transcripts in a given genomic annotation to the nearest downstream instance of a specified motif, within a maximum specified distance. It first scans the provided genome for occurrences of the motif to create a set of potential UTR extension points. Then, it adjusts the transcript models in the annotation, extending the 3' UTRs to these points where applicable.
#'
#' @param annotation_gr A `GRanges` object containing the genomic annotations, typically including transcripts for which 3' UTRs need to be extended. This object should be the result of importing genomic annotation data, such as from a GFF3 file.
#' @param genome A `BSgenome` object representing the genome to be scanned for the specified motif. This object is expected to encapsulate the entire genomic sequence for the organism of interest, facilitating the search for motif occurrences.
#' @param motif A character string representing the DNA sequence motif to search for within the genome as potential 3' UTR extension points. Default is "AATAAA", a common polyadenylation signal motif.
#' @param dist An integer specifying the maximum distance (in base pairs) to search downstream of each transcript's 3' end for the motif. The function will only consider motifs within this distance for extending the 3' UTR. Default value is 500 bp.
#'
#' @return The original `GRanges` object (annotation_gr) with modified transcript models, where 3' UTRs have been potentially extended to the nearest downstream motif within the specified distance. The function does not explicitly return the object; ensure to capture the result into a variable.
#'
#' @export
#'
#' @examples
#' # Assuming 'annotation_gr' is loaded from a GFF3 file and 'BSgenome.Oidioi.OIST.OKI2018.I69' is loaded
#' updated_annotations <- add_three_prime_utr_wrapper(annotation_gr, BSgenome.Oidioi.OIST.OKI2018.I69, motif = "AATAAA", dist = 500)
#' # 'updated_annotations' now contains transcripts with potentially extended 3' UTRs.
add_three_prime_utr_wrapper <- function(annotation_gr, genome, motif = "AATAAA", dist = 500){

  utr_motif_gr <- scan_motif(genome, motif)
  add_three_prime_utr(annotation_gr, utr_motif_gr)

}


