% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_3primeutr.R
\name{add_three_prime_utr}
\alias{add_three_prime_utr}
\title{Extend Transcripts to Nearest Downstream UTR Motif}
\usage{
add_three_prime_utr(annotation_gr, utr_motif_gr, dist = 500)
}
\arguments{
\item{annotation_gr}{A \code{GRanges} object containing genomic annotations, specifically transcripts that need to be adjusted. The annotations should include at least the basic features like sequence names, start and end positions, and strand information.}

\item{utr_motif_gr}{A \code{GRanges} object representing UTR motifs that serve as potential extension points for the 3' ends of transcripts. These motifs are used to identify the nearest downstream motif for each transcript.}

\item{dist}{An integer specifying the maximum distance (in base pairs) within which to search for a downstream UTR motif relative to the end of each transcript. The search distance is adjusted internally by subtracting one to ensure precise boundary matching. Default value is 500 bp.}
}
\value{
A modified \code{GRanges} object containing the original annotations with adjusted 3' ends of transcripts. Each transcript's 3' end is extended to the nearest downstream UTR motif within the specified maximum distance, when possible. The function also integrates new 3' UTR annotations into the dataset and performs necessary adjustments to related genomic features such as exons and genes.
}
\description{
This function adjusts the 3' ends of transcripts in an annotation set to extend to the nearest downstream UTR motif, based on a specified maximum distance. It is designed to refine gene models by leveraging external data on UTR motifs, potentially improving annotations for downstream analyses.
}
