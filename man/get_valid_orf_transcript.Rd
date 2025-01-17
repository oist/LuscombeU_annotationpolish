% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_valid_ORF_transcripts.R
\name{get_valid_orf_transcript}
\alias{get_valid_orf_transcript}
\title{Extract and Refine Transcripts with Valid ORFs}
\usage{
get_valid_orf_transcript(gr)
}
\arguments{
\item{gr}{A GenomicRanges object, usually from liftoff containing the genomic data with fields for \code{valid_ORF}, \code{transcript_id}, \code{type}, and \code{gene_id}. This should include annotations that distinguish different genomic features such as genes and transcripts.}
}
\value{
A GenomicRanges object that has been filtered and adjusted to include only transcripts with valid ORFs, along with their corresponding genes. The boundaries of genes are adjusted to encompass all associated transcripts.
}
\description{
This function processes genomic range data to select and refine transcripts that contain valid open reading frames (ORFs). It filters the genomic ranges to include only those transcripts marked with a valid ORF and ensures that the gene boundaries are adjusted accordingly.
}
\examples{
# Assume `gr` is a GenomicRanges object pre-loaded with necessary metadata
# Example usage:
refined_gr <- get_valid_orf_transcript(gr)

}
