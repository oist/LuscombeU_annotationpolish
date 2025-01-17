% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fix_liftoff_oikobase.R
\name{fix_liftoff_oikobase}
\alias{fix_liftoff_oikobase}
\title{Title: Fix Liftoff data from Oikobase gff}
\usage{
fix_liftoff_oikobase(gr)
}
\arguments{
\item{gr}{A GenomicRanges object or a tibble that represents genomic features with columns like type, Name, ID, and potentially others related to genomic annotations. This input is expected to be the result of a Liftoff process or similar annotation transfer tool.}
}
\value{
A GenomicRanges object with standardized and fixed IDs, and cleaned genomic feature annotations. It includes operations like renaming, ID conversion, and the establishment of gene-transcript relationships. The function also ensures the retention of important columns and adds a unique ID for each feature.
}
\description{
Description: This function processes a GenomicRanges object (or compatible tibble) to standardize and fix IDs for genes and transcripts based on specified naming conventions. It performs several operations including removing descriptions from names, changing transcript types, assigning gene and transcript IDs, converting transcript IDs to gene IDs (and vice versa), selecting important columns, and dealing with parent-child relationships in genomic data. The output is a cleaned and standardized GenomicRanges object ready for further analysis.
}
\examples{
liftoff <- system.file("extdata", "split_example_liftoff.gff3", package="annotationpolish") |> rtracklayer::import()
gr_cleaned <- fix_liftoff_oikobase(liftoff)
print(gr_cleaned)
}
