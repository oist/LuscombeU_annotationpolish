#' Extract First Parent from GRanges Object Metadata
#'
#' This function takes a `GRanges` object and processes the `Parent` metadata column. For each entry in the `Parent` column, it extracts the first element if the column contains a list. This is useful in scenarios where the `Parent` metadata is stored as lists of identifiers, and you are interested in simplifying the structure by retaining only the first parent identifier for each range.
#'
#' @param gr A `GRanges` object. The function expects this object to have a `Parent` metadata column, potentially containing list-like structures.
#'
#' @return A modified `GRanges` object with the `Parent` metadata column simplified to contain only the first element of any list present. If the original `Parent` column contains `NULL` or an empty list for a particular range, that entry is replaced with `NA`.
#'
#' @export
#'
#' @examples
#' # Assuming `gr` is a GRanges object with a complex 'Parent' metadata column
#' # Simplify the 'Parent' column to only keep the first parent identifier
#' simplified_gr <- unlist_parent(gr)
#'
#' # Example GRanges object creation (for users unfamiliar with GRanges):
#' # gr <- GRanges(seqnames = "chr1",
#' #               ranges = IRanges(start = c(100, 200), end = c(150, 250)),
#' #               strand = c("+", "-"),
#' #               Parent = list(c("gene1", "gene2"), "gene3"))
#' # This will modify the 'Parent' column of `gr` to have "gene1" for the first range
#' # and "gene3" for the second range.
#' @param gr
#'
#' @return
#' @export
#'
#' @examples
unlist_parent <- function(gr){
  mcols(gr)$Parent <- sapply(mcols(gr)$Parent, function(x) {
    if (is.null(x) || length(x) == 0) {
      return(NA)
    } else {
      return(x[[1]])
    }
  })

  gr
}


