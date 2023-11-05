qenv <- new.env()

#' Main function to perform calculation of relative expression of qpcr data
#' @param df data frame of qpcr data
#' @param hkg character vector of housekeeping genes
#' @param efficiency named list of efficiency values
#' @examples
#' qpcr(df, hkg = c("HKG"))
#' @export
qpcR <- function(df, hkg = NULL, reference = "control", groups = NULL) {
    assign("reference_group", reference, qenv)
    assign("groups", groups, qenv)
    prep_data <- prepare(df)
    requested_groups <- make_groups(prep_data, groups)
    clean_pairs <- sanitize(requested_groups, hkg)
    rel_expr <- pair_wise(clean_pairs, hkg)
    return(do.call("rbind", unname(rel_expr)))
}