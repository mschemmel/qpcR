qenv <- new.env()

#' main function to perform calculation of relative expression of qpcr data
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
    data_clean <- cleanup(requested_groups, hkg)
    rel_expr <- pair_wise(data_clean, hkg)
    return(do.call("rbind", unname(rel_expr)))
}