qenv <- new.env()

#' main function to perform calculation of mean relative expression of qpcr data
#' @param df data frame of qpcr data
#' @param hkg character vector of housekeeping genes
#' @param efficiency named list of efficiency values
#' @examples
#' qpcr(df, hkg = c("HKG"))
#' @export
qpcR <- function(df, hkg = NULL, efficiency = NULL, cmp = "control", groups = NULL) {
    assign("control_group", cmp, qenv)
    prep_data <- prepare(df)
    requested_groups <- make_groups(prep_data, groups)
    data_clean <- cleanup(requested_groups, hkg)
    e_val <- unlist(lapply(set_efficiency(unique(df$gene), efficiency), get_e))
    rel_expr <- pair_wise(data_clean, e_val, hkg)
    names(rel_expr) <- setdiff(unique(df$gene), hkg)
    return(rel_expr)
}