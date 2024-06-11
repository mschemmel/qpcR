qenv <- new.env()

#' Main function to perform calculation of relative expression of qpcr data
#' @param df data frame of qpcr data
#' @param hkg character vector of housekeeping genes
#' @param reference character vector of reference sample
#' @param groups character vector of conditions to group on
#' @examples
#' qpcr(df, hkg = c("HKG"))
#' @export
qpcR <- function(df, hkg = NULL, reference = "control", groups = NULL, aggregate = TRUE) {
    if(is.null(hkg)) stop("No housekeeping gene provided.")
    assign("reference_group", reference, qenv)
    assign("groups", groups, qenv)
    requested_groups <- make_groups(prepare(df), hkg, groups)
    rel_expr <- do.call("rbind", unname(pair_wise(requested_groups, hkg)))
    return(conflate(rel_expr, do = aggregate))
}

