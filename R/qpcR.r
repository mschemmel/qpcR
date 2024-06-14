qenv <- new.env()

#' Main function to perform calculation of relative expression of qpcr data
#' @param df data frame of qpcr data
#' @param hkg character vector of housekeeping genes
#' @param reference character vector of reference sample (default = control)
#' @param groups character vector of conditions to group on
#' @param aggregate boolean if output should be aggregated by groups or not (reference will equal 1, default = TRUE)
#' @param outlier boolean if outlier should be detected and removed before calculation of relative expression values (default = TRUE)
#' @param outlier.method character string to choose the method to remove outliers (default: interquartile)
#' @returns data frame of relative expression values compared by provided groups and housekeeping genes
#' @examples
#' qpcr(df, hkg = c("HKG"))
#' @export
qpcR <- function(df, hkg = NULL, reference = "control", groups = NULL, aggregate = TRUE, outlier = TRUE, outlier.method = "interquartile") {
    if(is.null(hkg)) stop("No housekeeping gene provided.")
    assign("reference_group", reference, qenv)
    assign("groups", groups, qenv)
    requested_groups <- make_groups(prepare(df), hkg, groups)
    requested_groups <- filter_outlier(requested_groups, method = outlier.method, do = outlier)
    rel_expr <- do.call("rbind", unname(pair_wise(requested_groups, hkg)))
    return(conflate(rel_expr, do = aggregate))
}

