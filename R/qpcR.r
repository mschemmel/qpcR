qenv <- new.env()

#' main function to perform calculation of mean relative expression of qpcr data
#' @param df data frame of qpcr data
#' @param hkg character vector of housekeeping genes
#' @param efficiency named list of efficiency values
#' @examples
#' qpcr(df, hkg = c("HKG"))
#' @export
qpcR <- function(df, hkg = NULL, efficiency = NULL, cmp = "control") {
    assign("control_group", cmp, qenv)

    dat_clean <- cleanup(prepare(df), hkg)
    # calculate E value
    e_val <- unlist(lapply(set_efficiency(unique(df$gene), efficiency), get_e))
    # get mean of control groups
    rel_expr <- pair_wise(dat_clean, e_val, hkg)
    names(rel_expr) <- setdiff(unique(df$gene), hkg)
    data.frame(rel_expr)
}