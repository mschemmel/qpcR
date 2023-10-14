
# main function
#' @export
qpcR <- function(df, hkg = NULL, efficiency = NULL) {
    dat_clean <- cleanup(prepare(df), hkg)
    # calculate E value
    e_val <- unlist(lapply(set_efficiency(unique(df$gene), efficiency), get_e))
    # get mean of control groups
    rel_expr <- pair_wise(dat_clean, e_val, hkg)
    names(rel_expr) <- setdiff(unique(df$gene), hkg)
    rel_expr
}