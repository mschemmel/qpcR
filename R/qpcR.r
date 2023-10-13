
# main function
#' @export
qpcR <- function(df, hkg = NULL, efficiency = NULL) {
    dat_clean <- cleanup(prepare(df), hkg)

    # calculate E value
    e_val <- unlist(lapply(set_efficiency(unique(df$gene), efficiency), get_e))
    # get mean of control groups
    lapply(dat_clean, function(x) {
        cmean <- control_mean(x)
        dcq <- delta_cq(x, cmean)
        # get treatment variables
        treats <- x[x$gene == hkg, ]$treatment
        # calculate ratio by mean ratio of all targets
        mean_ratio <- ratio_by_mean_ratio(dcq, e_val = e_val[names(dcq)], hkg_ = hkg, treatm = treats)
        # calculate mean relative expression
        mean_relative_expression(mean_ratio)
        
    })
}