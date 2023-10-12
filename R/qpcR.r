
# main function
#' @export
qpcR <- function(df, hkg = NULL, efficiency = NULL) {
    df <- df[c("gene", "treatment", "cq")]
    df$cq <- as.numeric(gsub(",", ".", df$cq))

    # calculate E value
    e_val <- unlist(lapply(set_efficiency(unique(df$gene), efficiency), get_e))

    # get treatment variables
    treats <- df[df$gene == hkg, ]$treatment

    # get mean of control groups
    cmean <- control_mean(df)
    dcq <- delta_cq(df, cmean)

    # calculate ratio by mean ratio of all targets
    mean_ratio <- ratio_by_mean_ratio(dcq, e_val, hkg_ = hkg, treatm = treats)

    # calculate mean relative expression
    return(mean_relative_expression(mean_ratio))
}