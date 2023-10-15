#' Remove NA sample wise as foundation of valid mean relative expression calculation
#' @param df initial data frame of expression data
#' @param hkg_ character string of provided housekeeping genes
#' @examples
#' cleanup(df, khg_ = c("HKG"))
cleanup <- function(df, hkg_) {
    lapply(setdiff(unique(df$gene), hkg_), function(x) {
        pair_comp <- df[df$gene %in% c(x, hkg_), ]
        subset(pair_comp[!(pair_comp$id %in% pair_comp[is.na(pair_comp$cq), ]$id), ], select = -c(brep, trep, id))
    })
}

#' get mean expression values of all control groups
#' @param df clean data frame of expression data
#' @examples
#' control_mean(df)
control_mean <- function(df) {
    contr <- get_control_group(df)
    return(tapply(contr$cq, contr$gene, mean))
}

#' calculate delta cq values per comparison
#' @param df clean data frame of expression data
#' @param contr_mean object of control_mean function
#' @examples
#' delta_cq(df, contr_mean)
delta_cq <- function(df, contr_mean) {
    dcq <- lapply(seq_along(contr_mean), function(x) {
        contr_mean[x] - df[df$gene == names(contr_mean)[x], ]$cq
    })
    names(dcq) <- names(contr_mean)
    return(dcq)
}

#' calculate the ratio compared to the mean ratio per gene
#' @param d_cq object of delta_cq
#' @param e_val named list of efficiency values
#' @param hkg_ character vector of housekeeping genes
#' @param treatm character vector of treatments
#' @examples
#' ratio_by_mean_ratio(d_cq, e_val, hkg_, treatm)
ratio_by_mean_ratio <- function(d_cq, e_val, hkg_, treatm) {
    target <- setdiff(names(e_val), hkg_)
    cmp <- e_val[[target]]^d_cq[[target]] / e_val[[hkg_]]^d_cq[[hkg_]]
    cpratio <- data.frame(treatment = treatm, cmp = cmp)
    rbyr <- cpratio$cmp / mean(get_control_group(cpratio)$cmp)
    return(data.frame(treatment = treatm, rbyr = as.numeric(rbyr)))
}

#' calculate mean relative expression
#' @param df data frame of treatment and corresponding ratio values
#' @examples
#' mean_relative_expression(df)
mean_relative_expression <- function(df) {
    return(tapply(df$rbyr, df$treatment, FUN = mean))
}

#' apply calculation for every hkg vs. target pair
#' @param data_list named list of pairs (data)
#' @param e_values E value pair of requested comparison
#' @param hkg character vector of housekeeping genes
#' @examples
#' pair_wise(data_list, e_values, hkg)
pair_wise <- function(data_list, e_values, hkg) {
    lapply(data_list, function(x) {
            cmean <- control_mean(x)
            dcq <- delta_cq(x, cmean)
            # get treatment variables
            treats <- x[x$gene == hkg, ]$treatment
            # calculate ratio by mean ratio of all targets
            mean_ratio <- ratio_by_mean_ratio(dcq, e_val = e_values[names(dcq)], hkg_ = hkg, treatm = treats)
            # calculate mean relative expression
            mean_relative_expression(mean_ratio)
    })
}
