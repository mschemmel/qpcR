cleanup <- function(df, hkg_) {
    lapply(setdiff(unique(df$gene), hkg_), function(x) {
        pair_comp <- df[df$gene %in% c(x, hkg_), ]
        subset(pair_comp[!(pair_comp$id %in% pair_comp[is.na(pair_comp$cq), ]$id), ], select = -c(brep, trep, id))
    })
}

control_mean <- function(df, treat = "control") {
    contr <- df[df$treatment == treat, ]
    return(tapply(contr$cq, contr$gene, mean))
}

delta_cq <- function(df, contr_mean) {
    dcq <- lapply(seq_along(contr_mean), function(x) {
        contr_mean[x] - df[df$gene == names(contr_mean)[x], ]$cq
    }) #TODO: preserve names during lapply
    names(dcq) <- names(contr_mean)
    return(dcq)
}

ratio_by_mean_ratio <- function(d_cq, e_val, hkg_, treatm) {
    target <- setdiff(names(e_val), hkg_)
    cmp <- e_val[[target]]^d_cq[[target]] / e_val[[hkg_]]^d_cq[[hkg_]]
    cpratiodata <- data.frame(treatment = treatm, cmp = cmp)
    cpratio_control <- cpratiodata[cpratiodata$treatment == "control", ][names(cpratiodata) != "treatment"]
    cpratiodata <- cpratiodata[names(cpratiodata) != "treatment"]
    mean_ratio <- lapply(cpratio_control, mean)
    rbyr <- data.frame(sapply(seq_along(mean_ratio), function(x) {
                cpratiodata[names(cpratiodata) == names(mean_ratio)[x]] / mean_ratio[x]
    }))
    return(cbind(treatment = treatm, rbyr))
}

#' @import stats
mean_relative_expression <- function(df) {
    mean_aggr <- stats::aggregate(df[, -1], by = list(df$treatment), FUN = mean, simplify = TRUE)
    names(mean_aggr)[names(mean_aggr) == "Group.1"] <- "treatment"
    names(mean_aggr)[names(mean_aggr) == "x"] <- "relative mean expression"
    return(mean_aggr)
}

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