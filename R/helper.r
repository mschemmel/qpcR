#' prepare input data
#' @param df input data.frame of qPCR data
#' @examples
#' prepare(df)
prepare <- function(df) {
    df <- df[c("gene", "treatment", "cq")]
    df$cq <- as.numeric(gsub(",", ".", df$cq))
    return(df)
}

#' calculate E value of standard samples
#' @param efficiency_of_standard efficiency values in percent (list, 0-100)
#' @examples
#' get_e(named_list)
get_e <- function(efficiency_of_standard) {
    return((efficiency_of_standard * 0.01) + 1)
}

set_efficiency <- function(genes, efficiency, default = 100) {
    no_known_efficiency <- setdiff(genes, names(efficiency))
    default_efficiencies <- NULL
    if (!identical(no_known_efficiency, character(0))) {
        default_efficiencies <- rep(default, length(setdiff(genes, names(efficiency))))
        names(default_efficiencies) <- no_known_efficiency
    }
    return(c(efficiency, default_efficiencies))
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
    targets <- setdiff(names(e_val), hkg_)
    cmp <- data.frame(sapply(targets, function(x) {
                e_val[[x]]^d_cq[[x]] / e_val[[hkg_]]^d_cq[[hkg_]]
    })) #TODO: preserve names during lapply

    cpratiodata <- cbind(treatment = treatm, cmp)
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
    return(mean_aggr)
}