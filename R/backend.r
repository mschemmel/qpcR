#' Create chunks of groups for calculation
#' @param df data frame of expression data
#' @param groups character vector of requested groups to compare
#' @examples
#' make_groups(df, groups = c("plate", "diet", "timepoint"))
make_groups <- function(df, groups) {
    return(split(df, as.list(df[groups]), drop = TRUE))
}

#' Remove NA sample wise as foundation of valid mean relative expression calculation
#' @param list_of_groups output of 'make_groups' containing a list of groups to compare expression
#' @param hkg_ character string of provided housekeeping genes
#' @examples
#' cleanup(df, khg_ = c("HKG"))
cleanup <- function(list_of_groups, hkg_) {
    unlist(lapply(list_of_groups, function(gr) {
        gr$id <- paste0(gr$treatment, gr$brep, gr$trep)
        lapply(setdiff(unique(gr$gene), hkg_), function(x) {
            pair_comp <- gr[gr$gene %in% c(x, hkg_), ]
            subset(pair_comp[!(pair_comp$id %in% pair_comp[is.na(pair_comp$cq), ]$id), ], select = -c(brep, trep, id))
        })
    }), recursive = FALSE)
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
ratio_by_mean_ratio <- function(df, d_cq, e_val, hkg_ = hkg) {
    target <- setdiff(names(e_val), hkg_)
    cpratio <- df[c("treatment", "gene", get("groups", qenv))]
    cpratio$cmp <- e_val[[target]]^d_cq[[target]] / e_val[[hkg_]]^d_cq[[hkg_]]
    cpratio$repxr <- as.numeric(cpratio$cmp / mean(get_control_group(cpratio)$cmp))
    return(cpratio[cpratio$gene != hkg_, ])
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
    return(lapply(data_list, function(x) {
                cmean <- control_mean(x)
                dcq <- delta_cq(x, cmean)
                ratio_by_mean_ratio(x, dcq, e_val = e_values[names(dcq)], hkg_ = hkg)
    }))
}
