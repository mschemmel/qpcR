#' Set efficiency of genes
#' @param genes character vector of gene names
#' @param efficiency named list of percentual efficiency of genes
#' @param default efficiency if not specified (default: 100)
#' @examples
#' set_efficiency(c("gene1", "gene2"), list("gene1" = 98, "gene2" = 95))
set_efficiency <- function(genes, efficiency, default = 100) {
    no_known_efficiency <- setdiff(genes, names(efficiency))
    default_efficiencies <- NULL
    if (!identical(no_known_efficiency, character(0))) {
        default_efficiencies <- rep(default, length(setdiff(genes, names(efficiency))))
        names(default_efficiencies) <- no_known_efficiency
    }
    return(c(efficiency, default_efficiencies))
}

#' Create chunks of groups for calculation
#' @param df data frame of expression data
#' @param groups character vector of requested groups to compare
#' @examples
#' make_groups(df, groups = c("plate", "diet", "timepoint"))
make_groups <- function(df, groups) {
    if (is.null(groups)) return(list(df))
    return(split(df, as.list(df[groups]), drop = TRUE))
}

#' Remove NA sample wise as foundation of valid mean relative expression calculation
#' @param list_of_groups output of 'make_groups' containing a list of groups to compare expression
#' @param hkg_ character string of provided housekeeping genes
#' @examples
#' sanitize(list_of_groups, khg_ = c("HKG"))
sanitize <- function(list_of_groups, hkg_) {
    unlist(lapply(list_of_groups, function(gr) {
        gr$id <- paste0(gr$treatment, gr$brep, gr$trep)
        lapply(setdiff(unique(gr$gene), hkg_), function(x) {
            pair_comp <- gr[gr$gene %in% c(x, hkg_), ]
            drop_columns(pair_comp[!(pair_comp$id %in% pair_comp[is.na(pair_comp$cq), ]$id), ], c("brep", "trep", "id"))
        })
    }), recursive = FALSE)
}

#' Get mean expression values of all control groups
#' @param df clean data frame of expression data
#' @examples
#' control_mean(df)
control_mean <- function(df) {
    contr <- get_reference_group(df)
    return(tapply(contr$cq, contr$gene, mean))
}

#' Calculate delta cq values per comparison
#' @param df clean data frame of expression data
#' @param contr_mean object of control_mean function
#' @examples
#' delta_cq(df, contr_mean)
delta_cq <- function(df, contr_mean) {
    return(setNames(lapply(seq_along(contr_mean), function(x) {
        contr_mean[x] - df[df$gene == names(contr_mean)[x], ]$cq
    }), names(contr_mean)))
}

#' Calculate the ratio compared to the mean ratio per gene
#' @param df data frame of requested groups
#' @param d_cq object of delta_cq
#' @param e_val named list of efficiency values
#' @param hkg_ character vector of housekeeping genes
#' @examples
#' ratio_by_mean_ratio(df, d_cq, e_val, hkg_, treatm)
ratio_by_mean_ratio <- function(df, d_cq, e_val, hkg) {
    target <- setdiff(names(e_val), hkg)
    cpratio <- df[c("treatment", "gene", get("groups", qenv))]
    cpratio$cmp <- e_val[[target]]^d_cq[[target]] / e_val[[hkg]]^d_cq[[hkg]]
    cpratio$rexpr <- as.numeric(cpratio$cmp / mean(get_reference_group(cpratio)$cmp))
    return(drop_columns(cpratio[cpratio$gene != hkg, ], "cmp"))
}

#' Calculate mean relative expression
#' @param df data frame of treatment and corresponding ratio values
#' @examples
#' mean_relative_expression(df)
mean_relative_expression <- function(df) {
    return(tapply(df$rbyr, df$treatment, FUN = mean))
}

#' Apply calculation for every hkg vs. target pair
#' @param data_list named list of pairs (data)
#' @param e_values E value pair of requested comparison
#' @param hkg character vector of housekeeping genes
#' @examples
#' pair_wise(data_list, e_values, hkg)
pair_wise <- function(data_list, hkg) {
    return(lapply(data_list, function(x) {
        dcq <- delta_cq(x, control_mean(x))
        e_values <- unlist(lapply(set_efficiency(unique(x$gene), lapply(split(x$efficiency, x$gene), unique)), get_e))
        ratio_by_mean_ratio(x, dcq, e_val = e_values[names(dcq)], hkg)
    }))
}
