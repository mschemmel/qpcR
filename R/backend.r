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
#' sanitize(list_of_groups, khg = c("HKG"))
sanitize <- function(list_of_groups, hkg) {
    unlist(lapply(list_of_groups, function(gr) {
        gr$id <- paste0(gr$treatment, gr$brep, gr$trep) #TODO: Allow user selection of unique ID
        lapply(setdiff(gr$gene, hkg), function(x) {
            pair_comp <- gr[gr$gene %in% c(x, hkg), ]
            drop_columns(pair_comp[!(pair_comp$id %in% pair_comp[is.na(pair_comp$cq), ]$id), ], c("brep", "trep", "id"))
        })
    }), recursive = FALSE)
}

#' Calculate delta cq values per comparison
#' @param df clean data frame of expression data
#' @param contr_mean object of control_mean function
#' @examples
#' delta_cq(df, contr_mean)
delta_cq <- function(df) {
    ref <- get_reference_group(df)
    ref_mean <- tapply(ref$cq, ref$gene, mean)
    df$ref_mean <- ref_mean[match(df$gene, names(ref_mean))]
    df$delta_cq <- df$ref_mean - df$cq
    return(df)
}

#' Calculate the ratio compared to the mean ratio per gene
#' @param df data frame of requested groups
#' @param hkg character vector of housekeeping genes
#' @examples
#' ratio_by_mean_ratio(df, hkg)
ratio_by_mean_ratio <- function(df, hkg) {
    df <- delta_cq(df)
    target_df <- df[df$gene == setdiff(df$gene, hkg), ]
    target_df$ratio <- target_df$E^target_df$delta_cq
    hkg_df <- df[df$gene %in% hkg, ]
    hkg_df$ratio <- hkg_df$E^hkg_df$delta_cq
    ratio <- apply(data.frame(do.call(cbind, lapply(split(hkg_df, hkg_df$gene), function(x) x$ratio))), 1, geometric_mean)
    target_df$cmp_ratio <- target_df$ratio / ratio
    target_df$rexpr <- as.numeric(target_df$cmp_ratio / mean(get_reference_group(target_df)$cmp_ratio))
    return(drop_columns(target_df, c("control_mean", "delta_cq", "cmp_ratio")))
}

#' Apply calculation for every HKG vs. target pair
#' @param pair named list of pairs (data)
#' @param hkg character vector of housekeeping genes
#' @examples
#' pair_wise(pair, hkg)
pair_wise <- function(pair, hkg) {
    return(lapply(pair, function(x) {
        ratio_by_mean_ratio(x, hkg)
    }))
}
