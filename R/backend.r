#' Create chunks of groups for calculation
#' @param df data frame of expression data
#' @param hkg_ character string of provided housekeeping genes
#' @param groups character vector of requested groups to compare (default: NULL)
#' @examples
#' make_groups(df, groups = c("plate", "diet", "timepoint"))
make_groups <- function(df, hkg, groups = NULL) {
    cols <- colnames(df)
    non_affected <- c("cq", "E", "efficiency")
    df[!(cols %in% non_affected)] <- lapply(df[!(cols %in% non_affected)], as.character)
    if (is.null(groups)) {
        return(lapply(setdiff(df$gene, hkg), function(x) { df[df$gene %in% c(x, hkg), ] }))
    }
    if (!all(groups %in% colnames(df))) {
        stop("Not all group(s) in data provided.")
    }
    return(split(df, as.list(df[groups]), drop = TRUE))
}

#' Remove NA sample wise as foundation of valid mean relative expression calculation
#' @param list_of_groups output of 'make_groups' containing a list of groups to compare expression
#' @param hkg_ character string of provided housekeeping genes
#' @examples
#' sanitize(list_of_groups, hkg = c("HKG"))
sanitize <- function(list_of_groups, hkg) {
    unlist(lapply(list_of_groups, function(gr) {
        if (any(is.na(gr$cq))) {
            gr$id <- paste0(gr$treatment, gr$brep, gr$trep) #TODO: Allow user selection of unique ID
            lapply(setdiff(gr$gene, hkg), function(x) {
                pair_comp <- gr[gr$gene %in% c(x, hkg), ]
                drop_columns(pair_comp[!(pair_comp$id %in% pair_comp[is.na(pair_comp$cq), ]$id), ], c("brep", "trep", "id"))
            })
        } else {
            list(gr)
        }
    }), recursive = FALSE)
}

#' Calculate delta cq values per comparison
#' @param df clean data of given group
#' @examples
#' delta_cq(df, contr_mean)
delta_cq <- function(df) {
    ref <- get_reference_group(df)
    ref_mean <- tapply(ref$cq, ref$gene, mean)
    df$ref_mean <- ref_mean[match(df$gene, names(ref_mean))]
    df$delta_cq <- df$ref_mean - df$cq
    return(df)
}

#' Calculate efficiency based ratio of delta cq value
#' @param df data frame to add ratio column to
#' @examples
#' ratio(df)
ratio <- function(df) {
    df$ratio <- df$E^df$delta_cq
    return(df)
}

#' Calculate the ratio compared to the mean ratio per gene
#' @param df data frame of requested groups
#' @param hkg character vector of housekeeping genes
#' @examples
#' ratio_by_mean_ratio(df, hkg)
ratio_by_mean_ratio <- function(df, hkg) {
    df <- delta_cq(df)
    target_df <- ratio(df[!(df$gene %in% hkg), ])
    hkg_df <- ratio(df[df$gene %in% hkg, ])
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

#' Summarize expression data
#' @param df data frame of relative expression data
#' @param do boolean if to aggregate the input data (default: TRUE)
#' @examples
#' conflate(df)
conflate <- function(df, do = TRUE) {
    if (!do) return(df)
    formula_string <- paste0("rexpr ~ ", paste(c("treatment", "gene", get("groups", qenv)), collapse = "+"))
    stats <- function(x) c("mean" = mean(x), "sd" = sd(x), "n" = length(x))
    return(do.call(data.frame, aggregate(as.formula(formula_string), df, stats)))
}

