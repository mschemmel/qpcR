#' Get all combinations of target and housekeeping genes as list
#' @param df data frame of every group
#' @param hkg_ character string of housekeeping gene(s)
#' @returns a list of all target x housekeeping genes combinations
#' @examples
#' comp_gene_pair(df, "Actin")
comp_gene_pair <- function(df, hkg) {
    return(lapply(setdiff(df$gene, hkg), function(x) { drop_columns(df[df$gene %in% c(x, hkg), ]) }))
}

#' Create chunks of groups for calculation
#' @param df data frame of expression data
#' @param hkg character string of provided housekeeping gene(s)
#' @param groups character vector of requested groups to compare (default: NULL)
#' @returns a list of all target x housekeeping x group combinations
#' @examples
#' make_groups(df, hkg = "HKG", groups = c("plate", "diet", "timepoint"))
make_groups <- function(df, hkg, groups = NULL) {
    if (!all(groups %in% colnames(df))) stop("Not all group(s) in data provided.")
    if (is.null(groups)) return(comp_gene_pair(df, hkg))
    pairs <- lapply(split(df, as.list(df[groups]), drop = TRUE), function(gr) {
                if (any(is.na(gr$cq))) gr <- filter_NA(gr)
                comp_gene_pair(gr, hkg)
            })
    return(unlist(pairs, recursive = FALSE))
}

#' Filter outlier of expression values, gene and treatment specific
#' @param df data frame of expression data
#' @param method method to apply to find outliers of numeric vector
#' @param do boolean if expression values should be filtered by outliers
#' @returns a data frame with removed outlier values based on a selected method ('interquartile', 'z-score', 'hampel')
#' @examples
#' filter_outlier(df, method = "interquartile", do = TRUE)
filter_outlier <- function(df, method, do) {
    if (!do) return(df)
    return(lapply(df, function(x) {
        outlier <- Reduce('+', cut_in_half(with(x, ave(cq, list(gene, treatment), FUN = function(x) outlier_method(x, method))), chunk_size = length(unique(x$gene)))) == 0
        return(x[outlier, ])
    }))
}

#' Calculate delta cq values per comparison
#' @param df clean data of given group
#' @returns a data frame with added columns 'ref_mean' and 'delta_cq' as prerequisite for relative expression calculation
#' @examples
#' delta_cq(df)
delta_cq <- function(df) {
    ref <- get_reference_group(df)
    ref_mean <- tapply(ref$cq, ref$gene, mean)
    df$ref_mean <- ref_mean[match(df$gene, names(ref_mean))]
    df$delta_cq <- df$ref_mean - df$cq
    return(df)
}

#' Calculate efficiency based ratio of delta cq value
#' @param df data frame to add ratio column to
#' @returns a data frame with added column 'ratio' containing efficiency based ratio values of delta cq values
#' @examples
#' ratio(df)
ratio <- function(df) {
    df$ratio <- df$E^df$delta_cq
    return(df)
}

#' Calculate the ratio compared to the mean ratio per gene
#' @param df data frame of requested groups
#' @param hkg character vector of housekeeping genes
#' @returns a data frame with calculted ratio of a gene compared to the mean ratio
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
#' @returns a lists of calculated ratios (data frame) for every housekeeping vs. target gene pair
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
#' @returns the final output of calculated relative expression values
#' @examples
#' conflate(df)
conflate <- function(df, do) {
    if (!do) return(df[order(df$treatment), ])
    formula_string <- paste0("rexpr ~ ", paste(c("treatment", "gene", get("groups", qenv)), collapse = "+"))
    stats <- function(x) c("mean" = mean(x), "sd" = sd(x), "se" = se(x), "n" = length(x))
    aggregate_out <- do.call(data.frame, aggregate(as.formula(formula_string), df, stats))
    return(aggregate_out[order(aggregate_out$treatment), ])
}