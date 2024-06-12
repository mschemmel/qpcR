#' Prepare input data
#' @param df input data.frame of qPCR data
#' @examples
#' prepare(df)
prepare <- function(df) {
    df$cq <- as.numeric(gsub(",", ".", df$cq))
    if (!("efficiency" %in% colnames(df))) df$efficiency <- 100
    df$efficiency[is.na(df$efficiency)] <- 100
    df$E <- get_E(as.numeric(gsub(",", ".", df$efficiency)))
    cols <- colnames(df)
    non_affected <- c("cq", "E", "efficiency")
    df[!(cols %in% non_affected)] <- lapply(df[!(cols %in% non_affected)], as.character)
    return(df)
}

#' Sanitize groups if NA in input data
#' @param df input data.frame containing NA
#' @examples
#' sanitize(df)
sanitize <- function(df) {
    if (!all(c("brep", "trep") %in% colnames(df))) stop("Found NA in 'cq' column. Column 'brep' or 'trep' are required, but not provided.")
    df$id <- paste0(df$treatment, df$brep, df$trep)
    df <- df[!(df$id %in% df[is.na(df$cq), ]$id), ]
    return(df)
}

#' Calculate E value of standard samples
#' @param efficiency numeric vector of efficiency values in percent (vector, 0-100)
#' @examples
#' get_E(c(87,67,98,78))
get_E <- function(efficiency) { return((efficiency * 0.01) + 1) }

#' Select control group of data
#' @param df data frame of provided genes
#' @examples
#' get_reference_group(df)
get_reference_group <- function(df) {
    return(df[df$treatment == base::get("reference_group", qenv), ])
}

#' Drop columns based on character vector
#' @param df input data frame
#' @param cols character vector of columns to drop
#' @examples
#' drop_columns(df, cols)
drop_columns <- function(df, cols = c("brep", "trep", "id")) {
    return(df[, names(df)[!(names(df) %in% cols)]])
}

#' Calculate geometric mean if multiple housekeeping genes are given
#' @param values vector of numeric values
#' @examples
#' geometric_mean(c(0.23,0.53,0.12))
geometric_mean <- function(values) return(exp(mean(log(na.omit(values)))))

#' Calculate standard error
#' @param x numeric vector of expression values
#' @examples
#' se(runif(10))
se <- function(x) return(sd(x) / sqrt(length(x)))

#' Check for outlier in expression values
#' @param x numeric vector of expression data
#' @examples
#' find_outlier(x)
interquartile_range <- function(x) return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))

#' Cut vector into half
#' @param x numeric vector
#' @param chunk_size number of elements per chunk (default: 2)
#' @examples
#' cut_in_half(c(1:10))
cut_in_half <- function(x, chunk_size = 2) {
    if (length(x) %% chunk_size != 0) stop("Chunks not equal in length.")
    return(split(x, ceiling(seq_along(x) / (length(x) / chunk_size))))
}