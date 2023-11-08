#' Prepare input data
#' @param df input data.frame of qPCR data
#' @examples
#' prepare(df)
prepare <- function(df) {
    df$cq <- as.numeric(gsub(",", ".", df$cq))
    if (!("efficiency" %in% colnames(df))) df$efficiency <- 100
    df$E <- get_e(df$eff)
    return(df)
}

#' Calculate E value of standard samples
#' @param efficiency numeric vector of efficiency values in percent (vector, 0-100)
#' @examples
#' get_e(c(87,67,98,78))
get_e <- function(efficiency) {
    efficiency[is.na(efficiency)] <- 100
    return(as.numeric((efficiency * 0.01) + 1))
}

#' Select control group of data
#' @param df data frame of provided genes
#' @examples
#' get_control_group(df)
get_reference_group <- function(df) {
    return(df[df$treatment == get("reference_group", qenv), ])
}

#' Drop columns based on character vector
#' @param df input data frame
#' @param cols character vector of columns to drop
#' @examples
#' drop_columns(df, cols)
drop_columns <- function(df, cols = NULL) {
    if (is.null(cols)) return(df)
    return(df[, -which(names(df) %in% cols)])
}

#' Calculate geometric mean if multiple housekeeping genes are given
#' @param values vector of numeric values
#' @examples
#' geometric_mean(c(0.23,0.53,0.12))
geometric_mean <- function(values) {
    return(exp(mean(log(na.omit(values)))))
}
