#' Sanitize groups if NA in input data
#' @param df input data.frame containing NA
#' @return data frame without NA values in cq column by removing sample specific values per group to ensure balanced number of observations
#' @keywords internal
filter_NA <- function(df) {
    assert(all(c("brep", "trep") %in% colnames(df)), "Found NA in 'cq' column. Column 'brep' or 'trep' are required, but not provided.")
    df$id <- paste0(df$treatment, df$brep, df$trep)
    df <- df[!(df$id %in% df[is.na(df$cq), ]$id), ]
    return(df)
}

#' Calculate E value of standard samples
#' @param efficiency numeric vector of efficiency values in percent (vector, 0-100)
#' @return numeric vector of E values based on provided efficiency
#' @keywords internal
get_E <- function(efficiency) { return((efficiency * 0.01) + 1) }

#' Select control group of data
#' @param df data frame of provided genes
#' @return data frame of selected reference group (treatment == reference)
#' @keywords internal
get_reference_group <- function(df) {
    return(df[df$treatment == base::get("reference_group", qenv), ])
}

#' Drop columns based on character vector
#' @param df input data frame
#' @param cols character vector of columns to drop
#' @return data frame with omitted columns
#' @keywords internal
drop_columns <- function(df, cols = c("brep", "trep", "id")) {
    return(df[, names(df)[!(names(df) %in% cols)]])
}

#' Calculate geometric mean if multiple housekeeping genes are given
#' @param values vector of numeric values
#' @return geometric mean of numeric input vector
#' @keywords internal
#' @importFrom stats na.omit
geometric_mean <- function(values) return(exp(mean(log(na.omit(values)))))

#' Calculate standard error
#' @param x numeric vector of expression values
#' @return standard error of numeric input vector
#' @keywords internal
se <- function(x) return(sd(x) / sqrt(length(x)))

#' Check for outliers in expression values
#' @param x numeric vector of expression data (cq column)
#' @param method method to apply to find outliers of numeric vector
#' @importFrom stats quantile IQR median mad
#' @return boolean vector of numeric input data if outlier were detected by choosed method
#' @keywords internal
outlier_method <- function(x, method) {
    assert(is.numeric(x), "Input vector for outlier detection is not numeric.")
    assert(method %in% c("interquartile", "z-score", "hampel"), "Outlier detection method not supported. Please choose 'interquartile', 'z-score' or 'hampel'.")
    switch(method,
           "interquartile" = x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x),
           "z-score" = abs((x - mean(x)) / sd(x)) > 3,
           "hampel" = x < (median(x) - 3 * mad(x, constant = 1)) | x > (median(x) + 3 * mad(x, constant = 1)))
}

#' Cut vector into half
#' @param x numeric vector
#' @param chunk_size number of elements per chunk (default: 2)
#' @return vector half in length of input vector
#' @keywords internal
cut_in_half <- function(x, chunk_size = 2) {
    assert(length(x) %% chunk_size == 0, "Chunks not equal in length.")
    return(split(x, ceiling(seq_along(x) / (length(x) / chunk_size))))
}

#' Simple function to assert given conditions
#' @param statement vector of comparison to proof
#' @param err_message string of error message when not all comparison are TRUE
#' @return character string used as error message or nothing if all TRUE
#' @keywords internal
assert <- function(statement, err_message = NULL) {
    if (!all(statement)) stop(err_message)
}