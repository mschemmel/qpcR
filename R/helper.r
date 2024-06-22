#' Prepare input data
#' @param df input data.frame of qPCR data
#' @returns clean data frame of input data (type conversion, adding default values, calculate E values of efficiency)
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
#' @returns data frame without NA values in cq column by removing sample specific values per group to ensure balanced number of observations
#' @examples
#' sanitize(df)
filter_NA <- function(df) {
    if (!all(c("brep", "trep") %in% colnames(df))) stop("Found NA in 'cq' column. Column 'brep' or 'trep' are required, but not provided.")
    df$id <- paste0(df$treatment, df$brep, df$trep)
    df <- df[!(df$id %in% df[is.na(df$cq), ]$id), ]
    return(df)
}

#' Calculate E value of standard samples
#' @param efficiency numeric vector of efficiency values in percent (vector, 0-100)
#' @returns numeric vector of E values based on provided efficiency
#' @examples
#' get_E(c(87,67,98,78))
get_E <- function(efficiency) { return((efficiency * 0.01) + 1) }

#' Select control group of data
#' @param df data frame of provided genes
#' @returns data frame of selected reference group (treatment == reference)
#' @examples
#' get_reference_group(df)
get_reference_group <- function(df) {
    return(df[df$treatment == base::get("reference_group", qenv), ])
}

#' Drop columns based on character vector
#' @param df input data frame
#' @param cols character vector of columns to drop
#' @returns data frame with omitted columns
#' @examples
#' drop_columns(df, cols)
drop_columns <- function(df, cols = c("brep", "trep", "id")) {
    return(df[, names(df)[!(names(df) %in% cols)]])
}

#' Calculate geometric mean if multiple housekeeping genes are given
#' @param values vector of numeric values
#' @returns geometric mean of numeric input vector
#' @examples
#' geometric_mean(c(0.23,0.53,0.12))
geometric_mean <- function(values) return(exp(mean(log(na.omit(values)))))

#' Calculate standard error
#' @param x numeric vector of expression values
#' @returns standard error of numeric input vector
#' @examples
#' se(runif(10))
se <- function(x) return(sd(x) / sqrt(length(x)))

#' Check for outliers in expression values
#' @param x numeric vector of expression data (cq column)
#' @param method method to apply to find outliers of numeric vector
#' @returns boolean vector of numeric input data if outlier were detected by choosed method
#' @examples
#' outlier_method(x)
outlier_method <- function(x, method) {
    if (!is.numeric(x)) stop("Input vector for outlier detection is not numeric.")
    if (!(method %in% c("interquartile", "z-score", "hampel"))) stop("Outlier detection method not supported. Please choose 'interquartile', 'z-score' or 'hampel'.")
    switch(method,
           "interquartile" = x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x),
           "z-score" = abs((x - mean(x)) / sd(x)) > 3,
           "hampel" = x < (median(x) - 3 * mad(x, constant = 1)) | x > (median(x) + 3 * mad(x, constant = 1)))
}

#' Cut vector into half
#' @param x numeric vector
#' @param chunk_size number of elements per chunk (default: 2)
#' @returns vector half in length of input vector
#' @examples
#' cut_in_half(c(1:10))
cut_in_half <- function(x, chunk_size = 2) {
    if (length(x) %% chunk_size != 0) stop("Chunks not equal in length.")
    return(split(x, ceiling(seq_along(x) / (length(x) / chunk_size))))
}

#' Detect variable of input data used as 'reference'
#' @param treatment vector of unique character strings found in input data ('treatment' column)
#' @param reference vector of single reference value
#' @returns character string used as 'reference' group for further calculations
#' @examples
#' detect_reference(c("control", "infected", "inoculated"), "control")
detect_reference <- function(treatment, reference) {
    if (is.null(reference)) {
        ref_defaults <- c("Control", "control", "Mock", "mock", "Kontrolle", "kontrolle", "C", "c", "K", "k")
        ref_ <- treatment[treatment %in% ref_defaults]
        if (identical(unique(ref_), character(0))) stop("No 'reference' variable provided and non auto-detected.")
        if (length(ref_) > 1) stop("Found more than one 'reference' variable in 'treatment' column.")
        return(ref_)
    }
    if (!(any(reference %in% unique(treatment)))) stop("Variable provided as 'reference' not present in input data.")
    return(reference)
}