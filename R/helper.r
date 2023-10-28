#' Prepare input data
#' @param df input data.frame of qPCR data
#' @examples
#' prepare(df)
prepare <- function(df) {
    df$cq <- as.numeric(gsub(",", ".", df$cq))
    if (!("efficiency" %in% colnames(df))) df$efficiency <- 100
    return(df)
}

#' Calculate E value of standard samples
#' @param efficiency_of_standard efficiency values in percent (list, 0-100)
#' @examples
#' get_e(named_list)
get_e <- function(efficiency_of_standard) {
    return((efficiency_of_standard * 0.01) + 1)
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