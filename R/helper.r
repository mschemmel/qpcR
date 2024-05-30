#' Prepare input data
#' @param df input data.frame of qPCR data
#' @examples
#' prepare(df)
prepare <- function(df) {
    df$cq <- as.numeric(gsub(",", ".", df$cq))
    if (!("efficiency" %in% colnames(df))) df$efficiency <- 100
    df$E <- get_e(as.numeric(gsub(",", ".", df$efficiency)))
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
     df$id <- paste0(df$treatment, df$brep, df$trep) #TODO: Allow user selection of unique ID
     df <- df[!(df$id %in% df[is.na(df$cq), ]$id), ]
     return(df)
}


#' Calculate E value of standard samples
#' @param efficiency numeric vector of efficiency values in percent (vector, 0-100)
#' @examples
#' get_e(c(87,67,98,78))
get_e <- function(efficiency) {
    efficiency[is.na(efficiency)] <- 100
    return((efficiency * 0.01) + 1)
}

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
geometric_mean <- function(values) {
    return(exp(mean(log(na.omit(values)))))
}
