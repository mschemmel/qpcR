#' prepare input data
#' @param df input data.frame of qPCR data
#' @examples
#' prepare(df)
prepare <- function(df) {
    df$cq <- as.numeric(gsub(",", ".", df$cq))
    return(df)
}

#' calculate E value of standard samples
#' @param efficiency_of_standard efficiency values in percent (list, 0-100)
#' @examples
#' get_e(named_list)
get_e <- function(efficiency_of_standard) {
    return((efficiency_of_standard * 0.01) + 1)
}

#' set efficiency of genes
#' @param genes character vector of gene names
#' @param efficiency named list of percentual efficiency of genes
#' @param default efficiency if not specified
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

#' select control group of data
#' @param df data frame of provided genes
#' @examples
#' get_control_group(df)
get_reference_group <- function(df) {
    return(df[df$treatment == get("control_group", qenv), ])
}

#' drop columns based on character vector
#' @param df input data frame
#' @param cols character vector of columns to drop
#' @examples
#' drop_columns(df, cols)
drop_columns <- function(df, cols) {
    return(df[, -which(names(df) %in% cols)])
}