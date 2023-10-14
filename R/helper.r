#' prepare input data
#' @param df input data.frame of qPCR data
#' @examples
#' prepare(df)
prepare <- function(df) {
    df <- df[c("gene", "treatment", "cq", "brep", "trep")]
    df$cq <- as.numeric(gsub(",", ".", df$cq))
    df$id <- paste0(df$treatment, df$brep, df$trep)
    return(df)
}

#' calculate E value of standard samples
#' @param efficiency_of_standard efficiency values in percent (list, 0-100)
#' @examples
#' get_e(named_list)
get_e <- function(efficiency_of_standard) {
    return((efficiency_of_standard * 0.01) + 1)
}

set_efficiency <- function(genes, efficiency, default = 100) {
    no_known_efficiency <- setdiff(genes, names(efficiency))
    default_efficiencies <- NULL
    if (!identical(no_known_efficiency, character(0))) {
        default_efficiencies <- rep(default, length(setdiff(genes, names(efficiency))))
        names(default_efficiencies) <- no_known_efficiency
    }
    return(c(efficiency, default_efficiencies))
}
