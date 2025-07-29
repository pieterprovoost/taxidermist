#' taxidermist
#'
#' @docType package
#' @name taxidermist
"_PACKAGE"

#' @import dplyr stringr purrr
NULL

#' Replace string 'NA' with actual NA in all character columns
#' 
#' @param df A data frame.
#' @return A data frame with 'NA' strings in character columns replaced by NA values.
replace_na_strings <- function(df) {
  df %>%
    dplyr::mutate(
      dplyr::across(
        where(is.character),
        ~dplyr::na_if(., "NA")
      )
    )
}

#' Parse taxonomy strings into a data frame
#'
#' @param input Character vector. Each element must contain 'tax=' followed by a comma-separated list of key:value pairs (e.g., d:Eukaryota,p:NA,...).
#' @return A data frame with columns for each taxonomic rank and one row per input string.
#' @export
parse_taxonomy_string <- function(input) {
  ranks_map = c("superkingdom" = "sk", "domain" = "d", "phylum" = "p", "class" = "c", "order" = "o", "family" = "f", "genus" = "g", "species" = "s")

  valid_pattern <- "^.*tax=([a-z]+:[^,]+)(,[a-z]+:[^,]+)*$"
  if (any(!grepl(valid_pattern, input))) {
    stop("Each input string must contain 'tax=' followed by a comma-separated list of key:value pairs (e.g., d:Eukaryota,p:NA,...)")
  }

  input %>%
    str_extract("(?<=tax=)[^\\s]+") %>% 
    str_split(",") %>% 
    map(str_split, ":") %>% 
    map(function(lst) {
      setNames(sapply(lst, `[`, 2), sapply(lst, `[`, 1))
    }) %>% 
    map(function(lst) {
      lst[setdiff(unname(ranks_map), names(lst))] <- NA
      lst[unname(ranks_map)]
    }) %>%
    bind_rows() %>% 
    rename(any_of(ranks_map)) %>%
    replace_na_strings() %>%
    mutate(scientificName = do.call(coalesce, across(all_of(rev(names(ranks_map)))))) %>%
    rowwise() %>%
    mutate(
      taxonRank = {
        vals <- c_across(all_of(names(ranks_map)))
        non_na <- which(!is.na(vals))
        if (length(non_na)) names(ranks_map)[max(non_na)] else NA_character_
      }
    ) %>%
    ungroup()
}

#' Add taxonomy to a data frame by parsing a taxonomy string column
#' 
#' @param df A data frame.
#' @param col The column containing taxonomy strings formatted as 'tax=' followed by a comma-separated list of key:value pairs (e.g., d:Eukaryota,p:NA,...).
#' @return A data frame.
#' @export
parse_taxonomy <- function(df, col) {
  col <- rlang::enquo(col)
  tax_df <- parse_taxonomy_string(pull(df, !!col))
  dplyr::bind_cols(df, tax_df)
}
