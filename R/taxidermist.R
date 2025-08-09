#' taxidermist
#'
#' @docType package
#' @name taxidermist
"_PACKAGE"

#' @import dplyr tidyr stringr purrr
NULL

ranks_map <- c("superkingdom" = "sk", "domain" = "d", "phylum" = "p", "class" = "c", "order" = "o", "family" = "f", "genus" = "g", "species" = "s")

#' Replace string 'NA' with actual NA in all character columns
#' 
#' @param df A data frame.
#' @return A data frame with 'NA' strings in character columns replaced by NA values.
#' @export
replace_na_strings <- function(df) {
  df %>%
    dplyr::mutate(
      dplyr::across(
        where(is.character),
        ~dplyr::na_if(., "NA")
      )
    )
}

#' Update scientificName and taxonRank based on rank columns
#' 
#' @param df A data frame.
#' @return A data frame.
#' @export
update_scientificname_rank <- function(df) {
  tax_columns <- intersect(names(ranks_map), names(df))
  df %>%
    mutate(scientificName = do.call(coalesce, across(all_of(rev(tax_columns))))) %>%
    rowwise() %>%
    mutate(
      taxonRank = {
        vals <- c_across(all_of(tax_columns))
        non_na <- which(!is.na(vals))
        if (length(non_na)) tax_columns[max(non_na)] else NA_character_
      }
    ) %>%
    ungroup()
}

#' Parse taxonomy strings into a data frame
#'
#' @param input Character vector. Each element must contain 'tax=' followed by a comma-separated list of key:value pairs (e.g., d:Eukaryota,p:NA,...).
#' @return A data frame with columns for each taxonomic rank and one row per input string.
#' @export
parse_taxonomy_string <- function(input) {
  valid_pattern <- "^.*;tax=([a-z]{1,2}:.*)(,[a-z]{1,2}:.*)+$"
  if (any(!grepl(valid_pattern, input))) {
    stop("Each input string must contain 'tax=' followed by a comma-separated list of key:value pairs (e.g., d:Eukaryota,p:NA,...)")
  }

#str_split(">HQ616763;tax=d:Eukaryota,p:Mollusca,c:Gastropoda,o:Nudibranchia,f:Aeolidiidae,g:Spurilla,s:Spurilla_neapolitana_(Delle_Chiaje,_1841)", "(?=[sk|d|p|c|o|f|g|s]:)")
# input <- ">HQ616763;tax=d:Eukaryota,p:Mollusca,c:Gastropoda,o:Nudibranchia,f:Aeolidiidae,g:Spurilla,s:Spurilla_neapolitana_(Delle_Chiaje,_1841)"
# input %>%
#   str_extract("(?<=tax=)[^\\s]+") %>% 
#   str_split(",(?=[sk|d|p|c|o|f|g|s]:)") 



  input %>%
    str_extract("(?<=tax=)[^\\s]+") %>% 
    str_split(",(?=[sk|d|p|c|o|f|g|s]:)") %>% 
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
    update_scientificname_rank()
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

#' Get distinct names from all rank columns.
#' 
#' @param df A data frame.
#' @return A data frame.
#' @export
get_distinct_names <- function(df) {
  tax_columns <- intersect(names(ranks_map), names(df))
  df %>%
    select(all_of(tax_columns)) %>%
    pivot_longer(everything(), values_to = "taxon") %>%
    filter(!is.na(taxon)) %>%
    distinct(taxon) %>%
    pull(taxon)
}

#' Parse names using gnparser.
#' 
#' @param input A vector of names.
#' @return A data frame.
gn_parse_names <- function(input) {
  parsed <- gnparser::parse(input)
  parsed$input <- input
  parsed
}

#' Parse names in taxonomy columns using gnparser, remove unparsable names, and update scientificName and taxonRank.
#' 
#' @param df A data frame.
#' @return A data frame.
#' @export
remove_unparsable_names <- function(df) {
  tax_columns <- intersect(names(ranks_map), names(df))

  unparsable_names <- df %>% 
    get_distinct_names() %>% 
    gn_parse_names() %>% 
    filter(parsed == FALSE | quality > 2) %>% 
    pull(input)

  df %>%
    mutate(across(all_of(tax_columns), ~ ifelse(. %in% unparsable_names, NA, .))) %>% 
    update_scientificname_rank()
}

#' Populate the species column using scientificName and taxonRank.
#' 
#' @param df A data frame.
#' @return A data frame.
#' @export
populate_species <- function(df) {
  df %>%
    mutate(species = if_else(taxonRank == "species", scientificName, NA_character_))
}
