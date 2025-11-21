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
    stop(glue("Each input string must contain 'tax=' followed by a comma-separated list of key:value pairs. Failed to parse {input}"))
  }
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
  parsed <- gnparser::parse(input, cores = 2)
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

#' Match scientificName with WoRMS and populate scientificName and scientificNameID.
#' Exact matches only, in case of multiple exact matches, the first one is used.
#' 
#' @param df A data frame.
#' @return A data frame.
#' @export
match_exact_worms <- function(df, accepted=TRUE, replace_taxonomy=TRUE, taxamatch=TRUE) {
  if (taxamatch) {
    fun <- worrms::wm_records_taxamatch
    batch_size <- 50
  } else {
    fun <- worrms::wm_records_names
    batch_size <- 50
  }
  unique_names <- na.omit(unique(df$scientificName))
  chunks <- split(unique_names, ceiling(seq_along(unique_names) / batch_size))
  valid_matches <- purrr::map(chunks, function(chunk) {
    worms_results <- tryCatch({
      fun(chunk, marine_only = FALSE)
    }, error = function(e) {
      message("WoRMS request failed, skipping batch - ", e)
      rep(list(data.frame()), length(chunk))
    })
    names(worms_results) <- chunk
    valid_ids <- purrr::imap(worms_results, function(result_set, input) {
      if (nrow(result_set) == 0) {
        data.frame(input = input, valid_AphiaID = NA_integer_)
      } else {
        result_set %>%
          filter(match_type == "exact") %>% 
          mutate(sort_key = status != "accepted") %>% 
          arrange(sort_key) %>% 
          select(-sort_key) %>% 
          slice(1) %>% 
          mutate(input = input) %>%
          mutate(across(!c(AphiaID, valid_AphiaID), as.character))
      }
    }) %>% 
      bind_rows() %>% 
      select(input, valid_AphiaID) %>% 
      filter(!is.na(valid_AphiaID)) 
    if (nrow(valid_ids) > 0) {
      worrms::wm_record(valid_ids$valid_AphiaID) %>% 
        mutate(input = valid_ids$input)
    } else {
      data.frame(input = character(0), scientificname = character(0), lsid = character(0), phylum = character(0), class = character(0), order = character(0), family = character(0), genus = character(0), rank = character(0))
    }
  }, .progress = TRUE) %>% 
    bind_rows()
  if (replace_taxonomy) {
    valid_matches <- valid_matches %>%
      select(input, scientificName = scientificname, scientificNameID = lsid, phylum, class, order, family, genus, taxonRank = rank) %>% 
      mutate(taxonRank = tolower(taxonRank))
  } else {
    valid_matches <- valid_matches %>%
      select(input, scientificName = scientificname, scientificNameID = lsid)
  }
  df %>% 
    left_join(valid_matches, by = c("scientificName" = "input"), suffix = c("_old", "")) %>% 
    select(-ends_with("_old"))
}

#' Populate taxonomy fields using a numeric AphiaID column.
#' 
#' @param df A data frame.
#' @return A data frame.
#' @export
populate_worms <- function(df, col="AphiaID", replace_taxonomy=TRUE) {
  batch_size <- 50
  ids <- df[,col]
  chunks <- split(ids, ceiling(seq_along(ids) / batch_size))
  valid_matches <- purrr::map(chunks, function(chunk) {
    worms_results <- tryCatch({
      worrms::wm_record(chunk, marine_only = FALSE)
    }, error = function(e) {
      message("WoRMS request failed, skipping batch - ", e)
      rep(list(data.frame()), length(chunk))
    }) %>% split(seq(nrow(.)))
    names(worms_results) <- chunk
    valid_ids <- purrr::imap(worms_results, function(result_set, input) {
      if (nrow(result_set) == 0) {
        data.frame(input = input, valid_AphiaID = NA_integer_)
      } else {
        result_set %>%
          filter(match_type == "exact") %>% 
          mutate(sort_key = status != "accepted") %>% 
          arrange(sort_key) %>% 
          select(-sort_key) %>% 
          slice(1) %>% 
          mutate(input = input) %>%
          mutate(across(!c(AphiaID, valid_AphiaID), as.character))
      }
    }) %>% 
      bind_rows() %>% 
      select(input, valid_AphiaID) %>% 
      filter(!is.na(valid_AphiaID)) 
    if (nrow(valid_ids) > 0) {
      worrms::wm_record(valid_ids$valid_AphiaID) %>% 
        mutate(input = valid_ids$input)
    } else {
      data.frame(input = character(0), scientificname = character(0), lsid = character(0), phylum = character(0), class = character(0), order = character(0), family = character(0), genus = character(0), rank = character(0))
    }
  }, .progress = TRUE) %>% 
    bind_rows()
  if (replace_taxonomy) {
    valid_matches <- valid_matches %>%
      select(input, scientificName = scientificname, scientificNameID = lsid, phylum, class, order, family, genus, taxonRank = rank) %>% 
      mutate(taxonRank = tolower(taxonRank))
  } else {
    valid_matches <- valid_matches %>%
      select(input, scientificName = scientificname, scientificNameID = lsid)
  }
  df %>% 
    left_join(valid_matches %>% mutate(input = as.numeric(input)), by = c("AphiaID" = "input"), suffix = c("_old", "")) %>% 
    select(-ends_with("_old"))
}

# replace_accepted <- function(df, replace_taxonomy=TRUE) {
#   unique_names <- na.omit(unique(df$scientificName))
#   chunks <- split(unique_names, ceiling(seq_along(unique_names) / batch_size))
#   valid_matches <- purrr::map(chunks, function(chunk) {
#     worms_results <- tryCatch({
#       fun(chunk, marine_only = FALSE)
#     }, error = function(e) {
#       message("WoRMS request failed, skipping batch - ", e)
#       rep(list(data.frame()), length(chunk))
#     })
#     names(worms_results) <- chunk
#     valid_ids <- purrr::imap(worms_results, function(result_set, input) {
#       if (nrow(result_set) == 0) {
#         data.frame(input = input, valid_AphiaID = NA_integer_)
#       } else {
#         result_set %>%
#           filter(match_type == "exact") %>% 
#           mutate(sort_key = status != "accepted") %>% 
#           arrange(sort_key) %>% 
#           select(-sort_key) %>% 
#           slice(1) %>% 
#           mutate(input = input) %>%
#           mutate(across(!c(AphiaID, valid_AphiaID), as.character))
#       }
#     }) %>% 
#       bind_rows() %>% 
#       select(input, valid_AphiaID) %>% 
#       filter(!is.na(valid_AphiaID)) 
#     if (nrow(valid_ids) > 0) {
#       worrms::wm_record(valid_ids$valid_AphiaID) %>% 
#         mutate(input = valid_ids$input)
#     } else {
#       data.frame(input = character(0), scientificname = character(0), lsid = character(0), phylum = character(0), class = character(0), order = character(0), family = character(0), genus = character(0), rank = character(0))
#     }
#   }, .progress = TRUE) %>% 
#     bind_rows()
#   if (replace_taxonomy) {
#     valid_matches <- valid_matches %>%
#       select(input, scientificName = scientificname, scientificNameID = lsid, phylum, class, order, family, genus, taxonRank = rank) %>% 
#       mutate(taxonRank = tolower(taxonRank))
#   } else {
#     valid_matches <- valid_matches %>%
#       select(input, scientificName = scientificname, scientificNameID = lsid)
#   }
#   df %>% 
#     left_join(valid_matches, by = c("scientificName" = "input"), suffix = c("_old", "")) %>% 
#     select(-ends_with("_old"))
# }
