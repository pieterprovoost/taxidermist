library(testthat)
library(taxidermist)

test_that("parse_taxonomy_string parses strings correctly", {
  input <- c(
    ">AF290071;tax=d:Eukaryota,p:NA,c:NA,o:NA,f:NA,g:NA,s:Uncultured_marine_alveolate_Group_II_DH147-EKD16",
    ">GQ892044;tax=sk:Eukaryota,p:Endomyxa,c:Phytomyxea,o:Plasmodiophorida,f:Plasmodiophoridae,g:Polymyxa,s:Polymyxa_graminis",
    ">PV277008;tax=d:Eukaryota,p:Annelida,c:Polychaeta,o:NA,f:Capitellidae,g:Notodasus,s:NA",
    ">XR_012195772;tax=d:NA,p:NA,c:NA,o:NA,f:NA,g:NA,s:NA"
  )
  df <- parse_taxonomy_string(input)
  expect_equal(nrow(df), 4)
  expect_true(all(c("superkingdom","domain","phylum","class","order","family","genus","species") %in% colnames(df)))
  expect_equal(df$domain[1], "Eukaryota")
  expect_equal(df$phylum[2], "Endomyxa")
  expect_equal(df$scientificName[3], "Notodasus")
  expect_equal(df$taxonRank[3], "genus")
  expect_true(is.na(df$phylum[1]))
  expect_true(is.na(df$superkingdom[1]))
})

test_that("parse_taxonomy_string converts 'NA' strings to NA", {
  input <- ">X;tax=d:NA,p:NA,c:NA,o:NA,f:NA,g:NA,s:NA"
  df <- parse_taxonomy_string(input)
  expect_true(all(is.na(df)))
})

test_that("parse_taxonomy_string errors on invalid input", {
  expect_error(parse_taxonomy_string(">X;d:Eukaryota,p:NA"))
  expect_error(parse_taxonomy_string(">X;tax=d:Eukaryota,pNA"))
  expect_error(parse_taxonomy_string(">X;tax=d:Eukaryota,p:NA,c:NA,"))
})

test_that("parse_taxonomy adds taxonomy columns correctly", {
  input <- data.frame(taxonomy = c(
    ">AF290071;tax=d:Eukaryota,p:NA,c:NA,o:NA,f:NA,g:NA,s:Uncultured_marine_alveolate_Group_II_DH147-EKD16",
    ">PV277008;tax=d:Eukaryota,p:Annelida,c:Polychaeta,o:NA,f:Capitellidae,g:Notodasus,s:NA",
    ">XR_012195772;tax=d:NA,p:NA,c:NA,o:NA,f:NA,g:NA,s:NA"
  ))
  df <- input %>% 
    parse_taxonomy(taxonomy)
  expect_equal(df$scientificName[1], "Uncultured_marine_alveolate_Group_II_DH147-EKD16")
  expect_equal(df$scientificName[2], "Notodasus")
  expect_true(is.na(df$scientificName[3]))
  expect_equal(df$domain[1], "Eukaryota")
  expect_equal(df$domain[2], "Eukaryota")
  expect_true(is.na(df$domain[3]))
})

test_that("remove_unparsable_names functions as expected", {
  input <- data.frame(taxonomy = c(
    ">AF290071;tax=d:Eukaryota,p:NA,c:NA,o:NA,f:NA,g:NA,s:Uncultured_marine_alveolate_Group_II_DH147-EKD16",
    ">PV277008;tax=d:Eukaryota,p:Annelida,c:Polychaeta,o:NA,f:Capitellidae,g:Notodasus,s:NA",
    ">XR_012195772;tax=d:NA,p:NA,c:NA,o:NA,f:NA,g:NA,s:NA"
  ))
  df <- parse_taxonomy(input, taxonomy) %>%
    remove_unparsable_names()
  expect_equal(df$scientificName[1], "Eukaryota")
  expect_equal(df$taxonRank[1], "domain")
  expect_equal(df$scientificName[2], "Notodasus")
  expect_equal(df$taxonRank[2], "genus")
  expect_true(is.na(df$scientificName[3]))
  expect_true(is.na(df$taxonRank[3]))
})

test_that("populate_species functions as expected", {
  input <- data.frame(
    scientificName = c("Abra alba", "Mola"),
    genus = c("Abra", "Mola"),
    taxonRank = c("species", "genus")
  )
  df <- populate_species(input)
  expect_equal(df$species[1], "Abra alba")
  expect_true(is.na(df$species[2]))
})
