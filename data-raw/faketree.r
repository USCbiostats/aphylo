faketree <- as.matrix(readr::read_csv2("data-raw/faketree.csv"))
devtools::use_data(faketree, overwrite = TRUE)
