fakeexperiment <- as.matrix(readr::read_csv2("data-raw/fakeexperiment.csv"))
devtools::use_data(fakeexperiment, overwrite = TRUE)

