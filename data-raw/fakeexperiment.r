rm(list = ls())

file.copy("data-raw/fakeexperiment.csv", "data/fakeexperiment.csv",
          overwrite = TRUE)