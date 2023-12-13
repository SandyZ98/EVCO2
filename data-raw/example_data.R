## code to prepare `example_data` dataset goes here

example_data <- readxl::read_excel('data-raw/Aranet4 0C1B8_2023-11-15T18_47_56-0500.xlsx')
usethis::use_data(example_data, overwrite = TRUE)
