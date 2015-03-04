library("testthat")
library("refundDevel")

Sys.setenv(NOT_CRAN = "true")
test_check("refundDevel")

