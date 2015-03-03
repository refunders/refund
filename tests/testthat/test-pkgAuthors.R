context("Package Authors in DESCRIPTION")
library(refundDevel)

test_that("Package Author Family Names are correct", {
    last.names <- unlist(eval(parse(text=packageDescription("refundDevel")[["Authors@R"]]))$family)
    contributors <- c("Crainiceanu", "Reiss", "Goldsmith", "Huang", "Huo", "Scheipl",
                      "Swihart", "Greven", "Harezlak", "Kundu", "Zhao", "McLean", "Xiao",
                      "Gellar")
    expect_identical(last.names, contributors)
})
