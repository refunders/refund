context("Package Authors in DESCRIPTION")
library(refundDevel)

test_that("Package Author Family Names are correct", {
    last.names <- unlist(eval(parse(text=packageDescription("refundDevel")[["Authors@R"]]))$family)
    contributors <- c("Huang", "Scheipl", "Goldsmith", "Gellar", "Harezlak", "McLean", "Swihart",
                      "Xiao", "Crainiceanu", "Reiss", "Chen", "Greven", "Huo", "Kundu", "Wrobel")
    expect_identical(last.names, contributors)
})
