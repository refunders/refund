context("Test mfpca")

test_that("all mfpca functions work on the DTI example", {
  skip_on_cran()

  data(DTI)
  DTI <- subset(DTI, Nscans < 6)  ## example where all subjects have 6 or fewer visits
  id <- DTI$ID
  Y <- DTI$cca

  mfpca.sc.DTI <- mfpca.sc(Y = Y, id = id, twoway = TRUE)
  mfpca.face.DTI <- mfpca.face(Y = Y, id = id, twoway = TRUE)

  expect_equal(dim(mfpca.face.DTI$Yhat)[1], dim(mfpca.sc.DTI$Yhat)[1])
  expect_equal(dim(mfpca.face.DTI$Yhat.subject)[1], dim(mfpca.sc.DTI$Yhat.subject)[1])
  expect_equal(dim(mfpca.face.DTI$scores$level1)[1], dim(mfpca.sc.DTI$scores$level1)[1])
  expect_equal(dim(mfpca.face.DTI$scores$level2)[1], dim(mfpca.sc.DTI$scores$level2)[1])
})

test_that("mfpca.face options work", {
  skip_on_cran()

  data(DTI)
  DTI <- subset(DTI, Nscans < 6)  ## example where all subjects have 6 or fewer visits
  id <- DTI$ID
  Y <- DTI$cca
  mfpca.base <- mfpca.face(Y = Y, id = id)

  # visit argument
  mfpca.visit <- mfpca.face(Y = Y, id = id, visit = DTI$visit)
  expect_equal(mfpca.base$npc$level1, mfpca.visit$npc$level1)

  # weight argument
  mfpca.weight <- mfpca.face(Y = Y, id = id, weight = "subj")
  expect_equal(dim(mfpca.base$scores$level1)[1], dim(mfpca.weight$scores$level1)[1])

  # pve argument
  mfpca.pve <- mfpca.face(Y = Y, id = id, pve = 0.95)
  expect_equal(dim(mfpca.base$scores$level1)[1], dim(mfpca.pve$scores$level1)[1])

  # npc argument
  mfpca.npc <- mfpca.face(Y = Y, id = id, npc = 5)
  expect_equal(mfpca.npc$npc$level1, 5)
  expect_equal(mfpca.npc$npc$level2, 5)
})


