context("Test hrcompet estimation and data return")

test_that("Wheter return gives desired data dimensions", {

  set.seed(1)
  data <- hrcomprisk::dat_ckid
  dat_final<-npcrest(data, exit, event, exposure=b1nb0, entry, eoi=2)
  expect_equal(nrow(dat_final$cuminc), 187)
  expect_equal(ncol(dat_final$cuminc), 11)

})
