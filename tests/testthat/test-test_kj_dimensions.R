test_that("emuFit stops when test_kj indices are out of bounds", {
  # 1. Load small example data
  data(wirbel_otu_small)
  data(wirbel_sample_small)
  
  # Case 1: When k is out of bounds. e.g. p = 2 (Group = CRC or CTR) so k = 3 is invalid
  bad_k = data.frame(k = 3, j = 1)
  
  expect_error(
    emuFit(formula = ~ Group,
           data    = wirbel_sample_small,
           Y       = wirbel_otu_small,
           test_kj = bad_k),
    "k exceeds the number of groups p"
  )
  
  # Case 2: Test j out of bounds (J = 47)
  bad_j = data.frame(k = 2, j = 48)
  expect_error(
    emuFit(formula = ~ Group,
           data    = wirbel_sample_small,
           Y       = wirbel_otu_small,
           test_kj = bad_j),
    "j exceeds the number of taxa J"
  )
  
  # Case 3: Bad string k
  expect_error(
    emuFit(formula = ~ Group,
           data    = wirbel_sample_small,
           Y       = wirbel_otu_small,
           test_kj = data.frame(k = "FakeVariable", j = 1)
  ),
  "Make sure that the values of `k` in `test_kj` are numeric or correspond to column names of the `X` matrix."
  )
  
  # Case 4: No issue case should still work
  expect_silent(
    emuFit(
      formula = ~ Group,
      data    = wirbel_sample_small,
      Y       = wirbel_otu_small,
      test_kj = data.frame(k = 2, j = 47)
    )
  )
  
  # Case 5: No issue if valid string
  expect_silent(
    emuFit(
      formula = ~ Group,
      data    = wirbel_sample_small,
      Y       = wirbel_otu_small,
      test_kj = data.frame(k = "GroupCTR", j = 47)
    )
  )
})