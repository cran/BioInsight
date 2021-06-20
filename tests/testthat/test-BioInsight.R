test_that("Verify Variable Type: Count as Integer",{
  sample_of_event = c(1L, 2L, 3L)
  expect_equal(typeof(sample_of_event), 'integer')
})

test_that("Verify Variable Type: Annotation Matrix",{
  gene = 'ABC1'
  gene_biotype = 'protein_coding'
  start_position = 1L
  end_position = 999L
  expect_equal(typeof(gene), 'character')
  expect_equal(typeof(gene_biotype), 'character')
  expect_equal(typeof(start_position), 'integer')
  expect_equal(typeof(end_position), 'integer')
})

test_that("Verify Variable Type: Groups as Factor",{
  groups = rep(as.factor(c("1","2")), each = 5)
  x = is.factor(groups)
  expect_true(x)
})
