#' Tests for Phase 1: Tiered Gene Coverage Analysis
#'
#' TDD tests for gene coverage calculation and tiered merging functions
#'
#' @author Expression Integration Pipeline
#' @date 2025

library(testthat)

# Source the Phase 1 module - handle different working directories
# When run via testthat::test_file(), working directory is project root
# When run directly, may need to adjust path
module_path <- "scripts/integrative_analysis/phase1_tiered_coverage/tiered_coverage.R"
if (!file.exists(module_path)) {
  # Try parent directory (if running from tests folder)
  module_path <- file.path("..", module_path)
}
if (!file.exists(module_path)) {
  stop("Cannot find tiered_coverage.R module")
}
source(module_path)

# ============================================================================
# Test Data Setup
# ============================================================================

#' Create mock expression data for testing
create_mock_datasets <- function() {
  # Dataset 1: genes A, B, C, D
  dataset1 <- data.frame(
    row.names = c("1", "2", "3", "4"),  # ENTREZID style
    Sample1 = c(5.1, 6.2, 7.3, 8.4),
    Sample2 = c(5.2, 6.3, 7.4, 8.5)
  )

  # Dataset 2: genes A, B, E (missing C, D)
  dataset2 <- data.frame(
    row.names = c("1", "2", "5"),
    Sample3 = c(5.3, 6.4, 9.1),
    Sample4 = c(5.4, 6.5, 9.2)
  )

  # Dataset 3: genes A, C, E, F (missing B, D)
  dataset3 <- data.frame(
    row.names = c("1", "3", "5", "6"),
    Sample5 = c(5.5, 7.5, 9.3, 10.1),
    Sample6 = c(5.6, 7.6, 9.4, 10.2)
  )

  # Dataset 4: genes A, B, C, D, E, F (all genes)
  dataset4 <- data.frame(
    row.names = c("1", "2", "3", "4", "5", "6"),
    Sample7 = c(5.7, 6.7, 7.7, 8.7, 9.5, 10.3),
    Sample8 = c(5.8, 6.8, 7.8, 8.8, 9.6, 10.4)
  )

  return(list(
    dataset1 = dataset1,
    dataset2 = dataset2,
    dataset3 = dataset3,
    dataset4 = dataset4
  ))
}

# Expected gene coverage:
# Gene 1 (A): in 4/4 = 100%
# Gene 2 (B): in 3/4 = 75%  (missing from dataset3)
# Gene 3 (C): in 3/4 = 75%  (missing from dataset2)
# Gene 4 (D): in 2/4 = 50%  (missing from dataset2, dataset3)
# Gene 5 (E): in 3/4 = 75%  (missing from dataset1)
# Gene 6 (F): in 2/4 = 50%  (missing from dataset1, dataset2)

# ============================================================================
# Tests for get_gene_presence_matrix()
# ============================================================================

test_that("get_gene_presence_matrix returns correct dimensions", {
  datasets <- create_mock_datasets()

  presence_matrix <- get_gene_presence_matrix(datasets)

  # Should have 6 unique genes (rows) and 4 datasets (cols)
  expect_equal(nrow(presence_matrix), 6)
  expect_equal(ncol(presence_matrix), 4)
})

test_that("get_gene_presence_matrix correctly identifies gene presence", {
  datasets <- create_mock_datasets()

  presence_matrix <- get_gene_presence_matrix(datasets)

  # Gene 1 should be in all datasets
  expect_equal(sum(presence_matrix["1", ]), 4)

  # Gene 2 should be in 3 datasets (not dataset3)
  expect_equal(sum(presence_matrix["2", ]), 3)
  expect_equal(presence_matrix["2", "dataset3"], 0)

  # Gene 6 should be in 2 datasets (dataset3, dataset4)
  expect_equal(sum(presence_matrix["6", ]), 2)
  expect_equal(presence_matrix["6", "dataset1"], 0)
  expect_equal(presence_matrix["6", "dataset2"], 0)
})

# ============================================================================
# Tests for calculate_gene_coverage()
# ============================================================================

test_that("calculate_gene_coverage returns correct percentages", {
  datasets <- create_mock_datasets()

  coverage <- calculate_gene_coverage(datasets)

  expect_equal(unname(coverage["1"]), 1.0)       # 4/4
  expect_equal(unname(coverage["2"]), 0.75)      # 3/4
  expect_equal(unname(coverage["3"]), 0.75)      # 3/4
  expect_equal(unname(coverage["4"]), 0.50)      # 2/4
  expect_equal(unname(coverage["5"]), 0.75)      # 3/4
  expect_equal(unname(coverage["6"]), 0.50)      # 2/4
})

test_that("calculate_gene_coverage handles single dataset", {
  datasets <- list(dataset1 = create_mock_datasets()$dataset1)

  coverage <- calculate_gene_coverage(datasets)

  # All genes should have 100% coverage in single dataset
  expect_true(all(coverage == 1.0))
})

test_that("calculate_gene_coverage handles empty list", {
  expect_error(calculate_gene_coverage(list()))
})

# ============================================================================
# Tests for filter_genes_by_coverage()
# ============================================================================

test_that("filter_genes_by_coverage at 100% returns only universal genes", {
  datasets <- create_mock_datasets()
  coverage <- calculate_gene_coverage(datasets)

  tier1_genes <- filter_genes_by_coverage(coverage, threshold = 1.0)

  expect_equal(length(tier1_genes), 1)
  expect_equal(tier1_genes, "1")
})

test_that("filter_genes_by_coverage at 75% returns correct genes", {
  datasets <- create_mock_datasets()
  coverage <- calculate_gene_coverage(datasets)

  tier2_genes <- filter_genes_by_coverage(coverage, threshold = 0.75)

  # Should include genes 1, 2, 3, 5 (all >= 75%)
  expect_equal(length(tier2_genes), 4)
  expect_true(all(c("1", "2", "3", "5") %in% tier2_genes))
  expect_false("4" %in% tier2_genes)
  expect_false("6" %in% tier2_genes)
})

test_that("filter_genes_by_coverage at 50% returns more genes", {
  datasets <- create_mock_datasets()
  coverage <- calculate_gene_coverage(datasets)

  tier3_genes <- filter_genes_by_coverage(coverage, threshold = 0.50)

  # Should include all 6 genes
  expect_equal(length(tier3_genes), 6)
})

test_that("filter_genes_by_coverage at 0% returns all genes", {
  datasets <- create_mock_datasets()
  coverage <- calculate_gene_coverage(datasets)

  all_genes <- filter_genes_by_coverage(coverage, threshold = 0.0)

  expect_equal(length(all_genes), 6)
})

# ============================================================================
# Tests for generate_coverage_report()
# ============================================================================

test_that("generate_coverage_report returns expected structure", {
  datasets <- create_mock_datasets()

  report <- generate_coverage_report(datasets)

  # Should contain required elements
  expect_true("presence_matrix" %in% names(report))
  expect_true("gene_coverage" %in% names(report))
  expect_true("summary" %in% names(report))
  expect_true("tier_counts" %in% names(report))
})

test_that("generate_coverage_report tier_counts are correct", {
  datasets <- create_mock_datasets()

  report <- generate_coverage_report(datasets)

  # Check tier counts
  expect_equal(report$tier_counts$tier1_100pct, 1)   # Only gene 1
  expect_equal(report$tier_counts$tier2_75pct, 4)   # Genes 1,2,3,5
  expect_equal(report$tier_counts$tier3_50pct, 6)   # All genes
})

# ============================================================================
# Tests for merge_expression_data_tiered()
# ============================================================================

test_that("merge_expression_data_tiered at 100% matches current behavior", {
  skip("Requires file-based test data - run with actual data")

  # This test compares tiered merge at 100% with original merge
  # Should produce identical results
})

test_that("merge_expression_data_tiered at 75% includes more genes", {
  skip("Requires file-based test data - run with actual data")

  # Tiered merge at 75% should include more genes than 100%
})

test_that("merge_expression_data_tiered fills NA for missing values", {
  skip("Requires file-based test data - run with actual data")

  # When a gene is missing from a dataset, should be NA
})

# ============================================================================
# Tests for merge_datasets_tiered() (in-memory version)
# ============================================================================

test_that("merge_datasets_tiered at 100% keeps only common genes", {
  datasets <- create_mock_datasets()

  merged <- merge_datasets_tiered(datasets, coverage_threshold = 1.0)

  # Should only have gene 1 (the only one in all 4 datasets)
  expect_equal(nrow(merged), 1)
  expect_equal(rownames(merged), "1")

  # Should have all 8 samples
  expect_equal(ncol(merged), 8)
})

test_that("merge_datasets_tiered at 75% includes genes in >= 3 datasets", {
  datasets <- create_mock_datasets()

  merged <- merge_datasets_tiered(datasets, coverage_threshold = 0.75)

  # Should have genes 1, 2, 3, 5
  expect_equal(nrow(merged), 4)
  expect_true(all(c("1", "2", "3", "5") %in% rownames(merged)))
})

test_that("merge_datasets_tiered at 50% includes all genes", {
  datasets <- create_mock_datasets()

  merged <- merge_datasets_tiered(datasets, coverage_threshold = 0.50)

  # Should have all 6 genes
  expect_equal(nrow(merged), 6)
})

test_that("merge_datasets_tiered produces NA for missing genes", {
  datasets <- create_mock_datasets()

  merged <- merge_datasets_tiered(datasets, coverage_threshold = 0.50)

  # Gene 6 is missing from datasets 1 and 2
  # Samples 1-4 are from datasets 1 and 2, should be NA for gene 6
  expect_true(is.na(merged["6", "Sample1"]))
  expect_true(is.na(merged["6", "Sample2"]))
  expect_true(is.na(merged["6", "Sample3"]))
  expect_true(is.na(merged["6", "Sample4"]))

  # Samples 5-8 are from datasets 3 and 4, should have values
  expect_false(is.na(merged["6", "Sample5"]))
  expect_false(is.na(merged["6", "Sample7"]))
})

test_that("merge_datasets_tiered preserves expression values correctly", {
  datasets <- create_mock_datasets()

  merged <- merge_datasets_tiered(datasets, coverage_threshold = 1.0)

  # Check that values are preserved for gene 1
  expect_equal(merged["1", "Sample1"], 5.1)
  expect_equal(merged["1", "Sample3"], 5.3)
  expect_equal(merged["1", "Sample5"], 5.5)
  expect_equal(merged["1", "Sample7"], 5.7)
})

# ============================================================================
# Integration Tests
# ============================================================================

test_that("tiered merge produces more genes than inner join", {
  datasets <- create_mock_datasets()

  # Inner join (current behavior) - only genes in ALL datasets
  inner_join_genes <- Reduce(intersect, lapply(datasets, rownames))

  # Tiered at 75%
  merged_75 <- merge_datasets_tiered(datasets, coverage_threshold = 0.75)

  # Tiered should have more genes

  expect_gt(nrow(merged_75), length(inner_join_genes))
})

test_that("gene counts decrease monotonically with higher thresholds", {
  datasets <- create_mock_datasets()

  merged_50 <- merge_datasets_tiered(datasets, coverage_threshold = 0.50)
  merged_75 <- merge_datasets_tiered(datasets, coverage_threshold = 0.75)
  merged_100 <- merge_datasets_tiered(datasets, coverage_threshold = 1.0)

  expect_gte(nrow(merged_50), nrow(merged_75))
  expect_gte(nrow(merged_75), nrow(merged_100))
})

# ============================================================================
# Run tests if script is executed directly
# ============================================================================

if (sys.nframe() == 0) {
  test_results <- test_file("tests/test_phase1_tiered_coverage.R")
  print(test_results)
}
