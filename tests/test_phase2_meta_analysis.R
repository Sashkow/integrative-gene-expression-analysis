#' Tests for Phase 2: Meta-Analysis Methods
#'
#' Run with: testthat::test_file("tests/test_phase2_meta_analysis.R")

library(testthat)

# Find and source the module
module_paths <- c(

  "scripts/integrative_analysis/phase2_meta_analysis/meta_analysis.R",
  "../scripts/integrative_analysis/phase2_meta_analysis/meta_analysis.R"
)

module_found <- FALSE
for (p in module_paths) {
  if (file.exists(p)) {
    source(p)
    module_found <- TRUE
    break
  }
}

if (!module_found) {
  stop("Could not find meta_analysis.R module")
}

# ============================================================================
# Test Data Setup
# ============================================================================

#' Create mock study data for testing
create_mock_studies <- function(n_studies = 3, n_genes = 100, n_samples = 10) {
  set.seed(42)
  studies <- list()

  genes <- paste0("gene", 1:n_genes)

  for (i in 1:n_studies) {
    study_id <- paste0("Study", i)
    samples <- paste0(study_id, "_s", 1:n_samples)

    # Create expression matrix with some differential expression
    expr <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 2),
                   nrow = n_genes, dimnames = list(genes, samples))

    # Add differential expression for first 20 genes
    n_control <- n_samples %/% 2
    n_case <- n_samples - n_control
    groups <- c(rep(0, n_control), rep(1, n_case))

    # Upregulate first 10 genes in cases
    expr[1:10, groups == 1] <- expr[1:10, groups == 1] + 2

    # Downregulate genes 11-20 in cases
    expr[11:20, groups == 1] <- expr[11:20, groups == 1] - 2

    studies[[study_id]] <- list(
      expr = expr,
      groups = groups,
      n_baseline = n_control,
      n_contrast = n_case
    )
  }

  return(studies)
}


# ============================================================================
# Tests for filter_balanced_studies
# ============================================================================

test_that("filter_balanced_studies keeps only balanced studies", {
  studies <- list(
    Study1 = list(n_baseline = 5, n_contrast = 5),
    Study2 = list(n_baseline = 0, n_contrast = 10),  # Unbalanced
    Study3 = list(n_baseline = 3, n_contrast = 4)
  )

  balanced <- filter_balanced_studies(studies, min_per_group = 2)

  expect_equal(length(balanced), 2)
  expect_true("Study1" %in% names(balanced))
  expect_true("Study3" %in% names(balanced))
  expect_false("Study2" %in% names(balanced))
})

test_that("filter_balanced_studies respects min_per_group", {
  studies <- list(
    Study1 = list(n_baseline = 5, n_contrast = 5),
    Study2 = list(n_baseline = 1, n_contrast = 10),  # Below threshold
    Study3 = list(n_baseline = 3, n_contrast = 2)
  )

  balanced <- filter_balanced_studies(studies, min_per_group = 3)

  expect_equal(length(balanced), 1)
  expect_true("Study1" %in% names(balanced))
})

test_that("filter_balanced_studies returns empty list when no studies qualify", {
  studies <- list(
    Study1 = list(n_baseline = 0, n_contrast = 5),
    Study2 = list(n_baseline = 5, n_contrast = 0)
  )

  balanced <- filter_balanced_studies(studies, min_per_group = 2)

  expect_equal(length(balanced), 0)
})


# ============================================================================
# Tests for run_limma_single_study
# ============================================================================

test_that("run_limma_single_study returns expected columns", {
  mock <- create_mock_studies(n_studies = 1, n_genes = 50, n_samples = 10)
  study <- mock[[1]]

  result <- run_limma_single_study(study$expr, study$groups)

  expect_true(is.data.frame(result))
  expect_true("logFC" %in% colnames(result))
  expect_true("P.Value" %in% colnames(result))
  expect_true("adj.P.Val" %in% colnames(result))
  expect_true("gene_id" %in% colnames(result))
  expect_equal(nrow(result), 50)
})

test_that("run_limma_single_study detects differential expression", {
  mock <- create_mock_studies(n_studies = 1, n_genes = 50, n_samples = 20)
  study <- mock[[1]]

  result <- run_limma_single_study(study$expr, study$groups)

  # First 10 genes should be upregulated (positive logFC)
  up_genes <- result[result$gene_id %in% paste0("gene", 1:10), ]
  expect_true(mean(up_genes$logFC) > 0)

  # Genes 11-20 should be downregulated (negative logFC)
  down_genes <- result[result$gene_id %in% paste0("gene", 11:20), ]
  expect_true(mean(down_genes$logFC) < 0)
})


# ============================================================================
# Tests for run_per_study_de
# ============================================================================

test_that("run_per_study_de runs on all studies", {
  mock <- create_mock_studies(n_studies = 3, n_genes = 30, n_samples = 8)

  results <- run_per_study_de(mock)

  expect_equal(length(results), 3)
  expect_true(all(c("Study1", "Study2", "Study3") %in% names(results)))

  for (study_id in names(results)) {
    expect_true("study" %in% colnames(results[[study_id]]))
    expect_equal(unique(results[[study_id]]$study), study_id)
  }
})


# ============================================================================
# Tests for prepare_dexma_input
# ============================================================================

test_that("prepare_dexma_input creates correct structure", {
  skip_if_not_installed("DExMA")

  mock <- create_mock_studies(n_studies = 2, n_genes = 50, n_samples = 10)

  dexma_input <- prepare_dexma_input(mock)

  expect_equal(length(dexma_input), 2)

  for (study_id in names(dexma_input)) {
    expect_true("GEX" %in% names(dexma_input[[study_id]]))
    expect_true("Sample_Pheno" %in% names(dexma_input[[study_id]]))
    expect_true("group" %in% colnames(dexma_input[[study_id]]$Sample_Pheno))

    # Check that pheno rownames match expression colnames
    expect_equal(
      rownames(dexma_input[[study_id]]$Sample_Pheno),
      colnames(dexma_input[[study_id]]$GEX)
    )
  }
})

test_that("prepare_dexma_input creates correct group factors", {
  skip_if_not_installed("DExMA")

  mock <- create_mock_studies(n_studies = 1, n_genes = 20, n_samples = 10)

  dexma_input <- prepare_dexma_input(mock)

  groups <- dexma_input$Study1$Sample_Pheno$group
  expect_true(is.factor(groups))
  expect_equal(levels(groups), c("control", "case"))
  expect_equal(sum(groups == "control"), mock$Study1$n_baseline)
  expect_equal(sum(groups == "case"), mock$Study1$n_contrast)
})


# ============================================================================
# Tests for prepare_rankprod_input
# ============================================================================

test_that("prepare_rankprod_input creates combined matrix", {
  skip_if_not_installed("RankProd")

  mock <- create_mock_studies(n_studies = 3, n_genes = 50, n_samples = 8)

  rp_input <- prepare_rankprod_input(mock, min_studies = 2)

  expect_true("expr" %in% names(rp_input))
  expect_true("groups" %in% names(rp_input))
  expect_true("origin" %in% names(rp_input))

  # Combined samples should equal sum of individual studies
  total_samples <- sum(sapply(mock, function(s) ncol(s$expr)))
  expect_equal(rp_input$n_samples, total_samples)

  # Origin should have correct number of unique values
  expect_equal(length(unique(rp_input$origin)), 3)
})

test_that("prepare_rankprod_input filters by min_studies", {
  skip_if_not_installed("RankProd")

  # Create studies with different gene sets
  studies <- list(
    Study1 = list(
      expr = matrix(1:20, nrow = 5, dimnames = list(c("A", "B", "C", "D", "E"), paste0("s1_", 1:4))),
      groups = c(0, 0, 1, 1),
      n_baseline = 2, n_contrast = 2
    ),
    Study2 = list(
      expr = matrix(1:16, nrow = 4, dimnames = list(c("A", "B", "C", "F"), paste0("s2_", 1:4))),
      groups = c(0, 0, 1, 1),
      n_baseline = 2, n_contrast = 2
    )
  )

  # Only genes A, B, C appear in both studies
  rp_input <- prepare_rankprod_input(studies, min_studies = 2)

  expect_equal(nrow(rp_input$expr), 3)  # A, B, C
  expect_true(all(rownames(rp_input$expr) %in% c("A", "B", "C")))
})


# ============================================================================
# Tests for run_metafor_meta (integration test)
# ============================================================================

test_that("run_metafor_meta runs without error", {
  skip_if_not_installed("metafor")

  mock <- create_mock_studies(n_studies = 3, n_genes = 30, n_samples = 10)

  result <- run_metafor_meta(mock, min_studies = 2)

  expect_true(is.data.frame(result))
  expect_true("gene_id" %in% colnames(result))
  expect_true("logFC" %in% colnames(result))
  expect_true("pvalue" %in% colnames(result))
  expect_true("fdr" %in% colnames(result))
  expect_true(nrow(result) > 0)
})

test_that("run_metafor_meta detects true DEGs", {
  skip_if_not_installed("metafor")

  # Use larger sample size for more power
  mock <- create_mock_studies(n_studies = 3, n_genes = 30, n_samples = 20)

  result <- run_metafor_meta(mock, min_studies = 2)

  # Check that upregulated genes have positive logFC
  up_results <- result[result$gene_id %in% paste0("gene", 1:10), ]
  expect_true(mean(up_results$logFC, na.rm = TRUE) > 0)

  # Check that downregulated genes have negative logFC
  down_results <- result[result$gene_id %in% paste0("gene", 11:20), ]
  expect_true(mean(down_results$logFC, na.rm = TRUE) < 0)
})


# ============================================================================
# Tests for run_dexma_meta (integration test)
# ============================================================================

test_that("run_dexma_meta runs without error", {
  skip_if_not_installed("DExMA")

  mock <- create_mock_studies(n_studies = 2, n_genes = 50, n_samples = 10)

  result <- run_dexma_meta(mock, effect_size = "SMD", missAllow = 0.3)

  expect_true(is.data.frame(result))
  expect_true("gene_id" %in% colnames(result))
  expect_true("logFC" %in% colnames(result))
  expect_true("fdr" %in% colnames(result))
  expect_equal(nrow(result), 50)
})

test_that("run_dexma_meta detects true DEGs", {
  skip_if_not_installed("DExMA")

  mock <- create_mock_studies(n_studies = 2, n_genes = 30, n_samples = 16)

  result <- run_dexma_meta(mock, effect_size = "SMD", missAllow = 0.3)

  # Check direction of effect
  up_results <- result[result$gene_id %in% paste0("gene", 1:10), ]
  expect_true(mean(up_results$logFC, na.rm = TRUE) > 0)

  down_results <- result[result$gene_id %in% paste0("gene", 11:20), ]
  expect_true(mean(down_results$logFC, na.rm = TRUE) < 0)
})


# ============================================================================
# Tests for run_rankprod_meta (integration test)
# ============================================================================

test_that("run_rankprod_meta runs without error", {
  skip_if_not_installed("RankProd")

  mock <- create_mock_studies(n_studies = 2, n_genes = 50, n_samples = 10)

  result <- run_rankprod_meta(mock, min_studies = 2, num_perm = 100, logged = TRUE)

  expect_true(is.data.frame(result))
  expect_true("gene_id" %in% colnames(result))
  expect_true("direction" %in% colnames(result))
  expect_true("pfp" %in% colnames(result))
  expect_equal(nrow(result), 50)
})

test_that("run_rankprod_meta identifies direction correctly", {
  skip_if_not_installed("RankProd")

  mock <- create_mock_studies(n_studies = 2, n_genes = 30, n_samples = 16)

  result <- run_rankprod_meta(mock, min_studies = 2, num_perm = 100, logged = TRUE)

  # Upregulated genes should mostly have direction "up"
  up_genes <- result[result$gene_id %in% paste0("gene", 1:10), ]
  expect_true(sum(up_genes$direction == "up") > sum(up_genes$direction == "down"))

  # Downregulated genes should mostly have direction "down"
  down_genes <- result[result$gene_id %in% paste0("gene", 11:20), ]
  expect_true(sum(down_genes$direction == "down") > sum(down_genes$direction == "up"))
})


# ============================================================================
# Tests for combine_meta_results
# ============================================================================

test_that("combine_meta_results merges results correctly", {
  dexma_results <- data.frame(
    gene_id = c("A", "B", "C", "D"),
    logFC = c(1.5, -1.2, 0.5, 0.1),
    fdr = c(0.01, 0.02, 0.1, 0.5)
  )

  rankprod_results <- data.frame(
    gene_id = c("A", "B", "C", "E"),
    direction = c("up", "down", "up", "down"),
    pfp = c(0.01, 0.03, 0.2, 0.04)
  )

  combined <- combine_meta_results(dexma_results, rankprod_results, fdr_threshold = 0.05)

  expect_true(is.data.frame(combined))
  expect_equal(nrow(combined), 5)  # A, B, C, D, E

  # Check categorization
  expect_equal(combined$category[combined$gene_id == "A"], "both")
  expect_equal(combined$category[combined$gene_id == "B"], "both")
  expect_equal(combined$category[combined$gene_id == "D"], "not_significant")
  expect_equal(combined$category[combined$gene_id == "E"], "rankprod_only")
})

test_that("combine_meta_results checks direction agreement", {
  dexma_results <- data.frame(
    gene_id = c("A", "B"),
    logFC = c(1.5, -1.2),
    fdr = c(0.01, 0.01)
  )

  rankprod_results <- data.frame(
    gene_id = c("A", "B"),
    direction = c("up", "down"),
    pfp = c(0.01, 0.01)
  )

  combined <- combine_meta_results(dexma_results, rankprod_results, fdr_threshold = 0.05)

  # Both should have direction agreement
  expect_true(all(combined$direction_agree[combined$category == "both"]))
})


# ============================================================================
# Tests for print_meta_summary
# ============================================================================

test_that("print_meta_summary runs without error", {
  results <- data.frame(
    gene_id = paste0("gene", 1:100),
    fdr = c(rep(0.001, 10), rep(0.02, 20), rep(0.08, 30), rep(0.5, 40))
  )

  # Should not throw error
  expect_output(print_meta_summary(results, "TestMethod", "fdr"))
})


# ============================================================================
# Run all tests
# ============================================================================

if (interactive()) {
  cat("\n========================================\n")
  cat("Running Phase 2 Meta-Analysis Tests\n")
  cat("========================================\n\n")
}
