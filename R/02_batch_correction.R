#' Batch Correction Module
#'
#' Functions for batch effect correction using ComBat
#'
#' @author Expression Integration Pipeline
#' @date 2025

#' Apply ComBat batch correction
#'
#' @param exprs Expression matrix (genes x samples)
#' @param pdata Phenodata with batch and covariate information
#' @param batch_col Column name for batch variable (default: "secondaryaccession")
#' @param mod_formula Formula for model matrix (default: ~trim_term)
#' @param par_prior Use parametric priors (default: TRUE)
#' @param prior_plots Generate prior plots (default: FALSE)
#' @return Batch-corrected expression matrix
#' @export
apply_combat_correction <- function(exprs,
                                     pdata,
                                     batch_col = "secondaryaccession",
                                     mod_formula = ~trim_term,
                                     par_prior = TRUE,
                                     prior_plots = FALSE) {

  library(sva)

  # Create batch variable
  batch <- as.factor(pdata[[batch_col]])

  message("Batch correction using '", batch_col, "': ", nlevels(batch), " batches")
  message("  Batch sizes: ", paste(table(batch), collapse = ", "))

  # Create model matrix for biological variables
  mod <- model.matrix(mod_formula, data = pdata)
  message("  Adjusting for ", ncol(mod), " covariate(s)")

  # Check for confounding
  check_confounding(pdata, batch_col, mod_formula)

  # Apply ComBat
  exprs_corrected <- ComBat(
    dat = as.matrix(exprs),
    batch = batch,
    mod = mod,
    par.prior = par_prior,
    prior.plots = prior_plots
  )

  message("Batch correction completed")

  return(exprs_corrected)
}


#' Check for batch-covariate confounding
#'
#' @param pdata Phenodata
#' @param batch_col Batch variable column name
#' @param mod_formula Model formula
#' @return Invisibly returns test results
#' @export
check_confounding <- function(pdata, batch_col, mod_formula) {

  # Extract variables from formula
  vars <- all.vars(mod_formula)

  # Test each variable against batch
  for (var in vars) {
    if (var %in% colnames(pdata)) {

      tbl <- table(pdata[[batch_col]], pdata[[var]])

      # Fisher's exact test for small tables, chi-square for larger
      if (prod(dim(tbl)) <= 100) {
        test <- fisher.test(tbl, simulate.p.value = TRUE)
        message("  Fisher test (", batch_col, " vs ", var, "): p = ",
                format.pval(test$p.value, digits = 3))
      } else {
        test <- chisq.test(tbl)
        message("  Chi-square test (", batch_col, " vs ", var, "): p = ",
                format.pval(test$p.value, digits = 3))
      }

      if (test$p.value < 0.05) {
        warning("Potential confounding between ", batch_col, " and ", var)
      }
    }
  }

  invisible(NULL)
}


#' Create trimester variable from gestational age
#'
#' @param pdata Phenodata
#' @param ga_col Gestational age column name
#' @return Phenodata with added trim_term column
#' @export
add_trimester_variable <- function(pdata, ga_col = "Gestational.Age.Category") {

  # If trim_term already exists, return as-is
  if ("trim_term" %in% colnames(pdata)) {
    message("trim_term column already exists")
    message("Trimester distribution:")
    print(table(pdata$trim_term))
    return(pdata)
  }

  # Use Gestational.Age.Category if available
  if ("Gestational.Age.Category" %in% colnames(pdata)) {
    pdata$trim_term <- pdata$Gestational.Age.Category

    # Standardize term names
    pdata$trim_term <- ifelse(pdata$trim_term %in% c("Early Preterm", "Late Preterm"),
                              "Preterm", pdata$trim_term)
    pdata$trim_term <- ifelse(pdata$trim_term == "Unknown Gestational Category",
                              NA, pdata$trim_term)

    message("Created trim_term from Gestational.Age.Category")
  }
  # Otherwise use numeric gestational age if available
  else if (ga_col %in% colnames(pdata) && is.numeric(pdata[[ga_col]])) {
    pdata$trim_term <- ifelse(pdata[[ga_col]] <= 12, "First Trimester", "")
    pdata$trim_term <- ifelse(pdata[[ga_col]] > 12 & pdata[[ga_col]] <= 27,
                              "Second Trimester", pdata$trim_term)
    pdata$trim_term <- ifelse(pdata[[ga_col]] > 27 & pdata[[ga_col]] <= 37,
                              "Preterm", pdata$trim_term)
    pdata$trim_term <- ifelse(pdata[[ga_col]] > 37, "Term", pdata$trim_term)

    message("Created trim_term from numeric gestational age")
  } else {
    stop("Cannot create trim_term: no suitable gestational age column found")
  }

  message("Trimester distribution:")
  print(table(pdata$trim_term, useNA = "always"))

  return(pdata)
}
