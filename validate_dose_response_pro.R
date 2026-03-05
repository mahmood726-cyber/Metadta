# ============================================================================
# Dose Response Pro v18.1 - R Validation Script
# ============================================================================
# Purpose: Validate the CLI implementation directly against R packages.
# Packages: dosresmeta, metafor, mvmeta
#
# Usage:
#   "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" tests/validate_dose_response_pro.R
#
# Outputs:
#   tests/r_validation_results.json
#   tests/r_validation_results.txt
# ============================================================================

required_packages <- c("dosresmeta", "metafor", "mvmeta", "jsonlite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(sprintf("Missing required packages: %s", paste(missing_packages, collapse = ", ")))
}
invisible(lapply(required_packages, function(pkg) library(pkg, character.only = TRUE)))

VALIDATION_BETA_TOLERANCE <- 0.02
VALIDATION_SE_TOLERANCE <- 0.05
VALIDATION_TAU2_TOLERANCE <- 0.02
MIN_POSITIVE <- 1e-12

# ============================================================================
# TEST DATASETS
# ============================================================================

create_test_data_1 <- function() {
  data.frame(
    study = rep(c("Study1", "Study2", "Study3"), each = 4),
    dose = rep(c(0, 1, 2, 3), 3),
    cases = c(45, 52, 61, 58, 38, 44, 51, NA, 52, 58, 63, 71),
    n = c(50000, 48000, 45000, 40000, 42000, 40000, 38000, NA, 55000, 52000, 48000, 45000),
    stringsAsFactors = FALSE
  )
}

create_test_data_2 <- function() {
  data.frame(
    study = rep(c("Study1", "Study2", "Study3", "Study4"), each = 5),
    dose = rep(c(0, 1, 2, 3, 4), 4),
    cases = c(
      30, 35, 45, 50, 55,
      25, 30, 38, 42, 48,
      28, 33, 40, 46, 52,
      22, 28, 35, 40, 45
    ),
    n = rep(40000, 20),
    stringsAsFactors = FALSE
  )
}

create_test_data_3 <- function() {
  data.frame(
    study = rep(c("Study1", "Study2", "Study3", "Study4", "Study5"), each = 3),
    dose = rep(c(0, 1, 2), 5),
    cases = c(
      20, 40, 60,
      15, 25, 35,
      25, 35, 50,
      30, 50, 75,
      18, 30, 45
    ),
    n = c(
      30000, 30000, 30000,
      25000, 25000, 25000,
      35000, 35000, 35000,
      40000, 40000, 40000,
      28000, 28000, 28000
    ),
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# PATH + CLI HELPERS
# ============================================================================

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) == 0) return(getwd())
  normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = FALSE)
}

get_project_root <- function() {
  normalizePath(file.path(get_script_dir(), ".."), winslash = "/", mustWork = TRUE)
}

discover_python <- function() {
  py <- Sys.which("python")
  if (!identical(py, "")) return(py)

  py3 <- Sys.which("python3")
  if (!identical(py3, "")) return(py3)

  stop("Python executable not found in PATH (python/python3).")
}

run_cli_analysis <- function(test_data, model = "quadratic") {
  project_root <- get_project_root()
  cli_script <- file.path(project_root, "dose-response-cli.py")
  if (!file.exists(cli_script)) {
    stop(sprintf("CLI script not found: %s", cli_script))
  }

  python_bin <- discover_python()
  input_csv <- tempfile(pattern = "dose_cli_input_", fileext = ".csv")
  output_json <- tempfile(pattern = "dose_cli_output_", fileext = ".json")
  on.exit(unlink(c(input_csv, output_json), force = TRUE), add = TRUE)

  cli_input <- data.frame(
    study_id = test_data$study,
    dose = test_data$dose,
    cases = test_data$cases,
    n = test_data$n,
    stringsAsFactors = FALSE
  )
  write.csv(cli_input, input_csv, row.names = FALSE, na = "")

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)

  args <- c(
    "dose-response-cli.py",
    "--input", input_csv,
    "--output", output_json,
    "--model", model
  )
  cmd_out <- system2(python_bin, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(cmd_out, "status")
  if (is.null(status)) status <- 0L

  if (!identical(status, 0L)) {
    stop(sprintf(
      "CLI execution failed (status=%d): %s",
      as.integer(status),
      paste(tail(cmd_out, 20), collapse = "\n")
    ))
  }
  if (!file.exists(output_json)) {
    stop("CLI did not produce expected output JSON.")
  }

  parsed <- jsonlite::fromJSON(output_json, simplifyVector = TRUE)
  if (is.null(parsed$beta) || length(parsed$beta) < 3) {
    stop("CLI output missing expected quadratic coefficients.")
  }
  if (is.null(parsed$se) || length(parsed$se) < 3) {
    stop("CLI output missing expected quadratic standard errors.")
  }
  parsed
}

# ============================================================================
# DOSRESMETA PREP
# ============================================================================

prepare_dosresmeta_data <- function(test_data) {
  clean <- test_data[complete.cases(test_data), c("study", "dose", "cases", "n")]
  clean <- clean[order(clean$study, clean$dose), ]
  per_study <- split(clean, clean$study)

  transformed <- lapply(per_study, function(df) {
    if (nrow(df) < 3) {
      stop(sprintf("Study '%s' needs at least 3 dose rows for quadratic fit", as.character(df$study[1])))
    }

    ref_idx <- which.min(df$dose)
    ref_cases <- df$cases[ref_idx]
    ref_n <- df$n[ref_idx]

    rr <- (df$cases / df$n) / (ref_cases / ref_n)
    logrr <- log(rr)

    se <- rep(NA_real_, nrow(df))
    for (i in seq_len(nrow(df))) {
      if (i != ref_idx) {
        se[i] <- sqrt(max((1 / df$cases[i]) - (1 / df$n[i]) + (1 / ref_cases) - (1 / ref_n), MIN_POSITIVE))
      }
    }

    df$logrr <- logrr
    df$se <- se
    df$type <- "ir"
    df
  })

  out <- do.call(rbind, transformed)
  rownames(out) <- NULL
  out
}

# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================

validate_dosresmeta <- function(test_name, test_data) {
  cat(sprintf("\n--- %s ---\n", test_name))

  tryCatch({
    cli_results <- run_cli_analysis(test_data, model = "quadratic")
    dr_data <- prepare_dosresmeta_data(test_data)

    dr_fit <- dosresmeta::dosresmeta(
      formula = logrr ~ dose + I(dose^2),
      id = study,
      type = type,
      se = se,
      cases = cases,
      n = n,
      data = dr_data
    )

    dr_summary <- summary(dr_fit)

    our_beta <- as.numeric(cli_results$beta[2:3])
    our_se <- as.numeric(cli_results$se[2:3])
    our_tau2 <- as.numeric(cli_results$tau2)
    r_beta <- as.numeric(stats::coef(dr_fit))
    r_se <- as.numeric(sqrt(diag(stats::vcov(dr_fit))))

    if (length(r_beta) != 2 || length(r_se) != 2) {
      stop("Unexpected dosresmeta output shape; expected 2 coefficients")
    }

    beta_diff <- our_beta - r_beta
    se_diff <- our_se - r_se

    comparison <- data.frame(
      parameter = c("linear", "quadratic"),
      our_beta = our_beta,
      r_beta = r_beta,
      beta_diff = beta_diff,
      abs_beta_diff = abs(beta_diff),
      our_se = our_se,
      r_se = r_se,
      se_diff = se_diff,
      abs_se_diff = abs(se_diff),
      stringsAsFactors = FALSE
    )

    beta_match <- all(comparison$abs_beta_diff <= VALIDATION_BETA_TOLERANCE)
    se_match <- all(comparison$abs_se_diff <= VALIDATION_SE_TOLERANCE)

    psi <- tryCatch(as.matrix(dr_summary$Psi), error = function(e) matrix(NA_real_, nrow = 1, ncol = 1))
    dr_tau2 <- if (all(!is.finite(diag(psi)))) NA_real_ else mean(diag(psi), na.rm = TRUE)
    tau2_diff <- if (is.finite(dr_tau2) && is.finite(our_tau2)) as.numeric(our_tau2 - dr_tau2) else NA_real_
    tau2_abs_diff <- if (is.finite(tau2_diff)) abs(tau2_diff) else NA_real_
    tau2_match <- isTRUE(is.finite(tau2_abs_diff) && tau2_abs_diff <= VALIDATION_TAU2_TOLERANCE)

    pass <- isTRUE(beta_match && se_match && tau2_match)

    cat(sprintf("  Beta match (tol %.3f): %s\n", VALIDATION_BETA_TOLERANCE, ifelse(beta_match, "PASS", "FAIL")))
    cat(sprintf("  SE match (tol %.3f): %s\n", VALIDATION_SE_TOLERANCE, ifelse(se_match, "PASS", "FAIL")))
    cat(sprintf("  Tau2 match (tol %.3f): %s\n", VALIDATION_TAU2_TOLERANCE, ifelse(tau2_match, "PASS", "FAIL")))

    list(
      name = test_name,
      pass = pass,
      fit_ok = TRUE,
      beta_match = isTRUE(beta_match),
      se_match = isTRUE(se_match),
      tau2_match = isTRUE(tau2_match),
      n_studies = length(unique(dr_data$study)),
      n_rows = nrow(dr_data),
      tolerance = list(
        beta = VALIDATION_BETA_TOLERANCE,
        se = VALIDATION_SE_TOLERANCE,
        tau2 = VALIDATION_TAU2_TOLERANCE
      ),
      tau2 = list(
        our = our_tau2,
        r = as.numeric(dr_tau2),
        diff = as.numeric(tau2_diff),
        abs_diff = as.numeric(tau2_abs_diff)
      ),
      comparison = comparison,
      cli_metadata = list(
        model = cli_results$model,
        n_studies = as.integer(cli_results$n_studies),
        estimation_mode = cli_results$estimation_mode
      ),
      message = "PASS requires joint agreement on beta, SE, and tau2."
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
    list(
      name = test_name,
      pass = FALSE,
      fit_ok = FALSE,
      beta_match = FALSE,
      se_match = FALSE,
      tau2_match = FALSE,
      error = conditionMessage(e)
    )
  })
}

run_all_tests <- function() {
  cat("============================================================\n")
  cat("Dose Response Pro v18.1 - R Validation Test Suite\n")
  cat("============================================================\n")

  test_definitions <- list(
    list(name = "Test 1: Simple Linear Trend", data = create_test_data_1()),
    list(name = "Test 2: Quadratic Trend", data = create_test_data_2()),
    list(name = "Test 3: High Heterogeneity", data = create_test_data_3())
  )

  test_results <- lapply(test_definitions, function(td) {
    validate_dosresmeta(td$name, td$data)
  })

  test_pass <- vapply(test_results, function(x) isTRUE(x$pass), logical(1))
  overall_pass <- all(test_pass)

  cat("\n=== VALIDATION SUMMARY ===\n")
  for (i in seq_along(test_definitions)) {
    label <- test_definitions[[i]]$name
    status <- ifelse(test_pass[i], "PASS", "FAIL")
    cat(sprintf("%s: %s\n", label, status))
  }
  cat(sprintf("Overall: %s\n", ifelse(overall_pass, "ALL TESTS PASSED", "SOME TESTS FAILED")))

  list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    r_version = R.version.string,
    package_versions = list(
      dosresmeta = as.character(packageVersion("dosresmeta")),
      metafor = as.character(packageVersion("metafor")),
      mvmeta = as.character(packageVersion("mvmeta"))
    ),
    overall_pass = overall_pass,
    summary = list(
      total_tests = length(test_results),
      passed_tests = sum(test_pass),
      failed_tests = sum(!test_pass)
    ),
    tolerances = list(
      beta = VALIDATION_BETA_TOLERANCE,
      se = VALIDATION_SE_TOLERANCE,
      tau2 = VALIDATION_TAU2_TOLERANCE
    ),
    tests = test_results
  )
}

render_validation_report <- function(results) {
  lines <- c(
    "Dose Response Pro v18.1 - R Validation Results",
    sprintf("Generated: %s", results$timestamp),
    sprintf("R: %s", results$r_version),
    sprintf("dosresmeta: %s | metafor: %s | mvmeta: %s",
      results$package_versions$dosresmeta,
      results$package_versions$metafor,
      results$package_versions$mvmeta
    ),
    sprintf("Overall: %s", ifelse(results$overall_pass, "PASS", "FAIL")),
    sprintf("Total tests: %d | Passed: %d | Failed: %d",
      results$summary$total_tests,
      results$summary$passed_tests,
      results$summary$failed_tests
    ),
    sprintf(
      "Tolerance(beta): %.3f | Tolerance(SE): %.3f | Tolerance(tau2): %.3f",
      results$tolerances$beta,
      results$tolerances$se,
      results$tolerances$tau2
    ),
    "",
    "Per-test details",
    paste(rep("=", 60), collapse = "")
  )

  for (test in results$tests) {
    status <- ifelse(isTRUE(test$pass), "PASS", "FAIL")
    lines <- c(lines, sprintf("[%s] %s", status, test$name))

    if (!is.null(test$error)) {
      lines <- c(lines, sprintf("  Error: %s", test$error), "")
      next
    }

    lines <- c(lines, sprintf("  Fit OK: %s", ifelse(isTRUE(test$fit_ok), "yes", "no")))
    lines <- c(lines, sprintf("  Beta match: %s", ifelse(isTRUE(test$beta_match), "yes", "no")))
    lines <- c(lines, sprintf("  SE match: %s", ifelse(isTRUE(test$se_match), "yes", "no")))
    lines <- c(lines, sprintf("  Tau2 match: %s", ifelse(isTRUE(test$tau2_match), "yes", "no")))

    cmp <- test$comparison
    for (i in seq_len(nrow(cmp))) {
      lines <- c(lines, sprintf(
        "  %s -> beta diff: %.6f | se diff: %.6f",
        cmp$parameter[i], cmp$beta_diff[i], cmp$se_diff[i]
      ))
    }

    if (!is.null(test$tau2)) {
      lines <- c(lines, sprintf(
        "  tau2(ours): %.6f | tau2(R): %.6f | diff: %.6f",
        test$tau2$our, test$tau2$r, test$tau2$diff
      ))
    }
    if (!is.null(test$cli_metadata)) {
      lines <- c(lines, sprintf(
        "  CLI model: %s | CLI studies: %d | mode: %s",
        test$cli_metadata$model,
        test$cli_metadata$n_studies,
        test$cli_metadata$estimation_mode
      ))
    }

    lines <- c(lines, "")
  }

  paste(lines, collapse = "\n")
}

write_validation_outputs <- function(results) {
  out_dir <- get_script_dir()
  json_path <- file.path(out_dir, "r_validation_results.json")
  txt_path <- file.path(out_dir, "r_validation_results.txt")

  json_text <- jsonlite::toJSON(results, pretty = TRUE, auto_unbox = TRUE, na = "null")
  writeLines(json_text, con = json_path, useBytes = TRUE)
  writeLines(render_validation_report(results), con = txt_path, useBytes = TRUE)

  list(json = json_path, txt = txt_path)
}

main <- function() {
  results <- run_all_tests()
  output_paths <- write_validation_outputs(results)

  cat(sprintf("\nJSON results saved: %s\n", output_paths$json))
  cat(sprintf("Text report saved: %s\n", output_paths$txt))

  if (isTRUE(results$overall_pass)) {
    quit(save = "no", status = 0)
  }

  quit(save = "no", status = 1)
}

if (sys.nframe() == 0) {
  main()
}
