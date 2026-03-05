#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  required_packages <- c("dosresmeta", "jsonlite")
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(sprintf("Missing required packages: %s", paste(missing, collapse = ", ")))
  }
  library(dosresmeta)
  library(jsonlite)
})

MIN_POSITIVE <- 1e-12

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript strict_r_batch_dosresmeta.R <input_json> <output_json>")
}

input_json <- args[[1]]
output_json <- args[[2]]

payload <- fromJSON(input_json, simplifyDataFrame = FALSE)
datasets <- payload$datasets

prepare_dosresmeta_data <- function(rows_df) {
  clean <- rows_df[complete.cases(rows_df), c("study", "dose", "cases", "n")]
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

fit_one <- function(dataset) {
  dataset_id <- dataset$dataset_id
  rows <- dataset$rows

  out <- list(
    dataset_id = dataset_id,
    ok = FALSE,
    runtime_ms = NA_real_,
    beta_linear = NA_real_,
    beta_quadratic = NA_real_,
    se_linear = NA_real_,
    se_quadratic = NA_real_,
    tau2 = NA_real_,
    n_studies = NA_integer_,
    n_rows = NA_integer_,
    error = NULL
  )

  tryCatch({
    if (is.data.frame(rows)) {
      rows_df <- rows
    } else {
      rows_df <- do.call(
        rbind,
        lapply(rows, function(r) {
          data.frame(
            study = as.character(r$study),
            dose = as.numeric(r$dose),
            cases = as.numeric(r$cases),
            n = as.numeric(r$n),
            stringsAsFactors = FALSE
          )
        })
      )
    }
    rows_df$dose <- as.numeric(rows_df$dose)
    rows_df$cases <- as.numeric(rows_df$cases)
    rows_df$n <- as.numeric(rows_df$n)
    rows_df <- rows_df[is.finite(rows_df$dose) & is.finite(rows_df$cases) & is.finite(rows_df$n), ]
    rows_df <- rows_df[rows_df$cases > 0 & rows_df$n > 0, ]

    start <- proc.time()[["elapsed"]]
    dr_data <- prepare_dosresmeta_data(rows_df)
    fit <- dosresmeta::dosresmeta(
      formula = logrr ~ dose + I(dose^2),
      id = study,
      type = type,
      se = se,
      cases = cases,
      n = n,
      data = dr_data
    )
    elapsed_ms <- (proc.time()[["elapsed"]] - start) * 1000

    beta <- as.numeric(stats::coef(fit))
    se <- as.numeric(sqrt(diag(stats::vcov(fit))))
    psi <- tryCatch(as.matrix(summary(fit)$Psi), error = function(e) matrix(NA_real_, nrow = 1, ncol = 1))
    tau2 <- if (all(!is.finite(diag(psi)))) NA_real_ else mean(diag(psi), na.rm = TRUE)

    out$ok <- TRUE
    out$runtime_ms <- as.numeric(elapsed_ms)
    out$beta_linear <- if (length(beta) >= 1) beta[[1]] else NA_real_
    out$beta_quadratic <- if (length(beta) >= 2) beta[[2]] else NA_real_
    out$se_linear <- if (length(se) >= 1) se[[1]] else NA_real_
    out$se_quadratic <- if (length(se) >= 2) se[[2]] else NA_real_
    out$tau2 <- as.numeric(tau2)
    out$n_studies <- length(unique(dr_data$study))
    out$n_rows <- nrow(dr_data)
    out
  }, error = function(e) {
    out$error <- conditionMessage(e)
    out
  })
}

results <- lapply(datasets, fit_one)

response <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  package_versions = list(
    dosresmeta = as.character(packageVersion("dosresmeta"))
  ),
  input_count = length(datasets),
  results = results
)

json_text <- jsonlite::toJSON(response, pretty = TRUE, auto_unbox = TRUE, na = "null")
writeLines(json_text, con = output_json, useBytes = TRUE)
