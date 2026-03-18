# print definition for class
#' @export
#' @importFrom crayon bold
print.dx_bin <- function(x, ...) {

  cat("\nDxPower - Confirmatory Diagnostic Accuracy Study \n")
  cat("--------------------------------------------------\n")

  cat("Design:           ", x$design, "\n")
  cat("Framework:        Co-primary (IUT)\n")
  # Alpha line
  if (x$design == "single-test design") {
    alpha_line <- paste0(x$alpha, " (", x$alternative, ")")
  } else {
    alpha_line <- as.character(x$alpha)
  }

  cat("Alpha level:      ", alpha_line, "\n")
  cat("Prevalence:       ", x$prev, "\n\n")

  if (!is.null(x$hypothesis)) {
    cat("Hypotheses:\n")
    cat("  Sensitivity :   ", x$hypothesis$sensitivity, "\n")
    cat("  Specificity :   ", x$hypothesis$specificity, "\n\n")
  }

  if (x$analysis == "samplesize") {

    cat(bold("Required sample size per group:\n"))
    cat("  Total:          ", bold(ceiling(x$N_tot)), "\n")
    cat("  Diseased:       ", bold(ceiling(x$N_se)), "\n")
    cat("  Non-diseased:   ", bold(ceiling(x$N_sp)), "\n\n")

    cat("Target overall power:  ", round(x$target_power, 3), "\n\n")

    cat("Allocated powers:\n")
    cat("  Sensitivity:     ", round(x$power_se, 3), "\n")
    cat("  Specificity:     ", round(x$power_sp, 3), "\n\n")

  } else if (x$analysis == "power") {

    cat(bold("Achieved overall power: "),
        round(x$overall_power, 3), "\n\n")

    cat("Endpoint-specific power:\n")
    cat("  Sensitivity:     ", round(x$power_se, 3), "\n")
    cat("  Specificity:     ", round(x$power_sp, 3), "\n\n")

    cat("Sample size per group:\n")
    cat("  Total:           ", ceiling(x$N_tot), "\n")
    cat("  Diseased:        ", ceiling(x$N_se), "\n")
    cat("  Non-diseased:    ", ceiling(x$N_sp), "\n\n")
  }

  invisible(x)
}
