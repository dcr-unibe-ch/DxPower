# print definition for class
#' @export
#' @importFrom crayon bold
print.dx_bin <- function(x, ...) {

  cat("\nDxPower - Confirmatory Diagnostic Accuracy Study \n")
  cat("--------------------------------------------------\n")
  cat("Framework: Co-primary",x$hypothesis, "testing (IUT)\n")
  cat("Alpha level:", x$alpha, "\n")
  cat("Alternative:", x$alternative, "\n")
  cat("Prevalence:", x$prev, "\n")
  cat("CI method:", x$ci_method, "\n\n")

  if (x$type == "samplesize") {

    cat(bold("Required sample size:\n"))
    cat("  Total:        ", bold(ceiling(x$N_tot)), "\n")
    cat("  Diseased:     ", bold(ceiling(x$N_se)), "\n")
    cat("  Non-diseased: ", bold(ceiling(x$N_sp)), "\n\n")

    cat("Target overall power: ",
        round(x$target_power, 3), "\n\n")

    cat("Allocated powers:\n")
    cat("  Sensitivity: ", round(x$pwr_se, 3), "\n")
    cat("  Specificity: ", round(x$pwr_sp, 3), "\n\n")



  } else if (x$type == "power") {

    cat(bold("Achieved overall power: ",
        round(x$overall, 3), "\n\n"))

    cat("Endpoint-specific power:\n")
    cat("  Sensitivity: ", round(x$pwr_se, 3), "\n")
    cat("  Specificity: ", round(x$pwr_sp, 3), "\n\n")

    cat("Sample size:\n")
    cat("  Total:        ", ceiling(x$N_tot), "\n")
    cat("  Diseased:     ", ceiling(x$N_se), "\n")
    cat("  Non-diseased: ", ceiling(x$N_sp), "\n\n")
  }

  invisible(x)
}
