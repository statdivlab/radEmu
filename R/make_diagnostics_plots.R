#' makes plots to investigate convergence of estimation under null hypothesis
#' 
#' @param diagnostic_df Dataframe with relevant information
#' 
#' @return Plots with diagnostic information. If the fisher scoring algorithm has been fit, this will be two
#' plots showing how the log likelihood and test statistic change over iterations. If the augmented lagrangian
#' algorithm has been fit, this will also include how the constraint gap and maximum changing element of B change
#' over iterations.
#' 
make_diagnostics_plots <- function(diagnostic_df) {
  p1 <- ggplot2::ggplot(diagnostic_df, ggplot2::aes(x = it, y = lik)) + 
    ggplot2::geom_line() + 
    ggplot2::geom_point() + 
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Iteration", y = "Log likelilhood")
  p2 <- ggplot2::ggplot(diagnostic_df, ggplot2::aes(x = it, y = test_stat)) + 
    ggplot2::geom_line() +
    ggplot2::geom_point() + 
    ggplot2::theme_bw() +
    ggplot2::ylim(0, NA) + 
    ggplot2::labs(x = "Iteration", y = "Test statistic")
  plots <- list(p1 = p1, p2 = p2)
  if ("max_abs_B" %in% names(diagnostic_df)) {
    p3 <- ggplot2::ggplot(diagnostic_df, ggplot2::aes(x = it, y = constraint_diff)) + 
      ggplot2::geom_line() + 
      ggplot2::geom_point() + 
      ggplot2::theme_bw() +
      ggplot2::ylim(0, NA) + 
      ggplot2::labs(x = "Iteration", y = "Constraint difference")
    p4 <- ggplot2::ggplot(diagnostic_df, ggplot2::aes(x = it, y = max_abs_B)) + 
      ggplot2::geom_line() + 
      ggplot2::geom_point() + 
      ggplot2::theme_bw() +
      ggplot2::ylim(0, NA) + 
      ggplot2::labs(x = "Iteration", y = "Max change to B")
    plots$p3 <- p3
    plots$p4 <- p4
  }
  return(plots)
}