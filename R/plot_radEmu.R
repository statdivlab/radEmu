#' Plotting function
#'
#' @param x Output from emuFit()
#' @param plot_key (Optional) Default \code{NULL}. List of named vectors containing names in the "covariate" column of the `coef` output of the radEmu model object. If you wish for multiple covariate values to be plotted on the same plot, then those variables should be included in the same named vector. By default, each column of the design matrix receives its own plot.
#' @param title (Optional). Default \code{NULL}. Character string. The main title for the graphic.
#' @param taxon_names (Optional). Default \code{NULL}. Data frame. If \code{NULL}, keep taxon names as listed in radEmu model. Otherwise, users can input a data frame with two columns: one labelled "category" with the same levels as in the radEmu output and another labelled "cat_small" with the preferred labels.
#' @param display_taxon_names (Optional). Default \code{TRUE}. Boolean. If \code{FALSE}, remove sample names from the plot.
#' @param data_only (Optional). Default \code{FALSE}. Boolean. If \code{TRUE}, only returns data frame.
#' @param ... There are no optional parameters at this time.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return Object of class \code{ggplot}. Plot of \code{radEmu} model fit with 95% confidence intervals.
#'
#' @examples
#' data(wirbel_sample)
#' data(wirbel_otu)
#' 
#' subset_studies <- which(wirbel_sample$Study %in% c("FR-CRC", "US-CRC", "AT-CRC"))
#' 
#' chosen_genera <- c("Eubacterium", "Faecalibacterium", "Fusobacterium", "Porphyromonas")
#' 
#' mOTU_names <- colnames(wirbel_otu)
#' mOTU_name_df <- data.frame(name = mOTU_names) %>% 
#'   dplyr::mutate(base_name = stringr::str_remove(mOTU_names, "unknown ") %>%
#'   stringr::str_remove("uncultured ")) %>%
#'   dplyr::mutate(genus_name = stringr::word(base_name, 1))
#' 
#' restricted_mOTU_names <- mOTU_name_df %>%
#'   dplyr::filter(genus_name %in% chosen_genera) %>%
#'   dplyr::pull(name)
#' 
#' small_Y <- wirbel_otu[subset_studies, restricted_mOTU_names]
#' category_to_rm <- which(colSums(small_Y) == 0)
#' 
#' small_sample <- wirbel_sample[subset_studies, ]
#' 
#' ch_fit <- emuFit(formula = ~ Group + Study, 
#'                  data = small_sample,
#'                  Y = small_Y,
#'                  run_score_tests = FALSE)
#' 
#' plot_key <- list(p1 = c("Control" = "GroupCTR"),
#'                  p2 = c("FR-Control" = "StudyFR-CRC",
#'                         "US-Control" = "StudyUS-CRC"))
#' 
#' out <- plot.radEmu(x = ch_fit,
#'                    plot_key = plot_key,
#'                    display_taxon_names = FALSE)
#' 
#' out$plots$p1
#' out$plots$p2
#' @export

plot.radEmu <- function(x,
                        plot_key = NULL,
                        title = NULL,
                        taxon_names = NULL,
                        display_taxon_names = TRUE,
                        data_only = FALSE, ...) {
  mod <- x
  
  # confirm that values in plot_key are actually in the coefficient table
  if (!all(unlist(plot_key) %in% unique(mod$coef$covariate))) {
    stop("At least one of the coefficient names included in the plot key is not in
         the covariate column of the coefficient table of the radEmu() output.")
  }
  
  # confirm that no coefficient is listed multiple times
  if (length(unlist(plot_key)) != length(unique(unlist(plot_key)))) {
    stop("One of the coefficient names is included in the plot key multiple times,
         however each coefficient can only be included in one plot.")
  }
  
  # determine which levels are not included in the plot key
  remaining_variables <- setdiff(unique(mod$coef$covariate), unlist(plot_key))

  # construct a list of variables that should be plotted together
  plot_key_default <- append(plot_key,
                             lapply(as.list(remaining_variables),
                                    function(vec){
                                      setNames(vec, vec)
                                    }))
  
  plot_key_default <- lapply(plot_key_default, function(vec){
    if (is.null(names(vec))) {
      setNames(vec, vec)
    } else {
      setNames(vec, ifelse(names(vec) != "", names(vec), vec))
    }
    
  })
  
  # for each covariate, determine which it should be included in
  covariate_groups <- sapply(plot_key_default, function(key){
    mod$coef$covariate %in% key
    })
  
  # add a variable to identify which plot each coefficient will be plotted on
  mod$coef$plot_key <- apply(covariate_groups, 1, which)
  
  # rename variable levels as given in the plot key
  new_variable_names <- names(unlist(unname(plot_key_default)))
  old_variable_names <- unname(unlist(unname(plot_key_default)))
  mod$coef <- mod$coef %>%
    dplyr::mutate(covariate = factor(covariate,
                                     levels = old_variable_names,
                                     labels = new_variable_names))
  
  # now separate output coefficients into a list for each plot
  coef_list <- mod$coef %>%
    split(.$plot_key)
  
  # match coefficient names using user-provided "taxon_names" data frame
  coef_list_renamed <- lapply(coef_list, function(coef_subset){
    if (!is.null(taxon_names)) {
      dplyr::left_join(coef_subset, taxon_names, by = "category") %>%
        dplyr::arrange(covariate, estimate) %>%
        dplyr::mutate(cat_small = reorder(cat_small, estimate))
    } else {
      coef_subset %>%
        dplyr::mutate(cat_small = category) %>%
        dplyr::arrange(covariate, estimate) %>%
        dplyr::mutate(cat_small = factor(cat_small,
                                         levels = unique(cat_small),
                                         labels = unique(cat_small)))
    }
  })
  
  # return plots
  if (!data_only) {
    
    p_list <- lapply(coef_list_renamed, function(coef_subset){
      p <- ggplot2::ggplot(coef_subset) +
        ggplot2::geom_point(ggplot2::aes(x = estimate,
                                         y = as.character(cat_small),
                                         color = covariate,
                                         group = covariate),
                            position = ggplot2::position_dodge2(width = 0.5),
                            size = 2) +
        ggplot2::geom_errorbar(ggplot2::aes(y = cat_small,
                                            xmin = lower,
                                            xmax = upper,
                                            color = covariate,
                                            group = covariate),
                               position = ggplot2::position_dodge2(width = 0.5),
                               width = 0.5) +
        ggplot2::geom_vline(xintercept = 0, alpha = 0.5) +
        ggplot2::theme_bw() +
        ggplot2::labs(title = title) +
        ggplot2::guides(color = ggplot2::guide_legend(title = "Comparison")) +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(y = "Category") + 
        ggplot2::labs(x = "Estimate")
      
      if (!display_taxon_names) {
        p <- p +
          ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                         axis.ticks.y = ggplot2::element_blank())
      }
      
      p
    })
    
    p_list <- setNames(p_list, paste0("p", 1:length(p_list)))
    
    return(list(plots = invisible(p_list),
                data = dplyr::bind_rows(coef_list_renamed)))
    
  } else {
    
    # return data
    return(list(plots = NULL,
                data = dplyr::bind_rows(coef_list_renamed)))
  }
}


