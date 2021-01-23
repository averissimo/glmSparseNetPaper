show.relative.risk.distribution <- function(fit.df, ix) {
  fit.df %>% as.tibble() %>%
    mutate(type.str = prettify.labels(type),
           type.model = prettify.labels(type, ret.value = 'model'),
           type.base  = paste0('base model: ', prettify.labels(type, 'string.base'))) %>%
    filter(set == ix) %>%
    ggplot() +
    geom_density(aes(relative.risk, color = type.model, fill = type.model), alpha = .1, adjust = 1/5) +
    theme_minimal() +
    #geom_vline(aes(xintercept = mean), linetype = 'dotted', color = '#999999') +
    #geom_text(aes(x = mean, y = -.5, label = sprintf(' Mean of relative risk %g', mean), hjust = 0), color = '#999999', inherit.aes = FALSE, check_overlap = TRUE) +
    facet_wrap( ~ type.base + set, ncol = 1) +
    scale_color_discrete('') +
    scale_fill_discrete('') +
    xlab('Relative risks') +
    ylab('Density') +
    ggtitle('Distribution of risk in Classic model -- {ix} set' %>% glue::glue()) +
    theme(legend.position = 'bottom') +
    xlim(0, 5)
}
