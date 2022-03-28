





figure_2_plot <- function(
  data,
  insertType = c('Insert_1', 'gDNA'),
  batch = 'b1',
  scales = 'free_x',
  hide.legend = FALSE
  ){

  cons.variants <- data %>%
    dplyr::mutate(
      subsitutions_raw = 100*(rawA+rawC+rawG+rawT)/rawDepth,
      subsitutions_cons = 100*(consA+consC+consG+consT)/consDepth,
      indels_raw = 100*(rawD+rawI)/rawDepth,
      indels_cons = 100*(consD+consI)/consDepth
    ) %>%
    dplyr::filter(
      InsertType %in% insertType,
      Batch %in% batch
    ) %>%
    tidyr::unite('Type', Manufacturer, Purification, sep = ' ') %>%
    dplyr::group_by(Position, Type, InsertType, Batch) %>%
    dplyr::summarise(
      `Substitutions (%)` = mean(subsitutions_cons),
      subs_std = sd(subsitutions_cons),
      `Deletions (%)` = mean(indels_cons),
      indels_std = sd(indels_cons)
    )

  cons.variants$Type <- forcats::fct_relevel(cons.variants$Type, 'control MCF7', after = Inf)

  indels_plot <- ggplot(
    data = cons.variants,
    mapping = aes(
      x = Position,
      y = `Deletions (%)`,
      fill = Type)
  ) +
    geom_bar(stat = 'identity') +
    geom_errorbar(
      mapping = aes(
        ymin=`Deletions (%)`,
        ymax=`Deletions (%)`+ indels_std
      ), size = .1, width=.4
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank(),
      legend.position = 'bottom',
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks.y = element_line()
    ) +
    guides(x = 'none') +
    facet_grid(Batch ~ Type, scales = scales, space = "free_x")
    #geom_hline(yintercept=0, linetype='solid', col = 'gray40', width = 0.1)

  subs_plot <- ggplot(
    data = cons.variants,
    mapping = aes(
      x = Position,
      y = `Substitutions (%)`,
      fill = Type)
  ) +
    geom_bar(stat = 'identity') +
    geom_errorbar(
      mapping = aes(
        ymin=`Substitutions (%)`,
        ymax=`Substitutions (%)`+ subs_std
      ), size = .1, width=0.4 ) +
    facet_grid(Batch ~ Type, scales = scales, space = "free_x") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank(),
      legend.position = 'bottom',
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks.y = element_line()
    ) +
    guides(x = 'none') +
    scale_y_reverse()
    #geom_hline(yintercept=0, linetype='solid', col = 'gray40', width = 0.1)

  if(hide.legend){
    indels_plot <- indels_plot + theme(legend.position = "none")
    subs_plot <- subs_plot + theme(legend.position = "none")
  }

  figure_2 <- grid.arrange(indels_plot,subs_plot,ncol = 1,nrow = 2)

  return(figure_2)
}
