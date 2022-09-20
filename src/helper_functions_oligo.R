#' Check and install packages
#' 
#' Checks if packages are available and install them if necessary. Does not
#' work for BioConductor packages (e.g. Rsamtools).
#' 
#' @export
#' 
#' @param requiredPackages vector containing the required packages
#'
checkPackages <- function(requiredPackages){
  for(package in requiredPackages){
    if(!require(package, character.only = TRUE)){
      install.packages(package, dependencies = TRUE)
      require(package, character.only = TRUE)
    }
  }
}


#' Count deletions in cigar string
#' 
#' @importFrom stringr str_extract_all str_remove str_count
#' 
#' @param cigar CIGAR string
#' @export
processCigar <- function(cigar_list, pattern = '\\d+D'){
  
  out <- c()
  for(j in 1:length(cigar_list)){
    cigar = cigar_list[j]
    
    # Get elements with digits before the deletion character "D"
    # Example: 
    # cigar = "32M1D37M1D29M"
    # cstring = "13" "1"
    d <- stringr::str_extract_all(string = cigar, pattern = pattern)[[1]] %>%
      stringr::str_remove("D")
    
    n_deleted_bases <- sum(as.numeric(d))
    
    # Prepare output variables
    out <- append(out, n_deleted_bases)
  }
  
  out[is.na(out)] <- 0
  
  # Return variables
  return(out)
}

#' Plot substitution and indel error per base
#' 
#' @importFrom dplyr mutate filter group_by summarise
#' @importFrom tidyr unite
#' 
#' @export
#' 
#' @param data
#' @param insertType
#' @param batch
#' @param scales
#' @param hide.legend
#' @param threshold
#'
perBaseErrorPlot <- function(
  data,
  insertType = c('Insert_1', 'gDNA'),
  batch = 'b1',
  scales = 'free_x',
  hide.legend = FALSE,
  threshold = 10
  ){

  cons.variants <- data %>%
    dplyr::mutate(
      subsitution_error = 100*(A+C+G+T)/Coverage,
      indel_error = 100*(D+I)/Coverage,
      total_error =  100*(A+C+G+T+D+I)/Coverage
    ) %>%
    dplyr::filter(
      `Consensus group size` == threshold,
      InsertType %in% insertType,
      Batch %in% batch
    ) %>%
    tidyr::unite('Type', Manufacturer, Purification, sep = ' ') %>%
    dplyr::group_by(Position, Type, InsertType, Batch) %>%
    dplyr::summarise(
      `Substitutions (%)` = mean(subsitution_error),
      subs_std = sd(subsitution_error),
      `Deletions (%)` = mean(indel_error),
      indels_std = sd(indel_error)
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

  perBaseErrorPlot <- grid.arrange(indels_plot,subs_plot,ncol = 1,nrow = 2)

  return(perBaseErrorPlot)
}


#' Amplicon ribbon plot
#'
#' @import dplyr
#' @import ggplot2
#'
#' @param data filtered consensus table object
#' @param manufacturers vector containing manufacturer to use
#' @param purifications vector containing purifications to use
#' @param insertTypes vector containing insertTypes to use
#' @param start_pos First position to plot
#' @param end_pos Last position to plot
#'
#' @export
#' 
ribbonPlot <- function(data, manufacturers, purifications, insertTypes, start_pos, end_pos){
  
  cons.variants <- data %>%
    dplyr::mutate(
      subsitution_error = 100*(A+C+G+T)/Coverage,
      indel_error = 100*(D+I)/Coverage,
      total_error =  100*(A+C+G+T+D+I)/Coverage
    ) %>%
    dplyr::filter(
      `Consensus group size` == consensus_cutoff,
      Manufacturer %in% manufacturers,
      Purification %in% purifications,
      InsertType %in% insertTypes,
      Batch %in% c('b1','b2', 'b3'),
      Position >= start_pos, 
      Position <= end_pos
    ) %>%
    dplyr::group_by(Position, Manufacturer,Purification, InsertType, Batch, Reference) %>%
    dplyr::summarise(
      `Consensus 10 error (%)` = mean(total_error),
      sd = sd(total_error)
    )
  
  cons.variants <- drop_na(cons.variants)
  
  labels <- cons.variants %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Position, Reference)
  
  amplicon_plot <- ggplot(
    data = cons.variants, 
    mapping = aes(
      x = Position,
      y = `Consensus 10 error (%)`,
      fill = Batch
    )
  ) + 
    geom_line(mapping = aes(col=Batch)) +
    geom_ribbon(
      mapping = aes(
        ymin=`Consensus 10 error (%)`-sd,
        ymax=`Consensus 10 error (%)`+sd
      ), alpha=0.5
    ) +
    theme_minimal() + 
    expand_limits(y = 1.5) +
    facet_wrap(Purification ~ Manufacturer,scales = 'free', ncol = 1,strip.position = 'right') +
    scale_x_continuous(breaks=labels$Position, labels=labels$Reference) +
    scale_fill_manual(values=c("gray55", '#f6be4c', '#469ed4')) +
    scale_color_manual(values=c('gray55', '#f6be4c', '#469ed4')) + 
    theme(
      axis.text.x = element_text(size = 6), 
      panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = 'black'),
      axis.ticks.x = element_line(colour = 'black'),
      axis.ticks.y = element_line(colour = 'black')
    ) +
    geom_hline(
      yintercept = 0,
      colour = 'gray55',
      linetype = 'dashed'
    ) 
  
  return(amplicon_plot)
}
