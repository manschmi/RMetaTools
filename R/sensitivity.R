#' Sensitivity
#'
#' Computes Sensitivity for Metagene Plots.
#'
#' @param df tidy data frame of values and conditions.
#' @param value_name column name of the column containing the values
#' @param ctrl_factor name of the factor containing the control
#' @param ctrl_name name of the control condition
#' @param pseudocount pseudocount value
#' @param pseudocount_method method used for pseudocounting (default is 'fill_with_min')
#'
#'
#' @details Loads a \code{tidy} data.frame. Containing 1 column for each factor and all values in 1 column. Need to provide column name of the column containing the control (ie 'siRNA') and the name of the control (ie 'siFFL'). Computes the sensitivity: sensitivity = (expression(kd) + pseudocount) / (expression(ctrl) + pseudocount)$.
#' Pseudocounting: Either adds pseudocount to all values, or if pseudocount_method == "fill_with_min" uses the min of all values for pseudocounting.
#'
#' @return A tidy data frame with all factors but containing sensitivity instead of values.
#'
#' @examples
#'  path_to_RNAseq <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/'
#'  m <- load_deeptoolsmatrix(paste0(path_to_RNAseq, 'snR_TES_joined.gz'), na.omit = FALSE, na.fill = 0) %>%  separate(sample_name, c('x', 'sname' ), sep='_') %>% mutate(sname = sub('rep-2plus', 'repplus', sname)) %>% separate(sname, c('siRNA', 'f', 'rep', 'bedstrand' ), sep='-') %>% dplyr::select(id, rel_pos, siRNA, value)
#'  head(m)
#'  sensitivity(m, ctrl_factor = 'siRNA', ctrl_value = 'siFFL',)
#'
#'  #for scaling:
#'  dplyr::group_by_(df, .dots=lapply(split_to, as.symbol)) %>%
#'  mutate(value = scale(value)) %>%
#'  ggplot(., aes(x=rel_pos, y=value, color=siRNA, linetype=rep)) + geom_line() + facet_grid(.~ab)
#'
#' @export
sensitivity <- function(df, value_name = 'value',
                        ctrl_factor, ctrl_name,
                        pseudocount=NA, pseudocounting_method='fill_with_min'){

  if (pseudocounting_method == 'fill_with_min') {
    min_val = min(df[,value_name])
    df[,value_name] <- df[,value_name] + min_val
  }

  ctrl <- df[df[,ctrl_factor]==ctrl_name,]
  kd <- df[df[,ctrl_factor]!=ctrl_name,]




  reader <- function(flist) {
    lapply(seq_along(flist), function(i) {read.table(paste0(path,flist[i]), header=F, col.names=c('rel_pos', 'value')) %>%
      mutate(sample = simple_name(flist[i])) %>%
      separate(sample, split_to, split_by)}) %>% bind_rows(.)
  }

  files <- dir(path)
  for(re in regex){
    files <- files[grep(re, files)]
  }

  if(stranded){
    plus_files <- files[grep('plus', files)]
    df_plus <- reader(plus_files)

    minus_files <- files[grep('minus', files)]
    df_minus <- reader(minus_files)

    total_cnt = plus_cnt + minus_cnt
    df <- dplyr::mutate(df_minus, rel_pos = -rel_pos) %>%
      dplyr::rename(minus_value=value) %>%
      dplyr::left_join(df_plus, .) %>%
      dplyr::mutate(value = (value*plus_cnt + minus_value*minus_cnt)/(total_cnt)) %>%
      dplyr::select(-minus_value)
  } else {
    df <- reader(files)
  }

  return(df)
}
