#' load_bwtool_agg
#'
#' Loads a bwtool aggregates file.
#'
#' @param path path to aggregate files.
#' @param regex regular expression or list of regex, for selection of the files.
#' @param stranded search for 'plus' and 'minus' files on top of regex.
#' @param junk_list a list of substrings to be removed before splitting to levels.
#' @param split_to split filenames into levels
#' @param split_by delimiter for splitting, default='_'
#' @param plus_cnt count of features used for aggregate on the plus strand (default=1)
#' @param minus_cnt count of features used for aggregate on the minus strand (default=1)
#'
#'
#' @details Loads a bwtool aggregate file and collects this in a \code{tidy} data.frame.
#'
#' @return A tidy data frame.
#'
#' @examples
#'  path <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/bw/aggregates/exon1/'
#'  regex <- c('rep[3,4]', '_exon1_less200_', '.bed_starts_aggregate.txt$')
#'  junk_list <- c('_exon1_less200', '.bed_starts_aggregate.txt$', '_N20', '_3D12', '_plus', '_minus')
#'  df <- load_bwtool_agg(path, regex, stranded=TRUE, junk_list)
#'  ggplot(df, aes(x=rel_pos, y=value, color=siRNA, linetype=rep)) + geom_line() + facet_grid(.~ab)
#'
#'  #for scaling:
#'  dplyr::group_by_(df, .dots=lapply(split_to, as.symbol)) %>%
#'  mutate(value = scale(value)) %>%
#'  ggplot(., aes(x=rel_pos, y=value, color=siRNA, linetype=rep)) + geom_line() + facet_grid(.~ab)
#'
#' @export
load_bwtool_agg <- function(path, regex,
                            stranded = TRUE,
                            junk_list = c(),
                            split_to=c('siRNA', 'ab', 'rep'), split_by='_',
                            plus_cnt=1, minus_cnt=1){

  simple_name <- function(name) {
    sname <- name
    for(f in junk_list) {
      sname <- sub(f, '', sname)
    }
    sname
  }

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



#' load_bwtool_summary
#'
#' Loads a bwtool summary file.
#'
#' @param path path to summary files.
#' @param regex regular expression or list of regex, for selection of the files.
#' @param header logical. bwtool summary done with -header? Default=TRUE
#' @param keep_bed logical. bwtool summary done with -keep_bed? Default=FALSE
#' @param with_quantiles logical. bwtool summary done with -with-quantiles? Default=FALSE
#' @param with_sum logical. bwtool summary done with -with-sum? Default=FALSE
#' @param with_sum_of_squares logical. bwtool summary done with -with-sum-of-squares? Default=FALSE
#' @param junk_list a list of substrings to be removed before splitting to levels.
#' @param stranded search for 'plus' and 'minus' files on top of regex.
#' @param split_to split filenames into levels
#' @param split_by delimiter for splitting, default='_'
#'
#'
#' @details Loads a bwtool summary file and collects this in a \code{tidy} data.frame.
#'
#' @return A tidy data frame.
#'
#' @examples
#' path <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/bw/summaries/'
#' regex <- c('*rep[1,2]*_refseq_genes_summary.txt', 'FFL|Z18')
#' junk_list <- c('_refseq_genes_summary.txt', '_N20', '_3D12')
#' df <- load_bwtool_summary(path, regex, junk_list=junk_list, with_quantiles = T, with_sum = T, with_sum_of_squares = T)
#' ggplot(df, aes(x=rep, y=log2(sum), fill=rep)) + geom_violin() + facet_grid(siRNA~ab)
#'
#'
#' @export
load_bwtool_summary <- function(path,
                                regex,
                                header = TRUE,
                                keep_bed = FALSE,
                                with_quantiles = FALSE,
                                with_sum = FALSE,
                                with_sum_of_squares = FALSE,
                                stranded = FALSE,
                            junk_list = c(),
                            split_to=c('siRNA', 'ab', 'rep'),
                            split_by='_'){

  simple_name <- function(name) {
    sname <- name
    for(f in junk_list) {
      sname <- sub(f, '', sname)
    }
    sname
  }


  if (keep_bed) {
    colnames <- c('chr', 'start', 'end', 'name', 'size', 'num_data', 'min', 'max', 'mean')
  } else {
    colnames <- c('chr', 'start', 'end', 'size', 'num_data', 'min', 'max', 'mean')
  }
  if (with_quantiles) {
    colnames <- c(colnames, 'first_10%','1st_quart','median','3rd_quart','last_10%')
  }
  if (with_quantiles) {
    colnames <- c(colnames, 'sum_of_squares')
  }
  if (with_sum) {
    colnames <- c(colnames, 'sum')
  }

  reader <- function(flist) {
    lapply(seq_along(flist), function(i) {read.table(paste0(path,flist[i]), header=header, col.names=colnames) %>%
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

    df <- bind_rows(df_plus, df_minus)
  } else {
    df <- reader(files)
  }


  return(df)
}
