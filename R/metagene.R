
#' Coerce metagene object to tibble
#'
#' Collects info from metagene object into a single tibble.
#'
#' @param m a metagene object
#' @param ... other arguments, currently ignored
#'
#' @details
#'
#' @return
#' a tibble with columns, rel_pos, value, sample, group, and name, where name is the id column of the metagene object regions, usually col4 of the bed file used to create the object.
#'
#' @examples
#'  fname <- system.file("extdata", "matrix_grouped.gz", package = "RMetaTools")
#'  m <- load_deeptoolsmatrix3_new(fname)
#'  to_tibble(m)
#' @export
to_tibble <- function(m, ...){
  lapply(seq_along(m$matrices), function(i) {
    as_tibble(m$matrices[[i]]) %>%
      set_colnames(m$rel_pos) %>%
      gather(rel_pos, value) %>%
      mutate(sample = m$samples[[i]],
             group = m$groups[[i]],
             name = rep(m$regions[[m$groups[[i]]]]$id,
                        length(m$rel_pos))
      )
  }) %>%
    bind_rows %>%
    dplyr::select(sample, group, name, rel_pos, value)
}





#' Select and rename metagene object sample or group names
#'
#' Select and/or rename sample and/or group names
#'
#' @param m a metagene object
#' @param select_samples regex or list of regex(es) for samples to select
#' @param remove_samples regex or list of regex(es) for samples to remove
#' @param select_groups regex or list of regex(es) for groups to remove
#' @param remove_groups regex or list of regex(es) for groups to remove
#' @param rename_samples a regexes and substitute string, default is substitute is empty
#' @param rename_groups a regexes and substitute string, default is substitute is empty
#' @param verbose report changes to stdout
#'
#' @details
#'
#' @return
#' a tibble with columns, rel_pos, value, sample, group, and name, where name is the id column of the metagene object regions, usually col4 of the bed file used to create the object.
#'
#' @examples
#'  fname <- system.file("extdata", "matrix_grouped.gz", package = "RMetaTools")
#'  m <- load_deeptoolsmatrix3_new(fname)
#'
#' @export
modify.metagene <- function(m,
                            select_samples = NULL,
                            remove_samples = NULL,
                            select_groups = NULL,
                            remove_groups = NULL,
                            rename_samples = NULL,
                            rename_groups = NULL,
                            verbose = FALSE){
  sel <- rep(TRUE, length(m$matrices))
  if (!missing(select_samples)) {
    sel <- sel & grepl(paste(select_samples, collapse='|'), m$samples)
  }

  if (!missing(remove_samples)) {
    sel <- sel & !grepl(paste(remove_samples, collapse='|'), m$samples)
  }

  if (!missing(select_groups)) {
    sel <- sel & grepl(paste(select_groups, collapse='|'), m$groups)
  }

  if (!missing(remove_groups)) {
    sel <- sel & !grepl(paste(remove_groups, collapse='|'), m$groups)
  }

  if(verbose){
    print(paste0('selecting ', sum(sel), ' of ', length(sel), ' matrices.'))
  }

  out <- m
  out$matrices <- out$matrices[sel]
  out$samples <- out$samples[sel]
  out$groups <- out$groups[sel]
  out$rel_pos <- out$rel_pos[sel]
  out$regions <- out$regions[out$groups]
  out$header <- lapply(out$header, function(h)
    {
      if(length(h) == length(sel)){
        h <- h[sel]
      }
    h
    })

  sel_idx <- which(sel)
  sample_sizes <- as.numeric(out$header$sample_boundaries[sel_idx+1])-as.numeric(out$header$sample_boundaries[sel_idx])
  out$header$sample_boundaries <- c(0, cumsum(sample_sizes))

  group_sizes <- as.numeric(out$header$group_boundaries[sel_idx+1])-as.numeric(out$header$group_boundaries[sel_idx])
  out$header$group_boundaries <- c(0, cumsum(group_sizes))

  if (!missing(rename_samples)) {
    if(verbose){
      print(paste0('sample names: ', paste(unique(out$samples))))
    }

    if(length(rename_samples) == 1){
      rename_samples <- c(rename_samples, '')
    }

    out$samples %<>% sub(rename_samples[1], rename_samples[2], .)

    if(verbose){
      print(paste0(' changed to: ', paste(unique(out$samples))))
    }
  }

  if (!missing(rename_groups)) {
    if(verbose){
      print(paste0('group names: ', paste(unique(out$groups))))
    }

    if(length(rename_groups) == 1){
      rename_groups <- c(rename_groups, '')
    }

    out$groups %<>% sub(rename_samples[1], rename_samples[2], .)
    names(out$rel_pos) %<>% sub(rename_samples[1], rename_samples[2], .)
    names(out$regions) %<>% sub(rename_samples[1], rename_samples[2], .)

    if(verbose){
      print(paste0(' changed to: ', paste(unique(out$groups))))
    }
  }

  out

}



#' Metagene Aggregate
#'
#' Metagene aggregate from an S3 metagene object.
#'
#' @param m the object to create the metagene aggregate for
#' @param aggregate_fun the function name to aggregate values such as 'events', 'sum', 'mean', 'ttest' (default is 'mean')
#' @param transform_fun a function to be applied to vlues before aggregate function is applied.
#' @param ... other arguments passed to fun
#'
#' @details Aggregates values.
#' transform_fun takes a matrix as input and returns a similar sized matrix. ie log2p1 <- function(m){log2(m+1)}
#' aggreagate_fun can be mean, sum, events or t.test or vector containing several of those.
#'
#' @return
#' a tibble with column, rel_pos, sample, group, and aggregate values where column name is mean, sum, events or estimate+conf.low+conf.high for ttest
#' @examples
#'  fname <- system.file("extdata", "matrix_grouped.gz", package = "RMetaTools")
#'  m <- load_deeptoolsmatrix3_new(fname)
#'  metagene_aggregate(m)
#'  metagene_aggregate(m, aggregate_fun = 'events')
#'  metagene_aggregate(m, aggregate_fun = 'ttest')
#'  log2p1 <- function(m){log2(m+1)}
#'  metagene_aggregate(m, aggregate_fun = 'mean', transform_fun = log2p1)
#'
#' @export
metagene_aggregate <- function(m, aggregate_fun = 'mean', transform_fun = NULL, ... ) {

  if (missing(transform_fun)) {
    mat_list <- m$matrices
  }else{
    mat_list <- lapply(m$matrices, transform_fun)
  }

  if( 'mean' %in% aggregate_fun ){
    agg_fun <- function(m){data.frame(mean = colMeans(m))}
  }else if( 'sum' %in% aggregate_fun ){
    agg_fun <- function(m){data.frame(mean = colSums(m))}
  }else if( 'events' %in% aggregate_fun ){
    agg_fun <- function(m){
      data.frame(events = apply(m,2,function(mi)sum(mi > 0)))
    }
  }else if( aggregate_fun == 'ttest' ){
    agg_fun <- function(m){
      meta_ttests <- apply(m, 2, function(mi) t.test(mi))
      tibble(
        mean =  lapply(meta_ttests, function(tti) tti$estimate) %>% simplify,
        conf.low =  lapply(meta_ttests, function(tti) tti$conf.int[1]) %>% simplify,
        conf.high = lapply(meta_ttests, function(tti) tti$conf.int[2]) %>% simplify
      )
    }
  }

  lapply(seq_along(mat_list), function(i) {
    tibble(
      rel_pos = m$rel_pos,
      sample = m$samples[[i]],
      group = m$groups[[i]]) %>%
      bind_cols(., agg_fun(mat_list[[i]]))
  }) %>%
    bind_rows
}


#' Metagene Ratio
#'
#' Computes Ratios between matrices to one control
#'
#' @param m the object to create the metagene aggregate for
#' @param samples sample names for one or more matrices from m to be used
#' @param control sample name for one of the matrices from m to be used as denominator
#' @param transform_fun a function to be applied to vlues before aggregate function is applied.
#' @param ... other arguments passed to fun
#'
#' @details Ratio of values.
#' transform_fun takes a matrix as input and returns a similar sized matrix. ie p1 <- function(m){m+1}
#' Computes the ratio to the control matrix.
#' Handles all groups if there is more than one group present.
#'
#' @return
#' a metagene object
#'
#' @examples
#'  fname <- system.file("extdata", "matrix_grouped.gz", package = "RMetaTools")
#'  m <- load_deeptoolsmatrix3_new(fname)
#'  p1 <- function(m){return(m+1)}
#'  ratios <- metagene_ratio(m, c('RRP40', 'ZCCHC8', 'ZFC3H1'), 'MTR4', transform_fun=p1)
#' @export
metagene_ratio <- function(m, samples, control, transform_fun = NULL, ... ) {
  ctrl_idx <- which(m$samples == control)
  samples_idx <- which(m$samples %in% samples)

  ctrl_mats  <- m$matrices[ctrl_idx]
  ctrl_groups <- m$groups[ctrl_idx]
  samples_mat <- m$matrices[samples_idx]
  samples_groups <- m$groups[samples_idx]
  sample_names <-

  if (!missing(transform_fun)) {
    ctrl_mats <- lapply(ctrl_mats, transform_fun)
    samples_mat <- lapply(samples_mat, transform_fun)
  }

  ratios_mat <- lapply(seq_along(samples_mat), function(i) samples_mat[[i]]/ctrl_mats[[which(ctrl_groups == samples_groups[i])]])

  names(ratios_mat) <- paste0(m$samples[samples_idx], '_rel_', control, '_', m$groups[samples_idx])

  s3 <- structure(list(), class = "metagene")
  s3$header <- m$header

  s3$matrices <- ratios_mat

  s3$samples <- paste0(m$samples[samples_idx], '_rel_', control)
  s3$groups <- m$groups[samples_idx]

  s3$rel_pos <- m$rel_pos

  s3$regions <- m$regions

  s3
}

# agg_fun <- function(mat){
#   apply(mat, 2, function(mi) tidy(t.test(mi))) %>%
#     bind_rows %>%
#     dplyr::select(estimate, conf.low, conf.high) %>%
#     dplyr::rename(mean=estimate)
# }
#   # tibble(
#   #   mean =  lapply(meta_ttests, function(tti) tti$estimate) %>% simplify,
#   #   conf.low =  lapply(meta_ttests, function(tti) tti$conf.int[1]) %>% simplify,
#   #   conf.high = lapply(meta_ttests, function(tti) tti$conf.int[2]) %>% simplify
#   # )
# }


#' plot_heatmap
#'
#' Plots tidy dataframe from a deeptools matrix file to a ggplot2 object.
#'
#' @param df the df to load.
#' @param sort_y how to sort y axis of the heatmap ('region_len', 'sum', 'mean', 'max' )
#' @param facets a string for facetting of the plot.
#'
#' @details Plots a df from a deeptools matrix tibble.
#'
#' @return An ggplot2 plot object.
#'
#' @examples
#'
#'  fname <- "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/pA_seq/deeptools_out/Evgenia_PASseq_all_snRNA_joined.gz"
#'  df <- load_deeptoolsmatrix(fname)
#'  head(df)
#'  plot_heatmap(df, sort_y='sum', facets='.~sample_name')
#'
#'  #can be further modified
#'  p <- plot_heatmap(df, sort_y='sum', facets='.~sample_name')
#'  p + ylab('')
#'  p + scale_fill_gradient2(low='white', high='darkred')
#'  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#'  p + scale_fill_gradient2(low='white', high='darkred', name ='log10(norm. reads)')
#'
#' @export
plot_heatmap <- function( df, sort_y='sum', facets=NULL ) {

  if ( missing(sort_y) ) {
    y_rank <- df %>% distinct(id) %>% mutate(sort_value = 1:nrow(.))
  } else if (sort_y == 'sum') {
    y_rank <- df %>% group_by(id) %>% summarize(sort_value = sum(value)) %>% select(id, sort_value)
  } else if (sort_y == 'mean') {
    y_rank <- df %>% group_by(id) %>% summarize(sort_value = mean(value)) %>% select(id, sort_value)
  } else if (sort_y == 'max') {
    y_rank <- df %>% group_by(id) %>% summarize(sort_value = max(value)) %>% select(id, sort_value)
  } else if (sort_y == 'region_length') {
    y_rank <- df %>% distinct(id, .keep_all=T) %>% mutate(sort_value = end-start)
  }

  y_rank %<>% mutate(yrank = rank(sort_value, ties.method='first'))

  p <- left_join(df,y_rank)  %>%
    ggplot(., aes(x=rel_pos, y=yrank, fill=log10(value))) +
    geom_raster(interpolate=FALSE) +
    theme_bw()

  if (facets) {
    p <- p +
      facet_grid(facets)
  }

  p
}



#' Sort metagene object regions
#'
#' Sort metagene object regions as sum per line
#'
#' @param m metagene object
#'
#' @details Sorts the matrix and regions.
#'
#' @return The sorted metagene object
#'
#' @examples
#'
#'  matrix_name <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/eRNA_tsspm5000_sorted_by_INT_binding_joined_sensitivity.gz'
#'  bed_name <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Shiekkatar_INT_ChIP/eRNAs_sorted_by_INTS_binding.bed'
#'
#'  m_sorted <- sort_matrix(matrix_name, bed_name)
#'
#'  write.table(m_sorted, file='/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/eRNA_tsspm5000_sorted_by_INT_binding_joined_sensitivity_sorted', row.names = FALSE, col.names = FALSE, quote = F, dec = '.')
#'
#'
#' @export
sort_matrix <- function( m, bed_file ) {

  groups <- unique(m$groups)
  grp_order <- lapply(groups, function(group) {
    lapply(m$matrices[m$groups == group], rowSums) %>%
      bind_rows %>%
      rowSums %>%
      order
  })

  sorted <- m
  sorted$matrices <- lapply(seq_along(m$matrices), function(i)
    m$matrices[[i]][grp_order[[m$groups[[i]]]]])

  sorted$regions <- lapply(seq_along(m$regions), function(i)
    m$regions[[i]][grp_order[[i]]])

  b <- read.table(bed_file)[,1:6]
  colnames(b) <- c('chr', 'start', 'end', 'id', 'name', 'strand')
  b %<>% as_tibble

  semi_join(m, b)
}




#' Plot Heatmap
#'
#' Heatmap of specified tidy metagene tbl.
#'
#' @param tbl tidy metagene tbl
#' @param color_by name of column from tbl to use as value (default='value')
#' @param do_interpolate interpolate geom_raster output (default=FALSE)
#'
#' @details uses ggplot2 to create plot of the heatmap
#'
#' @return ggplot2 object
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' anno=bedfile;anchor='center';upstream=1000;downstream=1000;window_size=10
#' collapse_fun = rowMeans;negate_neg_strand_values=FALSE
#' metamat <- metagene_matrix(bw_plus, bw_minus, bedfile, anchor, upstream, downstream, window_size)
#' plot_htmp(metamat, color_by='value', do_interpolate=FALSE)
#'
#' #log2 values incl pseudocount
#' metamat %<>% mutate(log2_value = log2(value+1))
#' plot_htmp(metamat, color_by='log2_value', do_interpolate=FALSE)
#'
#' @export
plot_htmp <- function(tbl, color_by='value', do_interpolate=FALSE){
  ggplot(tbl, aes_string(x='rel_pos', y='name', fill=color_by)) +
    geom_raster(interpolate=do_interpolate) +
    scale_fill_gradient(low='white', high='black')+
    scale_x_continuous(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}



#' Create Metagene Averages
#'
#' Per position bin averages of specified tidy metagene tbl.
#'
#' @param tbl tidy metagene tbl
#'
#' @details uses t.test to get mean ('estimate'), confidence interval etc for each position bin converted to columns using broom::tidy
#'
#' @return tbl
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' anno=bedfile;anchor='center';upstream=1000;downstream=1000;window_size=10
#' collapse_fun = rowMeans;negate_neg_strand_values=FALSE
#' metamat <- metagene_matrix(bw_plus, bw_minus, bedfile, anchor, upstream, downstream, window_size)
#' tidy_avgs <- meta_average(metamat)
#'
#' #alt:
#' tidy_log2_avgs <- meta_average(mutate(metamat, value = log2(value+1)))
#'
#' @export
meta_average <- function(tbl){
  tbl %>%
    group_by_at(vars(-value, -name)) %>%
    do(tidy(t.test(.$value)))
}



#' Plot Metagene Average Profile
#'
#' Plot Metagene Average Profile
#'
#' @param tbl tidy metagene tbl
#' @param color_by name of column from tbl to use for coloring of lines and ranges (default=NULL)
#'
#' @details plots line from ('estimate'), and range from conf.low and conf.high columns in tbl
#'
#' @return tbl
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' anno=bedfile;anchor='center';upstream=1000;downstream=1000;window_size=10
#' collapse_fun = rowMeans;negate_neg_strand_values=FALSE
#' metamat <- metagene_matrix(bw_plus, bw_minus, bedfile, anchor, upstream, downstream, window_size)
#' tidy_log2_avgs <- meta_average(mutate(metamat, value = log2(value+1)))
#' plot_profile(tidy_log2_avgs)
#'
#' @export
plot_profile <- function(meta_tbl, color_by=NULL){
  ggplot(meta_tbl) +
    geom_ribbon(aes_string(x='rel_pos', ymin='conf.low', ymax='conf.high', fill=color_by), alpha=.2) +
    geom_line(aes_string(x='rel_pos', y='estimate', color=color_by)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.grid = element_blank())
}
