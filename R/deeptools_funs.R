#' load_deeptoolsmatrix
#'
#' Loads a deeptools matrix file.
#'
#' @param file the file to load.
#' @param na.omit Skip na values in the matrix.
#' @param na.fill Fills NAs with this value zero.
#' @param as.matrix Return matrix as matrix (default is FALSE, returns tidy data.frame)
#'
#' @details Loads a gzipped deeptools matrix file and collects this in a \code{tidy} data.frame. If as.matrix is true, the header information is added as attribute.
#'
#' @return An object of class \code{deeptools_matrix} with provided means and covariances
#'   computed as 1/\code{alpha}.
#'
#' @examples
#'
#'  fname <- system.file("extdata", "matrix.gz", package = "RMetaTools")
#'  df <- load_deeptoolsmatrix(fname)
#'  head(df)
#'  df_withNA <- load_deeptoolsmatrix(fname, na.omit=FALSE)
#'  head(df_withNA)
#'  df_withNAas0 <- load_deeptoolsmatrix(fname, na.omit=FALSE, na.fill=0)
#'  head(df_withNAas0)
#'
#' @export
load_deeptoolsmatrix <- function(fname, na.omit=TRUE, na.fill=0, as.matrix=FALSE){
  header <- parse_matrix_header(fname)
  values <- read.table(gzfile(fname), header=F, skip=1)

  bed_info <- values[,1:6]
  msamples <- unlist(header['sample_labels'])
  sample_bounds <- as.numeric(unlist(header['sample_boundaries']))
  bin <- as.numeric(header['bin size'])
  body <- as.numeric(header['body'])
  up <- as.numeric(header['upstream'])
  dn <- as.numeric(header['downstream'])
  xlabels <- seq(-up,body+dn-1,bin)

  mgroups <- unlist(header['group_labels'])
  group_bounds <- as.numeric(unlist(header['group_boundaries']))

  if (as.matrix) {
    attr(values, 'header') <- header
    bednames <- c('chr', 'start', 'end', 'id', 'name', 'strand')

    name_vec <- c()
    for(i in seq_along(msamples)) {
      name <- msamples[i]
      left <- sample_bounds[i] + 7
      right <- sample_bounds[i+1] + 6
      name_vec <- c(name_vec, paste0(name, xlabels))
    }
    colnames(values) <- c(bednames, name_vec)
    return(values)
  }

  dfs <- list()
  for(i in seq_along(msamples)) {
    name <- msamples[i]
    left <- sample_bounds[i] + 7
    right <- sample_bounds[i+1] + 6
    s_values <- cbind(bed_info, values[,left:right])
    colnames(s_values) <- c('chr', 'start', 'end', 'id', 'name', 'strand', xlabels)

    s_values$sample_name <- name
    s_values$group <- mgroups[1]
    for(g in seq_along(mgroups)) {
      upper <- group_bounds[g] + 1
      lower <- group_bounds[g+1]
      s_values$group[upper:lower] <- mgroups[g]
    }
    dfs[[name]] <- gather(s_values, rel_pos, value, -sample_name, -group, -c(1:6)) %>%
      mutate(rel_pos = as.numeric(rel_pos))
  }

  d <- bind_rows(dfs)
  if(na.omit) {
    d <- filter(d, !is.na(value))
  } else if (!missing(na.fill)) {
    fill_val <- as.numeric(na.fill)
    d <- mutate(d, value = ifelse(is.na(value), fill_val, value))
  }

  as_data_frame(d)
}


#' load_deeptoolsmatrix3
#'
#' Loads a deeptools v3 matrix file.
#'
#' @param fname the file to load.
#' @param na.omit Skip na values in the matrix. Default is TRUE.
#' @param na.fill Fills NAs with this value zero. Default is 0.
#' @param as.matrix Return matrix as matrix (default is FALSE, ie returns tidy data.frame instead)
#' @param as.S3 Return S3 class metagene object, see details for more info
#'
#' @details Loads a gzipped deeptools matrix file and collects this in a
#'   \code{tidy} data.frame. If as.matrix is true, the header information is
#'   added as attribute.
#'   as.S3 returns a metagene object, with slots: $header, $regions, $samples,
#'   $groups and $matrices. Where $header contains the header information as a
#'   named vector. Matrices is a list of matrices where each matrix is specific
#'   to 1 sample and 1 group. $samples are the samplenames for each matrix.
#'   $groups are the groups for each matrix.
#'
#' @return An object of class \code{deeptools_matrix} with provided means and covariances
#'   computed as 1/\code{alpha}.
#'
#' @examples
#'
#'  fname <- system.file("extdata", "deeptools3_matrix.gz", package = "RMetaTools")
#'  df <- load_deeptoolsmatrix3(fname)
#'  head(df)
#'  df_withNA <- load_deeptoolsmatrix3(fname, na.omit=FALSE)
#'  head(df_withNA)
#'  df_withNAas0 <- load_deeptoolsmatrix3(fname, na.omit=FALSE, na.fill=0)
#'  head(df_withNAas0)
#'  s3obj <- load_deeptoolsmatrix3(fname, na.omit=FALSE, na.fill=0, as.S3=TRUE)
#'  head(s3obj$regions)
#'
#' @export
load_deeptoolsmatrix3 <- function(fname, na.omit=TRUE, na.fill=0, as.matrix=FALSE, as.S3=FALSE){
  header <- parse_matrix_header(fname)
  values <- read.table(gzfile(fname), header=F, skip=1)

  bed_info <- values[,1:6]
  msamples <- unlist(header[['sample_labels']])
  sample_bounds <- as.numeric(unlist(header[['sample_boundaries']]))
  bin <- as.numeric(header[['bin size']])
  body <- as.numeric(header[['body']])
  up <- as.numeric(header[['upstream']])
  dn <- as.numeric(header[['downstream']])
  unscaled5 <- as.numeric(header[['unscaled 5 prime']])
  unscaled3 <- as.numeric(header[['unscaled 3 prime']])
  xlabels <- lapply(seq_along(bin), function(i) seq(-up[i],unscaled5[i]+body[i]+unscaled3[i]+dn[i]-1,bin[i]))

  mgroups <- unlist(header[['group_labels']])
  group_bounds <- as.numeric(unlist(header[['group_boundaries']]))

  if (as.matrix) {
    attr(values, 'header') <- header
    bednames <- c('chr', 'start', 'end', 'id', 'name', 'strand')

    name_vec <- c()
    for(i in seq_along(msamples)) {
      name <- msamples[i]
      left <- sample_bounds[i] + 7
      right <- sample_bounds[i+1] + 6
      name_vec <- c(name_vec, paste0(name, '@', as.character(xlabels[[1]])))
    }
    colnames(values) <- c(bednames, name_vec)
    return(values)
  }

  if (as.S3) {
    s3 <- structure(list(), class = "metagene")
    s3$header <- header
    colnames(bed_info) <- c('chr', 'start', 'end', 'id', 'name', 'strand')
    s3$regions <- tbl_df(bed_info)

    indeces <- expand.grid(seq_along(msamples), seq_along(mgroups))

    mat_ranges <- data.frame(t(data.frame(sample_start_idx = sample_bounds[indeces[,1]] + 7,
                                          sample_end_idx = sample_bounds[indeces[,1]+1] + 6,
                                          group_start_idx = group_bounds[indeces[,2]]+1,
                                          group_end_idx = group_bounds[indeces[,2]+1])))

    s3$matrices <- lapply(mat_ranges, function(x) {
      values[x[3]:x[4], x[1]:x[2]]
    })

    s3$samples <- rep(msamples, length(mgroups))
    s3$groups <- rep(mgroups, each=length(msamples))
    s3$rel_pos <- seq(-up,unscaled5+body+unscaled3+dn-1,bin)

    return(s3)
  }

  dfs <- list()
  for(i in seq_along(msamples)) {
    name <- msamples[i]
    left <- sample_bounds[i] + 7
    right <- sample_bounds[i+1] + 6
    s_values <- cbind(bed_info, values[,left:right])
    colnames(s_values) <- c('chr', 'start', 'end', 'id', 'name', 'strand', xlabels[[i]])

    s_values$sample_name <- name
    s_values$group <- mgroups[1]
    for(g in seq_along(mgroups)) {
      upper <- group_bounds[g] + 1
      lower <- group_bounds[g+1]
      s_values$group[upper:lower] <- mgroups[g]
    }
    dfs[[name]] <- gather(s_values, rel_pos, value, -sample_name, -group, -c(1:6)) %>%
      mutate(rel_pos = as.numeric(rel_pos))
  }

  d <- bind_rows(dfs)
  if(na.omit) {
    d <- filter(d, !is.na(value))
  } else if (!missing(na.fill)) {
    fill_val <- as.numeric(na.fill)
    d <- mutate(d, value = ifelse(is.na(value), fill_val, value))
  }

  as_data_frame(d)
}




#' parse_matrix_header
#'
#' Parses a deeptools matrix file header.
#' @param file the file to get the header from.
#'
#' @details Loads a gzipped deeptools matrix file and collects this in a \code{tidy} data.frame.
parse_matrix_header <- function(fname) {

  header <-  scan(gzfile(fname), nlines=1, skip=0, what='character', sep='#') %>% sub('@\\{', '', .) %>% sub('\\}', '', .) %>% sub('"', '', .)

  ar_starts <- gregexpr('\\[', header)[[1]]
  ar_ends <- gregexpr('\\]', header)[[1]]

  if ( length(ar_starts) != length(ar_ends) ) {
    stop("problem in parsing header information in matrix")
  }

  header <- sub(';', '_', header) #to allow hack and use ';' as separator inside json arrays
  for ( i in seq_along(ar_starts) ) {
    header <- paste0(substring(header, 1, ar_starts[i]-1), gsub(',', ';', substring(header, ar_starts[i], ar_ends[i])), substring(header, ar_ends[i]+1, nchar(header)))
  }

  chunks <- unlist(strsplit(header, ','))

  tag_ids <- sapply(chunks, function(m) strsplit(m, ':')[[1]][1])
  tag_vals <- lapply(chunks, function(m) strsplit(m, ':')[[1]][2])

  names(tag_vals) <- tag_ids

  tag_vals <- lapply(tag_vals, function(m) sub('\\[', '', m))
  tag_vals <- lapply(tag_vals, function(m) sub('\\]', '', m))
  tag_vals <- lapply(tag_vals, function(m) strsplit(m, ';')[[1]])

  tag_vals
}




#' metagene
#'
#' Metagene aggregate from an S3 metagene object.
#'
#' @param obj the object to create the metagene aggregate for
#' @param aggregate_fun the function name to aggregate values such as 'events', 'sum', 'mean', 'ttest' (default is 'mean')
#' @param transform_fun a function to be applied to vlues before aggregate function is applied.
#' @param ... other arguments passed to fun
#'
#' @details Aggregates values.
#' transform_fun takes a matrix as input and returns a similar sized matrix. ie log2p1 <- function(m){log2(m+1)}
#' aggreagate_fun can be mean, sum, events or t.test or vector containing several of those.
#'
#' @return
#' a tbl_df with column, rel_pos, sample, group, and aggregate values where column name is mean, sum, events or estimate+conf.low+conf.high for ttest
#' @examples
#'  fname <- system.file("extdata", "matrix_grouped.gz", package = "RMetaTools")
#'  df <- load_deeptoolsmatrix3(fname, as.S3=TRUE)
#'  meta_df <- metagene_aggregate(df)
#'  meta_df <- metagene_aggregate(df, aggregate_fun = 'events')
#'  ggplot(meta_df, aes(x=rel_pos, y=events, color=sample)) + geom_line() + facet_grid(group~.)
#'  meta_df <- metagene_aggregate(df, aggregate_fun = 'ttest')
#'  log2p1 <- function(m){log2(m+1)}
#'  meta_df <- metagene_aggregate(df, aggregate_fun = 'ttest', transform_fun = log2p1)
#'
#' @export
metagene_aggregate <- function( obj, aggregate_fun = 'mean', transform_fun = NULL, ... ) {

  if (missing(transform_fun)) {
    mat_list <- obj$matrices
  }else{
    mat_list <- lapply(obj$matrices, transform_fun)
  }

  meta_df = tbl_df(data.frame(
    rel_pos = rep(s3$rel_pos, length(s3$matrices)),
    sample = rep(s3$samples, each=length(s3$rel_pos)),
    group = rep(s3$groups, each=length(s3$rel_pos))))

  if( 'mean' %in% aggregate_fun ){
    meta_df$mean = lapply(mat_list, colMeans) %>% simplify
  }else if( 'sum' %in% aggregate_fun ){
    meta_df$sum = lapply(mat_list, colSums) %>% simplify
  }else if( 'events' %in% aggregate_fun ){
    meta_df$events = lapply(mat_list, function(m) apply(m,2,function(mi) sum(mi > 0))) %>% simplify
  }else if( 'ttest' %in% aggregate_fun ){
    meta_ttests <- lapply(mat_list, function(mi) apply(mi,2,t.test))

    meta_df$estimate = lapply(meta_ttests, function(tt) lapply(tt, function(tti) tti$estimate) %>% simplify) %>% simplify
    meta_df$conf.low =  lapply(meta_ttests, function(tt) lapply(tt, function(tti) tti$conf.int[1]) %>% simplify) %>% simplify
    meta_df$conf.high =  lapply(meta_ttests, function(tt) lapply(tt, function(tti) tti$conf.int[2]) %>% simplify) %>% simplify

  }

  meta_df
}




#' plot_heatmap
#'
#' Plots tidy dataframe from a deeptools matrix file to a ggplot2 object.
#'
#' @param df the df to load.
#' @param sort_y how to sort y axis of the heatmap ('region_len', 'sum', 'mean', 'max' )
#' @param facets a string for facetting of the plot.
#'
#' @details Plots a df from a deeptools matrix tbl_df.
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



#' Sort Matrix
#'
#' Sorts a deeptools matrix file using a bed file.
#'
#' @param matrix_file filename for the matrix (gzipped).
#' @param bed_file filename for the bed file used for sorting
#'
#' @details Sorts the rows of the matrix without header.
#'
#' @return The sorted matrix.
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
sort_matrix <- function( matrix_file, bed_file ) {

  m <- load_deeptoolsmatrix(matrix_file, as.matrix = TRUE) %>%
    tbl_df

  b <- read.table(bed_file)[,1:6]
  colnames(b) <- c('chr', 'start', 'end', 'id', 'name', 'strand')
  b %<>% tbl_df

  semi_join(m, b)
}




# ```{r}
# mod <- function(mats, fun){
#   lapply(mats, fun)
# }
# library(broom)
# tmat <- list(matrix(1:10,ncol=2), matrix(11:20,ncol=2))
# names(tmat) <- c('s1', 's2')
# add1_and_log2 <- function(mat, pseudocount=1){ log2(mat+pseudocount) }
# ttcols <- function(mat){bind_rows(apply(mat,2,function(col) tidy(t.test(col))))}
# tmat2 <- mod(tmat, add1_and_log2)
# tmat_tt <- mod(tmat2, ttcols)
#
# tmat3 <- imap(tmat, ~ log2(.x+1))
# names(tmat3)
#
# tmat3 <- imap(tmat, ~ add1_and_log2(., pseudocount=1))
#
# transform_fun <- add1_and_log2
# tmat3 <- imap(tmat, ~ transform_fun(., pseudocount=1))
#
# ```

#'  fname <- system.file("extdata", "deeptools3_matrix.gz", package = "RMetaTools")
#'  df <- load_deeptoolsmatrix3(fname)
#'  head(df)
#'  df_withNA <- load_deeptoolsmatrix3(fname, na.omit=FALSE)
#'  head(df_withNA)
#'  df_withNAas0 <- load_deeptoolsmatrix3(fname, na.omit=FALSE, na.fill=0)
#'  head(df_withNAas0)
#'  s3obj <- load_deeptoolsmatrix3(fname, na.omit=FALSE, na.fill=0, as.S3=TRUE)
#'  head(df_withNAas0)
# fname <- '/Users/au220280/ms_tools/MS_Metagene_Tools/test/matrix_grouped.gz'
#
# header <- parse_matrix_header(fname)
# values <- read.table(gzfile(fname), header=F, skip=1)
#
# bed_info <- values[,1:6]
# msamples <- unlist(header[['sample_labels']])
# sample_bounds <- as.numeric(unlist(header[['sample_boundaries']]))
# bin <- as.numeric(header[['bin size']])
# body <- as.numeric(header[['body']])
# up <- as.numeric(header[['upstream']])
# unscaled5 <- as.numeric(header[['unscaled 5 prime']])
# unscaled3 <- as.numeric(header[['unscaled 3 prime']])
# dn <- as.numeric(header[['downstream']])
# xlabels <- lapply(seq_along(bin), function(i) seq(-up[i],body[i]+dn[i]-1,bin[i]))
#
# mgroups <- unlist(header[['group_labels']])
# group_bounds <- as.numeric(unlist(header[['group_boundaries']]))
#
#
# indeces <- expand.grid(seq_along(msamples), seq_along(mgroups))
#
# mat_ranges <- data.frame(t(data.frame(sample_start_idx = sample_bounds[indeces[,1]] + 7,
#                      sample_end_idx = sample_bounds[indeces[,1]+1] + 6,
#                      group_start_idx = group_bounds[indeces[,2]]+1,
#                      group_end_idx = group_bounds[indeces[,2]+1])))
#
# s3$matrices <- lapply(mat_ranges, function(x) {
#   values[x[3]:x[4], x[1]:x[2]]
# })
#
# s3$samples <- rep(msamples, length(mgroups))
# s3$groups <- rep(mgroups, each=length(msamples))
# s3$rel_pos <- seq(-up,unscaled5+body+unscaled3+dn-1,bin)
#
# aggfun <- colSums
#
#
#
# ggplot(meta_ev, aes(x=rel_pos, y=events, color=sample)) +
#   geom_line() +
#   facet_wrap(~group)
