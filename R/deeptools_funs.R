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
#'  df_withNAas0 <- load_deeptoolsmatrix3(fname, na.fill=0)
#'  head(df_withNAas0)
#'  s3obj <- load_deeptoolsmatrix3(fname, na.fill=0, as.S3=TRUE)
#'  head(s3obj$regions)
#'  lapply(s3obj$matrices, head)
#'  s3obj$samples
#'  s3obj$groups
#'  lapply(s3obj$rel_pos, head)
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


  if ( !missing(na.fill) ) {
     values_df <- as.data.frame(values[,7:ncol(values)])
     values_df[is.na(values_df)] <- na.fill
     values[,7:ncol(values)] <- as.matrix(values_df)
   }

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
    s3$regions <- as_tibble(bed_info)

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
    s3$rel_pos <- rep(xlabels, length(mgroups))

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

  as_data_frame(d)
}


#' load_deeptoolsmatrix3_new
#'
#' Loads a deeptools v3 matrix file new version.
#'
#' @param fname the file to load.
#' @param na.fill Fills NAs with this value zero. Default is NA.
#'
#' @details Always returns and S3 object. Loads a gzipped deeptools matrix file and collects this in a
#'   metagene object, with slots: $header, $regions, $samples,
#'   $groups and $matrices. Where $header contains the header information as a
#'   named vector. Matrices is a list of matrices where each matrix is specific
#'   to 1 sample and 1 group. $samples are the samplenames for each matrix.
#'   $groups are the groups for each matrix. $rel_pos is the list of relative positions to anchorpoint for each group. Names are the group names. $regions same as rel_pos are the region names (ie col4 of bed file).
#'
#' @return An object of class \code{metagene}.
#'
#' @examples
#'
#'  fname <- system.file("extdata", "deeptools3_matrix.gz", package = "RMetaTools")
#'  m <- load_deeptoolsmatrix3_new(fname)
#'  m$header
#'  m$samples
#'  m$groups
#'  lapply(m$rel_pos, head)
#'  lapply(m$regions, head)
#'  lapply(m$matrices, dim)
#'
#' @export
load_deeptoolsmatrix3_new <- function(fname, na.fill){
  header <- parse_matrix_header(fname)
  values <- read.table(gzfile(fname), header=F, skip=1, stringsAsFactors = F)

  bed_info <- as_tibble(values[,1:6])
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


  if ( !missing(na.fill) ) {
    values_df <- as.data.frame(values[,7:ncol(values)])
    values_df[is.na(values_df)] <- na.fill
    values[,7:ncol(values)] <- as.matrix(values_df)
  }

  s3 <- structure(list(), class = "metagene")
  s3$header <- header

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

  s3$rel_pos <- xlabels[[1]]

  colnames(bed_info) <- c('chr', 'start', 'end', 'id', 'score', 'strand')
  s3$regions <- lapply(seq_along(mgroups), function(i) {
    bed_info[(group_bounds[i]+1):(group_bounds[i+1]),]
  })
  names(s3$regions) <- mgroups

  return(s3)

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





