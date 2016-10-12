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
#'  fname <- "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/pA_seq/Evgenia_PASseq_all_snRNA_joined.gz"
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
  up <- as.numeric(header['upstream'])
  dn <- as.numeric(header['downstream'])
  if (header['ref point'] != 'null') {
    xlabels <- seq(-up,dn-1,bin)
  } else {
    xlabels <- seq(-up,1000+dn-1,bin)
  }

  mgroups <- unlist(header['group_labels'])
  group_bounds <- as.numeric(unlist(header['group_boundaries']))

  if (as.matrix) {
    attr(values, 'header') <- header
    colnames(values)[1:6] <- c('chr', 'start', 'end', 'id', 'name', 'strand')

    name_vec <- c()
    for(i in seq_along(msamples)) {
      name <- msamples[i]
      left <- sample_bounds[i] + 7
      right <- sample_bounds[i+1] + 6
      name_vec <- c(name_vec, rep(name, right-left))
    }
    values$sample <- rep(msamples,)
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

  tbl_df(d)
}


#' parse_matrix_header
#'
#' Parses a deeptools matrix file header.
#' @param file the file to get the header from.
#'
#' @details Loads a gzipped deeptools matrix file and collects this in a \code{tidy} data.frame.
parse_matrix_header <- function(fname) {

  header <-  scan(gzfile(fname), nlines=1, skip=0, what='character', sep='#') %>% sub('@\\{', ',', .) %>% sub('\\}', '', .) %>% sub('"', '', .)

  chunks <- strsplit(header, ',')

  meta <- c()
  keep_pasting <- FALSE
  for (c in chunks[[1]][2:length(chunks[[1]])]) {
    if (keep_pasting){
      meta[length(meta)] <- paste(meta[length(meta)], c, sep=';')
      if (grepl('\\]',c)) {
        keep_pasting <- FALSE
      }
      next()
    }
    if (grepl('\\[',c) && !grepl('\\]',c)) {
      keep_pasting <- TRUE
    }

    meta <- c(meta, c)
  }

  tag_id <- sapply(meta, function(m) strsplit(m, ':')[[1]][1])
  tag_vals <- lapply(meta, function(m) strsplit(m, ':')[[1]][2])

  names(tag_vals) <- tag_id

  tag_vals <- lapply(tag_vals, function(m) sub('\\[', '', m))
  tag_vals <- lapply(tag_vals, function(m) sub('\\]', '', m))
  tag_vals <- lapply(tag_vals, function(m) strsplit(m, ';')[[1]])

  tag_vals
}
