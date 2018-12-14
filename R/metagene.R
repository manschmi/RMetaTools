#' Metagene Regions
#'
#' Creates regions for metagene analysis.
#'
#' @param bed filename for bed file.
#' @param anchor anchor point TSS, TES or center (default='TSS')
#' @param upstream bp upstream of anchor (default=1000)
#' @param downstream bp downstream of anchor (default=1000)
#'
#' @details Loads the bed file, extracts the position of the anchor point and creates ranges from x bp upstream to x bp downstream.
#'
#' @return GRanges with regions.
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000)
#'
#' @export
meta_regions <- function(bed,
                        anchor = 'TSS',
                        upstream = 1000,
                        downstream = 1000) {

  anno <- rtracklayer::import(bed)
  anno$original_start <- start(anno)
  anno$original_end <- end(anno)
  if( !('name' %in% colnames(mcols(anno))) ) {
    anno$name <- paste0('region', seq(length(anno)))
  }

  strands <- as.vector(strand(anno))
  anchor_position <- case_when(
    (anchor == 'TSS' & strands == '-') ~ end(anno),
    anchor == 'TSS' ~ start(anno),
    anchor == 'TES' & strands == '-' ~ start(anno),
    anchor == 'TES' ~ end(anno),
    anchor == 'center' ~ as.integer(start(anno) + (end(anno)-start(anno))/2)
  )

  start(anno) <- ifelse(strand(anno) == '+', anchor_position - upstream, anchor_position - downstream)
  end(anno) <- ifelse(strand(anno) == '+', anchor_position + downstream - 1, anchor_position + upstream - 1) #0 counts as first position

  return(anno)
}


#' Single Strand Metagene Matrix
#'
#' Creates metagene matrix for single bigwig file.
#'
#' @param bw filename for bigwig file.
#' @param anno GRanges object of regions to query.
#' @param upstream bp upstream of anchor (default=1000)
#' @param downstream bp downstream of anchor (default=1000)
#' @param strand strand of anno to use (default='+')
#'
#' @details Loads the bed file, extracts the position of the anchor point and creates ranges from x bp upstream to x bp downstream.
#'
#' @return matrix with columns from most upstream to most downstream and rows are the individual regions from specified strand. rownames are set as name column from anno. colnames are set to relative position of each nt
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000)
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' mat <- get_single_strand_matrix(bw_plus, regions, 1000, 5000, '+')
#'
#' @export
get_single_strand_matrix <- function(bw, anno, upstream=1000, downstream=1000, strand="+") {

  ##handling specific strand and order to ensure correct output naming
  seqlevels(anno) <- sort(seqlevels(anno))
  anno <- sort(anno)

  ##ensure to query only chrs present in bigwig,
  ## ie otherwise import(bw, ...) will return in error
  bw_chrs <- seqlevels(BigWigFile(bw))
  contained_anno_rows <- which(seqnames(anno) %in% bw_chrs)

  if (length(contained_anno_rows) != length(anno)) {
    warning(paste0('ignoring chromosomes from annotation file not present in the bigwig file for strand: ', strand,
                   ': ', paste0(seqlevels(anno[-contained_anno_rows]))))
  }

  mat <- matrix(NA, nrow = length(anno), ncol = upstream+downstream)
  matc <- as.matrix(rtracklayer::import(bw,
                                       which=anno[contained_anno_rows],
                                       as='NumericList'))
  mat[contained_anno_rows,] <- matc

  rownames(mat) <- paste(seqnames(anno),
                         anno$original_start,
                         anno$original_end,
                         anno$name,
                         strand(anno),
                         sep = '\t')

  colnames(mat) <- seq(-upstream, downstream-1)

  mat
}



#' Collapse Matrix To Windows
#'
#' Averages Matrix Columns into Windows using an Averaging Function
#'
#' @param mat matrix
#' @param collapse_fun function for averaging values inside windows (Default = rowMeans, see Details).
#' @param window_size size for binning in bp (default=10, ie average over 10nt)
#' @param ... additional arguments passed to collapse_fun
#'
#' @details for a matrix with ncol=100 and window_size=10, will return a matrix with ncol=10. collapse_fun is a function that functions as ie rowMeans (the default) or rowSums.
#'
#' @return matrix with ncol(input matrix)/window_size columns
#'
#' @examples
#'
#' mat <- matrix(1:1000, nrow=10, ncol=100)
#' matc <- collapse_to_window_size(mat, rowMeans, 10)
#' dim(matc)
#'
#' log2RowMeans <- function(mat, pseudocount=1){apply(mat,1,function(row) mean(log2(row+pseudocount)))}
#' matc <- collapse_to_window_size(mat, log2RowMeans, 10, pseudocount = .1)
#' dim(matc)
#'
#' @export
collapse_to_window_size <- function(mat, collapse_fun = rowMeans, window_size = 10, ...) {
  if (window_size > 1) {
    size <- ncol(mat)
    wins <- seq(1,size, window_size)

    matc <- lapply(wins, function(i) collapse_fun(mat[,i:(i+window_size-1)], ...)) %>%
      bind_cols %>%
      as.matrix

    colnames(matc) <- sapply(wins, function(i) colnames(mat)[i])
    rownames(matc) <- rownames(mat)

    mat <- matc
  }

  mat
}


#' Metagene Matrix
#'
#' Creates metagene matrix for both strands using 1 or 2 bigwig files.
#'
#' @param bw_plus filename for plus strand bigwig file.
#' @param bw_minus filename for minus strand bigwig file.
#' @param anno GRanges object of regions to query.
#' @param upstream bp upstream of anchor (default=1000)
#' @param downstream bp downstream of anchor (default=1000)
#' @param window_size size for binning in bp (default=1, ie no binning)
#' @param collapse_fun function for averaging values inside windows (Default = rowMeans, see Details).
#' @param collapse_fun_args additional arguments passed to collapse_fun (default=NULL)
#' @param negate_neg_strand_values negate values from minus strand bigwig (default=FALSE).
#' @details Loads the bed file, extracts the position of the anchor point and creates ranges from x bp upstream to x bp downstream using collapse_fun average over window_size regions. See \code{\link{collapse_to_window_size}} for more info about averaging over window size. negate_neg_strand_values can be set to TRUE to deal with minus strand bigwigs which contain all negative values as used by some labs.
#'
#' @return matrix with columns from most upstream to most downstream and rows are the individual regions. Rownames are the Name column from the bed file.
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000)
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 5000, 1)
#'
#' #if bigwigs are not stranded, ie from ChIP experiment both plus and minus are the same file
#' bw <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_hg38.bw", package = "RMetaTools")
#' mat <- get_matrix(bw, bw, regions, 1000, 5000, 50)
#'
#' @export
get_matrix <- function(bw_plus, bw_minus, anno, upstream=1000, downstream=1000, window_size = 1,
                       collapse_fun = rowMeans, collapse_fun_args = NULL,
                       negate_neg_strand_values=FALSE) {
  smat <- get_single_strand_matrix(bw_plus, anno[strand(anno)=="+"], upstream, downstream)
  asmat <- get_single_strand_matrix(bw_minus, anno[strand(anno)=="-"], upstream, downstream)
  asmat <- asmat[,ncol(asmat):1]
  if (negate_neg_strand_values) {
    asmat <- -asmat
  }
  mat <- rbind(smat, asmat)
  if(window_size > 1) {
    mat <- collapse_to_window_size(mat, collapse_fun, window_size, collapse_fun_args)
  }
  mat
}



#' Tidy Metagene Matrix
#'
#' Converts metagene matrix into a tidy tbl.
#'
#' @param mat metagene matrix created from get_matrix
#'
#' @details Converts metagene matrix to tbl using rownames and colnames information to assign a column *name* and *rel_pos* respectively.
#'
#' @return matrix with columns from most upstream to most downstream and rows are the individual regions.
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000)
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 1000, 50)
#' tidy_meta <- mat_to_tbl(mat)
#'
#' @export
mat_to_tbl <- function(mat) {

  nm <- mat %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = 'name') %>%
    tidyr::gather(., relpos, value, -name) %>%
    mutate(relpos = as.numeric(relpos)) %>%
    tbl_df

  nm
}



#' Metagene Matrix
#'
#' Computes a metagene matrix from scratch.
#'
#' @param bw_plus bigwig filename for plus strand signal
#' @param bw_minus bigwig filename for minus strand signal
#' @param anno Bed filename for annotations to use OR a GRanges object (see Details).
#' @param anchor 'TSS', 'TES' or 'center.'
#' @param upstream bp upstream of anchor (default=1000)
#' @param downstream bp downstream of anchor (default=1000)
#' @param window_size size for binning in bp (default=1, ie no binning)
#' @param collapse_fun function for averaging values inside windows (Default = rowMeans, see Details).
#' @param collapse_fun_args additional arguments passed to collapse_fun (default=NULL)
#' @param negate_neg_strand_values negate values from minus strand bigwig (default=FALSE).
#'
#' @details Loads the bed file, extracts the position of the anchor point and creates ranges from x bp upstream to x bp downstream. Queries from both bigwigs and the correct strand. Summarizes signal of window_size bins and returns a tidy data.frame. negate_neg_strand_values can be set to TRUE to deal with minus strand bigwigs which contain all negative values as used by some labs.
#' When using a GRanges object for anno, this assumes the bigwigs will be queries from start to end! All ranges must be exactly same size!
#' @return Tidy data frame with columns: gene (name, ie col4 from *bed* file), rel_pos (position in bp relative to *anchor*), value (averaged value from bigwig file).
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' metamat <- metagene_matrix(bw_plus, bw_minus, bedfile, 'center', 1000, 1000, 50)
#'
#' @export
metagene_matrix <- function(bw_plus, bw_minus, anno, upstream=1000, downstream=1000, window_size = 1,
                            collapse_fun = rowMeans, collapse_fun_args = NULL,
                            negate_neg_strand_values=FALSE) {
  if( class(anno) == 'GRanges' ) {
    regions <- anno
  } else if ( is.character(class(anno)) ) {
    tryCatch(
      regions <- RMetaTools::meta_regions(anno, anchor, upstream, downstream),
      error = function(c) {
        c$message <- paste0("Error while trying to create metagene regions from bed file: ", anno, '\n', c$message)
        stop(c)
      }
    )
  } else {
    stop(paste0('failed to interpret regions', anno))
  }

  tryCatch(
    mat <- RMetaTools::get_matrix(bw_plus, bw_minus, regions, upstream, downstream, window_size,
                                  collapse_fun, collapse_fun_args, negate_neg_strand_values),
    error = function(c) {
      c$message <- paste0("Error while trying to create  metagene matrix.\n", c$message)
      stop(c)
    }
  )

  tryCatch(
    {
      mat <- RMetaTools::mat_to_tbl(mat)
    },
    error = function(c) {
      c$message <- paste0("Error while converting raw metagene matrix to tbl.\n", c$message)
      stop(c)
    }
  )

  mat
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
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000)
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 1000, 50)
#' tidy_meta <- mat_to_tbl(mat)
#' plot_htmp(tidy_meta, color_by='value', do_interpolate=FALSE)
#'
#' #log2 values incl pseudocount
#' tidy_meta %<>% mutate(log2_value = log2(value+1))
#' plot_htmp(tidy_meta, color_by='log2_value', do_interpolate=FALSE)
#'
#' @export
plot_htmp <- function(tbl, color_by='value', do_interpolate=FALSE){
  ggplot(tbl, aes_string(x='relpos', y='name', fill=color_by)) +
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
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000)
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 1000, 50)
#' tidy_meta <- mat_to_tbl(mat)
#' tidy_avgs <- meta_average(tidy_meta)
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
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000)
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 1000, 50)
#' tidy_meta <- mat_to_tbl(mat)
#' tidy_avgs <- meta_average(tidy_meta)
#' plot_profile(tidy_avgs)
#'
#' @export
plot_profile <- function(meta_tbl, color_by=NULL){
  ggplot(meta_tbl) +
    geom_ribbon(aes_string(x='relpos', ymin='conf.low', ymax='conf.high', fill=color_by), alpha=.2) +
    geom_line(aes_string(x='relpos', y='estimate', color=color_by)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(panel.grid = element_blank())
}
