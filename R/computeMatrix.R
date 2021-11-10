#' Metagene Regions
#'
#' Creates regions for metagene analysis.
#'
#' @param bed filename for bed file.
#' @param anchor anchor point TSS, TES, center or NULL (default='TSS')
#' @param upstream bp upstream of anchor (default=1000)
#' @param downstream bp downstream of anchor (default=1000)
#'
#' @details Loads the bed file. If anchor is NULL, returns region from TSS-upstream to TES+downstream. If anchor, extracts the position of the anchor point and creates ranges from x bp upstream to x bp downstream of the anchor point.
#'
#' @return GRanges with regions.
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' meta_regions(bedfile, 'TSS', 1000, 5000)
#'
#' @export
meta_regions <- function(bed,
                        anchor = 'TSS',
                        upstream = 1000,
                        downstream = 1000) {

  if( class(bed) == 'GRanges' ) {
     anno <- bed
  } else if ( class(bed) == 'character' ) {
    tryCatch(
      anno <- rtracklayer::import(bed),
      error = function(c) {
        c$message <- paste0("Error while trying to create metagene regions from bed file: ", anno, '\n', c$message)
        stop(c)
      }
    )
  } else {
    stop(paste0('failed to interpret regions', anno))
  }

  anno$original_start <- start(anno)
  anno$original_end <- end(anno)
  if( !('name' %in% colnames(mcols(anno))) ) {
    anno$name <- paste0('region', seq(length(anno)))
  }

  strands <- as.vector(strand(anno))

  if(is.null(anchor)){
    start(anno) <- ifelse(strand(anno) == '+', start(anno) - upstream, start(anno) - downstream)
    end(anno) <- ifelse(strand(anno) == '+', end(anno) + downstream, end(anno) + upstream) #0 counts as first position
  }else{
    anchor_position <- case_when(
      (anchor == 'TSS' & strands == '-') ~ end(anno),
      anchor == 'TSS' ~ start(anno),
      anchor == 'TES' & strands == '-' ~ start(anno),
      anchor == 'TES' ~ end(anno),
      anchor == 'center' ~ as.integer(start(anno) + (end(anno)-start(anno))/2)
    )
    start(anno) <- ifelse(strand(anno) == '+', anchor_position - upstream, anchor_position - downstream)
    end(anno) <- ifelse(strand(anno) == '+', anchor_position + downstream, anchor_position + upstream) #0 counts as first position
  }

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
#'
#' @details Loads the bed file, extracts the position of the anchor point and creates ranges from x bp upstream to x bp downstream. I upstream and downstream are both NULL, collects signal within region specified in anno.
#'
#' @return matrix with columns from most upstream to most downstream and rows are the individual regions from specified strand. rownames are set as name column from anno. colnames are set to relative position of each nt
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000)
#' bw <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' mat <- get_single_strand_matrix(bw=bw, anno=regions, upstream=1000, downstream=5000)
#'
#' @export
get_single_strand_matrix <- function(bw, anno, upstream=1000, downstream=1000) {

  ##handling specific strand and order to ensure correct output naming
  seqlevels(anno) <- sort(seqlevels(anno))
  anno <- sort(anno)

  ##handling uneven region sizes without defined window size
  if(is.null(upstream) & is.null(downstream)) {
    message('using scale region mode, upstream and downstream arguments are ignored')
    upstream <- 0
    downstream <- max(width(anno))
  }

  ##ensure to query only chrs present in bigwig,
  ## ie otherwise import(bw, ...) will return in error
  bw_chrs <- seqlevels(BigWigFile(bw))
  contained_anno_rows <- which(as.character(seqnames(anno)) %in% bw_chrs)
  if (length(contained_anno_rows) != length(anno)) {
    warning(paste0('ignoring chromosomes from annotation file not present in the bigwig file for strand: ', strand,
                   ': ', paste0(seqlevels(anno[-contained_anno_rows]))))
  }
  mat <- matrix(NA, nrow = length(anno), ncol = upstream+downstream+1)
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
  colnames(mat) <- seq(-upstream, downstream)

  mat
}



#' Collapse Matrix To Windows by Window Size
#'
#' Averages Matrix Columns into Windows using an Averaging Function amd a defined window size
#'
#' @param mat matrix
#' @param collapse_fun function for averaging values inside windows (Default = rowMeans, see Details).
#' @param window_size size for binning in bp (default=10, ie average over 10nt)
#'
#' @details for a matrix with ncol=100 and window_size=10, will return a matrix with ncol=10. collapse_fun is a function that functions as ie rowMeans (the default) or rowSums.
#'
#' @return matrix with ncol(input matrix)/window_size columns
#'
#' @examples
#'
#' mat <- matrix(1:1000, nrow=10, ncol=100)
#' matc <- collapse_to_n_windows(anno, mat, mean, 10)
#' dim(matc)
#'
#' log2RowMeans <- function(mat, pseudocount=1){apply(mat,1,function(row) mean(log2(row+pseudocount)))}
#' matc <- collapse_to_window_size(mat, log2RowMeans, 10)
#' dim(matc)
#'
#' @export
collapse_to_window_size <- function(mat, collapse_fun = rowMeans, window_size = 10) {
  if (window_size > 1) {
    size <- ncol(mat)
    if( size %% window_size != 0) {
      wins <- seq(1,window_size*round(size/window_size,0), window_size)
    }else{
      wins <- seq(1,size, window_size)
    }

    matc <- lapply(wins, function(i) collapse_fun(mat[,i:(i+window_size-1)])) %>%
      bind_cols %>%
      as.matrix

    colnames(matc) <- sapply(wins, function(i) colnames(mat)[i])
    rownames(matc) <- rownames(mat)

    mat <- matc
  }

  mat
}


#' Collapse Matrix To Windows by Window Number
#'
#' Averages Matrix Columns into Windows using an Averaging Function amd a defined number of windows
#'
#' @param mat matrix
#' @param collapse_fun function for averaging values inside windows (Default = rowMeans, see Details).
#' @param n_windows number of windows (default=10, ie collapse to 10 values eah)
#'
#' @details for a matrix with ncol=100 or 50 and n_windows=10, will return a matrix with ncol=10 in both cases. collapse_fun is a function that functions as ie rowMeans (the default) or rowSums.
#'
#' @return matrix with n_windows columns
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "snRNA_MS31_hg38.bed", package = "RMetaTools")
#' regions <- rtracklayer::import(bedfile)
#' mat <- get_single_strand_matrix(bw=bw, anno=regions, upstream=NULL, downstream=NULL)
#' dim(mat)
#' matc <- collapse_to_n_windows(anno, mat, mean, 10)
#' dim(matc)
#'
#' @export
collapse_to_n_windows <- function(anno, mat, collapse_fun = mean, n_windows = 10) {

  matc <- t(sapply(seq_along(anno), function(i) {
      windows <- as.integer(seq(1,ncol(mat)+1,length.out = n_windows+1))

      window_starts <- windows[1:length(windows)-1]
      window_ends <- windows[2:length(windows)]-1
      wins <- data.frame(s=window_starts, e=window_ends)

      apply(wins, 1, function(win) collapse_fun(mat[i,win['s']:win['e']]))
    }))

  colnames(matc) <- 1:n_windows
  rownames(matc) <- rownames(mat)

  matc
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
#' @param collapse_to_win_size size for binning in bp (default=1, ie no binning)
#' @param collapse_to_n_wins number of bins to collapse to (default=NULL, ie no binning)
#' @param collapse_avg_fun function for averaging values inside windows (Default = rowMeans, see Details).
#' @param negate_neg_strand_values negate values from minus strand bigwig (default=FALSE).
#' @details Loads the bed file. Anchored mode requires a anchor, upstream and downstream, extracts the position of the anchor point and
#'   creates ranges from x bp upstream to x bp downstream. I anchor, upstream and downstream are NULL, collects values over specified entire regions.
#'   If collapse_to_win_size is set to a integer > 0, uses collapse_avg_fun (rowMeans, rowSums or similar) to overage signals within this range See \code{\link{collapse_to_window_size}} for more info about averaging over window size.
#'   If collapse_to_n_wins is not NULL, ups: this requires collapse_avg_fun to be set to mean, sum or similar (ie instead of the default rowMeans) to overage signals into x equally sized bins. Only used if collapse_to_win_size is not > 1. See \code{\link{collapse_to_n_windows}} for more info about averaging into defined number of windows.
#'   Negate_neg_strand_values can be set to TRUE to deal with minus strand bigwigs which contain all
#'   negative values as used by some labs.
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
#' dim(mat)
#' mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 5000, 50)
#' dim(mat)
#' mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 5000, collapse_to_n_wins=5, collapse_avg_fun=mean)
#' dim(mat)
#'
#' #if bigwigs are not stranded, ie from ChIP experiment both plus and minus are the same file
#' bw <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_hg38.bw", package = "RMetaTools")
#' mat <- get_matrix(bw, bw, regions, 1000, 5000, 50)
#'
#'
#' @export
get_matrix <- function(bw_plus, bw_minus, anno,
                       upstream=1000, downstream=1000,
                       collapse_to_win_size = NULL,
                       collapse_to_n_wins = NULL,
                       collapse_avg_fun = rowMeans,
                       negate_neg_strand_values=FALSE) {
  smat <- get_single_strand_matrix(bw_plus, anno[strand(anno)=="+"], upstream, downstream)
  asmat <- get_single_strand_matrix(bw_minus, anno[strand(anno)=="-"], upstream, downstream)
  asmat <- asmat[,ncol(asmat):1]
  if (negate_neg_strand_values) {
    asmat <- -asmat
  }
  if(ncol(smat) > ncol(asmat)){
    asmat <- cbind(asmat, matrix(NA, nrow=nrow(asmat), ncol=(ncol(smat)-ncol(asmat)) ))
  }else if(ncol(smat) < ncol(asmat)){
    smat <- cbind(smat, matrix(NA, nrow=nrow(smat), ncol=(ncol(asmat)-ncol(smat)) ))
  }
  mat <- rbind(smat, asmat)
  if(!is.null(collapse_to_win_size)) {
    mat <- collapse_to_window_size(mat, collapse_avg_fun, collapse_to_win_size)
  }else if(!is.null(collapse_to_n_wins)){
    mat <- collapse_to_n_windows(anno, mat, collapse_avg_fun, collapse_to_n_wins)
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
#' mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 5000, 50)
#' tidy_meta <- mat_to_tbl(mat)
#'
#' @export
mat_to_tbl <- function(mat) {

  nm <- mat %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = 'name') %>%
    tidyr::gather(., rel_pos, value, -name) %>%
    mutate(rel_pos = as.numeric(rel_pos)) %>%
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
#' @param anchor 'TSS', 'TES', 'center' or NULL (default='TSS')
#' @param upstream bp upstream of anchor (default=1000)
#' @param downstream bp downstream of anchor (default=1000)
#' @param collapse_to_win_size size for binning in bp (default=1, ie no binning)
#' @param collapse_to_n_wins number of bins to collapse to (default=NULL, ie no binning)
#' @param collapse_avg_fun function for averaging values inside windows (Default = rowMeans, see Details).
#' @param negate_neg_strand_values negate values from minus strand bigwig (default=FALSE).
#'
#' @details Loads the bed file. If anchor is given, extracts the position of the anchor point and creates ranges from x bp upstream to x bp downstream. If anchor is NULL, collects signal in specified region from TSS-upstream to TES+downstream. Queries from both bigwigs and the correct strand. Summarizes signal of window_size bins and returns a tidy data.frame. negate_neg_strand_values can be set to TRUE to deal with minus strand bigwigs which contain all negative values as used by some labs.
#' When using a GRanges object for anno, this assumes the bigwigs will be queries from start to end! All ranges must be exactly same size!
#'
#' @return Tidy data frame with columns: gene (name, ie col4 from *bed* file), rel_pos (position in bp relative to *anchor*), value (averaged value from bigwig file).
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' anno=bedfile;anchor='center';upstream=1000;downstream=1000;window_size=10
#' collapse_fun = rowMeans;negate_neg_strand_values=FALSE
#' metamat <- metagene_matrix(bw_plus, bw_minus, bedfile, anchor, upstream, downstream, window_size)
#' metamat <- metagene_matrix(bw_plus, bw_minus, bedfile, anchor, upstream, downstream, collapse_to_n_wins=5, collapse_avg_fun=mean)
#' bedfile <- system.file("extdata", "snRNA_MS31_hg38.bed", package = "RMetaTools")
#' metamat <- metagene_matrix(bw_plus, bw_minus, bedfile, anchor, upstream, downstream, window_size)
#' metamat <- metagene_matrix(bw_plus, bw_minus, bedfile, anchor=NULL, 10, 10, collapse_to_n_wins=5, collapse_avg_fun=mean)
#'
#' # for pre-loaded regions (useful for ie modifying regions by hand)
#' # this most useful
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000)
#' start(regions) <- start(regions) + 1000
#' end(regions) <- end(regions) - 1000
#' regions <- regions[width(regions)>0,]
#' metamat <- metagene_matrix(bw_plus, bw_minus, regions, anchor=NULL, upstream=0, downstream=0, collapse_to_n_wins=5, collapse_avg_fun=mean)
#' metamat <- metagene_matrix(bw_plus, bw_minus, regions, anchor='TSS', upstream=100, downstream=100, collapse_to_win_size=50, collapse_avg_fun=rowMeans)
#' metamat <- metagene_matrix(bw_plus, bw_minus, bedfile, anchor='TSS', upstream=100, downstream=100, collapse_to_win_size=50, collapse_avg_fun=rowMeans)
#'
#' @export
metagene_matrix <- function(bw_plus, bw_minus, anno,
                            anchor='TSS', upstream=1000, downstream=1000,
                            collapse_to_win_size = NULL,
                            collapse_to_n_wins = NULL,
                            collapse_avg_fun = rowMeans,
                            negate_neg_strand_values=FALSE) {

  meta_args <- c(anchor, upstream, downstream)
  names(meta_args) <- c('anchor', 'upstream', 'downstream')

  if (!is.null(collapse_to_win_size)) {meta_args['collapse_to_win_size'] = collapse_to_win_size}
  if (!is.null(collapse_to_n_wins)) {meta_args['collapse_to_n_wins'] = collapse_to_n_wins}


  tryCatch(
      regions <- meta_regions(anno, anchor, upstream, downstream),
      error = function(c) {
        c$message <- paste0("Error while trying to create metagene regions \n", c$message)
        stop(c)
      })

  if(is.null(anchor)){
    upstream <- NULL
    downstream <- NULL
  }

  tryCatch(
    mat <- get_matrix(bw_plus, bw_minus, regions,
                      upstream, downstream,
                      collapse_to_win_size, collapse_to_n_wins,
                      collapse_avg_fun, negate_neg_strand_values),
    error = function(c) {
      c$message <- paste0("Error while trying to create  metagene matrix.\n", c$message)
      stop(c)
    }
  )

  tryCatch(
    {
      mat <- mat_to_tbl(mat)
    },
    error = function(c) {
      c$message <- paste0("Error while converting raw metagene matrix to tbl.\n", c$message)
      stop(c)
    }
  )

  for(i in seq_along(meta_args)){
    attr(mat, names(meta_args)[i]) <- meta_args[i]
  }

  mat
}



#' Compute Metagene Values For Region Around A Reference-point
#'
#' Computes a metagene object from scratch, similar to metagene_matrix but as more compact form.
#'
#' @param bw_plus bigwig filename for plus strand signal
#' @param bw_minus bigwig filename for minus strand signal
#' @param anno Bed filename for annotations to use OR a GRanges object (see Details).
#' @param anchor 'TSS', 'TES', 'center' or NULL (default='TSS')
#' @param upstream bp upstream of anchor (default=1000)
#' @param downstream bp downstream of anchor (default=1000)
#' @param bin_size size for binning in bp (default=1, ie no binning)
#' @param collapse_avg_fun function for averaging values inside windows (Default = rowMeans, see Details).
#' @param negate_neg_strand_values negate values from minus strand bigwig (default=FALSE).
#'
#' @details Loads the bed file. Extracts the position of the anchor point and
#'   creates ranges from x bp upstream to x bp downstream.
#'   Queries from both bigwigs and the correct strand. Summarizes signal of
#'   bin_size bins.
#'   negate_neg_strand_values can be set to TRUE to deal with minus strand bigwigs which contain all
#'   negative values as used by some labs. When using a GRanges object for anno,
#'   this assumes the bigwigs will be queries from start to end! All ranges must
#'   be exactly same size! anchor, upstream and downstream info is still reguired.
#'
#' @return Metagene object
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' anno=bedfile;anchor='center';upstream=1000;downstream=1000;window_size=10
#' collapse_fun = rowMeans;negate_neg_strand_values=FALSE
#' m <- metagene_rel_anchor(bw_plus, bw_minus, bedfile, anchor, upstream, downstream, window_size)
#' m
#' @export
metagene_rel_anchor <- function(bw_plus,
                     bw_minus,
                     anno,
                     anchor='TSS', upstream=1000, downstream=1000,
                     bin_size = 100,
                     collapse_avg_fun = rowMeans,
                     negate_neg_strand_values=FALSE,
                     sample_name = NULL,
                     group_name = NULL) {

  s3 <- structure(list(), class = "metagene")
  s3$header <- c(anchor, upstream, downstream, bin_size)
  names(s3$header) <- c('ref point', 'upstream', 'downstream', 'bin size')

  tryCatch(
    regions <- meta_regions(anno, anchor, upstream, downstream),
    error = function(c) {
      c$message <- paste0("Error while trying to create metagene regions \n", c$message)
      stop(c)
    })

  if(is.null(sample_name)){
    sample_name <- sub('plus', '', bw_plus)
  }

  if(is.null(group_name)){
    group_name <-'genes'
  }

  tryCatch(
    mat <- get_matrix(bw_plus, bw_minus, regions,
                      upstream, downstream,
                      collapse_to_win_size = bin_size,
                      collapse_to_n_wins = NULL,
                      collapse_avg_fun, negate_neg_strand_values),
    error = function(c) {
      c$message <- paste0("Error while trying to create  metagene matrix.\n", c$message)
      stop(c)
    }
  )

  tryCatch(
    {
      s3$samples <- sample_name
      s3$groups <- group_name

      s3$rel_pos <- colnames(mat)
      colnames(mat) <- NULL

      mat_region_ids <- sapply(rownames(mat), function(s) str_split(s, '\t', simplify = T)[4])
      names(mat_region_ids) <- NULL

      s3$regions <- list(
        as_tibble(
          regions[order(factor(regions$name, levels = mat_region_ids)),]) %>%
          dplyr::rename(id=name))

      names(s3$regions) <- group_name
      rownames(mat) <- NULL

      s3$matrices <- list(mat)
    },
    error = function(c) {
      c$message <- paste0("Error while converting raw metagene matrix to tbl.\n", c$message)
      stop(c)
    }
  )

  s3
}
