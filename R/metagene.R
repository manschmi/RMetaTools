#' Metagene Regions
#'
#' Creates regions for metagene analysis.
#'
#' @param bed filename for bed file.
#' @param anchor anchor point TSS, TES or center (default='TSS')
#' @param upstream bp upstream of anchor (default=1000)
#' @param downstream bp downstream of anchor (default=1000)
#' @param window_size size for binning in bp (default=1, ie no binning)
#'
#' @details Loads the bed file, extracts the position of the anchor point and creates ranges from x bp upstream to x bp downstream.
#'
#' @return GRanges with regions.
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000, 50)
#'
#' @export
meta_regions <- function(bed,
                        anchor = 'TSS',
                        upstream = 1000,
                        downstream = 1000,
                        window_size = 1) {
  print('loading file')
  anno <- rtracklayer::import(bed)

  print('get range')
  if(anchor == 'TSS'){
    ranges(anno) <- start(ranges(anno))
  }else if(anchor == 'TES'){
    ranges(anno) <- end(ranges(anno))
  }else if(anchor == 'center'){
    ends <- end(ranges(anno))
    starts <- start(ranges(anno))
    ranges(anno) <- starts + (end-starts)/2
  }else{
    stop('unknown anchor point, use TSS, TES or center')
  }
  start(ranges(anno)) <- start(ranges(anno)) - upstream
  end(ranges(anno)) <- end(ranges(anno)) + downstream - 1 #0 counts as first position

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
#' @param window_size size for binning in bp (default=1, ie no binning)
#' @param strand strand of anno to use (default='+')
#'
#' @details Loads the bed file, extracts the position of the anchor point and creates ranges from x bp upstream to x bp downstream.
#'
#' @return matrix with columns from most upstream to most downstream and rows are the individual regions from specified strand.
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000, 50)
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' mat <- get_single_strand_matrix(bw_plus, regions, 1000, 1000, 50, '+')
#'
#' @export
get_single_strand_matrix <- function(bw, anno, upstream=1000, downstream=1000, window=1, strand="+") {
  mat <- as.matrix(rtracklayer::import(bw, which=anno[strand(anno)==strand], as='NumericList'))
  size <- upstream+downstream-1
  if(strand=="+"){
    wins <- seq(1,size, window)
  }else if(strand=="-"){
    wins <- rev(seq(1,size, window))
  }
  lapply(wins, function(i) rowMeans(mat[,i:(i+window-1)])) %>% bind_cols
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
#' @param strand strand of anno to use (default='+')
#'
#' @details Loads the bed file, extracts the position of the anchor point and creates ranges from x bp upstream to x bp downstream.
#'
#' @return matrix with columns from most upstream to most downstream and rows are the individual regions.
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000, 50)
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 1000, 50)
#'
#' #if bigwigs are not stranded, ie from ChIP experiment both plus and minus are the same file
#' bw <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_hg38.bw", package = "RMetaTools")
#' mat <- get_matrix(bw, bw, regions, 1000, 1000, 50)
#'
#' @export
get_matrix <- function(bw_plus, bw_minus, anno, upstream=1000, downstream=1000, window_size = 1) {
  smat <- get_single_strand_matrix(bw_plus, anno, upstream, downstream, window_size, strand="+")
  asmat <- get_single_strand_matrix(bw_minus, anno, upstream, downstream, window_size, strand="-")
  bind_rows(smat, asmat)
}



#' Tidy Metagene Matrix
#'
#' Converts metagene matrix into a tidy tbl.
#'
#' @param mat metagene matrix created from get_matrix
#' @param anno GRanges object of regions to query.
#' @param anchor 'TSS', 'TES' or 'center.'
#' @param upstream bp upstream of anchor (default=1000)
#' @param downstream bp downstream of anchor (default=1000)
#' @param window_size size for binning in bp (default=1, ie no binning)
#' @param strand strand of anno to use (default='+')
#'
#' @details Loads the bed file, extracts the position of the anchor point and creates ranges from x bp upstream to x bp downstream.
#'
#' @return matrix with columns from most upstream to most downstream and rows are the individual regions.
#'
#' @examples
#'
#' bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000, 50)
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 1000, 50)
#' tidy_meta <- mat_to_tbl(mat)
#'
#' @export
mat_to_tbl <- function(mat, anno, anchor, upstream, downstream, window_size) {
  nanno <- length(anno)
  size <- upstream + downstream
  nbins <- size/window_size
  nm <- gather(mat, relpos, value)
  nm$gene <- rep(as.character(seq(1:nanno)),nbins)
  nm$relpos <- rep(seq(-upstream,(downstream-1),window_size),each=nanno)
  nm
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
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000, 50)
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
  ggplot(tbl, aes_string(x='relpos', y='gene', fill=color_by)) +
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
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000, 50)
#' bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
#' bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
#' mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 1000, 50)
#' tidy_meta <- mat_to_tbl(mat)
#' tidy_avgs <- meta_average(tidy_meta)
#'
#' @export
meta_average <- function(tbl){
  tbl %>%
    group_by_at(vars(-value, -gene)) %>%
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
#' regions <- meta_regions(bedfile, 'TSS', 1000, 5000, 50)
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



