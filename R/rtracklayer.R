#' LoadBigWigSet
#'
#' Loads a specified range from specified bigwigs files.
#'
#' @param file files to load in dataframe with metadata and factors for grouping.
#' @param region named list with a chr, start, end and strand slot
#'
#' @details Loads the values from the region from the file(s) and collects everything in a tidy dataframe.
#'   Files are a dataframe or named list containing at a minimum these slots: 'path' 'filename' and named category
#'    levels
#'
#' @return A tidy data.frame with columns for chr start end strand and all levels specified the files.
#'
#' @examples
#'
#' path <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/bw'
#' files <- create_bw_file_set(path, c('rep3', '.bw$'), c('_N20', '_3D12'), c('siRNA', 'ab', 'rep'), '_')
#' region <- list(chr='chr1', start=1000000, end=1001000, strand='+')
#' set <- load_bw_set(files, region)
#' ggplot(set, aes(x=starts, y=scores, color=rep)) + geom_line() + facet_grid(siRNA~ab)
#'
#' @export
load_bw_set <- function(files, region){
  chunks <- GRanges(region$chr, IRanges(region$start, region$end))


  d <- apply(files, 1, function(s) {
    import(paste(s['path'], s['file'], sep='/'), which = chunks)
  })

  len <- sapply(d, length)

  df <- data.frame(starts = as.numeric(sapply(d, function(x) start(ranges(x)))),
             scores = as.numeric(sapply(d, score)))

  flevels <- names(files)[!(names(files) %in% c('path','file'))]
  for(fl in flevels) {
    df[,fl] <- rep(files[,fl],len)
  }

  df
}



#' CreateBigWigSampleSet
#'
#' Creates a dataframe holding information on bigwig files.
#'
#' @param path path where the .bw files are located.
#' @param regex a regular expression (or list of) for selection files in path, default = '.bw'
#' @param junk_list a list of substrings to be removed before splitting to levels.
#' @param split_to split filenames into levels
#' @param split_by delimiter for splitting, default='_'
#'
#'
#' @details Creates a data frame with path, filename. Uses the substrings in junk to clean-up the filename and the  simplified filename is used with split_to create levels using tidyr::separate function and split_by is the separator used. By default regex is always also removed from filename during simplification.
#'
#' @return A tidy data.frame with columns for path filename and all levels.
#'
#' @examples
#'
#' path <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/bw'
#' create_bw_file_set(path, c('_rep[3,4]', '.bw$'), c('_N20', '_3D12'), c('siRNA', 'ab', 'rep'), '_')
#'
#' @export
create_bw_file_set <- function(path, regex, junk_list, split_to, split_by){

  simple_name <- function(name) {
    sname <- name
    for(re in regex){
      sname <- sub(re, '', name)
    }
    for(f in junk_list) {
      sname <- sub(f, '', sname)
    }
    sname
  }

  bw_files <- dir(path)
  for(re in regex){
    bw_files <- bw_files[grep(re, bw_files)]
  }

  data.frame(path=rep(path, length(bw_files)),
             file = bw_files,
             sample = sapply(bw_files, simple_name)) %>%
  tidyr::separate(sample, split_to, split_by)
}



#' getRegionForGene
#'
#' gets region for a specified gene.
#'
#' @param anno a data.frame holding the annotations
#' @param id id of the gene. Default =''.
#' @param gene_name alternative if id is not provided. Default =''.
#' @param anchor region around gene, TSS or TES. Default ='gene'.
#' @param left_offset extend by this upstream start or anchor
#' @param right_offset extend by this downstream end or anchor
#'
#' @details Retrieves the information in a named vector.
#'
#' @return A named vector.
#'
#' @examples
#'
#' anno_gencode <- read.table('/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/bw/gencode_v23_hg19_genes.bed', col.names = c('chr', 'start', 'end', 'id', 'type', 'strand', 'name'))
#' region <- get_region_for_gene(anno_gencode, gene_name='GAPDH', anchor='TSS', 2000, 5000)
#' files <- create_bw_file_set(path, c('[y,I]','rep[3,4]', '.bw$'), c('_N20', '_3D12'), c('siRNA', 'ab', 'rep'), '_')

#' set <- load_bw_set(files, region)
#' ggplot(set, aes(x=starts, y=scores, color=rep)) + geom_line() + facet_grid(siRNA~ab)
#'
#' region <- get_region_for_gene(anno_gencode, gene_name='HIST1H4C', anchor='gene', left_offset=2000, right_offset=5000)
#' files <- create_bw_file_set(path, c('[y,I]','rep[3,4]', '.bw$'), c('_N20', '_3D12'), c('siRNA', 'ab', 'rep'), '_')

#' set <- load_bw_set(files, region)
#' ggplot(set, aes(x=starts, y=scores, color=rep)) + geom_line() + facet_grid(siRNA~ab)
#' ggplot(set, aes(x=starts, y=scores, color=siRNA)) + geom_line() + facet_grid(.~ab)
#'
#' scaled_set <- set %>% group_by(siRNA, ab, rep) %>% mutate(scores=scale(scores))
#' ggplot(scaled_set, aes(x=starts, y=scores, color=siRNA)) + geom_line() + facet_grid(.~ab)
#'
#' @export
get_region_for_gene <- function(anno, gene_id='', gene_name='', anchor='gene', left_offset=0, right_offset=0){

  if (id != '') {
    region <- dplyr::filter(anno, id==gene_id)
  } else if (gene_name != '') {
    region <- dplyr::filter(anno, name==gene_name)
  } else {
    print('No id or gene_name provided')
    return()
  }

  print(nrow(region))

  if (anchor == 'TSS') {
    region$end = region$start
  } else if (anchor == 'TES') {
    region$start = region$end
  }

  if (region$strand == '+') {
    region$start = region$start - left_offset
    region$end = region$end + right_offset
  } else {
    region$start = region$start - right_offset
    region$end = region$end + left_offset
  }

  region
}
