#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

##this is a script that parses command line options and creates a metagene matrix that is saved in an .RData file

## argument list


if (args[1] == '-h' | args[1] == '--help') {
   cat(
'This is a script that parses command line options and creates a metagene matrix and
saves the matrix as .RData file based using RMetaTools.
(ie This script depends on RMetaTools and its dependencies are installed in working environment!)

Required Arguments:
  -region=100:TSS:100/10  .... ie 100bp upstream to 100bp downstream of the TSS as 10bp windows
  -bwpath=bw_path ...comma-separated list of bigwig files
  -plus_bw_regex=KD1*plus*.bw,KD2*plus*.bw
  -bed=bed_path/anno.bed ... path to bedfile
  -out=my_RData_name.RData

Optional Arguments:
  -minus_bw_regex=KD1*plus*.bw,KD2*plus*.bw ... optional, bw_plus will be used for minus strand if missing.
  -collapse_fun=rowSums ... optional, how are values from single positions averaged into windows (default=rowMeans).
  -negate_neg_strand_values ... optional, use when minus strand bigwigs values are negatives.

Example Usage:
  Rscript metagene.R -region=1000:TSS:4000/250 -bwpath=~/data/myproject/bigwigs -plus_bw_regex=_plus.bw -minus_bw_regex=_minus.bw -bed=~/annotations/genes.bed -out=~/results/projectbw_rel_genes_1000TSS4000by250.RData')
    quit()
} else {
  args=sapply(args, function(arg) strsplit(arg, split = '='))
  argnames = lapply(args, function(arg) gsub('^-', '', arg[1]))
  argvalues = sapply(args, function(arg) arg[2])
  names(argvalues) <- argnames
  args <- argvalues
}

if( !('region' %in% names(args)) ) {
  stop('required arg -region not found')
}
region <- args['region']
window_size <- as.integer(sub('.*/', '', region))
region <- strsplit(sub('/.*', '', region), ':')[[1]]
upstream <- as.integer(region[1])
anchor <- region[2]
downstream <- as.integer(region[3])

if( !('bed' %in% names(args)) ) {
  stop('required arg -bed not found')
}
bed_fname <- args['bed']
if( dir.exists(bed_fname) ){
  stop(paste0('bedfilename ', bed_fname, ' not found!'))
}

if( !('bwpath' %in% names(args)) ) {
  stop('required arg -bwpath not found')
}
bwpath <- args['bwpath']

if( !('plus_bw_regex' %in% names(args)) ) {
  stop('required arg -plus_bw_regex not found')
}
plus_bw_regex <- args['plus_bw_regex']
plus_bws <- dir(bwpath)[grepl(plus_bw_regex, dir(bwpath))]

if( 'minus_bw_regex' %in% names(args) ) {
  minus_bw_regex <- args['minus_bw_regex']
  minus_bws <- sub(plus_bw_regex, minus_bw_regex, plus_bws)
}else{
  message('minus bws not found, assuming unstranded data from plus_bw_regex files')
  minus_bws <- plus_bws
}

if( !('out' %in% names(args)) ) {
  stop('required arg -out not found')
}
out_fname <- args['out']

collapse_fun <- rowMeans
if( 'collapse_fun' %in% names(args) ) {
  collapse_fun <- args['collapse_fun']
  print(paste0('UPS: using alternative window average fun: ', collapse_fun))
}

negate_neg_strand_values <- FALSE
if( 'negate_neg_strand_values' %in% names(args) ) {
  print('UPS: negating neg strand values!')
  negate_neg_strand_values <- TRUE
}

print(paste0('rel: ', anchor, '   upstream: ', upstream, '   downstream: ', downstream, '   win: ', window_size))

tryCatch(suppressPackageStartupMessages(library(RMetaTools)),
         error=function(e) print("Error while loading required package RMetaTools or some if its dependencies!!"))

print('preparing regions')
regions <- RMetaTools::meta_regions(bed_fname, anchor, upstream, downstream)


print('starting to compute individual matrices')
mat_list <- list()

for(i in seq_along(plus_bws)){
  cat(paste0('   for bigwigs: ', bwpath, ' > ', plus_bws[i], '\t', minus_bws[i], '\n'))
  mat_list[[i]] <- RMetaTools::metagene_matrix(bw_plus=paste0(bwpath, plus_bws[i]),
                                          bw_minus=paste0(bwpath, minus_bws[i]),
                                          regions, anchor, upstream, downstream, window_size,
                                          collapse_fun, negate_neg_strand_values) %>%
    mutate(fname = plus_bws[i])
}

print(paste0('saving combined matrix to ', out_fname))
mat <- bind_rows(mat_list)
save(mat, file=out_fname)
