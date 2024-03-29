% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtracklayer.R
\name{load_bw_set}
\alias{load_bw_set}
\title{LoadBigWigSet}
\usage{
load_bw_set(files, region)
}
\arguments{
\item{region}{named list with a chr, start, end and strand slot}

\item{file}{files to load in dataframe with metadata and factors for grouping.}
}
\value{
A tidy data.frame with columns for chr start end strand and all levels specified the files.
}
\description{
Loads a specified range from specified bigwigs files.
}
\details{
Loads the values from the region from the file(s) and collects everything in a tidy dataframe.
  Files are a dataframe or named list containing at a minimum these slots: 'path' 'filename' and named category
   levels
}
\examples{

path <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/bw'
files <- create_bw_file_set(path, c('rep3', '.bw$'), c('_N20', '_3D12'), c('siRNA', 'ab', 'rep'), '_')
region <- list(chr='chr1', start=1000000, end=1001000, strand='+')
set <- load_bw_set(files, region)
ggplot(set, aes(x=starts, y=scores, color=rep)) + geom_line() + facet_grid(siRNA~ab)

}
