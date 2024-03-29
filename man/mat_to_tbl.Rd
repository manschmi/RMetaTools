% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeMatrix.R
\name{mat_to_tbl}
\alias{mat_to_tbl}
\title{Tidy Metagene Matrix}
\usage{
mat_to_tbl(mat)
}
\arguments{
\item{mat}{metagene matrix created from get_matrix}
}
\value{
matrix with columns from most upstream to most downstream and rows are the individual regions.
}
\description{
Converts metagene matrix into a tidy tbl.
}
\details{
Converts metagene matrix to tbl using rownames and colnames information to assign a column *name* and *rel_pos* respectively.
}
\examples{

bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
regions <- meta_regions(bedfile, 'TSS', 1000, 5000)
bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
mat <- get_matrix(bw_plus, bw_minus, regions, 1000, 5000, 50)
tidy_meta <- mat_to_tbl(mat)

}
