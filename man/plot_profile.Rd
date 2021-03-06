% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metagene.R
\name{plot_profile}
\alias{plot_profile}
\title{Plot Metagene Average Profile}
\usage{
plot_profile(meta_tbl, color_by = NULL)
}
\arguments{
\item{color_by}{name of column from tbl to use for coloring of lines and ranges (default=NULL)}

\item{tbl}{tidy metagene tbl}
}
\value{
tbl
}
\description{
Plot Metagene Average Profile
}
\details{
plots line from ('estimate'), and range from conf.low and conf.high columns in tbl
}
\examples{

bedfile <- system.file("extdata", "Chen_PROMPT_TSSs_liftedTohg38.bed", package = "RMetaTools")
bw_plus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_plus_hg38.bw", package = "RMetaTools")
bw_minus <- system.file("extdata", "GSM1573841_mNET_8WG16_siLuc_minus_hg38.bw", package = "RMetaTools")
anno=bedfile;anchor='center';upstream=1000;downstream=1000;window_size=10
collapse_fun = rowMeans;negate_neg_strand_values=FALSE
metamat <- metagene_matrix(bw_plus, bw_minus, bedfile, anchor, upstream, downstream, window_size)
tidy_log2_avgs <- meta_average(mutate(metamat, value = log2(value+1)))
plot_profile(tidy_log2_avgs)

}
