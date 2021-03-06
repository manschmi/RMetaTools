% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deeptools_funs.R
\name{sort_matrix}
\alias{sort_matrix}
\title{Sort Matrix}
\usage{
sort_matrix(matrix_file, bed_file)
}
\arguments{
\item{matrix_file}{filename for the matrix (gzipped).}

\item{bed_file}{filename for the bed file used for sorting}
}
\value{
The sorted matrix.
}
\description{
Sorts a deeptools matrix file using a bed file.
}
\details{
Sorts the rows of the matrix without header.
}
\examples{

 matrix_name <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/eRNA_tsspm5000_sorted_by_INT_binding_joined_sensitivity.gz'
 bed_name <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Shiekkatar_INT_ChIP/eRNAs_sorted_by_INTS_binding.bed'

 m_sorted <- sort_matrix(matrix_name, bed_name)

 write.table(m_sorted, file='/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/eRNA_tsspm5000_sorted_by_INT_binding_joined_sensitivity_sorted', row.names = FALSE, col.names = FALSE, quote = F, dec = '.')


}
