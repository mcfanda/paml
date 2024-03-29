% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{as.data.frame.pamlpowermix}
\alias{as.data.frame.pamlpowermix}
\title{Print results of paml power analysis}
\usage{
\method{as.data.frame}{pamlpowermix}(obj)
}
\arguments{
\item{obj}{of class \code{\link[paml]{pamlpowermix}}}
}
\value{
an object of class \code{\link[base]{data.frame}}
}
\description{
Prints in a nice tabular shape the results of
\code{\link[paml]{modelpower}}) applied to a mixed model.
}
\examples{
 methods(class="pamlpowermix")
}
\author{
Marcello Gallucci
}
