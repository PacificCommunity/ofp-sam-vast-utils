% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thicken.poly.r
\name{thicken.poly}
\alias{thicken.poly}
\title{Function to make polygon "thicker", useful when transforming to a different crs}
\usage{
thicken.poly(coords, n)
}
\arguments{
\item{coords}{matrix of x and y coordinates}

\item{n}{number of points to "thicken" each segment by}
}
\value{
Matrix of two columns: x and y
}
\description{
Function to make polygon "thicker", useful when transforming to a different crs
}
