% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setters.R
\name{setters}
\alias{setters}
\alias{.setParameters}
\alias{.setParameters,HiCaptuRe,list-method}
\alias{.setByBaits}
\alias{.setByBaits,HiCaptuRe,list-method}
\alias{.setByRegions}
\alias{.setByRegions,HiCaptuRe,list-method}
\title{Functions to set data held in a HiCaptuRe object.}
\usage{
.setParameters(x, y)

\S4method{.setParameters}{HiCaptuRe,list}(x, y)

.setByBaits(x, y)

\S4method{.setByBaits}{HiCaptuRe,list}(x, y)

.setByRegions(x, y)

\S4method{.setByRegions}{HiCaptuRe,list}(x, y)
}
\arguments{
\item{x}{A HiCaptuRe object}

\item{y}{Data to add (parameters, tibble with bait-centric information or GRanges with region-centric)}
}
\value{
A HiCaptuRe object
}
\description{
Use these functions to set data stored in each of the slots of a
HiCapture object.
}
\keyword{internal}
