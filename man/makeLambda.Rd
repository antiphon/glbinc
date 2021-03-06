% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make-lambda-vector.R
\name{makeLambda}
\alias{makeLambda}
\title{Create the lambda vector}
\usage{
makeLambda(nlambda = 100, lambda.min = 0.001, ..., lambda.max,
  log.scale = TRUE)
}
\arguments{
\item{nlambda}{default 100}

\item{lambda.min}{default 0.001}

\item{...}{passed on to max_penalty to compute the max lambda if lambda.max is missing}

\item{lambda.max}{no default}

\item{logscale}{use log-scale equidistant? Default: true}
}
\description{
Create the lambda vector
}
