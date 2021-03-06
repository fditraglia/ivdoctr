% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rho_bounds.R
\name{candidate3}
\alias{candidate3}
\title{Evaluates the edge where r_TstarU is on the boundary.}
\usage{
candidate3(r_TstarU_lower, r_TstarU_upper, k_lower, k_upper, obs)
}
\arguments{
\item{r_TstarU_lower}{Vector of lower bounds of endogeneity}

\item{r_TstarU_upper}{Vector of upper bounds of endogeneity}

\item{k_lower}{Vector of lower bounds on measurement error}

\item{k_upper}{Vector of upper bounds on measurement error}

\item{obs}{Observables generated by get_observables}
}
\value{
List containing vector of lower bounds and vector of upper bounds of r_uz
}
\description{
Evaluates the edge where r_TstarU is on the boundary.
}
