% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{infer_RHS}
\alias{infer_RHS}
\alias{formals_from_expr}
\title{Utilities for formulas and graphics arguments}
\usage{
infer_RHS(ex)

formals_from_expr(ex, others = character(0))
}
\arguments{
\item{ex}{An expression as in \code{quote(x^2)} or the left side of a tilde expression.}

\item{others}{Character string(s) with names of additional arguments to be
included in the formals}
}
\description{
\code{infer_RHS} turns a one-sided formula into a two-sided formula suitable for \code{makeFun()}.
\code{formals_from_expression} creates a "formals" list for creating a function.
The list will have arguments in the canonical order.
}
