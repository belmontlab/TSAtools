% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bw_smthr_zoo.R
\name{pred_r_squared}
\alias{pred_r_squared}
\title{Predictive R-squared}
\usage{
pred_r_squared(linear.model)
}
\arguments{
\item{linear.model}{A linear regression model (class 'lm'). Required.

Use anova() to get the sum of squares for the linear model
Calculate the total sum of squares}
}
\description{
returns the predictive r-squared. Requires the function PRESS(), which returns
             the PRESS statistic.
}
\author{
Thomas Hopper
}
