% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPlotdata.R
\name{getPlotdata}
\alias{getPlotdata}
\title{Utility function: extracting data from a \code{\link{dc_CA}} object for 
plotting a single axis by your own code or \code{\link{plot.dcca}}.}
\usage{
getPlotdata(
  x,
  axis = 1,
  envfactor = NULL,
  traitfactor = NULL,
  newnames = NULL,
  facet = TRUE,
  remove_centroids = FALSE
)
}
\arguments{
\item{x}{results from \code{\link{dc_CA}} of class \code{dcca}.}

\item{axis}{the axis number to get (default 1).}

\item{envfactor}{name of row factor to display as color and lines in the CWM
plot (default \code{NULL}). The default extracts the factor from the 
environmental model. If set to \code{NA}, no additional coloring and lines
are displayed in \code{\link{plot.dcca}}. The parameter sets the 
\code{groups} variable in the \code{CWM_SNC} data frame of the return 
value/in the plot.}

\item{traitfactor}{name of column factor to display as color and lines in
the SNC plot (default \code{NULL}). The default extracts the factor from
the trait model. If set to \code{NA}, no additional coloring and lines are
displayed in \code{\link{plot.dcca}}. The parameter sets the \code{groups} 
variable in the \code{CWM_SNC} data frame of the return value/in the plot.}

\item{newnames}{a list  with two elements: names for traits and for 
environmental variables, default \code{NULL} for names derived from the 
result of \code{\link{scores.dcca}} with \code{tidy = TRUE}.}

\item{facet}{logical. Default \code{TRUE} for CWMs and SNCs plots in 
separate panels. If \code{FALSE}, this parameter changes the position of 
the environmental centroid names (from left to right).}

\item{remove_centroids}{logical to remove any centroids from the plot data 
(default \code{FALSE}). Can be a two-vector, \emph{e.g.} 
\code{c(TRUE, FALSE)} to remove only the environmental centroids.}
}
\value{
A list with three components
\describe{
\item{CWM_SNC}{a data.frame containing plot data}
\item{trait_env_scores}{a vector of scores per trait/environment}
\item{newNameList}{a vector of new names to be used in the plot}
}
}
\description{
\code{getPlotdata} extracts data from a \code{\link{dc_CA}} object for 
plotting the CWMs and SNCs of a single axis.
}
\examples{
data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

# must delete "Sites" from response matrix or data frame
Y <- dune_trait_env$comm[, -1] # must delete "Sites"

out <- dc_CA(formulaEnv = ~ A1 + Moist + Use + Manure + Mag,
                   formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                   response = Y,
                   dataEnv = dune_trait_env$envir,
                   dataTraits = dune_trait_env$traits,
                   verbose = FALSE)
dat <- getPlotdata(out)
names(dat)
names(dat$CWM_SNC)
levels(dat$CWM_SNC$groups)

plot(out)
}
