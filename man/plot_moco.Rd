% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_moco.R
\name{plot_moco}
\alias{plot_moco}
\title{Plot Motion-Controlled Mean Functional Connectivity and Associations}
\usage{
plot_moco(
  result,
  wb_path,
  template_path,
  parcellation = "Yeo_7",
  plots_path,
  zlim_est = c(-0.3, 0.3),
  zlim_association = c(-0.1, 0.1)
)
}
\arguments{
\item{result}{A list containing the results of motion-controlled functional connectivity analysis.}

\item{wb_path}{The path to the workbench binary for ciftiTools.}

\item{template_path}{The path to the template file for the visualization.}

\item{parcellation}{The parcellation used in the analysis.}

\item{plots_path}{The path where the plots will be saved.}

\item{zlim_est}{The z-axis limits for functional connectivity plots.}

\item{zlim_association}{The z-axis limits for association plots.}
}
\value{
This function does not return any value. It generates and saves plots based on the provided parameters.
}
\description{
This function plots motion-controlled mean functional connectivity example results using the provided data and parameters.
}
