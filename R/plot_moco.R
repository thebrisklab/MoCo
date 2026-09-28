#' Plot parcel-wise MoCo estimates on a cortical surface
#' 
#' Maps group-specific estimates and significant adjusted contrasts to a CIFTI
#' surface parcellation. This optional helper requires the `ciftiTools` package
#' and Connectome Workbench.
#' 
#' @param result A result returned by [moco()].
#' @param wb_path Path to the Connectome Workbench executable.
#' @param template_path Path to a template CIFTI file.
#' @param parcellation A parcellation recognized by
#'   [ciftiTools::load_parc()]. The network-collapsing logic currently supports
#'   `"Yeo_7"`.
#' @param plots_path Directory in which PNG files are written.
#' @param zlim_est Numeric length-two color limits for group-specific means.
#' @param zlim_association Numeric length-two color limits for contrasts.
#' @param fwer_index Row of `result$significant_regions` to plot when several
#'   FWER levels were requested.
#' @return Invisibly returns the three output file paths.
#' 
#' @export

plot_moco <- function(result, wb_path, template_path, 
                      parcellation = "Yeo_7", 
                      plots_path,
                      zlim_est = c(-0.3, 0.3), 
                      zlim_association = c(-0.1, 0.1),
                      fwer_index = 1L){
  
  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("plot_moco() requires the optional package 'ciftiTools'.")
  }
  dir.create(plots_path, recursive = TRUE, showWarnings = FALSE)

  ciftiTools::ciftiTools.setOption('wb_path', wb_path) 
  
  # read in the parcellation
  parc <- ciftiTools::load_parc(parcellation)
  # convert the parcellation to a vector
  parc_vec <- c(as.matrix(parc)) 
  # change the index of components to network
  parcel_label <- rownames(parc$meta$cifti$labels$parcels)[-1]
  
  if(parcellation == "Yeo_7"){
    parc_vec[parc_vec %in% which(grepl("Vis", parcel_label))] <- 1
    parc_vec[parc_vec %in% which(grepl("SomMot", parcel_label))] <- 2
    parc_vec[parc_vec %in% which(grepl("DorsAttn", parcel_label))] <- 3
    parc_vec[parc_vec %in% which(grepl("Limbic", parcel_label))] <- 4
    parc_vec[parc_vec %in% which(grepl("SalVentAttn", parcel_label))] <- 5
    parc_vec[parc_vec %in% which(grepl("Cont", parcel_label))] <- 6
    parc_vec[parc_vec %in% which(grepl("Default", parcel_label))] <- 7
  }
  
  # retrive results
  est_A0 <- as.matrix(result$est)[1,]
  est_A1 <- as.matrix(result$est)[2,]
  adj_association <- result$adj_association
  significant_regions <- result$significant_regions
  if (is.matrix(significant_regions)) {
    if (fwer_index < 1L || fwer_index > nrow(significant_regions)) {
      stop("fwer_index is outside the available significance rows.")
    }
    significant_regions <- significant_regions[fwer_index, ]
  }
  
  # read in the template xii file
  xii <- ciftiTools::read_xifti(template_path)
  # replace the medial wall vertices in xii with NA
  xii <- ciftiTools::move_from_mwall(xii, NA)
  # visualize the results by creating a new xifti file
  # convert parc_vec (the parcel key of each vertex) 
  # to xii_seed (seed correlation of the parcel each voxel belongs to)
  xii_A0 <- c(NA, est_A0)[parc_vec + 1]
  xii_A1 <- c(NA, est_A1)[parc_vec + 1]
  xii_sig <- c(NA, ifelse(significant_regions, adj_association, NA))[parc_vec + 1]
  # initialize a one-column xifti and replace its data 
  xii1 <- ciftiTools::select_xifti(xii, 1)
  xii_A0 <- ciftiTools::newdata_xifti(xii1, xii_A0)
  xii_A1 <- ciftiTools::newdata_xifti(xii1, xii_A1)
  xii_sig <- ciftiTools::newdata_xifti(xii1, xii_sig)
  
  # plot 
  output_files <- file.path(
    plots_path,
    c("MoCo_FC_A0.png", "MoCo_FC_A1.png", "Sig_regions.png")
  )
  plot(xii_A0, # title = "Motion-controlled Mean Functional Connectivity A0", 
       fname = output_files[1L], zlim = zlim_est)
  plot(xii_A1, # title = "Motion-controlled Mean Functional Connectivity A1", 
       fname = output_files[2L], zlim = zlim_est)
  plot(xii_sig, # title = "Significant Motion-controlled Association", 
       fname = output_files[3L], zlim = zlim_association)

  invisible(output_files)
}
