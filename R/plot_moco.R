#' Plot Motion-Controlled Mean Functional Connectivity and Associations
#' 
#' This function plots motion-controlled mean functional connectivity example results using the provided data and parameters.
#' 
#' @param result A list containing the results of motion-controlled functional connectivity analysis.
#' @param wb_path The path to the workbench binary for ciftiTools.
#' @param template_path The path to the template file for the visualization.
#' @param parcellation The parcellation used in the analysis.
#' @param plots_path The path where the plots will be saved.
#' @param zlim_est The z-axis limits for functional connectivity plots.
#' @param zlim_association The z-axis limits for association plots.
#' @return This function does not return any value. It generates and saves plots based on the provided parameters.
#' 
#' @export

plot_moco <- function(result, wb_path, template_path, 
                      parcellation = "Yeo_7", 
                      plots_path,
                      zlim_est = c(-0.3, 0.3), 
                      zlim_association = c(-0.1, 0.1)){
  
  # set the path for ciftiTools workbench
  ciftiTools.setOption('wb_path', wb_path) 
  
  # read in the parcellation
  parc <- load_parc(parcellation)
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
  
  # read in the template xii file
  xii <- read_xifti(template_path)
  # replace the medial wall vertices in xii with NA
  xii <- move_from_mwall(xii, NA)
  # visualize the results by creating a new xifti file
  # convert parc_vec (the parcel key of each vertex) 
  # to xii_seed (seed correlation of the parcel each voxel belongs to)
  xii_A0 <- c(NA, est_A0)[parc_vec + 1]
  xii_A1 <- c(NA, est_A1)[parc_vec + 1]
  xii_sig <- c(NA, ifelse(significant_regions, adj_association, NA))[parc_vec + 1]
  # initialize a one-column xifti and replace its data 
  xii1 <- select_xifti(xii, 1)
  xii_A0 <- newdata_xifti(xii1, xii_A0)
  xii_A1 <- newdata_xifti(xii1, xii_A1)
  xii_sig <- newdata_xifti(xii1, xii_sig)
  
  # plot 
  setwd(plots_path)
  plot(xii_A0, # title = "Motion-controlled Mean Functional Connectivity A0", 
       fname = "MoCo_FC_A0.png", zlim = zlim_est)
  plot(xii_A1, # title = "Motion-controlled Mean Functional Connectivity A1", 
       fname = "MoCo_FC_A1.png", zlim = zlim_est)
  plot(xii_sig, # title = "Significant Motion-controlled Association", 
       fname = "Sig_regions.png", zlim = zlim_association)
}


