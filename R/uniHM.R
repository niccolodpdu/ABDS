#' @name redblue
#' @title Generate a colormap with shades of red and blue colors
#' @description The colormap is used to create a gradient from the blue color to the red color. 
#' The function takes an optional argument m which specifies the number of colors in the colormap.
#'
#' @param m an integer, number of colors in the colormap
#' @return a matrix with three columns (r, g, b) representing the color map
#' @export
redblue <- function(m = nrow(colormap())) {
            if (m %% 2 == 0) {
                        # From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
                        m1 <- m * 0.5
                        r <- (0:(m1 - 1))/max(m1 - 1, 1)
                        g <- r
                        r <- c(r, rep(1, m1))
                        g <- c(g, rev(g))
                        b <- rev(r)
            } else {
                        # From [0 0 1] to [1 1 1] to [1 0 0];
                        m1 <- floor(m * 0.5)
                        r <- (0:(m1 - 1))/max(m1, 1)
                        g <- r
                        r <- c(r, rep(1, m1 + 1))
                        g <- c(g, 1, rev(g))
                        b <- rev(r)
            }
            cbind(r, g, b)
}


#' @name heatmap_data_modification
#' @title Function to modify data as per new heatmap design
#' @description
#' This function modifies the input data for visualization in a heatmap by performing scaling, log transformation and centering.
#' 
#' @param tot_sg1 A numeric matrix containing the input data to be modified.
#'
#' @return A modified numeric matrix containing the scaled, log-transformed, and centered data.
#' @export
heatmap_data_modification <- function(tot_sg1) {
            # Calculate row sums and factor for scaling
            tot_sg_sum <- rowSums(tot_sg1)
            factor <- 1000000 / tot_sg_sum
            
            # Scale data and set values less than 1 to 1
            tot_sg_scale <- tot_sg1 * factor
            tot_sg_scale[tot_sg_scale < 1] <- 1
            
            # Log transform scaled data and calculate mean and standard deviation
            tot_sg_scale_log <- log(tot_sg_scale)
            mean_gb <- mean(as.matrix(tot_sg_scale_log))
            std_gb <- sd(as.matrix(tot_sg_scale_log, na.rm = TRUE))
            
            # Center data around mean and standard deviation
            tot_sg_scale_log_center <- (tot_sg_scale_log - mean_gb) / std_gb
            
            # Return modified data
            return(tot_sg_scale_log_center)
}


