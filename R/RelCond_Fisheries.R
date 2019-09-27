#' Relative condition factor
#'
#' Compute relative condition factor, used to compare condition of individuals of different lengths in fisheries. Kn= W/W', with W' predicted weight.  
#' 
#' @param pop_ID Optional ID 
#' @param pop_length Vector of length (in cm)
#' @param pop_weight Vector of weigth (in g)
#' @param size_limits Vector (length(size_limits) == 2) of min and max size of individuals included in the calculation. Default=NULL to use all individuals.
#' @param outliers Index of any point that should be discarded from computation
#' 
#' @export
#' @examples Compute_FultonK(pop_length = FSAdata::YPerchGL$fl/10, pop_weight = FSAdata::YPerchGL$w, bin_limits = c(0,7.5,10,17.5))
#'
#' @keywords fisheries
#' @keywords relative condition
#'
#'

RelCond_Fisheries <- function(pop_ID=NULL, pop_length, pop_weight, size_limits=NULL, outliers = NULL) {
  # Create list to save outputs
  list_out <- list()
  
  # Make sure we're working with numeric variables
  pop_weight = as.numeric(paste(pop_weight))
  pop_length = as.numeric(paste(pop_length))
  
  
  # Edit size_limits if entered in a wrong format
  if(length(size_limits)%%2 != 0) size_limits=NULL
  if(is.null(size_limits)) size_limits = c(min(pop_length, na.rm=T), max(pop_length, na.rm=T))
  size_limits <- size_limits[order(size_limits, decreasing = F)]
  which_size_limits <- which(pop_length>=size_limits[1] & pop_length<=size_limits[2])
  
  # warning message if wrong length of input vector
  if (length(pop_weight) != length(pop_length)) stop("\n You do not have the same number of observations for weight and length. Review your data and run the function again.")
  if (!is.null(pop_ID) | length(pop_weight) != length(pop_ID)) pop_ID=1:length(pop_weight)
  
  
  # Remove outliers, then create output dataframe
  if(!is.null(outliers) & is.numeric(outliers)) {
    pop_ID     <- pop_ID[-outliers]
    pop_length <- pop_length[-outliers]
    pop_weight <- pop_weight[-outliers]
  }
  
  pop_ID     <- pop_ID[which_size_limits]
  pop_length <- pop_length[which_size_limits]
  pop_weight <- pop_weight[which_size_limits]
  
  list_out$Data <- data.frame(pop_ID, pop_length,pop_weight)
  
  
  # Compute optimal weight for population
  lm_temp = lm(log10(pop_weight)~log10(pop_length))
  a_temp = 10^(lm_temp$coefficients[1])
  b_temp = lm_temp$coefficients[2]
  
  # Save output parameters
  list_out$a <- a_temp
  list_out$b <- b_temp
  
  # Compute Kn= W/W'
  list_out$Data$Wprim <- a_temp*list_out$Data$pop_length^b_temp
  list_out$Data$Kn <- list_out$Data$pop_weight / list_out$Data$Wprim
  
  # create plot from which a and b parameter were extracted
  p0 <- ggplot(list_out$Data, aes(pop_length,pop_weight)) + 
    geom_point() +
    labs(title="Weight:length relationship",
         x ="Length (cm))", y = "Weight (g)") +
    stat_smooth() + 
    theme_bw()
  
  # create plot from which a and b parameter were extracted
  p1 <- ggplot(list_out$Data, aes(log10(pop_length),log10(pop_weight))) + 
    geom_point() +
    labs(title="weight:length relationship (log10 scale)",
         subtitle = paste0("For this sample, a= ",signif(a_temp,3), ", b= ", signif(b_temp,3)),
         x ="log10(length)", y = "log10(weight)") +
    stat_smooth(method="lm") + 
    theme_bw()
  
  # create plot of Kn vs. length
  p2 <- ggplot(list_out$Data, aes(pop_length,Kn)) + 
    geom_point() +
    labs(title="Relative condition factor Kn",
         subtitle = bquote(K[n]*"=W/W', with W'=a"*L^b),
         x ="Length (cm)", y = "Kn") +
    stat_smooth(method="lm", col="coral") +
    geom_hline(yintercept=1) + 
    theme_bw()
  
  list_out$plot_WL <- p0
  list_out$plot_WL_log <- p1
  list_out$plot_KnL <- p2
  
  return(list_out)
}
