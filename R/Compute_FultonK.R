#' Fulton's condition factor
#'
#' Compute Fulton's condition factor, fisheries. K= W/L^3 * constant.
#'
#' @param pop_ID Optional ID
#' @param pop_length Vector of length (in cm)
#' @param pop_weight Vector of weigth (in g)
#' @param bin_limits Vector of numbers for bins creation. If left NULL, use the min and max reported length
#' @param outliers Index of any point that should be discarded from computation
#'
#' @export
#' @examples Compute_FultonK(pop_length = FSAdata::YPerchTL$length/10, pop_weight = FSAdata::YPerchTL$weight, bin_limits = c(0,7.5,10,17.5))
#'
#' @keywords fisheries
#' @keywords Fulton's condition factor
#'

Compute_FultonK <- function(pop_ID=NULL, pop_length, pop_weight, bin_limits=NULL, outliers = NULL) {
  # packages needed
  require(ggplot2)
  require(wesanderson)

  # Create list to save outputs
  list_out <- list()
  if(is.null(bin_limits)) bin_limits=c(min(pop_length,na.rm=T), max(pop_length,na.rm=T))

  # Color output
  colrs <- wes_palette("Darjeeling1", length(bin_limits)-1, type="continuous")

  # Make sure we're working with numeric variables
  pop_weight = as.numeric(paste(pop_weight))
  pop_length = as.numeric(paste(pop_length))

  # Order bins
  bin_limits <- bin_limits[order(bin_limits, decreasing = F)]

  # warning message if wrong length of input vector
  if (length(pop_weight) != length(pop_length)) stop("\n You do not have the same number of observations for weight and length. Review your data and run the function again.")
  if (is.null(pop_ID) | length(pop_weight) != length(pop_ID)) pop_ID=1:length(pop_weight)

  # Save input data
  list_out$Data <- data.frame(pop_ID, pop_length,pop_weight)

  # create empty column in the dataset for the k of the size class
  list_out$Data$size_group <- rep(NA, nrow(list_out$Data))
  list_out$Data$FultonK <- rep(NA, nrow(list_out$Data))
  list_out$Data$delta_to_K <- rep(NA, nrow(list_out$Data))

  # create empty matrix to save output
  param <- as.data.frame(matrix(rep(NA,(length(bin_limits)-1)*6), ncol=6))
  colnames(param) = c("Min_size", "Max_size", "a", "b", "meanK", "slopeK")

  # create empty plot to see the points for different bins
  p1 <- ggplot(list_out$Data, aes(pop_length,pop_weight)) + labs(title=paste0("Weight:length relationship (n=", length(bin_limits)-1, ")"),
                                                                 x ="Length (mm)", y = "Weigth (g)") + theme_bw()
  p2 <- ggplot(list_out$Data, aes(pop_length,pop_weight)) + labs(title=paste0("Weight:length relationship (log-transformed) (n=", length(bin_limits)-1, ")"),
                                                                 x ="log(length)", y = "log(weight)") + theme_bw()
  p3 <- ggplot() + labs(title=paste0("Fulton's factor K (n=", length(bin_limits)-1, ")"),
                        x ="Length (cm)", y = "K") + theme_bw()


  # Beginning of the loop to calculate Fulton's K per size class
  # We're going to generate a and b for each size class, then calculate the average condition from that.
  for (i in 2:length(bin_limits)) {
    which_in_bin <- which(pop_length>bin_limits[i-1] & pop_length<=bin_limits[i])
    if (any(which(is.na(pop_length)) %in% which_in_bin)) which_in_bin <- which_in_bin[-which(which_in_bin %in% which(is.na(pop_length)))]
    if (any(which(is.na(pop_weight)) %in% which_in_bin)) which_in_bin <- which_in_bin[-which(which_in_bin %in% which(is.na(pop_weight)))]
    if (any(outliers %in% which_in_bin)) which_in_bin <- which_in_bin[-which(which_in_bin %in% outliers)]
    if (length(which_in_bin)>1) {
      lm_temp = lm(log10(pop_weight[which_in_bin])~log10(pop_length[which_in_bin]))
      a_temp = 10^(lm_temp$coefficients[1])
      b_temp = lm_temp$coefficients[2]

      # Now that we have and b, compute the fulton factor
      list_out$Data[which_in_bin,"size_group"] <- LETTERS[i-1]
      list_out$Data[which_in_bin,"FultonK"]    <- (list_out$Data$pop_weight/((list_out$Data$pop_length/10)^3) * a_temp)[which_in_bin]
      # Normalizing
      list_out$Data[which_in_bin,"FultonK"] <- (list_out$Data[which_in_bin,"FultonK"]-min(list_out$Data[which_in_bin,"FultonK"]))/(max(list_out$Data[which_in_bin,"FultonK"])-min(list_out$Data[which_in_bin,"FultonK"]))+0.5
      K_temp <- mean(list_out$Data[which_in_bin,"FultonK"])
      # Delta to mean K
      list_out$Data[which_in_bin,"delta_to_K"] <- list_out$Data[which_in_bin,"FultonK"] - K_temp

      # Add points to plot
      p1 <- p1 + geom_point(data =list_out$Data[which_in_bin,],mapping = aes(x=pop_length, y=pop_weight), colour=colrs[i-1])
      p2 <- p2 + geom_point(data =list_out$Data[which_in_bin,],mapping = aes(x=log10(pop_length), y=log10(pop_weight)), colour=colrs[i-1])
      p3 <- p3 + geom_point(data =list_out$Data[which_in_bin,],mapping = aes(x=pop_length, y=FultonK), colour=colrs[i-1], alpha=.4) + stat_smooth(data =list_out$Data[which_in_bin,],aes(x=pop_length, y=FultonK), method="lm", se = FALSE, colour="white", lwd=1.5) + stat_smooth(data =list_out$Data[which_in_bin,],aes(x=pop_length, y=FultonK), method="lm", se = FALSE, colour=colrs[i-1])

      # Save output
      lm_temp <- lm(FultonK ~ pop_length,data=list_out$Data[which_in_bin,])
      param[i-1,] <- c(bin_limits[i-1], bin_limits[i], a_temp, b_temp, K_temp, lm_temp$coefficients[2])
    } else { param[i-1,] <- c(bin_limits[i-1], bin_limits[i], NA,NA,NA,NA) }


  }

  # Save number of outputs
  list_out$Growth_parameters <- param[!is.na(param[,3]),]

  aplot <- ggplot(param, aes((param$Min_size+param$Max_size)/2, a)) + geom_point() + stat_smooth() +xlab("Average size (cm)")
  bplot <- ggplot(param, aes((param$Min_size+param$Max_size)/2, b)) + geom_point() + stat_smooth() +xlab("Average size (cm)")
  list_out$plot_param_a <- aplot
  list_out$plot_param_b <- bplot

  list_out$plot_WL <- p1
  list_out$plot_WL_log <- p2
  list_out$plot_KL <- p3

  return(list_out)
}
