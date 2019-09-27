#' Relative weight
#'
#' Compute relative weight, used to compare a population to an ideal population.  Wr = W/Ws * 100
#'
#' @param pop_ID Optional ID
#' @param pop_length Vector of length (in cm)
#' @param pop_weight Vector of weigth (in g)
#' @param a 'a' parameter reported for this species (e.g., data from fishbase)
#' @param b 'b' parameter reported for this species (e.g., data from fishbase)
#' @param size_limits Vector (length(size_limits) == 2) of min and max size of individuals included in the calculation. Default=NULL to use all individuals.
#' @param outliers Index of any point that should be discarded from computation
#' @param source_data Optional. Create here a character string with a reference, or the link to the webpage from which a and b parameters were scrapped, e.g., fishbase
#'
#' @export
#' @examples RelWeight_Fisheries(pop_length = FSAdata::YPerchTL$length/10, pop_weight = FSAdata::YPerchTL$weight, a = 0.01230, b = 3.04, source_data = "https://www.fishbase.in/summary/Perca-flavescens.html")
#'
#' @keywords fisheries
#' @keywords relative condition
#'
#'

RelWeight_Fisheries <- function(pop_ID=NULL, pop_length, pop_weight, a, b, size_limits=NULL, outliers = NULL, source_data=NULL) {
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
  if(is.null(a)|is.null(b)|!is.numeric(a)|!is.numeric(b)) stop("\n Make sure you entered numeric values for a and b. The function needs a and b values to compute relative weight. Find origin of a and b in Murphy et al (1990, 1991), and parameters on fishbase or other online repositories.")


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


  # Save output parameters
  list_out$a <- a
  list_out$b <- b
  list_out$source_a_b_param <- source_data

  # Compute Wr=W/Ws*100
  list_out$Data$Ws <- a*(list_out$Data$pop_length^b)
  list_out$Data$Wr <- list_out$Data$pop_weight / list_out$Data$Ws #*100
  # no longer multiplying by 100 because it brings my ratio to 100? Shoudn't it be close to 1?

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
         subtitle = paste0("For this sample, a= ",signif(a,3), ", b= ", signif(b,3)),
         x ="log10(length)", y = "log10(weight)") +
    stat_smooth(method="lm") +
    theme_bw()

  # create plot of Wr vs. length
  p2 <- ggplot(list_out$Data, aes(pop_length,Wr)) +
    geom_point() +
    stat_smooth(method="lm") +
    labs(title="Relative weight Wr",
         subtitle = bquote(W[r]*"=W/"*W[s]*", with "*W[s]*"=a"*L^b*",a and b calculated with RPL method (Murphy et al, 1991)"),
         x ="Length (cm)", y = "Wr") +
    geom_hline(yintercept=1) +
    theme_bw()

  # Save plots in output
  list_out$plot_WL <- p0
  list_out$plot_WL_log <- p1
  list_out$plot_WrL <- p2

  return(list_out)
}
