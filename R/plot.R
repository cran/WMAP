#' Plot method for objects of class 'causal_estimates'
#'
#' @param x An object of class 'causal_estimates'.
#' @param ... Additional arguments including:
#' \describe{
#'   \item{y_limit}{The y-axis range. Default is \code{c(0,50)}.}
#'   \item{color}{The boxplot color. Default is "red".}
#' }
#' @return A boxplot of percent sample ESS for a specific weighting method (\code{FLEXOR}, \code{IC}, or \code{IGO})
#'
#' @importFrom ggplot2 aes theme element_text ggtitle xlab ylab ylim ggplot geom_boxplot
#' @importFrom grDevices graphics.off
#' @importFrom forcats fct_inorder

#' @export
#'
plot.causal_estimates = function(x, ...){

  args = list(...)
  # Set default values
  y_limit = if ("y_limit" %in% names(args)) args$y_limit else c(0, 50)
  color = if ("color" %in% names(args)) args$color else "red"


  ESS.mt = x$collatedESS

  # ESS.mt is a B by 1 matrix
  # It contains the ESS of each bootstrap sample

  df = data.frame(as.vector(ESS.mt), method=rep(x$method, each=B))
  names(df) = c("ESS","method")

  plot1 = ggplot(df, aes(x=method, y=ESS)) +
    geom_boxplot(fill = color) +
    ylab("Percent ESS") +
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) +
    aes(x = fct_inorder(method)) + xlab("Method") +
    ggtitle("Effective sample sizes") + theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) +
    ylim(y_limit[1], y_limit[2]) #+ theme(legend.position = "none")
  # color palettes given by scale_fill_brewer or scale_colour_hue
  # See http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #plot1 <- plot1 + scale_fill_brewer(palette="Set1") #scale_colour_hue(c=10, l=100)

  graphics.off()

  plot1
}

