#' Boxplot of percent ESS
#'
#' Plot method for objects of class 'causal_estimates' to generate a boxplot of percent sample effective sample size (ESS) for a specific weighting method using bootstrap samples.
#'
#' @param x An object of class 'causal_estimates', the number of bootstrap samples \code{B} must be greater than 0.
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
#'
#' @examples
#' data(demo)
#' set.seed(1)
#' output2 <- causal.estimate(S, Z, X, Y, B = 5, method = "IC", naturalGroupProp)
#' plot(output2)
#'
#' @export
#'
plot.causal_estimates = function(x, ...){
  if(is.null(x$collatedESS)){
    stop("Error: To draw a boxplot, the number of bootstrap samples 'B' must be greater than 0.")
  }

  args = list(...)
  # Set default values
  y_limit = if ("y_limit" %in% names(args)) args$y_limit else c(0, 50)
  color = if ("color" %in% names(args)) args$color else "red"

  ESS.mt = x$collatedESS

  # ESS.mt is a B by 1 matrix
  # It contains the ESS of each bootstrap sample
  B = length(ESS.mt)

  df = data.frame(as.vector(ESS.mt), method=rep(x$method, each=B))
  names(df) = c("ESS","method")

  plot1 = ggplot(df, aes(x=method, y=ESS)) +
    geom_boxplot(fill = color) +
    ylab("Percent ESS") +
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) +
    aes(x = fct_inorder(method)) + xlab("Method") +
    ggtitle("Effective sample sizes") + theme(plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) +
    ylim(y_limit[1], y_limit[2]) #+ theme(legend.position = "none")


  #graphics.off()

  plot1
}

