#' @name    biplots
#' @aliases biplots
#' @title   Biplots
#' @description biplots
#'
#'
#' @param model Output from `bayes_ammi()`. This should contain the results of the Bayesian AMMI model, including all sampled iterations.
#' @param burnin Numeric. Percentage of iterations to discard as burn-in to avoid the effects of random initializations during sampling. For example, `burnin = 0.1` removes the first 10\% of iterations.
#' @param thin Numeric. Proportion of sampled iterations to retain for analysis. For example, `thin = 0.2` keeps 20\% of the iterations, selecting 1 out of every 5 iterations.
#' @param pb Numeric. Significance levels for the contours in the plot. Smaller values of `pb` result in wider contours, while higher values create smaller, more specific contours.
#' @param plot_stable Logical. If `TRUE`, stable instances are highlighted in the output plot.
#' @param plot_unstable Logical. If `TRUE`, unstable instances are highlighted in the output plot.
#' @param ncolors Integer. Specifies the number of distinct colors to use in the plot. Adjust this to control the visual differentiation of elements in the plot.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{plot}{A plot displaying the contours and final biplot values.}
#'   \item{contour_data}{A `data.frame` containing the data used to create the contours.}
#'   \item{biplot_data}{A `data.frame` containing the data used to recreate the final biplot values.}
#' }
#'
#' @author
#' \enumerate{
#'     \item Julian Garcia Abadillo Velasco (\email{garciaabadillo.j@ufl.edu})
#'     \item Diego Jarquin (\email{diego.jarquin@gmail.com})
#'    }
#'
#' @references
#'  Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J.M., Viele, K., Liu, G., and Cornelius, P.L. (2011)
#'  Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model
#'  \emph{Crop Science}, \strong{51}, 1458–1469.
#'  (\href{https://acsess.onlinelibrary.wiley.com/doi/abs/10.2135/cropsci2010.06.0343}{doi: 10.2135/cropsci2010.06.0343})
#'
#'
#' @import ggplot2
#' @import lme4
#' @import mvtnorm
#' @import rlang
#' @import rstiefel
#' @import scales
#' @import ks
#' @import purrr
#' @import ggrepel
#' @import tidyr
#' @import tibble
#' @import tmvtnorm
#' @importFrom  dplyr select group_by count ungroup mutate n summarise bind_rows
#' @importFrom magrittr %>%
#' @importFrom MASS Null mvrnorm
#' @importFrom stats rexp rgamma rnorm runif
#'
#' @export
#'
#' @examples
#'
#' data(Maiz)
#'
#' fm1 <-
#'   bayes_ammi(
#'     .data = Maiz,
#'     .y = y,
#'     .gen = entry,
#'     .env = site,
#'     .rep = rep,
#'     .nIter = 200
#'   )
#' names(fm1)
#' fm1$mu1
#' fm1$tau1
#' fm1$tao1
#' fm1$delta1
#' fm1$lambdas1
#' fm1$alphas1
#' fm1$gammas1
#'
#' library(ggplot2)
#'
#' output_05 <- biplots(model = fm1, plot_stable = TRUE, plot_unstable = TRUE, pb = 0.05)
#' output_05$biplot
#'
#' output_95 <- biplots(model = fm1, plot_stable = TRUE, plot_unstable = TRUE, pb = 0.95)
#' output_95$biplot
#'
#'
# if(getRversion() >= "2.15.1"){
#   utils::globalVariables(
#     c(
#       "alphas0"
#       , "c0"
#       , "delta0"
#       , "gammas0"
#       , "lambdas0"
#       , "ge_means0"
#       , "mu0"
#       , "n0"
#       , "tao0"
#       , "tau0"
#     )
#   )
# }
#'
biplots <- function(model, burnin = 0.3, thin = 0.2, pb = 0.05, plot_stable = TRUE, plot_unstable = TRUE, ncolors = 5) {
  UseMethod("biplots")
}


#' @export
#' @rdname biplots

biplots.default <-
  function(model, # bayes_ammi object
           burnin = 0.3, # percentage of iterations to burn
           thin = 0.2, # percentage of interations to sample
           pb = 0.05, # significance level for contorns
           plot_stable = TRUE, # plot stable instances
           plot_unstable = TRUE, # plot unstable instances
           ncolors = 5 # how many different colors for plotting instances?
  ) {
    # process alphas1
    pre_alpha <- model$alphas1
    sub_alphas <- list()

    k <- ncol(pre_alpha) - 2
    r <- length(unique(pre_alpha[, 2]))
    for (i in 1:k + 2) {
      subdf <- as.data.frame(pre_alpha[, c(1, 2, i)]) %>%
        pivot_wider(names_from = 2, values_from = 3)
      colnames(subdf)[1:r + 1] <- sprintf("alpha%s_%s", colnames(subdf)[1:r + 1], i - 2)
      sub_alphas[[paste0(i)]] <- subdf[, -1]
    }

    alpha <- do.call("cbind", sub_alphas)
    colnames(alpha) <- do.call("rbind", strsplit(colnames(alpha), "\\."))[, 2]
    niter <- nrow(pre_alpha) / r
    burnin <- ceiling(niter * burnin)
    thin <- 1 / thin
    ind <- burnin + seq(1, niter - burnin, by = thin)
    alpha <- alpha[ind, ]

    # process gammas
    pre_gamma <- model$gammas1
    sub_gammas <- list()

    # k = ncol(pre_alpha) - 2 #no need to compute k again because should be the same
    c <- length(unique(pre_gamma[, 2]))
    for (i in 1:k + 2) {
      subdf <- as.data.frame(pre_gamma[, c(1, 2, i)]) %>%
        pivot_wider(names_from = 2, values_from = 3)
      colnames(subdf)[1:c + 1] <- sprintf("alpha%s_%s", colnames(subdf)[1:c + 1], i - 2)
      sub_gammas[[paste0(i)]] <- subdf[, -1]
    }

    gamma <- do.call("cbind", sub_gammas)
    colnames(gamma) <- do.call("rbind", strsplit(colnames(gamma), "\\."))[, 2]
    # niter = nrow(pre_gamma)/c # again no need to compute, should be the same
    gamma <- gamma[ind, ]

    svd.mcmc <- cbind(alpha, gamma)
    nsim <- nrow(svd.mcmc)
    s <- (niter - burnin) / thin

    cont <- function(x, y, pb, lod, cex = cex) {
      # d <- cbind(svd.mcmc[,73],svd.mcmc[,85]) %>%
      d <- cbind(x, y) %>%
        magrittr::set_colnames(c("x", "y")) %>%
        as_tibble()

      kd <- ks::kde(d, compute.cont = TRUE, h = 0.2)

      ## extract results
      get_contour <- function(kd_out = kd, prob = pb) {
        # print(prob)
        # prob = sprintf("%g%%", round(prob * 100, 2))
        contour_95 <- with(kd_out, contourLines(
          x = eval.points[[1]], y = eval.points[[2]],
          z = estimate, levels = cont[prob]
        )[[1]])
        as_tibble(contour_95) %>%
          mutate(prob = prob)
      }

      # dat_out <- map_dfr(c("5%"), ~get_contour(kd, .)) %>%
      dat_out <- map_dfr(c(paste(as.character(pb * 100), "%", sep = "")), ~ get_contour(kd, .)) %>%
        group_by(prob) %>%
        mutate(n_val = 1:n()) %>%
        ungroup()

      return(dat_out)
    }

    # save all conts for ggplot
    all_cont <- list()

    for (i in 1:c) {
      data_output <- cont(x = svd.mcmc[, 2 * r + i], y = svd.mcmc[, 2 * r + c + i], pb, lod = i, cex = 0.5)
      tmp <- data_output
      tmp$byplot <- "Environment"
      tmp$ID <- i
      # print(i)
      all_cont[[paste0("E", i)]] <- tmp
    }

    for (i in 1:r) {
      data_output <- cont(x = svd.mcmc[, i], y = svd.mcmc[, r + i], pb, lod = i, cex = 0.5)
      tmp <- data_output
      tmp$byplot <- "Genotype"
      tmp$ID <- i
      # print(i)
      all_cont[[paste0("G", i)]] <- tmp
    }

    all_cont <- do.call("rbind", all_cont)

    summary <- all_cont %>%
      group_by(byplot, ID) %>%
      summarise(x = mean(x), y = mean(y))

    label_stable <- function(data, summary, colx, coly, radius = 0) {
      data$stable <- FALSE
      summary$stable <- FALSE
      for (i in 1:nrow(summary)) {
        c.byplot <- summary$byplot[i]
        c.ID <- summary$ID[i]
        tmp <- filter(data, byplot == c.byplot, ID == c.ID)
        if (min(tmp[, colx]) <= radius & max(tmp[, colx]) >= radius & min(tmp[, coly]) <= radius & max(tmp[, coly]) >= radius) {
          data[data$byplot == c.byplot & data$ID == c.ID, "stable"] <- TRUE
          summary[i, "stable"] <- TRUE
        }
      }
      return(list(data = data, summary = summary))
    }

    # classify stable and non stable lines
    tmp <- label_stable(all_cont, summary, 2, 3)
    data.stable <- tmp$data
    summary.stable <- tmp$summary
    data.stable$name <- data.stable$ID
    summary.stable$name <- summary.stable$ID

    # define base plot limits
    x_lim <- max(abs(min(data.stable$x)), abs(max(data.stable$x)))
    y_lim <- max(abs(min(data.stable$y)), abs(max(data.stable$y)))
    c <- 1.01
    x_lim <- c * x_lim
    y_lim <- c * y_lim

    base_plot <- data.stable %>%
      ggplot(aes(x = x, y = y, group = ID, fill = as.factor(ID %% ncolors))) +
      xlim(c(-x_lim, x_lim)) +
      facet_grid(~byplot) +
      theme_classic() +
      geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
      geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
      xlab("") +
      ylab("") +
      theme(legend.position = "none") +
      theme(
        panel.background = element_rect(color = "black", fill = "transparent"),
        axis.line = element_blank(),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white", face = "bold", size = 10)
      )

    if (plot_stable) {
      base_plot <- base_plot +
        geom_polygon(
          data = filter(data.stable, stable),
          alpha = 0.2, color = "black"
        ) +
        geom_point(
          data = filter(summary.stable, stable),
          aes(x = x, y = y, fill = as.factor(ID %% ncolors)), shape = 21, color = "white", size = 2
        ) +
        geom_label_repel(
          data = filter(summary.stable, stable),
          aes(x = x, y = y, fill = as.factor(ID %% ncolors), label = name), color = "white",
          max.overlaps = nrow(summary)
        )
    }

    if (plot_unstable) {
      base_plot <- base_plot +
        geom_polygon(
          data = filter(data.stable, !stable),
          alpha = 0.2, color = "black"
        ) +
        geom_point(
          data = filter(summary.stable, !stable),
          aes(x = x, y = y, fill = as.factor(ID %% ncolors)), shape = 21, color = "white", size = 2
        ) +
        geom_label_repel(
          data = filter(summary.stable, !stable),
          aes(x = x, y = y, fill = as.factor(ID %% ncolors), label = name), color = "white",
          max.overlaps = nrow(summary)
        )
    }
    return(list(contorns = data.stable, biplot_data = summary.stable, biplot = base_plot))
  }
