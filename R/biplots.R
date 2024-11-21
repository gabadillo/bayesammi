library(bayesammi)
library(dplyr)
library(tidyr)
library(ks)
library(magrittr)
library(purrr)
library(ggplot2)
library(ggrepel)

biplots = function(model,                  # bayes_ammi object
                   burnin = 0.3,           # percentage of iterations to burn
                   thin = 0.2,             # percentage of interations to sample
                   pb = 0.05,              # significance level for contorns
                   plot_stable = TRUE ,    # plot stable instances
                   plot_unstable = TRUE,   # plot unstable instances
                   ncolors = 5){           # how many different colors for plotting instances?

  # process alphas1
  pre_alpha = model$alphas1
  sub_alphas = list()

  k = ncol(pre_alpha) - 2
  r = length(unique(pre_alpha[,2]))
  for (i in 1:k+2){
    subdf = as.data.frame(pre_alpha[,c(1,2,i)]) %>%
      pivot_wider(names_from = 2, values_from =  3)
    colnames(subdf)[1:r+1] = sprintf('alpha%s_%s', colnames(subdf)[1:r+1], i-2)
    sub_alphas[[paste0(i)]] = subdf[,-1]
  }

  alpha = do.call('cbind', sub_alphas)
  colnames(alpha) = do.call('rbind',strsplit(colnames(alpha), '\\.'))[,2]
  niter = nrow(pre_alpha)/r
  burnin = ceiling(niter*burnin)
  thin = 1/thin
  ind <- burnin + seq(1,niter-burnin,by=thin)
  alpha = alpha[ind,]

  # process gammas
  pre_gamma = model$gammas1
  sub_gammas = list()

  # k = ncol(pre_alpha) - 2 #no need to compute k again because should be the same
  c = length(unique(pre_gamma[,2]))
  for (i in 1:k+2){
    subdf = as.data.frame(pre_gamma[,c(1,2,i)]) %>%
      pivot_wider(names_from = 2, values_from =  3)
    colnames(subdf)[1:c+1] = sprintf('alpha%s_%s', colnames(subdf)[1:c+1], i-2)
    sub_gammas[[paste0(i)]] = subdf[,-1]
  }

  gamma = do.call('cbind', sub_gammas)
  colnames(gamma) = do.call('rbind',strsplit(colnames(gamma), '\\.'))[,2]
  #niter = nrow(pre_gamma)/c # again no need to compute, should be the same
  gamma = gamma[ind,]

  svd.mcmc = cbind(alpha,gamma)
  nsim <- nrow(svd.mcmc)
  s <- (niter-burnin)/thin

  cont <- function(x,y,pb,lod,cex=cex){
    #d <- cbind(svd.mcmc[,73],svd.mcmc[,85]) %>%
    d <- cbind(x,y) %>%
      magrittr::set_colnames(c("x", "y")) %>%
      as_tibble()

    kd <- ks::kde(d, compute.cont=TRUE, h=0.2)

    ## extract results
    get_contour <- function(kd_out=kd, prob=pb) {
      #print(prob)
      #prob = sprintf("%g%%", round(prob * 100, 2))
      contour_95 <- with(kd_out, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                              z=estimate, levels=cont[prob])[[1]])
      as_tibble(contour_95) %>%
        mutate(prob = prob)
    }

    #dat_out <- map_dfr(c("5%"), ~get_contour(kd, .)) %>%
    dat_out <- map_dfr(c(paste(as.character(pb*100),"%",sep='')), ~get_contour(kd, .)) %>%
      group_by(prob) %>%
      mutate(n_val = 1:n()) %>%
      ungroup()

    return(dat_out)
  }

  # save all conts for ggplot
  all_cont = list()

  for(i in 1:c){
    data_output <- cont(x=svd.mcmc[,2*r+i],y=svd.mcmc[,2*r+c+i],pb,lod=i,cex=0.5)
    tmp = data_output
    tmp$byplot = "Environment"
    tmp$ID = i
    #print(i)
    all_cont[[paste0("E",i)]] = tmp
  }

  for(i in 1:r){
    data_output <- cont(x=svd.mcmc[,i],y=svd.mcmc[,r+i],pb,lod=i,cex=0.5)
    tmp = data_output
    tmp$byplot = "Genotype"
    tmp$ID = i
    #print(i)
    all_cont[[paste0("G",i)]] = tmp
  }

  all_cont = do.call("rbind",all_cont)

  summary = all_cont %>%
    group_by(byplot,ID) %>%
    summarise(x=mean(x),y=mean(y))

  label_stable = function(data,summary,colx,coly,radius=0){
    data$stable = FALSE
    summary$stable = FALSE
    for (i in 1:nrow(summary)){
      c.byplot = summary$byplot[i]
      c.ID = summary$ID[i]
      tmp = filter(data,byplot==c.byplot,ID==c.ID)
      if (min(tmp[,colx])<=radius & max(tmp[,colx])>=radius & min(tmp[,coly])<=radius & max(tmp[,coly])>=radius){
        data[data$byplot==c.byplot & data$ID==c.ID,"stable"] = TRUE
        summary[i,"stable"] = TRUE
      }
    }
    return(list(data=data,summary=summary))
  }

  # classify stable and non stable lines
  tmp = label_stable(all_cont,summary,2,3)
  data.stable = tmp$data
  summary.stable = tmp$summary
  data.stable$name = data.stable$ID
  summary.stable$name = summary.stable$ID

  # define base plot limits
  x_lim = max(abs(min(data.stable$x)), abs(max(data.stable$x)))
  y_lim = max(abs(min(data.stable$y)), abs(max(data.stable$y)))
  c = 1.01
  x_lim = c*x_lim
  y_lim = c*y_lim

  base_plot = data.stable %>%
    ggplot(aes(x=x,y=y,group=ID,fill=as.factor(ID%%ncolors)))+
    xlim(c(-x_lim, x_lim))+
    facet_grid(~byplot)+
    theme_classic()+
    geom_hline(yintercept = 0,color="gray",linetype="dashed")+
    geom_vline(xintercept = 0,color="gray",linetype="dashed")+
    xlab("")+ylab("")+
    theme(legend.position = "none")+
    theme(panel.background =element_rect(color="black",fill="transparent"),
          axis.line = element_blank(),
          strip.background = element_rect(fill="black"),
          strip.text = element_text(color="white",face="bold",size=10))

  if (plot_stable){
    base_plot = base_plot +
      geom_polygon(data = filter(data.stable, stable),
                   alpha=0.2,color="black")+
      geom_point(data = filter(summary.stable, stable),
                 aes(x=x,y=y,fill=as.factor(ID%%ncolors)),shape=21,color="white",size=2)+
      geom_label_repel(data = filter(summary.stable,stable),
                       aes(x=x,y=y,fill=as.factor(ID%%ncolors),label=name),color="white",
                       max.overlaps = nrow(summary))
  }

  if (plot_unstable){
    base_plot = base_plot +
      geom_polygon(data = filter(data.stable, !stable),
                   alpha=0.2,color="black")+
      geom_point(data = filter(summary.stable, !stable),
                 aes(x=x,y=y,fill=as.factor(ID%%ncolors)),shape=21,color="white",size=2)+
      geom_label_repel(data = filter(summary.stable,!stable),
                       aes(x=x,y=y,fill=as.factor(ID%%ncolors),label=name),color="white",
                       max.overlaps = nrow(summary))
  }
  return(list(contorns = data.stable, biplot_data = summary.stable, biplot = base_plot))
}


data(Maiz)
fm1 <- bayes_ammi( .data = Maiz , .y = y , .gen = entry , .env = site , .rep = rep , .nIter = 200)

output_05 = biplots(fm1, plot_stable = TRUE, plot_unstable = TRUE, pb = 0.05)
output_05$biplot

output_95 = biplots(fm1, plot_stable = TRUE, plot_unstable = TRUE, pb = 0.95)
output_95$biplot



