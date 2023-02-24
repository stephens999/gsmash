# input: fit, topics, grouping

# poisson2multinom
#
library(fastTopics)
library(ggplot2)
structure_plot_general = function(Lhat,Fhat,grouping,title=NULL,
                                  print_plot=TRUE,
                                  seed=12345,
                                  n_samples = 2000,gap=40,LD = TRUE,
                                  remove_l0f0 = TRUE,
                                  topic_model=FALSE,
                                  colors = c('#a6cee3',
                                    '#1f78b4',
                                    '#b2df8a',
                                    '#33a02c',
                                    '#fb9a99',
                                    '#e31a1c',
                                    '#fdbf6f',
                                    '#ff7f00',
                                    '#cab2d6',
                                    '#6a3d9a',
                                    '#ffff99',
                                    '#b15928')){
  set.seed(seed)
  #s       <- apply(Lhat,2,max)
  #Lhat    <-	t(t(Lhat) / s)
  if(!topic_model){
    ldf = my_ldf(Lhat,Fhat)
  }else{
    ldf = list(l=Lhat,d=1)
  }

  if(LD){
    Lhat = ldf$l%*%diag(ldf$d)
  }else{
    Lhat = ldf$l
  }
  if(remove_l0f0){
    Lhat = Lhat[,-c(1,2),drop=F]
  }
  Fhat = matrix(1,nrow=3,ncol=ncol(Lhat))
  if(is.null(colnames(Lhat))){
    colnames(Lhat) <- paste0("k",1:ncol(Lhat))
  }
  fit_list     <- list(L = Lhat,F = Fhat)
  class(fit_list) <- c("multinom_topic_model_fit", "list")
  p <- structure_plot(fit_list,grouping = grouping,
                      n = n_samples,gap = gap,colors=colors,verbose=F) +
    labs(y = "loading",color = "dim",fill = "dim") + ggtitle(title)
  if(print_plot){
    print(p)
  }
  return(p)
}

my_ldf = function(Lhat,Fhat){
  dl = apply(Lhat,2,norm,type='2')
  df = apply(Fhat,2,norm,type='2')
  return(list(l = apply(Lhat,2,function(z){z/norm(z,'2')}),
              f = apply(Fhat,2,function(z){z/norm(z,'2')}),
              d = dl*df))
}
