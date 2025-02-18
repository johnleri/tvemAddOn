organize_plot_data = function(model.fit.list, time.variable = "time_grid", intercept = F){

  observed.data = model.fit.list[["observed_model"]]
  plot.data = do.call(rbind,(pool.variate.estimates(observed.data, return.term = "beta")))

  plot.data$variable = sapply(strsplit(rownames(plot.data), "\\."), '[[', 1)
  plot.data$time = rep(observed.data[[1]][[time.variable]], times = nrow(plot.data)/100)

  plot.data$OR = exp(plot.data$qbar)
  plot.data$OR_CI_ll = exp(plot.data$qbar - plot.data$t * 1.96)
  plot.data$OR_CI_ul = exp(plot.data$qbar + plot.data$t * 1.96)

  if(intercept == F){
    plot.data = plot.data[which(plot.data$variable != "(Intercept)"),]
  }

  return(plot.data)
}
