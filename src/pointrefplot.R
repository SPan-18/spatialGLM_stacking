# rm(list = ls())
# library(ggplot2)
# library(latex2exp)
# library(MBA)
# library(fields)


col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow", "red"))
col.RdBu <- colorRampPalette(c("blue", "white", "red"))
# col.br <- colorRampPalette(c('#ca0020','#f4a582','#ffffff','#bababa','#404040'))

pointref_plot=function(tab, colname, 
                x_lim = c(0,1), y_lim = c(0,1), 
                title = NULL,
                legend_title = "z"){
  
  surf <- mba.surf(tab[,c("s1","s2",colname)], 
                   no.X = 200, no.Y = 200, h = 8, m = 1, n = 1, 
                   extend=FALSE)$xyz.est
  
  surf_df <- data.frame(expand.grid(surf$x, surf$y), z = as.vector(surf$z))
  surf_df <- na.omit(surf_df)
  names(surf_df) <- c("x", "y", "z")
  
  plot <- ggplot(surf_df, aes(x = x, y = y)) +
    geom_raster(aes(fill = z)) +
    # scale_fill_gradientn(colours = col.br(100)) +
    # scale_fill_gradientn(colours = col.RdBu(100)) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1) +
    xlab("Easting") +
    ylab("Northing") +
    labs(fill = legend_title) +
    theme_bw() +
    theme(legend.title = element_text(size = 10, hjust = 0.25),
          # axis.title.x = element_text(family = "mono"),
          # axis.title.y = element_text(family = "mono"),
      # axis.line = element_line(color='black'),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # panel.border = element_blank(),
      aspect.ratio = 1)
  
  if(!is.null(title)){
    plot <- plot + labs(title = title)  +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  plot
}

# leg_title <- TeX('$z(s)$')
# pointref_plot(simdat,"z", legend_title = leg_title)
