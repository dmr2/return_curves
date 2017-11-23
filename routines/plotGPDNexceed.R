insert_minor <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )}


plotGPDNExceed <- function( df, historical){
  
  library(ggplot2)
  library(scales)
  
  print(historical$name)
  
  xmax <- ceiling(tail(df$height,n=1000)[which.min(abs(tail(df$freq,n=1000)-1e-4))])
  #xmax <- ceiling(df$height[which.min(abs(df$freq-1e-4))])
  if (xmax <= 3){
    xmax <- 4
  }
  
  lty=c("solid","dotted","dotted","dotted","solid","solid","solid","dashed","solid","dashed")
  colors = c("grey","#1f78b4","#33a02c","#ff7f00","#1f78b4","#33a02c","#ff7f00","grey","grey","grey")
  lwd = c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,.5,.5,.5)
  breaks=c("Historic Flood Return Curve","N+SL50 1.5°C","N+SL50 2.0°C","N+SL50 2.5°C",
           "Ne 1.5°C","Ne 2.0°C","Ne 2.5°C")
           
  
  p <- ggplot(df)  + 
  
    geom_line(data=df, aes(x=height,y=freq,colour=name,lty=name,size=name))+
    geom_point(data=historical, aes(x=z,y=freq,shape=group),size=5,color="black",
               stroke=1.5,alpha=0.5) +
    scale_size_manual("",values=lwd,guide=FALSE) +
    scale_linetype_manual("",breaks=breaks,values=lty) +
    scale_shape_manual("",values=1) +
    scale_color_manual("",breaks=breaks,values=colors) +
    guides(colour = guide_legend(override.aes = list(lwd=1.5))) +
    
    theme_bw(base_size = 27) + annotation_logticks(sides = "lr", size=1,
    short = unit(0.3, "cm"), mid = unit(0.3, "cm"), long = unit(0.5, "cm")) + 
    scale_x_continuous(breaks=seq(0,xmax,by=.5),limits=c(0,xmax),
                       labels=head(insert_minor(seq(0,xmax,by=1),1),-1),expand = c(0,0)) + 
    scale_y_log10(breaks=c(10,1,0.1,0.01,0.001,0.0001),
                  labels=trans_format("log10", scales::math_format(10^.x)),expand = c(0,0)) +
    coord_cartesian(ylim=c(1e-4, 10)) +
    
    theme(legend.justification = c(1, 1), legend.position = c(1, -.3), 
          legend.title = element_blank(),legend.key = element_blank(),
          legend.direction="horizontal",legend.text=element_text(size=18),
          panel.grid.minor = element_blank(), aspect.ratio=.5,
          panel.border = element_rect(linetype = "solid", colour = "black", size=2),
          axis.line = element_line(colour = 'black', size = 5),
          axis.ticks = element_line(colour = "black", size = 1),
          plot.margin=unit(c(1,1,5,0.5),"cm"),
          legend.key.size = unit(2, 'lines')) +
    labs(title=paste(df$loc," (",df$year[1],")",sep=""), 
         x="Flood Height (m)",y="Expected Events per Year") 
  
   fil = paste("flood_return_curves_",df$basin[1],"_",df$id[1],"_",df$year[1],".pdf",sep="")
   ggsave(fil,p,width=11, height=8.5)
  # Plot in Feet
  
  
}