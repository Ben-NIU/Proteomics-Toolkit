
## This is the function where we complete a ggplot first,then convert it to plotly,
## and modify the hoverinfo format.
PLOTLY<-function(data){
  library(plotly)
  pdf(NULL)
A<-ggplot(data, aes(x=mz, xend=mz,y=0, yend=percent)) + geom_segment(size=1,colour="brown3") + labs(y="Relative Intensity (%)", x="m/z")+ scale_x_continuous(breaks=pretty(data$mz,n=6))  +  theme(panel.border=element_rect(fill=NA, size=1, linetype="solid", color="black")) + coord_cartesian(ylim=c(0,105)) + scale_y_continuous(expand=c(0,0)) + theme(panel.background=element_blank(), axis.title.y=element_text(size=15, color="#006666"), axis.title.x=element_text(size=14, color="#006666")) + xlim(min(data$mz)-(max(data$mz)-min(data$mz))*1/8, max(data$mz)+(max(data$mz)-min(data$mz))*1/8)
P<-plotly_build(A)
P$data[[1]]$text<-paste("m/z: ", P$data[[1]]$x,"<br>int.: ", P$data[[1]]$y, sep="") 
P
}