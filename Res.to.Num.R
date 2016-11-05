Res.to.Num<-function(res){
	if(length(res)==0){
		return(NULL)
		} else {
  as.numeric(substr(res, 2, length(strsplit(res, "")[[1]])))}
  
}