j<-function(x, RES){
  SS<-sum(strsplit(x, split="")[[1]] %in% RES)
  if(SS>0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
  