
DIGEST<-function(sequence, sites="RK", term="C",missed=0, thresh=4){
  split_sq<-strsplit(sequence,"")[[1]]
  split_st<-strsplit(sites,"")[[1]]
  id0<-as.numeric()
  if(term=="C"){
  for(i in 1:length(split_sq)){
    if(split_sq[i] %in% split_st & split_sq[i+1] !="P" & i<length(split_sq)) {
      id<-i
      id0<-c(id0,id)}
  } ## After the for loop, this id0 contains the site numbers for all the specified sites (RK with the next residue not being P, for example).
  
  ## The following part is to calculate the start and stop sites. 
  start<-1
  for(t in id0){
    start0<-t+1
    start<-c(start,start0)}
  stop<-sort(c(length(split_sq),id0)) ## The start and stop here are for the 0 miscleavage.
  }
  
if(term=="N"){
  for(i in 1:length(split_sq)){
    if(split_sq[i] %in% split_st & i<length(split_sq)) {
      id<-i
      id0<-c(id0,id)}
  }
  start<-sort(c(1,id0))
  stop<-sort(c(id0-1, length(split_sq)))
}

  ## Below consider the missed cleavage, calculate the new start and stop number.
  mis<-0:missed
  pp_list0<-data.frame("peptide"=as.character(), "start"=as.numeric(),"stop"=as.numeric(),"miss cleav."=as.numeric(), check.names = FALSE)
  
  for(x in mis){
    a<-length(start)-x
    b<-1+x
    start_final<-start[1:a]
    stop_final<-stop[b:length(stop)]
    ## Generate the dataframe that contains peptides starting and ending number, for each miscleavage count. And row bind them all together.      
    pp_list<-data.frame("peptide"=rep("X",length(start_final)),"start"=start_final, "stop"=stop_final, "miss cleav."=x, check.names = FALSE)
    pp_list0<-rbind(pp_list0,pp_list)}
  
  pp<-as.character()
  for(n in 1:nrow(pp_list0)){
    pp0<-paste(split_sq[pp_list0$start[n]:pp_list0$stop[n]], collapse="")
    pp<-c(pp,pp0)}
  pp_list0$peptide<-pp  ## This pp_list0 corresponds to the digestion with the correct miscleavage number.
  pp_list0$Length<-pp_list0$stop-pp_list0$start+1
  RS<-subset(pp_list0, Length>thresh)
  RS.name<-RS$peptide
  RS<-RS[,-1]
  row.names(RS)<-RS.name
  return(RS)
}
  
  