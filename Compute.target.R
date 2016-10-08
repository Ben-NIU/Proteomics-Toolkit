Compute.target<-function(DIGest, CS.from, CS.to, MOD, min,max){
  sta<-DIGest$start
  stp<-DIGest$stop
  MD.list<-data.frame("md"=c("none","Oxidation","diOxidation","triOxidation","Carbonyl","Deamidation(N)","Phosphorylation(S,T,Y)"),"add.mass"=c(0,15.99491,31.98983,47.98474, 13.97926, 0.9840156,79.96633))
  pp<-row.names(DIGest)
  cp<-lapply(pp, ConvertPeptide, IAA=TRUE)
  MD<-MD.list[MD.list$md %in% MOD,]
  line<-NULL
  for(m in MD$add.mass){
  	for(i in CS.from:CS.to){
  		mz<-round((unlist(lapply(cp, MonoisotopicMass, charge=0)) + m + i*1.007276466)/i,4)
  		line0<-data.frame("m/z"=mz, "charge"=paste("+ ", i, sep=""), "Mod"=as.character(MD[which(MD$add.mass==m),]$md), check.names=FALSE)
  		line<-rbind(line,line0)}
  }
  line$Start<-sta
  line$Stop<-stp
  line$Sequence<-pp
  judge.deami<-unlist(lapply(line$Sequence, function(x){"N" %in% strsplit(x, "")[[1]]}))
  judge.phosph<-unlist(lapply(line$Sequence, function(x){"Y" %in% strsplit(x, "")[[1]]})) | unlist(lapply(line$Sequence, function(x){"S" %in% strsplit(x, "")[[1]]})) | unlist(lapply(line$Sequence, function(x){"T" %in% strsplit(x, "")[[1]]}))
  line$Judge.deami<-judge.deami
  line$Judge.phosph<-judge.phosph
  line<-subset(line, !(Mod=="Deamidation(N)" & Judge.deami==FALSE))
  line<-subset(line, !(Mod=="Phosphorylation(S,T,Y)" & Judge.phosph==FALSE))
  line<-subset(line, line[,1]>=min & line[,1]<=max)
  return(line)
}
  