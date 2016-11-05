Compute2.protein<-function(seq,carbam, resi1, ptm1, resi2, ptm2, resi3, ptm3, ptm4, ptm.form, x,y ){
  sq0<-ConvertPeptide(seq, IAA=FALSE)
  sq<-ConvertPeptide(seq, IAA=carbam)
  fm<-paste(paste(c("C","H","N","O","S"),c(sq0$C,sq0$H,sq0$N,sq0$O,sq0$S), sep=""), collapse = " ")
  PTM<-data.frame("C"=c(0,0,0,0,0,2,0,0,0,0,0,2,3,0,0),"H"=c(0,0,0,0,1,2,1,-1,3,-3,-2,4,6,0,-2),"N"=c(0,0,0,0,0,0,1,-1,1,-1,0,0,0,0,0),"O"=c(0,1,2,3,3,1,-1,1,0,0,-1,0,0,0,1),"S"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), "P"=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),"Br"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0))
  row.names(PTM)<-c("none","monoOxi","diOxi","triOxi","Phospho","Acetyl","Amidated","Deamidated","Ammonium-add","Ammonium-loss","Dehydrated","Dimethyl","Trimethyl","Bromo","Carbonyl")
  P1<-PTM[which(row.names(PTM)==ptm1),]
  P2<-PTM[which(row.names(PTM)==ptm2),]
  P3<-PTM[which(row.names(PTM)==ptm3),]
  P4<-PTM[row.names(PTM) %in% ptm4,] ## PTM4 corresponds to the extra PTM, and ptm.form corresponds to the input formula of PTM.
  
  LS<-data.frame(C=sq$C+resi1*P1$C +resi2*P2$C + resi3*P3$C + sum(P4$C) + ptm.form$C, H=sq$H+resi1*P1$H +resi2*P2$H + resi3*P3$H + sum(P4$H) + ptm.form$H,N=sq$N+resi1*P1$N +resi2*P2$N + resi3*P3$N + sum(P4$N)+ ptm.form$N,O=sq$O+resi1*P1$O +resi2*P2$O +resi3*P3$O + sum(P4$O)+ ptm.form$O,S=sq$S+resi1*P1$S +resi2*P2$S + resi3*P3$S + sum(P4$S)+ ptm.form$S,P=resi1*P1$P +resi2*P2$P + resi3*P3$P + sum(P4$P)+ ptm.form$P,Br=resi2*P2$Br+resi3*P3$Br + ptm.form$Br, Cl=ptm.form$Cl, F=ptm.form$F, Si=ptm.form$Si)
  
  MM<-(MonoisotopicMass(formula=list(C=sq$C+resi1*P1$C +resi2*P2$C + resi3*P3$C + sum(P4$C), H=sq$H+resi1*P1$H +resi2*P2$H + resi3*P3$H + sum(P4$H),N=sq$N+resi1*P1$N +resi2*P2$N + resi3*P3$N + sum(P4$N),O=sq$O+resi1*P1$O +resi2*P2$O +resi3*P3$O + sum(P4$O),S=sq$S+resi1*P1$S +resi2*P2$S + resi3*P3$S + sum(P4$S),P=resi1*P1$P +resi2*P2$P + resi3*P3$P + sum(P4$P),Br=resi2*P2$Br+resi3*P3$Br)) + MonoisotopicMass(formula=ptm.form) + (x:y)*1.007276466)/(x:y)
  
  tb<-data.frame("m/z"=MM, check.names=FALSE)
  row.names(tb)<-paste("+", x:y, sep="")
  list("TB"=t(tb),"Original.F"=fm,"ls"=LS)
}