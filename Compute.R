library(OrgMassSpecR)
Compute<-function(Seq, Carbam, resi1, ptm1, resi2, ptm2, resi3, ptm3, ptm4){
  Sq<-ConvertPeptide(Seq, IAA=Carbam)
  MM<-(MonoisotopicMass(formula=list(C=Sq$C, H=Sq$H, N=Sq$N, O=Sq$O, S=Sq$S)) + resi1*ptm1 + resi2*ptm2 + resi3*ptm3 + sum(as.numeric(ptm4)) + (1:6)*1.007276466)/(1:6)
  tb<-data.frame("+1 charge"=MM[1],"+2 charge"=MM[2],"+3 charge"=MM[3],"+4 charge"=MM[4],"+5 charge"=MM[5],"+6 charge"=MM[6], check.names = FALSE)
  row.names(tb)<-"m/z"
  return(tb)
}