Fixed.Mod<-function(peptides, ptm1, res1, ptm2, res2, ptm3, res3, ptmfml, resfml){ ## the ptm1, res1, ..., ptm3, res3, correspond to those on the sidebar. And the ptmfml and resfml correspond to the user-input PTM in the sidebar.
  source("Res.to.Num.R")
  source("form.R")
  PTM<-data.frame("C"=c(0,0,0,0,0,2,0,0,0,0,0,2,3,0,0),"H"=c(0,0,0,0,1,2,1,-1,3,-3,-2,4,6,0,-2),"N"=c(0,0,0,0,0,0,1,-1,1,-1,0,0,0,0,0),"O"=c(0,1,2,3,3,1,-1,1,0,0,-1,0,0,0,1),"S"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), "P"=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),"Br"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),"Cl"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),"Si"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),"F"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  if(substr(resfml,1,1)!="C"){
    PTM.add<-form(ptmfml)
    } else {
      PTM.add<-list("C"=form(ptmfml)$C-2, "H"=form(ptmfml)$H-3, "N"=form(ptmfml)$N-1, "O"=form(ptmfml)$O-1, "S"=form(ptmfml)$S, "P"=form(ptmfml)$P, "Br"=form(ptmfml)$Br, "Cl"=form(ptmfml)$Cl, "Si"=form(ptmfml)$Si, "F"=form(ptmfml)$F)}
  
  PTM<-rbind(PTM, PTM.add)
  row.names(PTM)<-c("none","monoOxi","diOxi","triOxi","Phospho","Acetyl","Amidated","Deamidated","Ammonium-add","Ammonium-loss","Dehydrated","Dimethyl","Trimethyl","Bromo","Carbonyl",ptmfml)
  
  al<-data.frame("RES"=c(Res.to.Num(res1), Res.to.Num(res2), Res.to.Num(res3),Res.to.Num(resfml)), "PTM"=c(rep(ptm1, length(res1)), rep(ptm2, length(res2)), rep(ptm3, length(res3)),rep(ptmfml, length(resfml))))
   
  peptides$Mods<-rep("",nrow(peptides))
  for(i in 1:nrow(peptides)){
    id<-which(al$RES>=peptides$start[i] & al$RES<=peptides$stop[i])
    if(length(id)==0){
      next} else {
        mods<-as.character(al$PTM[id])
        
        combine.mods<-PTM[row.names(PTM) %in% mods,]
        cm<-0
        for(s in 1:nrow(combine.mods)){
          cm0<-combine.mods[s,]*sum(mods %in% row.names(combine.mods)[s])
          cm<-cm+cm0}
        
        combine.mods<-paste(paste(names(cm), cm, sep=""), collapse = " ")
        peptides$Mods[i]<-combine.mods}}
  ## after this for loop, we add new column named "Mods", which essentially is the add-up of all the elemental compositions of PTMs for each peptide.
  ConvertPeptide2<-function(F){ ## F is any peptide sequence
    x<-ConvertPeptide(F, IAA=TRUE) ## x is the result from ConvertPeptide() function
    x<-list("C"=x$C, "H"=x$H, "N"=x$N,"O"=x$O,"S"=x$S,"P"=0,"Br"=0,"Cl"=0,"Si"=0,"F"=0)
    return(x)}
  ## The ConvertPeptide2 function simply makes the original ConvertPeptide function display more details, ouput the enlongated elemental list.
  peptides$Overall<-rep("", nrow(peptides))
  for(h in 1:nrow(peptides)){
    x<-unlist(ConvertPeptide2(row.names(peptides)[h])) + unlist(form(peptides$Mods[h]))
    x<-paste(paste(names(x), x, sep=""), collapse = " ")
    peptides$Overall[h]<-x}
  return(lapply(peptides$Overall, form))
}