library(shiny)
library(OrgMassSpecR)
source("ResCount.R")
source("Compute2.R")
source("Compute2.protein.R")
source("form.R")
source("PLOTLY.R")
source("DIGEST.R")
source("Compute.target.R")
shinyServer(function(input, output) {
######  
  seq.pp<-reactive({ 
    validate(
      need(input$seq.pp !="", "")
      )
    gsub(" ","", input$seq.pp) })
  
  seq<-reactive({ 
    validate(
      need(input$seq !="", "")
    )
    gsub(" ","", input$seq) })
######  
  IAA.pp<- reactive({switch(input$Cys.pp,
    "Yes" = TRUE, "No"=FALSE) })
  IAA<- reactive({switch(input$Cys,
     "Yes" = TRUE, "No"=FALSE) })
######
  Cov<-function(Seq){
    sp<-strsplit(Seq, "")[[1]]
    paste(sp, 1:length(sp),sep="")
  }
######
  output$Resi1.pp<-renderUI({
    selectInput("Resi1.pp", label=div("Residue(s)",style="font-family:'futura';font-size:10pt"), choices = Cov(seq.pp()), multiple = TRUE)
  })
  output$PTM1.pp<-renderUI({
    selectInput("PTM1.pp", label=div("PTM #1",style="font-family:'futura';font-size:10pt"), choices =list("none","monoOxi", "diOxi","triOxi"), multiple = FALSE, width = "80%" )
  })
  output$Resi2.pp<-renderUI({
    selectInput("Resi2.pp", label=div("Residue(s)",style="font-family:'futura';font-size:10pt"), choices = Cov(seq.pp()), multiple = TRUE)
  })
  output$PTM2.pp<-renderUI({
    selectInput("PTM2.pp", label=div("PTM #2",style="font-family:'futura';font-size:10pt"), choices =list("none", "Deamidated","Amidated", "Bromo","Phospho","Dimethyl","Trimethyl","Carbonyl"), multiple = FALSE, width="80%" )
  })
  output$Resi3.pp<-renderUI({
    selectInput("Resi3.pp", label=div("Residue(s)",style="font-family:'futura';font-size:10pt"), choices = Cov(seq.pp()), multiple = TRUE)
  })
  output$PTM3.pp<-renderUI({
    selectInput("PTM3.pp", label=div("PTM #3",style="font-family:'futura';font-size:10pt"), choices =list("none", "Deamidated","Amidated", "Bromo","Phospho","Dimethyl","Trimethyl","Carbonyl"), multiple = FALSE, width="80%" )
  })
  
  output$Resi1<-renderUI({
    selectInput("Resi1", label=div("Residue(s)",style="font-family:'futura';font-size:10pt"), choices = Cov(seq()), multiple = TRUE)
  })
  output$PTM1<-renderUI({
    selectInput("PTM1", label=div("PTM #1",style="font-family:'futura';font-size:10pt"), choices =list("none","monoOxi", "diOxi","triOxi"), multiple = FALSE, width = "80%" )
  })
  output$Resi2<-renderUI({
    selectInput("Resi2", label=div("Residue(s)",style="font-family:'futura';font-size:10pt"), choices = Cov(seq()), multiple = TRUE)
  })
  output$PTM2<-renderUI({
    selectInput("PTM2", label=div("PTM #2",style="font-family:'futura';font-size:10pt"), choices =list("none", "Deamidated","Amidated", "Bromo","Phospho","Dimethyl","Trimethyl","Carbonyl"), multiple = FALSE, width="80%" )
  })
  output$Resi3<-renderUI({
    selectInput("Resi3", label=div("Residue(s)",style="font-family:'futura';font-size:10pt"), choices = Cov(seq()), multiple = TRUE)
  })
  output$PTM3<-renderUI({
    selectInput("PTM3", label=div("PTM #3",style="font-family:'futura';font-size:10pt"), choices =list("none", "Deamidated","Amidated", "Bromo","Phospho","Dimethyl","Trimethyl","Carbonyl"), multiple = FALSE, width="80%" )
  })
######
  PF.pp<-reactive({form(input$ptmFormula.pp)})
  mz.pp<-reactive({ Compute2(seq.pp(), IAA.pp(), resi1=length(input$Resi1.pp), ptm1=input$PTM1.pp, resi2=length(input$Resi2.pp), ptm2=input$PTM2.pp, resi3=length(input$Resi3.pp), ptm3=input$PTM3.pp, ptm4=input$chkbx.pp, ptm.form=PF.pp()) })
  ct.pp<-reactive({ length(strsplit(seq.pp(), "")[[1]])})
  
  Fnow.pp<-reactive({paste(names(mz.pp()$ls)[mz.pp()$ls>0], mz.pp()$ls[mz.pp()$ls>0], sep="") })
  Itd.pp<-eventReactive(input$sim.pp,{ IsotopicDistribution(formula=list(C=mz.pp()$ls$C, H=mz.pp()$ls$H+input$ecs.pp, N=mz.pp()$ls$N, O=mz.pp()$ls$O, S=mz.pp()$ls$S, P=mz.pp()$ls$P, Br=mz.pp()$ls$Br), charge=input$ecs.pp)})
  
  
  PF<-reactive({form(input$ptmFormula)})
  mz<-reactive({ Compute2.protein(seq(), IAA(), resi1=length(input$Resi1), ptm1=input$PTM1, resi2=length(input$Resi2), ptm2=input$PTM2, resi3=length(input$Resi3), ptm3=input$PTM3, ptm4=input$chkbx, ptm.form=PF(), input$protein.fromcs, input$protein.tocs) })
  ct<-reactive({ length(strsplit(seq(), "")[[1]])})

  Fnow<-reactive({paste(names(mz()$ls)[mz()$ls>0], mz()$ls[mz()$ls>0], sep="") })
  Itd<-eventReactive(input$sim,{ IsotopicDistribution(formula=list(C=mz()$ls$C, H=mz()$ls$H+input$ecs, N=mz()$ls$N, O=mz()$ls$O, S=mz()$ls$S, P=mz()$ls$P, Br=mz()$ls$Br), charge=input$ecs)})
######
   
## The following three lines serve to reshape the plotting data (simulated MS),
 ## so that the peaks to be shown on the plot will occupy the 6/8 of total m/z range,
 ## the far-left and far-right side will be left blank, for better view.
  
######
  d.before.pp<-reactive({data.frame(mz=min(Itd.pp()$mz)-(max(Itd.pp()$mz)-min(Itd.pp()$mz))*1/8,intensity=0, percent=0) })
  d.after.pp<-reactive({ data.frame(mz=max(Itd.pp()$mz)+(max(Itd.pp()$mz)-min(Itd.pp()$mz))*1/8, intensity=0, percent=0) })
  d.pp<-reactive({rbind(d.before.pp(), Itd.pp(), d.after.pp()) })
  
  
  d.before<-reactive({data.frame(mz=min(Itd()$mz)-(max(Itd()$mz)-min(Itd()$mz))*1/8,intensity=0, percent=0) })
  d.after<-reactive({ data.frame(mz=max(Itd()$mz)+(max(Itd()$mz)-min(Itd()$mz))*1/8, intensity=0, percent=0) })
  d<-reactive({rbind(d.before(), Itd(), d.after()) })
######
  output$text1.pp<-renderText({
    seq.pp()
  })
  output$text2.pp<-renderText({
    mz.pp()$'Original.F'
  })
  output$text3.pp<-renderText({
    ct.pp()})
  output$tb1.pp<-renderTable({
    mz.pp()$TB
  }, digits=4, align=rep("c",7))
  
  observeEvent(input$act.pp, {
    output$tb2.pp<-DT::renderDataTable({
      ResCount(seq.pp())}, options=list(
        lengthMenu=list(c(5,10,-1),c("5","15","All")),
        pageLength=5, searching=FALSE, paging=FALSE
      )) })
  
  output$text4.pp<-renderText({
    Fnow.pp()})
  output$Idplot.pp<-renderPlotly({
    PLOTLY(d.pp())
  })  
  
 NCOL<-observe({2+as.numeric(input$protein.tocs)-as.numeric(input$protein.fromcs)})
  output$text1<-renderText({
    seq()
    })
  output$text2<-renderText({
    mz()$'Original.F'
    })
  output$text3<-renderText({
    ct()})
 observeEvent(as.numeric(input$protein.tocs)-as.numeric(input$protein.fromcs),{
   output$tb1<-renderTable({
   mz()$TB
    }, digits=4, align=rep("c",2+as.numeric(input$protein.tocs)-as.numeric(input$protein.fromcs))) })
  
  observeEvent(input$act, {
    output$tb2<-DT::renderDataTable({
    ResCount(seq())}, options=list(
      lengthMenu=list(c(5,10,-1),c("5","15","All")),
      pageLength=5, searching=FALSE, paging=FALSE
    )) })
  
  output$text4<-renderText({
    Fnow()})
  output$Idplot<-renderPlotly({
    PLOTLY(d())
    })
######
  
Insilico<-reactive({ DIGEST(seq(), sites=input$cleavSite, term=input$side ,missed=input$mc, thresh=input$thresh) })
  
 observeEvent(input$insilico, {
    output$tb3<-DT::renderDataTable({
      isolate(Insilico())}, options=list(
        lengthMenu=list(c(10,25,-1),c("10","25","All")),
        pageLength=25)) })
  
  observeEvent(input$targeted, {
   output$tb4<-DT::renderDataTable({
    isolate(Compute.target(Insilico(),input$csfrom, input$csto,input$target.ptm, input$lowmz, input$highmz, input$NM, input$Fml, input$allRes)[,c(-7,-8, -9)])}, options=list(
      lengthMenu=list(c(50,100,-1),c("50","100","All")),
      pageLength=100)) })
 output$dld.target<-downloadHandler(filename=function(){paste("Targeted-MS-",Sys.Date(),".csv",sep="")}, content=function(file){write.csv(Compute.target(Insilico(),input$csfrom, input$csto,input$target.ptm, input$lowmz, input$highmz,input$NM, input$Fml, input$allRes)[,c(-7,-8,-9)], file)})   
    
  
})
