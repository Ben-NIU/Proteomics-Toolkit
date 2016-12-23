library(shiny)
library(OrgMassSpecR)

source("ResCount.R")
source("Compute2.R")
source("Compute2.protein.R")
source("form.R")
source("PLOTLY.R")
source("DIGEST.R")
source("Fixed.Mod.R")
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
  output$ud<-renderUI({
    selectInput("ud", label=div(p(em(span("Residue(s)"))),style="font-family:'marker felt';color:purple; font-size:12pt"), choices = Cov(seq()), multiple = TRUE)
  })
######
  PF.pp<-reactive({form(input$ptmFormula.pp)})
  mz.pp<-reactive({ Compute2(seq.pp(), IAA.pp(), resi1=length(input$Resi1.pp), ptm1=input$PTM1.pp, resi2=length(input$Resi2.pp), ptm2=input$PTM2.pp, resi3=length(input$Resi3.pp), ptm3=input$PTM3.pp, ptm4=input$chkbx.pp, ptm.form=PF.pp()) })
  ct.pp<-reactive({ length(strsplit(seq.pp(), "")[[1]])})
  ## this fmz() function takes the ls from mz, and convert the ls to list with the correct number of H.
  fmz<-function(x, charge){
    y<-x$ls
    y$H<-y$H+charge
    as.list(y)}
  
  Fnow.pp<-reactive({paste(names(mz.pp()$ls)[mz.pp()$ls>0], mz.pp()$ls[mz.pp()$ls>0], sep="") })
  fj.pp<-reactive({fmz(mz.pp(), input$ecs.pp)})
  Itd.pp<-eventReactive(input$sim.pp,{ IsotopicDistribution(as.list(fj.pp()), charge=input$ecs.pp)})
  
  
  PF<-reactive({form(input$ptmFormula)})
  mz<-reactive({ Compute2.protein(seq(), IAA(), resi1=length(input$Resi1), ptm1=input$PTM1, resi2=length(input$Resi2), ptm2=input$PTM2, resi3=length(input$Resi3), ptm3=input$PTM3, ptm4=input$chkbx, ptm.form=PF(), input$protein.fromcs, input$protein.tocs) })
  ct<-reactive({ length(strsplit(seq(), "")[[1]])})

  Fnow<-reactive({paste(names(mz()$ls)[mz()$ls>0], mz()$ls[mz()$ls>0], sep="") })
 
  fj<-reactive({fmz(mz(), input$ecs)})
  Itd<-eventReactive(input$sim,{ IsotopicDistribution(as.list(fj()), charge=input$ecs)})

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
  }, digits=4, align="c")
  
  observeEvent(input$act.pp, {
    output$tb2.pp<-DT::renderDataTable({
      ResCount(seq.pp())}, options=list(
        lengthMenu=list(c(5,10,-1),c("5","15","All")),
        pageLength=5, searching=FALSE, paging=FALSE
      )) })
  
  output$text4.pp<-renderText({
    Fnow.pp()})
  output$Idplot.pp<-renderPlotly({
    PLOTLY(Itd.pp())
  })  
  
  output$text1<-renderText({
    seq()
    })
  output$text2<-renderText({
    mz()$'Original.F'
    })
  output$text3<-renderText({
    ct()})

   output$tb1<-renderTable({
   mz()$TB
    }, digits=4, align="c")
  
  observeEvent(input$act, {
    output$tb2<-DT::renderDataTable({
    ResCount(seq())}, options=list(
      lengthMenu=list(c(5,10,-1),c("5","15","All")),
      pageLength=5, searching=FALSE, paging=FALSE
    )) })
  
  output$text4<-renderText({
    Fnow()})
  output$Idplot<-renderPlotly({
    PLOTLY(Itd())
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
    isolate(Compute.target(Insilico(),input$csfrom, input$csto,input$target.ptm, input$lowmz, input$highmz, input$NM, input$Fml, input$allRes, input$PTM1, input$Resi1, input$PTM2, input$Resi2, input$PTM3, input$Resi3, input$ptmFormula, input$ud)[,c(-7,-8, -9)])}, options=list(
      lengthMenu=list(c(50,100,-1),c("50","100","All")),
      pageLength=100)) })
 output$dld.target<-downloadHandler(filename=function(){paste("Targeted-MS-",Sys.Date(),".csv",sep="")}, content=function(file){write.csv(Compute.target(Insilico(),input$csfrom, input$csto,input$target.ptm, input$lowmz, input$highmz,input$NM, input$Fml, input$allRes,input$PTM1, input$Resi1, input$PTM2, input$Resi2, input$PTM3, input$Resi3, input$ptmFormula, input$ud)[,c(-7,-8,-9)], file)})   
    
  
})
