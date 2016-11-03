library(shiny)
library(DT)
library(plotly)
shinyUI(fluidPage(

  navbarPage(title=span(strong("Proteomics Toolkit"), style="font-family:'rockwell extra bold';color:slateblue; font-size:20pt"),
    ## The 1st tabPanel is for peptide analysis, see below.
      tabPanel(title=span(strong("Peptide Analysis"),style="font-family:'luminari';color:darkgreen;font-size:13pt"),icon=icon("leaf",lib="font-awesome"),
        sidebarLayout(
          sidebarPanel(
            textInput("seq.pp", label=div(em("Input peptide sequence"),style="font-family:'marker felt';color:darkgreen; font-size:12pt"), placeholder = "for example: YGGFL"),
            radioButtons("Cys.pp", label=div(em("Carbamiodomethyl on Cys ?"),style="font-family:'marker felt';color:darkgreen; font-size:12pt"), choices = list("Yes","No"), selected = "Yes"),
            textInput("ptmFormula.pp", label = div(p(em(span("PTM in formula",style="font-family:'marker felt';color:darkgreen; font-size:12pt"))),  p(em(span('Note:',style="font-family:'marker felt';color:darkred; font-size:10pt"), span('A "SPACE" is required between every element! The allowed elements are',style="font-family:'marker felt';color:darkgreen; font-size:10pt"), span('C, H, N, O, S, P, Br, Cl, F, Si.', style="font-family:'marker felt';color:darkred; font-size:10pt")))), placeholder = "for example: C5 H8 O2 N Br2"),
            hr(),
            div(p(strong(em("Residue-specific PTMs"))), style="font-family:'marker felt';color:darkgreen; font-size:12pt"),
            fluidRow(column(7,uiOutput("PTM1.pp")),
                     column(4, uiOutput("Resi1.pp"))),
            fluidRow(column(7, uiOutput("PTM2.pp")),
                     column(4, uiOutput("Resi2.pp"))),
            fluidRow(column(7, uiOutput("PTM3.pp")),
                     column(4, uiOutput("Resi3.pp"))),
            checkboxGroupInput("chkbx.pp",label=div(strong(em("extra PTMs")),style="font-family:'marker felt';color:darkgreen; font-size:12pt"),choices=list("none","Acetyl", "Ammonium-add","Ammonium-loss","Dehydrated"), selected = "none")
        
               
          ),
        
        mainPanel(
          tabsetPanel(type="tabs", 
                tabPanel(title=span(strong("Peptide Computations"), style="font-family:'constantia';color:darkgreen;font-size:10pt"),icon = icon("calculator",lib = "font-awesome"),
                        br(),
                fluidRow(
                           column(5, div(h4(em(p("Peptide Sequence"))), style="font-family:'hannotate tc';color:darkgreen")),
                           column(5, div(h4(em(p("Elemental Composition"))), style="font-family:'hannotate tc';color:darkgreen")),
                           column(2, div(h4(em(p("Res.#"))), style="font-family:'hannotate tc';color:darkgreen"))
                          ),
                fluidRow(
                          column(5, verbatimTextOutput("text1.pp")),
                          column(5, verbatimTextOutput("text2.pp")),
                          column(2, verbatimTextOutput("text3.pp"))
                     ),
                  br(),
                 tableOutput("tb1.pp"),
                 hr(),
             fluidRow(
                column(3, actionButton("act.pp", label=span(strong("AA Composition"), style="font-family:'wawati tc';color:darkgreen; font-size:12pt"), icon = icon("cubes",lib="font-awesome"))),  
                column(8,DT::dataTableOutput("tb2.pp"), offset = 3)
             )
          ),       
          tabPanel(span(strong("MS Spectral Simulation"),style="font-family:'constantia';color:darkgreen;font-size:10pt"), icon=icon("laptop",lib = "font-awesome"),
                   br(),
                   div(h4(em(p("Formula (with PTMs)"))), style="font-family:'hannotate tc';color:darkgreen"),
                   fluidRow(
                     column(6, verbatimTextOutput("text4.pp"))
                   ),
                   hr(),
                   fluidRow(
                     column(5, div(h4(p(em("Charge State to Simulate: "))), style="font-family:'hannotate tc';color:darkgreen")),
                     column(2, numericInput("ecs.pp", label=NULL, value=1,width='80%', min=1, step=1)),
                     column(2, actionButton("sim.pp", label = span(strong("Simulate"), style="font-family:'wawati tc';color:darkgreen; font-size:12pt"), icon=icon("bar-chart",lib="font-awesome")))
                   ),
                   hr(),
                   plotlyOutput("Idplot.pp")
          )                                   
                                            
                                                                
  ))
               )
        ),
  ## The 2nd tabPanel is for the protein analysis, see below.
      tabPanel(title=span(strong("Protein Analysis"), style="font-family:'luminari';color:purple; font-size:13pt"), icon = icon("heartbeat",lib="font-awesome"),
  sidebarLayout(
    sidebarPanel(
      textInput("seq", label=div(em("Input protein sequence"),style="font-family:'marker felt';color:purple; font-size:12pt")),
      radioButtons("Cys", label=div(em("Carbamiodomethyl on Cys ?"),style="font-family:'marker felt';color:purple; font-size:12pt"), choices = list("Yes","No"), selected = "Yes"),
      textInput("ptmFormula", label = div(p(em(span("PTM in formula",style="font-family:'marker felt';color:purple; font-size:12pt")))), placeholder = "for example: C5 H10 S F2"),
      hr(),
      div(p(strong(em("Residue-specific PTMs"))), style="font-family:'marker felt';color:purple; font-size:12pt"),
      fluidRow(column(7,uiOutput("PTM1")),
               column(4, uiOutput("Resi1"))),
      fluidRow(column(7, uiOutput("PTM2")),
               column(4, uiOutput("Resi2"))),
      fluidRow(column(7, uiOutput("PTM3")),
               column(4, uiOutput("Resi3"))),
      checkboxGroupInput("chkbx",label=div(strong(em("extra PTMs")),style="font-family:'marker felt';color:purple; font-size:12pt"),choices=list("none","Acetyl", "Ammonium-add","Ammonium-loss","Dehydrated"), selected = "none")
    ),
    
    mainPanel(
      tabsetPanel(type="tabs",
          tabPanel(span(strong("Computations"), style="font-family:'constantia';color:purple;font-size:10pt"),icon = icon("calculator",lib = "font-awesome"),
         br(),
      fluidRow(
        column(5, div(h4(em(p("Protein Sequence"))), style="font-family:'hannotate tc';color:purple")),
        column(5, div(h4(em(p("Elemental Composition"))), style="font-family:'hannotate tc';color:purple")),
        column(2, div(h4(em(p("Res.#"))), style="font-family:'hannotate tc';color:purple"))
        ),
      fluidRow(
        column(5, verbatimTextOutput("text1")),
        column(5, verbatimTextOutput("text2")),
        column(2, verbatimTextOutput("text3"))
        ),
      br(),
      fluidRow(
        column(3, textInput("protein.fromcs", label=span(h4(em(p("from charge"))), style="font-family:'hannotate tc';color:purple;font-size:10pt"), value=5)),
        column(3, textInput("protein.tocs", label=span(h4(em(p("to charge"))), style="font-family:'hannotate tc';color:purple;font-size:10pt"), value=10))
      ),
      tableOutput("tb1"),
      hr(),
      fluidRow(
      column(3, actionButton("act", label=span(strong("AA Composition"), style="font-family:'wawati tc';color:purple; font-size:12pt"), icon = icon("cubes",lib="font-awesome"))),  
      column(8,DT::dataTableOutput("tb2"), offset = 3)
      )
          ),
         tabPanel(span(strong("MS Spectral Simulation"),style="font-family:'constantia';color:purple;font-size:10pt"),icon = icon("laptop",lib = "font-awesome"),
                  br(),
               div(h4(em(p("Formula (with PTMs)"))), style="font-family:'hannotate tc';color:purple"),
               fluidRow(
                 column(6, verbatimTextOutput("text4"))
                 ),
               hr(),
              fluidRow(
              column(6, div(h4(p(em("Charge State to Simulate: "))), style="font-family:'hannotate tc';color:purple")),
              column(2, numericInput("ecs", label=NULL, value=1,width='80%', min=1, step=1)),
              column(2, actionButton("sim", label = span(strong("Simulate"), style="font-family:'wawati tc';color:purple; font-size:12pt"), icon=icon("bar-chart",lib="font-awesome")))
              ),
              hr(),
              plotlyOutput("Idplot")
                  ),
         tabPanel(span(strong(em("in silico")), strong("Digestion"),style="font-family:'constantia';color:purple;font-size:10pt"),icon = icon("align-right",lib = "font-awesome"),
                  br(),
                  fluidRow(
                    column(3, textInput("cleavSite", label=div(h4(em("Cleavage Sites")), style="font-family:'hannotate tc';color:purple;font-size:12pt"), placeholder = "e.g., RK")),
                    column(2, selectInput("side", label=div(h4(em("Side")), style="font-family:'hannotate tc';color:purple;font-size:12pt"), choices=c("C","N"), selected="C", selectize = TRUE)),
                    column(2, numericInput("mc", label = div(h4(em("Missed")), style="font-family:'hannotate tc';color:purple;font-size:12pt"), value=1,min=0, step=1, width = "80%")),
                    column(3, numericInput("thresh", label=div(h4(em("Threshold")), style="font-family:'hannotate tc';color:purple; font-size:12pt"), value=4,min=0, step=1,width = "50%"))
                  ),
                  actionButton("insilico", label=span(strong("in silico Digest"), style="font-family:'wawati tc';color:purple; font-size:12pt"), icon = icon("sort-amount-desc", lib="font-awesome")),
                  hr(),
                  DT::dataTableOutput("tb3")
                  
                  ),
         tabPanel(span(strong("Targeted m/z"), style="font-family:'constantia';color:purple;font-size:10pt"),icon = icon("map-marker",lib = "font-awesome"),
                  br(),
                  fluidRow(
                    column(3, div(h4(em("Charge from")), style="font-family:'hannotate tc';color:purple; font-size:12pt", align="center")),
                    column(2, numericInput("csfrom", label=NULL,value=1, min=1, step=1, width = "80%")),
                    column(1, div(h4(em("to")), style="font-family:'hannotate tc';color:purple; font-size:12pt")),
                    column(2, numericInput("csto", label=NULL, value=5, min=2, step=1, width="80%"))
                    ),
                  fluidRow(
                    column(6, selectInput("target.ptm", label = div(h4(em("PTMs to include")),style="font-family:'hannotate tc';color:purple; font-size:14pt"), choices=list("none","Oxidation","diOxidation","triOxidation","Carbonyl","Deamidation(N)","Phosphorylation(S,T,Y)"),selected="none",selectize = TRUE, multiple = TRUE)),
                    column(3, numericInput("lowmz", label=span(h4(em("Lowest m/z")), style="font-family:'hannotate tc';color:purple; font-size:12pt"), value=150, min=50, max=600, step=1)),
                    column(3, numericInput("highmz", label=span(h4(em("Highest m/z")), style="font-family:'hannotate tc';color:purple; font-size:12pt"), value=4000, min=3000, max=5000, step=1))
                  ),
                  hr(),
                  div(h4(strong(em("User-defined PTM*"))), style="font-family:'hannotate tc';color:purple; font-size:12pt"),
                  fluidRow(
                  	column(3, textInput("NM", label=div(em("PTM name"), style="font-family:'hannotate tc';color:purple; font-size:12pt"), placeholder="e.g., FMe")),
                  	column(5, textInput("Fml", label=div(em("PTM formula"), style="font-family:'hannotate tc';color:purple; font-size:12pt"), placeholder="e.g., C16 H26")),
                  	column(3, selectInput("allRes", label=div(em("PTM targets"), style="font-family:'hannotate tc';color:purple;font-size:12pt"), choices=strsplit("ACDEFGHIKLMNPQRSTVWY", split="")[[1]], multiple=TRUE))
                  	),
             
                  fluidRow(
                    column(4,actionButton("targeted", label=span(strong("Show Targeted List"), style="font-family:'wawati tc';color:purple; font-size:12pt"),icon=icon("tags",lib="font-awesome"))),
                    column(4, downloadButton("dld.target", label=span(strong(em("Download Targeted List!")), style="font-family:'wawati tc';color:purple; font-size:12pt")))
                  ),
                  hr(),
                  DT::dataTableOutput("tb4")
               )
                  
                  
                  
      )
    )
  )
      )
  
  )
))
      
      
      
