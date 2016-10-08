library(shiny)
library(DT)

shinyUI(fluidPage(
  titlePanel(
  fluidRow(column(8, div(p(strong("Proteomics Toolkit")),style="font-family:'kokonor';color:black; font-size:30pt"),
                    h4(em(strong("from", span("Gross Lab - Ben NIU", style="font-family:'gabriola';color:blue; font-size:15pt"))))),
             column(4, img(src="wustl_name.png", width=210, height=70, align="right"))
  )),
  sidebarLayout(
    sidebarPanel(
      textInput("seq", label=div(em("Input peptide/protein sequence"),style="font-family:'marker felt';color:darkgreen; font-size:12pt"), placeholder = "for example: YGGFL"),
      radioButtons("Cys", label=div(em("Carbamiodomethyl on Cys ?"),style="font-family:'marker felt';color:darkgreen; font-size:12pt"), choices = list("Yes","No"), selected = "Yes"),
      textInput("ptmFormula", label = div(p(em(span("PTM in formula",style="font-family:'marker felt';color:darkgreen; font-size:12pt"))),  p('Note: A "SPACE" is required between every element! The allowed elements are', span('C, H, N, O, S, P, Br, Cl, F, Si.', style="color:red"))), placeholder = "e.g. C5 H10 O2 N Br2"),
      hr(),
      div(p(strong(em("Residue-specific PTMs"))), style="font-family:'marker felt';color:darkgreen; font-size:12pt"),
      fluidRow(column(7,uiOutput("PTM1")),
               column(4, uiOutput("Resi1"))),
      fluidRow(column(7, uiOutput("PTM2")),
               column(4, uiOutput("Resi2"))),
      fluidRow(column(7, uiOutput("PTM3")),
               column(4, uiOutput("Resi3"))),
      checkboxGroupInput("chkbx",label=div(strong(em("extra PTMs")),style="font-family:'marker felt';color:darkgreen; font-size:12pt"),choices=list("none","Acetyl", "Ammonium-add","Ammonium-loss","Dehydrated"), selected = "none")
    ),
    
    mainPanel(
      tabsetPanel(type="tabs",
          tabPanel(h5("Statistics"),
         br(),
      fluidRow(
        column(5, div(h4(strong(p("Sequence"))), style="font-family:'chalkboard se';color:darkblue")),
        column(5, div(h4(strong(p("Formula"))), style="font-family:'chalkboard se';color:darkblue")),
        column(2, div(h4(strong(p("Res.#"))), style="font-family:'chalkboard se';color:darkblue"))
        ),
      fluidRow(
        column(5, verbatimTextOutput("text1")),
        column(5, verbatimTextOutput("text2")),
        column(2, verbatimTextOutput("text3"))
        ),
      br(),
      tableOutput("tb1"),
      hr(),
      fluidRow(
      column(3, actionButton("act", label=div(strong("Composition"), style="font-family:'calibri';color:#3399FF; font-size:12pt"))),  
      column(8,DT::dataTableOutput("tb2"))
      )
          ),
         tabPanel(h5("Isotopic distribution"),
                  br(),
               div(h4(strong(p("Formula with PTMs"))), style="font-family:'chalkboard se';color:darkblue"),
               fluidRow(
                 column(6, verbatimTextOutput("text4"))
                 ),
               hr(),
              fluidRow(
              column(6, div(h4(p(strong("Charge State to Simulate: "))), style="font-family:'chalkboard se';color:darkblue")),
              column(2, numericInput("ecs", label=NULL, value=1,width='80%', min=1, step=1)),
              column(2, actionButton("sim", label = div(strong("Simulate"), style="font-family:'calibri';color:#3399FF; font-size:12pt")))
              ),
              div(h5(p(em("Allowed elements for simulation are C,H,N,O,S,P,Br,Cl,F,Si."))), style="font-family:'comic sans ms';color:purple;font-size:12pt"),
              hr(),
              plotlyOutput("Idplot")
                  ),
         tabPanel(h5(strong(em("in silico")), "Digestion"),
                  br(),
                  fluidRow(
                    column(3, textInput("cleavSite", label=div(strong("Cleavage Sites"), style="font-family:'chalkboard se';color:darkblue;font-size:12pt"), placeholder = "e.g., RK")),
                    column(2, textInput("side", label=div(strong("Side"), style="font-family:'chalkboard se';color:darkblue;font-size:12pt"), placeholder = "C or N")),
                    column(2, numericInput("mc", label = div(strong("Missed"), style="font-family:'chalkboard se';color:darkblue;font-size:12pt"), value=1,min=0, step=1, width = "80%")),
                    column(3, numericInput("thresh", label=div(strong("Threshold"), style="font-family:'chalkboard se';color:darkblue; font-size:12pt"), value=4,min=0, step=1,width = "50%"))
                  ),
                  actionButton("insilico", label=div(strong("in silico Digest"), style="font-family:'calibri';color:#3399FF; font-size:12pt")),
                  hr(),
                  DT::dataTableOutput("tb3")
                  
                  ),
         tabPanel(h5("Targeted m/z"),
                  br(),
                  fluidRow(
                    column(3, div(h4(strong("Charge from")), style="font-family:'chalkboard se';color:darkblue; font-size:12pt", align="center")),
                    column(2, numericInput("csfrom", label=NULL,value=1, min=1, step=1, width = "80%")),
                    column(1, div(h4(strong("to")), style="font-family:'chalkboard se';color:darkblue; font-size:12pt")),
                    column(2, numericInput("csto", label=NULL, value=5, min=2, step=1, width="80%")),
                    column(4, checkboxGroupInput("target.ptm", label = div(h4(strong("PTMs to include")),style="font-family:'chalkboard se';color:darkblue; font-size:12pt"), choices=list("none","Oxidation","diOxidation","triOxidation","Carbonyl","Deamidation(N)","Phosphorylation(S,T,Y)"),selected="none"))
                  ),
                  fluidRow(
                    column(4, numericInput("lowmz", label=span(strong("Lower Shreshold", style="font-family:'chalkboard se';color:darkgreen; font-size:12pt")), value=150, min=50, max=600, step=1)),
                    column(4, numericInput("highmz", label=span(strong("Higher Shreshold", style="font-family:'chalkboard se';color:darkgreen; font-size:12pt")), value=4000, min=3000, max=5000, step=1))
                    ),
                  fluidRow(
                    column(4,actionButton("targeted", label=div(strong("Show Targeted List"), style="font-family:'calibri';color:#3399FF; font-size:12pt"))),
                    column(4, downloadButton("dld.target", label=span(em("Download Targeted List!"), style="font-family:'calibri';color:#3399FF; font-size:12pt")))
                  ),
                  hr(),
                  DT::dataTableOutput("tb4")
               )
                  
                  
                  
      )
    )
  )
))
      
      
      
