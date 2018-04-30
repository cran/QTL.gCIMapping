ui <-tagList(
  navbarPage(
    "",id = "tabs",
    tabPanel(strong("QTL.gCIMapping"),
             
             h2("QTL.gCIMapping (QTL genome-wide Composite Interval Mapping)",align="center"),
             column(3,
                    br(),
                    h4(strong("Coding criteria")),
                    offset=4
             ),
             column(4,
                    br(),
                    h4(strong("Dataset example")),
                    offset=1
             ),
             column(4,     
                    tableOutput("codeexample"),
                    offset = 3
             ),
             column(4,      
                    tableOutput("dataexample")
             ),
             column(12, 
                    br(),
                    h4(strong("Reference"),align="center")
             ),
             
             
             
             column(7,
                    
                    h4("1. Wang Shi-Bo, Wen Yang-Jun, Ren Wen-Long, Ni Yuan-Li, Zhang Jin, Feng Jian-Ying, Zhang Yuan-Ming*.  
                  Mapping small-effect and linked quantitative trait loci for complex traits in backcross or DH populations via
                  a multi-locus GWAS methodology. Scientific Reports 2016,6:29951."),
                    h4("2. Wen Yang-Jun, Zhang Ya-Wen, Zhang Jin, Feng Jian-Ying, Jim M. Dunwell, Zhang Yuan-Ming*. An
                   efficient multi-locus mixed model framework for the detection of small and linked QTLs in F2. Submitted"),
                    br(),
                    br(),
                    h4("Authors: Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, Zhang Yuan-Ming"),
                    h4("Maintainer: Zhang Yuan-Ming (soyzhang at mail.hzau.edu.cn)"), 
                    h4("QTL.gCIMapping version 2.0, Realeased April 2018"),
                    offset = 3)
    ),
    
    
    tabPanel(strong("Start"),
             
             titlePanel("QTL.gCIMapping (QTL genome-wide Composite Interval Mapping)"),
             sidebarLayout(
               sidebarPanel(
                 radioButtons("dataformat", "Please select data format", choices = c("GCIM","WinQTLCart","QTLIciMapping"),selected="GCIM"),
                 fileInput("fileDataset", "Input dataset",multiple = TRUE),
                 conditionalPanel("input.dataformat == 'QTLIciMapping'",
                                  fileInput("fileCov", "Input covariate file",multiple = TRUE)
                 ),
                 selectInput("fileType","Show dataset:",choices=c("Genotype","Phenotype","Linkage map","Covariate")),
                 br(),
                 br(),
                 br(),
                 actionButton("parSe", label = "Parameter Settings",width=250,style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 br(),
                 br(),
                 br(),
                 actionButton("drawPl", label = "Figure",width = 250,style="color: #fff; background-color: #337ab7; border-color: #2e6da4"), 
                 br(),
                 br(),
                 br(),
                 actionButton("manl", label = "User manual",width=250,style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 uiOutput("manll")
               ),
               mainPanel(
                 tabsetPanel(id="inTabset",
                             tabPanel("Dataset",value = "DA",
                                      conditionalPanel("input.fileType == 'Genotype'",
                                                       #h3("Genotype"),
                                                       dataTableOutput("genTable")),
                                      conditionalPanel( "input.fileType == 'Phenotype'",
                                                        #h3("Phenotype"),
                                                        dataTableOutput("pheTable")),
                                      conditionalPanel( "input.fileType == 'Linkage map'",
                                                        #h3("Linkage map"),
                                                        dataTableOutput("mapTable")),
                                      conditionalPanel( "input.fileType == 'Covariate'",
                                                        #h3("Covariate"),
                                                        dataTableOutput("covTable")),
                                      conditionalPanel( "input.dataformat == 'WinQTLCart'",
                                                        #h3("Dataset"),
                                                        verbatimTextOutput("TableWIN"))
                             ),
                             tabPanel("Parameter Settings",value = "PA",
                                      fluidRow(  
                                        br(), 
                                        br(),
                                        column(4,  
                                               radioButtons("Pop", "Please select Population", choices = c("DH","RIL","BC1","BC2","F2"),selected="DH",inline = "TRUE")
                                        ),
                                        column(4,  
                                               radioButtons("model", "Please select Model", choices = c("Random model","Fixed model"),inline = "TRUE")
                                        ),
                                        column(4,  
                                               textInput("Walk", "Walk Speed for Genome-wide Scanning (cM):",value="1")
                                        ),
                                        column(4,
                                               textInput("Crilod", "Critical LOD score",value="2.5")
                                        ),
                                        column(4,
                                               br(),
                                               radioButtons("likelihood","Likelihood function (only for F2):",choices=c("REML","ML"),inline = "TRUE")
                                        ),
                                        column(4,
                                               br(),
                                               radioButtons("Rqtl","Completing CIM in the neighborhood (only for F2)",choices=c("TRUE","FALSE"),selected="FALSE",inline = "TRUE")
                                        ),
                                        column(4,
                                               br(),
                                               radioButtons("drawplot","Draw plot or not",choices=c("TRUE","FALSE"),inline = "TRUE")
                                        ),
                                        
                                        column(4,
                                               br(),      
                                               radioButtons("resolution","Resolution of plot",choices=c("Low","High"),inline = "TRUE")
                                        ),
                                        
                                        column(4,
                                               br(),       
                                               radioButtons("Plotformat","Plot format",choices=c("*.png","*.tiff","*.jpeg","*.pdf"),inline = "TRUE")
                                        ),
                                        
                                        column(4,
                                               h5(strong("Select trait ID")),
                                               textInput("Trait1","From",value="1")
                                        ),
                                        column(4,
                                               br(),
                                               br(),
                                               textInput("Trait2","To",value="1")
                                        ),
                                        
                                        column(4,
                                               br(),
                                               br(),
                                               textInput("SavePath", "Save path",value = "C:/Users/Administrator/Desktop")
                                        ),
                                        
                                        
                                        column(12, 
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               actionButton("run", label = "Run",width=280, icon("paper-plane"), 
                                                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                        ),
                                        column(12,
                                               br(),
                                               br(),
                                               uiOutput("result")
                                        )
                                      )
                             ),
                             
                             tabPanel("Figure",value = "FIG",
                                      h3("Genome-wide composite interval mapping (GCIM) figure"),
                                      fluidRow(
                                        column(12,
                                               radioButtons("Mose", "", c("Parameter Settings", "Draw plot"),inline = TRUE)
                                        ),
                                        
                                        conditionalPanel("input.Mose == 'Parameter Settings'",
                                                         column(12,
                                                                br(),
                                                                radioButtons("Resolution2", "Select resolution of plot", c("General resolution", "High resolution", "Set by yourself"),inline = TRUE)
                                                         ),
                                                         column(3,
                                                                textInput("widthGCIM","Width (px):",value="1500")
                                                         ),
                                                         column(3,
                                                                textInput("heightGCIM","Height (px):",value="600")
                                                         ),
                                                         column(3,
                                                                textInput("pointGCIM","Word resolution (1/72 inch, ppi):",value="12")
                                                         ),
                                                         column(3,
                                                                textInput("ppiGCIM","Figure resolution (ppi):",value="72")
                                                         ),
                                                         column(3,
                                                                textInput("leGCIM","Legend and tick marks:",value="1.0")
                                                         ),
                                                         column(3,
                                                                textInput("maGCIM","LOD line size:",value="1.0")
                                                         ),
                                                         column(3,
                                                                textInput("baGCIM","Size for -log10(P) curve:",value="0.5")
                                                         ),
                                                         column(3,
                                                                textInput("marGCIM","Margin space:",value="1.5")
                                                         ),
                                                         column(3,
                                                                textInput("axisGCIM","Space between tick marks and axis:",value="1.0")
                                                         ),
                                                         column(3,
                                                                textInput("logGCIM","Times for max{-log10(P)}:",value="1.5")
                                                         ),
                                                         
                                                         column(6,
                                                                textInput("lodGCIM","Critical LOD score:",value="2.5")
                                                         ),
                                                         
                                                         column(3,
                                                                selectInput("co1GCIM","LOD line color:",choices=c("red","black","blue","yellow","green","pink","purple","gray50","brown"))
                                                         ),
                                                         column(3,
                                                             
                                                                selectInput("co2GCIM","-log10(P) curve color1:",choices=c("gray50","black","blue","yellow","green","pink","purple","red","brown"))
                                                         ),
                                                         column(3,
                                                                
                                                                selectInput("co3GCIM","-log10(P) curve color2 (only for F2):",choices=c("green","black","blue","yellow","red","pink","purple","gray50","brown"))
                                                         )
                                        ),
                                        conditionalPanel("input.Mose == 'Draw plot'",
                                                         
                                                         column(6,
                                                                br(),
                                                           
                                                                selectInput("plotpopGCIM", "Select population",choices = c("DH BC1 BC2 RIL","F2"))   
                                                         ),                       
                                                         column(6,
                                                                br(),
                                                             
                                                                
                                                                fileInput("fileplotGCIM", "Input file to draw plot",multiple = TRUE,accept = ".xlsx")
                                                         ),   
                                                         
                                                         column(12,
                                                                radioButtons("Plotformat2","Plot format",choices=c("*.png","*.tiff","*.jpeg","*.pdf"),inline = "TRUE")
                                                         ),
                                                         column(12,
                                                                downloadButton("downloadplotGCIM", "Download plot",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                                         ),
                                                         column(12,
                                                                plotOutput("plotGCIM")
                                                         )
                                        )
                                      )
                             )
                 )
               )
             )
    )
    )
)