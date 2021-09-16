#
#
# ngsAnalyser is a Electron-Shiny application. 
# version: 1.1.4
#

packages <- c("shiny", "shinyWidgets", "shinyjs")

suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))

#load necessary packages

# options(shiny.maxRequestSize = 100*30*1024^2) 

eSPANwrapper <- source("eSPANwrapper.R")
ChIPwrapper <- source("ChIPwrapper.R")

# Define UI ----

# Edit this name if desired when starting a new app


VIRTUALENV_NAME = 'minipython'


# ------------------------- Settings (Do not edit) -------------------------- #

if (Sys.info()[['user']] == 'shiny'){

  # Running on shinyapps.io
  Sys.setenv(PYTHON_PATH = 'python3')
  Sys.setenv(VIRTUALENV_NAME = VIRTUALENV_NAME) # Installs into default shiny virtualenvs dir
  Sys.setenv(RETICULATE_PYTHON = paste0('/home/shiny/.virtualenvs/', VIRTUALENV_NAME, '/bin/python'))

} else if (Sys.info()[['user']] == 'rstudio-connect'){

  # Running on remote server
  Sys.setenv(PYTHON_PATH = '/opt/python/3.7.7/bin/python3')
  Sys.setenv(VIRTUALENV_NAME = paste0(VIRTUALENV_NAME, '/')) # include '/' => installs into rstudio-connect/apps/
  Sys.setenv(RETICULATE_PYTHON = paste0(VIRTUALENV_NAME, '/bin/python'))

} else {

  # Running locally
  options(shiny.port = 7450)
  Sys.setenv(PYTHON_PATH = 'python3')
  Sys.setenv(VIRTUALENV_NAME = VIRTUALENV_NAME) # exclude '/' => installs into ~/.virtualenvs/
  # RETICULATE_PYTHON is not required locally, RStudio infers it based on the ~/.virtualenvs path
 }


ui <- navbarPage(
  
    title = "ngsAnalyser", 
    
    
    ChIP_page <- tabPanel(
      
      setBackgroundColor(
        color = c("#2D2929", "#756E6E"),
        gradient = "radial",
        direction = c("top", "left")
      ),
      
      tags$head(tags$style(".shiny-notification {position: fixed; top: 77% ;left: 54%; right:1%; color: blue;font-size: 17px;font-style: italic;")),
      
      title = "ChIP Analysis",
      
                     fluidRow(
                       
                       column(4, offset = -1,
                              
                                sidebarPanel(style = "max-height: calc(100vh - 18rem); overflow-y: auto;",
                                             
                                  textInput("ChIPtitle", "Sample Name", placeholder = "Pol2/Pol2-HU/Pol2-30min"), 
                                  
                                  radioButtons("SamplType", "Sample Type",
                                               choiceNames = list(
                                                 "ChIP", "BrDU"
                                               ),
                                               choiceValues = list(
                                                 "ChIP", "BrDU"
                                               ),
                                               inline = TRUE),
                                  
                                  radioButtons("seqtype", "Sequencing Strategy",
                                               choiceNames = list(
                                                 "Paired-end", "Single-end" 
                                               ),
                                               choiceValues = list(
                                                 "paired", "single"
                                               ),
                                               inline = F),
                                  
                                  ##
                                  
                                  radioButtons("cw", "Coverage Window",
                                               choiceNames = list(
                                                 "Sliding", "Non-Overlapping"
                                               ),
                                               choiceValues = list(
                                                 "YES", "NO"
                                               ),
                                               inline = F),
                                  
                                  textInput("ChIPbin", "Bin (default recommended)", placeholder = "defaults to MFS"),
                                  
                                  textInput("ChIPsteps", "Slider steps (sliding windows)", placeholder = "10"),
                                  
                                  ##
                                  
                                  radioButtons("ChIPprofile", "Plot IP Profile",
                                               choiceNames = list(
                                                 "Genomewide", "Fixed-coordinates", "Both"
                                               ),
                                               choiceValues = list(
                                                 "Global", "Local", "Both"
                                               ),
                                               inline = F),
                                  
                                  textInput("textCoords", "Coordinates (chr:start-end)", placeholder = "4:50000-150000"),
                                  
                                  ##
                                  
                                  radioButtons("peaklist", "Plot Average Profiles At",
                                               choiceNames = list(
                                                 "Fired Origins", "Genomewide Peaks", "Peaks Outside Fired Origins", "TSS (not available!)"
                                               ),
                                               choiceValues = list(
                                                 "RO", "ALL", "ALL-RO", "TSS"
                                               ),
                                               inline = F),
                                  
                                  
                                  textInput("ChIPaveWin", "Range (upstream:downstream)", placeholder = "e.g. 3000:3000"),
                                  
                                  radioButtons("Norm", "Normalise against Input",
                                               choiceNames = list(
                                                 "Yes", "No"
                                               ),
                                               choiceValues = list(
                                                 "YES", "NO"
                                               ),
                                               inline = T),
                                  
                                  
                                  width = 20
                                ),
                                mainPanel()
                              
                              
                              ),
                       
                       column(8, offset = -1, style='padding-left:0px; padding-top:0px;',
                            
                                sidebarPanel(style="max-height: 113px;",
                                  
                                  tabsetPanel(
                                    
                                    tabPanel("Input_R1", 
                                             tags$style("li a {font-weight: bold;}"),
                                             textInput("text9", h1(), value = NULL, placeholder = "/Enter/full/path/to/fastq.gz ...")),
                                    
                                    tabPanel("Input_R2", textInput("text10", h1(),  value = NULL, placeholder = "/Enter/full/path/to/fastq.gz ...")),
                                    
                                    tabPanel("IP_R1", textInput("text11", h1(),  value = NULL, placeholder = "/Enter/full/path/to/fastq.gz ...")),
                                    
                                    tabPanel("IP_R2", textInput("text12", h1(),  value = NULL, placeholder = "/Enter/full/path/to/fastq.gz ..."))
                                    
                                    
                                  ),
                                  width = 30
                                ),
                                mainPanel( 
                                  shinyjs::useShinyjs(),
                                  textOutput("ChIP")
                                )
                              
                              
                              )
                       
                       
                   ),
      
      sidebarLayout(position = "left",
                    
                    sidebarPanel(style = "max-height:10px; left: 0px; right: 0px; margin-bottom: 0px",
                                 
                                 fluidRow(
                                   column(4, align="left", offset = 0,
                                          actionButton("go", "Run Analysis", icon=icon("play-circle"), width = 250, 
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4; font-weight: bold;")),
                                   
                                   column(4, align="center", offset = 0,
                                          useShinyjs(),
                                          extendShinyjs(text = "shinyjs.refresh = function() { {history.go(0); }", functions = c("refresh")),
                                          actionButton("cRefresh", "Refresh App", icon=icon("refresh"), width = 250, 
                                                       style="color: #fff; background-color: #33cc33; border-color: #33cc33; font-weight: bold;")),
                                   
                                   column(4, align="right", offset = 0,
                                          useShinyjs(),
                                          extendShinyjs(text = "shinyjs.closeWindow = function() { window.close(); }", functions = c("closeWindow")),
                                          actionButton("Cclose", "Quit App", icon = icon("sign-out-alt"), width = 250, 
                                                       style="color: #fff; background-color: #ff0000; border-color: #ff0000; font-weight: bold;"))
                                 ),
                                 width = 0
                    ),
                    mainPanel( )
      ),
      
      tags$style("#text9 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
      
      tags$style("#text10 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
      
      tags$style("#text11 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
      
      tags$style("#text12 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
  
      tags$style("#go {
                    margin-top: -23px;
                    margin-left: -11px;
                    }"),
      
      tags$style("#Cclose {
                    margin-top: -23px;
                    margin-left: -11px;
                    }"),
      
      tags$style("#cRefresh {
                    margin-top: -23px;
                   margin-left: -11px;
                    }"),
      
      tags$style("#ChIP {font-size:20px;
               display:block;
               color:red;
               text-align:left;
               position:absolute;
               width: 100%;
               left:50px;
               top:20px;}")
      
    ),
    
    ####
    
    eSPAN_page <- tabPanel(
      
      setBackgroundColor(
        color = c("#2D2929", "#756E6E"),
        gradient = "radial",
        direction = c("top", "left")
      ),
      
      tags$head(tags$style(".shiny-notification {position: fixed; top: 77% ;left: 54%; right:1%; color: blue;font-size: 17px;font-style: italic;")),
      
      title = "eSPAN Analysis",
      
      fluidRow(
        
        column(4, offset = -1, 
               
               sidebarPanel(style = "max-height: calc(100vh - 18rem); overflow-y: auto;",
                 
                 textInput("eSPANtitle", "Sample Name", placeholder = "Pol2/Pol2-HU/Pol2-30min"), 
                 
                 radioButtons("ew", "Coverage Window",
                              choiceNames = list(
                                "Sliding", "Non-Overlapping"
                              ),
                              choiceValues = list(
                                "YES", "NO"
                              ),
                              inline = F),
                 
                 textInput("eBin", "Bin (default recommended)", placeholder = "defaults to mean fragment size"),
                 
                 textInput("eStep", "Slider steps (sliding windows)", placeholder = "10"),
                 
                 ##
                 
                 radioButtons("eProfile", "Plot IP Profile",
                              choiceNames = list(
                                "Genomewide", "Fixed-coordinates", "Both"
                              ),
                              choiceValues = list(
                                "Global", "Local", "Both"
                              ),
                              inline = F),
                 
                 textInput("textCoords", "Coordinates (chr:start-end)", placeholder = "4:50000-150000"),
                 
                 ##
                 
                 radioButtons("eResults", "Plot Results",
                              choiceNames = list(
                                HTML("enrichment average"), 
                                HTML("watson/crick average"),
                                HTML("Strand-bias at ROs"),
                                HTML("all")
                              ),
                              choiceValues = list(
                                "Enrichment", "wat-cri", "bias-ROs", "all"
                              ),
                              inline = F,
                              selected = "all"),
                 
                 
  
                 width = 20
               ),
               mainPanel()
               
               
        ),
        
        column(8, offset = -1, style='padding-left:0px; padding-top:0px;',
               
               sidebarPanel(style = "max-height:203px;",
                 
                 tabsetPanel(
                   
                   tabPanel("Input_R1", 
                            tags$style("li a {font-weight: bold;}"),
                            textInput("text1", h1(), placeholder = "/Enter/full/path/to/fastq/files...")),
                   
                   tabPanel("Input_R2", textInput("text2", h1(), placeholder = "/Enter/full/path/to/fastq/files...")),
                   
                   tabPanel("ChIP_R1", textInput("text3", h1(), placeholder = "/Enter/full/path/to/fastq/files...")),
                   
                   tabPanel("ChIP_R2", textInput("text4", h1(), placeholder = "/Enter/full/path/to/fastq/files..."))
                   
                   
                 ),
                 
                 tabsetPanel(

                   tabPanel("BrDU_R1",
                            tags$style("li a {font-weight: bold;}"),
                            textInput("text5", h1(), placeholder = "/Enter/full/path/to/fastq/files...")),

                   tabPanel("BrDU_R2", textInput("text6", h1(), placeholder = "/Enter/full/path/to/fastq/files...")),

                   tabPanel("eSPAN_R1", textInput("text7", h1(), placeholder = "/Enter/full/path/to/fastq/files...")),

                   tabPanel("eSPAN_R2", textInput("text8", h1(), placeholder = "/Enter/full/path/to/fastq/files..."))


                 ),
                 
                 width = 30
               ),
               mainPanel( 
                 shinyjs::useShinyjs(),
                 textOutput("eSPAN")
               )
               
               
        )
        
        
      ),
      
      sidebarLayout(position = "left",
                    
                    sidebarPanel(style = "max-height:10px; margin-bottom: 0px;",
                                 
                                 fluidRow(
                                   column(4, align="left", offset = 0,
                                          actionButton("run", "Run Analysis", icon=icon("play-circle"), width = 250, 
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4; font-weight: bold;")),
                                   
                                   column(4, align="center", offset = 0,
                                          useShinyjs(),
                                          extendShinyjs(text = "shinyjs.refresh = function() { {history.go(0); }", functions = c("refresh")),
                                          actionButton("eRefresh", "Refresh App", icon=icon("refresh"), width = 250, 
                                                       style="color: #fff; background-color: #33cc33; border-color: #33cc33; font-weight: bold;")),
                                   
                                   column(4, align="right", offset = 0,
                                          useShinyjs(),
                                          extendShinyjs(text = "shinyjs.closeWindow = function() { window.close(); }", functions = c("closeWindow")),
                                          actionButton("Eclose", "Quit App", icon = icon("sign-out-alt"), width = 250, 
                                                       style="color: #fff; background-color: #ff0000; border-color: #ff0000; font-weight: bold;"))
                                 ),
                                 width = 0
                    ),
                    mainPanel( )
      ),
      
      tags$style("#text1 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
      
      tags$style("#text2 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
      
      tags$style("#text3 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
      
      tags$style("#text4 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
      
      tags$style("#text5 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
      
      tags$style("#text6 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
      
      tags$style("#text7 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
      
      tags$style("#text8 {
                    white-space: pre-wrap;
                    margin-top: -30px;
                    }"),
    
      tags$style("#run {
                    margin-top: -23px;
                    margin-left: -11px;
                    }"),
      
      tags$style("#Eclose {
                    margin-top: -23px;
                    margin-left: -11px;
                    }"),
      
      tags$style("#eRefresh {
                    margin-top: -23px;
                    margin-left: -11px;
                    }"),
      
      tags$style("#eSPAN {font-size:20px;
               display:block;
               color:red;
               text-align:left;
               position:absolute;
               width: 100%;
               left:50px;
               top:50px;}")
      
    ),
    
    
    
    # footer
  
    tags$footer(HTML("<br>
                    <br>
                    <!-- Footer -->
                           <footer class='page-footer font-large indigo'>
                           <!-- Copyright -->
                           <div style='left:0px;right:0px;bottom:10px;position:fixed;cursor:inherit;'>
                           <div class='footer-copyright text-center py-2'>ngsAnalyser-1.1.4 © 2021 :
                           <a href='https://www.cib.csic.es/research/cellular-and-molecular-biology/dna-replication-and-genome-integrity' target='_blank'> DNA Replication & Genome Integrity Lab</a>
                           </div>
                           <!-- Copyright -->
                           </footer>
                           <!-- Footer -->")),
    
    
    inverse = T,
    
    fluid = T,
    
    collapsible = T
    
  ) 
  
  
# Define any Python packages needed for the app here:
PYTHON_DEPENDENCIES = c('pip', 'MACS2', 'numpy')

# Define server logic ----
server <- function(input, output, session) {
  
  # ------------------ App virtualenv setup (Do not edit) ------------------- #
  
  virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
  python_path = Sys.getenv('PYTHON_PATH')
  
  # Create virtual env and install dependencies
  reticulate::virtualenv_create(envname = virtualenv_dir, python = python_path)
  reticulate::virtualenv_install(virtualenv_dir, packages = PYTHON_DEPENDENCIES, ignore_installed=TRUE)
  reticulate::use_virtualenv(virtualenv_dir, required = T)
  
  Sys.setenv(
    PATH = paste("~/.virtualenvs/minipython/bin", Sys.getenv("PATH"), sep = .Platform$path.sep),
    VIRTUAL_ENV = "~/.virtualenvs/minipython/venv"
  )
  Sys.unsetenv("PYTHONHOME") # works whether previously set or not
  
  observeEvent(input$run, {
    withCallingHandlers({
      shinyjs::html("eSPAN", "")
      
      withProgress(message = 'Analysis in progress... (this will take time)',
                   detail = 'Go get some coffee...', value = 0, {
                     
                     eSPANwrapper(
                                  Input_R1 = input$text1, 
                                  Input_R2 = input$text2,
                                  ChIP_R1 = input$text3, 
                                  ChIP_R2 = input$text4,
                                  BrDU_R1 = input$text5, 
                                  BrDU_R2 = input$text6,
                                  eSPAN_R1 = input$text7, 
                                  eSPAN_R2 = input$text8,
                                  
                                  ExpTitle = input$eSPANtitle,
                                  slidingWindow = input$ew, 
                                  bin = input$eBin,
                                  stepSize = input$eStep, 
                                  PlotIPprofile = input$eProfile,
                                  PlotAveProfile = input$eResults
                     
                     )
                   })
      
      
    },
    message = function(m) {
      shinyjs::html(id = "eSPAN", html = paste0(m$message,  sep="<br/>"), add = T)
    })
   
  })
  
  observeEvent(input$go, {
    withCallingHandlers({
      shinyjs::html("ChIP", "")
      
      withProgress(message = 'Analysis in progress... (this will take time)',
                   detail = 'Go get some coffee...', value = 0, {
                     ChIPwrapper(
                       Input_R1 = input$text9, 
                       Input_R2 = input$text10, 
                       IP_R1 = input$text11, 
                       IP_R2 = input$text12,
                       
                       ExpTitle = input$ChIPtitle,
                       SamplType = input$SamplType,
                       SeqType = input$seqtype,
                       
                       slidingWindow = input$cw, 
                       
                       bin = input$ChIPbin,
                       stepSize = input$ChIPsteps, 
                       
                       PlotIPprofile = input$ChIPprofile,
                       ChromCoords = input$textCoords, 
                       
                       peakSet = input$peaklist,
                       AveragingPan = input$ChIPaveWin, 
                       Normalise = input$Norm
                       
                     )
                   })
      
    },
    message = function(m) {
      shinyjs::html(id = "ChIP", html = paste0(m$message,  sep="<br/>"), add = T)
    })
  })
  
  observeEvent(c(input$Cclose, input$Eclose), {
    js$closeWindow()
    stopApp()
  }, ignoreNULL = T, ignoreInit = T)
  
  observeEvent(c(input$cRefresh,input$eRefresh), {
    js$refresh()
  }, ignoreNULL = T, ignoreInit = T)
  
}

# Run the app ----
shinyApp(ui = ui, server = server)


