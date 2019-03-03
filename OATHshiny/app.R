#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
conf=read.table("OATH.conf", as.is = T)
gear="java -jar gear.jar"
library(shiny)
options(shiny.maxRequestSize=conf[1,2]*1024^2, shiny.launch.browser=T)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("OATH"),
   hr(),
   fluidRow(
     column(6,
            fileInput('gfile_input',
                      paste0('3 source files (.bim, .bed, .fam) [< ', conf[1,2],' MB]'), 
                      multiple = TRUE,
                      accept = c("bed", "fam", "bim")
            )
     )
   ),

   fluidRow(
     column(6,
            fileInput('pfile_input',
                      paste0('Phenotype file')
            )
     ),
     column(6,
            textInput('pIdx', "Phenotype number", value=1)
     )
   ),

   fluidRow(
     column(6,
            fileInput('cfile_input', 
                      paste0('Covariate file')
            )
     ),
     column(6,
            textInput('cIdx', "Covariate number [1, 2, 3]", value=1)
     )
   ),
   
      # Show a plot of the generated distribution
   actionButton("run", "Run OATH", icon("refresh")),

   hr(),

   tabPanel('Pvalue-plot',
            tabPanel('EigenGWAS visualization',
                     sidebarPanel(
                       textInput("snpName", "SNP Input"),
                       selectInput('threshold',
                                   'p-value cutoff',
                                   choices = c(0.1, 0.05, 0.01, 0.005, 0.001), 
                                   selected = 0.05
                       ),
                       actionButton("Barplot", "Bar plot", icon("refresh"))
                     ),
                     mainPanel(plotOutput('barPlot'))
            )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   currentFile <-reactive({
     print("Reading files...")
     FileLoad=0
     str=""
     if(length(which(grepl("*.bed", input$gfile_input$name)))  != 1) {
       str=paste(str, "No bed file found.") 
     } else {
       FileLoad=FileLoad+1
     }
     
     if(length(which(grepl("*.bim", input$gfile_input$name)))  != 1) {
       str=paste(str, "\nNo bim file found.")
     } else {
       FileLoad=FileLoad+1
     }
     
     if(length(which(grepl("*.fam", input$gfile_input$name)))  != 1) {
       str=paste(str, "\nNo fam file found.")
     } else {
       FileLoad=FileLoad+1
     }
     
     if (FileLoad < 3) {
       showNotification(str, duration = 5, type = "error")
       return()
     } else if (FileLoad > 3) {
       showNotification("More than 3 files selected", duration = 5, type="error")
     }
     
     idx=grep(".bed$", input$gfile_input$datapath)
     if (length(idx)==1) {
       rt=substr(input$gfile_input$datapath[idx], 1, nchar(input$gfile_input$datapath[idx])-4)
     }
     for (i in 1:3) {
       if (i != idx) {
         f1 = input$gfile_input$datapath[i]
         tl = substr(f1, nchar(f1)-2, nchar(f1))
         file.symlink(f1, paste0(rt, ".", tl))
       }
     }
     
     froot = substr(input$gfile_input$datapath[idx], 1, nchar(input$gfile_input$datapath[idx])-4)
     return(froot)
   })

   currentCmd <- reactive({
     print("Reading files...")
     FileLoad=0
     str=""
     if(length(which(grepl("*.bed", input$gfile_input$name)))  != 1) {
       str=paste(str, "No bed file found.") 
     } else {
       FileLoad=FileLoad+1
     }
     
     if(length(which(grepl("*.bim", input$gfile_input$name)))  != 1) {
       str=paste(str, "\nNo bim file found.")
     } else {
       FileLoad=FileLoad+1
     }
     
     if(length(which(grepl("*.fam", input$gfile_input$name)))  != 1) {
       str=paste(str, "\nNo fam file found.")
     } else {
       FileLoad=FileLoad+1
     }
     
     if (FileLoad < 3) {
       showNotification(str, duration = 5, type = "error")
       return()
     } else if (FileLoad > 3) {
       showNotification("More than 3 files selected", duration = 5, type="error")
     }
     
     idx=grep(".bed$", input$gfile_input$datapath)
     if (length(idx)==1) {
       rt=substr(input$gfile_input$datapath[idx], 1, nchar(input$gfile_input$datapath[idx])-4)
     }
     for (i in 1:3) {
       if (i != idx) {
         f1 = input$gfile_input$datapath[i]
         tl = substr(f1, nchar(f1)-2, nchar(f1))
         file.symlink(f1, paste0(rt, ".", tl))
       }
     }
     
     froot = substr(input$gfile_input$datapath[idx], 1, nchar(input$gfile_input$datapath[idx])-4)

     gCmd=paste(gear, "oath-bus --bfile", froot, 
                "--pheno", input$pfile_input$datapath,
                "--mpheno", input$pIdx,
                "--covar", input$cfile_input$datapath,
                "--covar-number", input$cIdx,
                "--out", froot)
     return(gCmd)
   })
   
   observeEvent(input$run, {
     withProgress(message="OATH:", value=0, {
        incProgress(1/2, detail = paste0(" Running oath ... "))
       
        gCmd=currentCmd()
        if(!is.null(gCmd)) {
          print(gCmd)
          system(gCmd)
        }
        incProgress(1/2, detail = paste0(" Done ... "))
     })
   })
   
   observeEvent(input$Barplot, {
     output$barPlot <- renderPlot ({
       srcFile = currentFile()
       print(srcFile)
       info=read.table(paste0(srcFile, ".oath.info"), as.is = T, header = T)
       snpName=input$snpName
       print(snpName)
       sIdx=which(info$SNP==snpName)
       print(sIdx)
       if(length(sIdx) == 1) {
         mod=read.table(paste0(srcFile, ".oath.mod"), as.is = T, header = T)
         betaG=as.matrix(read.table(paste0(srcFile,".oath.beta"), as.is = T, header = T))
         seG=as.matrix(read.table(paste0(srcFile,".oath.se"), as.is = T, header = T))
         pMatG=-log10((1-pnorm(abs(betaG/seG)))*2)
         st=order(pMatG[sIdx,])
         barplot(pMatG[sIdx,st], border = F)
       } else {
         print("no good!")
      }
   })
  })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

