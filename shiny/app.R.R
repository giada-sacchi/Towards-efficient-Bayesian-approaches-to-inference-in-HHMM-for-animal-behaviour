##### GLOBAL
library(shiny)
library(shinythemes)
library(ggplot2)
library(forecast)
library(gridExtra)
library(grid)

# Create an array; each dimension corresponds to a different dataset
itns <- array(NA, dim=c(10000,51,8))

itns[,,1] <- readRDS("itns_MH.rds")
itns[,,2] <- cbind(itns[,1:30,1],t(readRDS("itns_block1.rds")))
itns[,,3] <- cbind(itns[,1:30,1],t(readRDS("itns_block2.rds")))
itns[,,4] <- cbind(itns[,1:30,1],t(readRDS("itns_MH1.rds")))
itns[,,5] <- cbind(itns[,1:30,1],t(readRDS("itns_MH2.rds")))
itns[,,6] <- t(readRDS("itns_pt4.rds"))
itns[,,7] <- t(readRDS("itns_pt7.rds"))
itns[,,8] <- t(readRDS("itns_MH_tot.rds"))

##### UI #####
ui <- fluidPage(theme = shinytheme("sandstone"),
                navbarPage(strong("BAYESIAN ESTIMATION"),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    radioButtons("algo",h4("Metropolis-Hastings Algorithm",style="color:lightseagreen"),
                                 c("Single Updates" = 1,
                                   "Single Updates (Complete Parameter Set)" = 8,
                                   "Block Updates 1" = 2,
                                   "Block Updates 2" = 3,
                                   "Single Updates 1" = 4,
                                   "Single Updates 2" = 5,
                                   "Parallel Tempering with 4 chains" = 6,
                                   "Parallel Tempering with 7 chains" = 7),
                                 selected = 1),
                    br(),
                    p(strong("Description of the selected algorithm"),style="color:lightseagreen"),
                    textOutput("descr")),
                  mainPanel(
                    tabsetPanel(type="tabs",
                                
                                tabPanel("TAB 1: Covariates",
                                         fluidRow(position="top",
                                                  column(6,
                                                         br(),
                                                         selectInput("cov",p("Covariate",style="color:cornflowerblue"),
                                                                     choices = list(
                                                                       "Dive Duration" = 31,
                                                                       "Maximum Depth" = 37,
                                                                       "Dive Wiggliness" = 43),
                                                                     selected = 37)),
                                                  column(6,
                                                         br(),
                                                         selectInput("par",p("Parameter",style="color:cornflowerblue"),
                                                                     choices = list("Mean" = 0,
                                                                                    "Standard Deviation" = 3,
                                                                                    "Proportion of 0's" = 6),
                                                                     selected = 3))),
                                         br(),
                                         p(strong("Production State 1"),style="color:lightseagreen"),
                                         plotOutput("plot1",height="auto"),
                                         br(),
                                         p(strong("Production State 2"),style="color:lightseagreen"),
                                         plotOutput("plot2",height="auto"),
                                         br(),
                                         p(strong("Production State 3"),style="color:lightseagreen"),
                                         plotOutput("plot3",height="auto")),
                                
                                tabPanel("TAB 2: Transition Probabilities",
                                         
                                         fluidRow(position="top",
                                                  column(6,
                                                         br(),
                                                         selectInput("lev",p("Level",style="color:cornflowerblue"),
                                                                     choices = list(
                                                                       "Lower Level" = 1,
                                                                       "Upper Level 1" = 7,
                                                                       "Upper Level 2" = 10),
                                                                     selected = 1)),
                                                  column(6,
                                                         br(),
                                                         selectInput("type",p("Parameter Type",style="color:cornflowerblue"),
                                                                     choices = list("Transition Probability Matrix" = 2,
                                                                                    "Stationary Distribution" = 0),
                                                                     selected = 0))),
                                         conditionalPanel(condition="(input.lev==10 & input.type==0)|(input.lev==7 & input.type==0)",
                                                          plotOutput("plotTPM1",height="auto")),
                                         conditionalPanel(condition="(input.lev==10 & input.type==2)|(input.lev==7 & input.type==2)",
                                                          plotOutput("plotTPM2",height="auto")),
                                         conditionalPanel(condition="(input$lev==1 & input$type==0)", 
                                                          plotOutput("plotTPM0",height="auto")),
                                         conditionalPanel(condition="(input$lev==1 & input$type==2)", 
                                                          plotOutput("plotTPM3",height="auto")))
                    )))))



##### SERVER #####
server <- function(input, output) {
  
  ### descriptions
  state <- reactiveValues()
  
  observe({
    state$x <- input$algo
    state$y <- ifelse(state$x  == 1, "Original single-update MH algorithm. Executed keeping fixed the transition probability matrices and the corresponding stationary distributions.",
                      ifelse(state$x  == 2, "Block-update MH algorithm. The simultaneous update of the means of each variable at every level is followed by the update of the corresponding standard deviation. The proportions of null values for the Dive Wiggliness are studied separately, by simultaneous update.Executed keeping fixed the transition probability matrices and the corresponding stationary distributions.", 
                             ifelse(state$x  == 3, "Block-update MH algorithm. The means of all the covariates are simultaneously updated level-wise, followed by the update of the corresponding standard deviations. The proportions of null values for the Dive Wiggliness are studied separately, by simultaneous update. Executed keeping fixed the transition probability matrices and the corresponding stationary distributions.", 
                                    ifelse(state$x  == 4, "Single-update MH algorithm [Alternative of the aformentioned Block Updates 1]. First, the update of the means for each variable at every level is followed by the update of the corresponding standard deviation. The proportions of null values for the Dive Wiggliness are studied at the end. Executed keeping fixed the transition probability matrices and the corresponding stationary distributions.", 
                                           ifelse(state$x  == 5, "Single-update MH algorithm [Alternative of the aformentioned Block Updates 2]. First, the means of all the covariates are updated, followed by the update of the standard deviations of all the covariates. The proportions of null values for the Dive Wiggliness are studied at the end. Executed keeping fixed the transition probability matrices and the corresponding stationary distributions.", 
                                                  ifelse(state$x  == 6, "MH algorithm involving parallel tempering with 4 different temperature chains. Executed keeping fixed the transition probability matrices and the corresponding stationary distributions.", 
                                                         ifelse(state$x  == 7, "MH algorithm involving parallel tempering with 7 different temperature chains. Executed keeping fixed the transition probability matrices and the corresponding stationary distributions.", 
                                                                "Original single-update MH algorithm for all the parameters, including the transition probability matrices and the corresponding stationary distributions.")))))))
  })
  
  output$descr <- renderText({ state$y })
  
  
  ### warning
  observe(if (input$par==6 & (input$cov==31 | input$cov==37)) {
    showNotification("Parameter not available for the selected covariate.", type="message", duration=5)
  } )
  
  ### plots
  lab <- reactiveValues()
  
  observe({
    lab$cov <- input$cov
    lab$y <- ifelse(lab$cov == 31, "Dive Duration (in seconds)",
                    ifelse(lab$cov == 37, "Maximum Depth (in metres)", 
                           ifelse(lab$cov == 43, "Dive Wiggliness (in metres)")))
  })
  
  output$plot1 <- renderPlot({
    
    if (!(input$par==6 & (input$cov==31 | input$cov==37))) {
      
      p1 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$par) + as.numeric(input$cov)),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab(lab$y) +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p2 <- ggAcf(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p3 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]), 
                   aes(x=itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab(lab$y) +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p <- grid.arrange(p1, p2, p3, ncol=3)
      
      p
      
    }
  }, height=200)
  
  output$plot2 <- renderPlot({
    
    if (!(input$par==6 & (input$cov==31 | input$cov==37))) {
      
      p3 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab(lab$y) +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p4 <- ggAcf(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p5 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]), 
                   aes(x=itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab(lab$y) +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p <- grid.arrange(p3, p4, p5, ncol=3)
      
      p
      
    }
  }, height=200)
  
  output$plot3 <- renderPlot({
    
    if (!(input$par==6 & (input$cov==31 | input$cov==37))) {
      
      p6 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+2),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$par) + as.numeric(input$cov)+2),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab(lab$y) +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p7 <- ggAcf(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+2),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p8 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+2),as.numeric(input$algo)]), 
                   aes(x=itns[,(as.numeric(input$par) + as.numeric(input$cov)+2),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab(lab$y) +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p <- grid.arrange(p6, p7, p8, ncol=3)
      
      p
      
    }
  }, height=200)
   
  observe(if (input$algo!=8) {
    showNotification("Transition proabilities (TAB 2) not available for the selected algorithm. Corresponding fixed values will be displayed.", type="warning", duration=5)
  } )
  
  output$plotTPM0 <- renderPlot({
    
    if (input$lev==1 & input$type==0){
      
      p1 <- ggplot(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("1st Entry") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p2 <- ggAcf(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p25 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)),as.numeric(input$algo)]), 
                    aes(x=itns[,(as.numeric(input$par) + as.numeric(input$cov)),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("1st Entry") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p3 <- ggplot(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+1),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$lev) + as.numeric(input$type)+1),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("2nd Entry") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p4 <- ggAcf(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+1),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p45 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]), 
                   aes(x=itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("2nd Entry") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p <- grid.arrange(p1, p2, p25, p3, p4, p45, ncol=3)
      p 
      
    }
    
  }, height=400)
  
  output$plotTPM3 <- renderPlot({
    
    if (input$lev==1 & input$type==2){
      
      p1 <- ggplot(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("1st Entry - 1st Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p2 <- ggAcf(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p25 <- ggplot(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)]), 
                   aes(x=itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("1st Entry - 1st Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p3 <- ggplot(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+1),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$lev) + as.numeric(input$type)+1),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("1st Entry - 2nd Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p4 <- ggAcf(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+1),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p45 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]), 
                   aes(x=itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("1st Entry - 2nd Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p5 <- ggplot(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+2),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$lev) + as.numeric(input$type)+2),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("2nd Entry - 1st Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p6 <- ggAcf(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+2),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p65 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+2),as.numeric(input$algo)]), 
                   aes(x=itns[,(as.numeric(input$par) + as.numeric(input$cov)+2),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("2nd Entry - 1st Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p7 <- ggplot(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+3),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$lev) + as.numeric(input$type)+3),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("2nd Entry - 2nd Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p8 <- ggAcf(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+3),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p85 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+3),as.numeric(input$algo)]), 
                   aes(x=itns[,(as.numeric(input$par) + as.numeric(input$cov)+3),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("2nd Entry - 2nd Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p <- grid.arrange(p1, p2, p25, p3, p4, p45, p5, p6, p65, p7, p8, p85, ncol=3)
      p 
      
    }
  }, height=800)
  
  output$plotTPM1 <- renderPlot({
    
    if ((input$lev==7|input$lev==10) & input$type==0){
      
      p1 <- ggplot(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("1st Entry") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p2 <- ggAcf(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p25 <- ggplot(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)]), 
                   aes(x=itns[,(as.numeric(input$lev) + as.numeric(input$type)),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("1st Entry") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p3 <- ggplot(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+1),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$lev) + as.numeric(input$type)+1),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("2nd Entry") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p4 <- ggAcf(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+1),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p45 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]), 
                   aes(x=itns[,(as.numeric(input$par) + as.numeric(input$cov)+1),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("2nd Entry") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p5 <- ggplot(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+2),as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,(as.numeric(input$lev) + as.numeric(input$type)+2),as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("3rd Entry") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p6 <- ggAcf(as.data.frame(itns[,(as.numeric(input$lev) + as.numeric(input$type)+2),as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p65 <- ggplot(as.data.frame(itns[,(as.numeric(input$par) + as.numeric(input$cov)+2),as.numeric(input$algo)]), 
                   aes(x=itns[,(as.numeric(input$par) + as.numeric(input$cov)+2),as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("3rd Entry") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p <- grid.arrange(p1, p2, p25, p3, p4, p45, p5, p6, p65, ncol=3)
      p 
      
    }
    
  },height=600)
  
  output$plotTPM2 <- renderPlot({
    
    if (input$lev==7 & input$type==2){

        p1 <- ggplot(as.data.frame(itns[,13,as.numeric(input$algo)]), aes(x=1:10000)) +
          geom_line(aes(y=as.vector(itns[,13,as.numeric(input$algo)]))) + 
          xlab("Iteration") +
          ggtitle("Trace Plot") +
          ylab("1st Entry - 1st Row") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        q1 <- ggAcf(as.data.frame(itns[,13,as.numeric(input$algo)])) +
          ggtitle("ACF Plot") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))

        p25 <- ggplot(as.data.frame(itns[,13,as.numeric(input$algo)]), 
                      aes(x=itns[,13,as.numeric(input$algo)]))+
          ggtitle("Density Plot") +
          xlab("1st Entry - 1st Row") +
          ylab("Density") +
          geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p2 <- ggplot(as.data.frame(itns[,14,as.numeric(input$algo)]), aes(x=1:10000)) +
          geom_line(aes(y=as.vector(itns[,14,as.numeric(input$algo)]))) + 
          xlab("Iteration") +
          ggtitle("Trace Plot") +
          ylab("1st Entry - 2nd Row") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        q2 <- ggAcf(as.data.frame(itns[,14,as.numeric(input$algo)])) +
          ggtitle("ACF Plot") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p45 <- ggplot(as.data.frame(itns[,14,as.numeric(input$algo)]), 
                      aes(x=itns[,14,as.numeric(input$algo)]))+
          ggtitle("Density Plot") +
          xlab("1st Entry - 2nd Row") +
          ylab("Density") +
          geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p3 <- ggplot(as.data.frame(itns[,15,as.numeric(input$algo)]), aes(x=1:10000)) +
          geom_line(aes(y=as.vector(itns[,15,as.numeric(input$algo)]))) + 
          xlab("Iteration") +
          ggtitle("Trace Plot") +
          ylab("1st Entry - 3rd Row") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        q3 <- ggAcf(as.data.frame(itns[,15,as.numeric(input$algo)])) +
          ggtitle("ACF Plot") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p65 <- ggplot(as.data.frame(itns[,15,as.numeric(input$algo)]), 
                      aes(x=itns[,15,as.numeric(input$algo)]))+
          ggtitle("Density Plot") +
          xlab("1st Entry - 3rd Row") +
          ylab("Density") +
          geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p4 <- ggplot(as.data.frame(itns[,16,as.numeric(input$algo)]), aes(x=1:10000)) +
          geom_line(aes(y=as.vector(itns[,16,as.numeric(input$algo)]))) + 
          xlab("Iteration") +
          ggtitle("Trace Plot") +
          ylab("2nd Entry - 1st Row") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        q4 <- ggAcf(as.data.frame(itns[,16,as.numeric(input$algo)])) +
          ggtitle("ACF Plot") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p85 <- ggplot(as.data.frame(itns[,16,as.numeric(input$algo)]), 
                      aes(x=itns[,1,as.numeric(input$algo)]))+
          ggtitle("Density Plot") +
          xlab("2nd Entry - 1st Row") +
          ylab("Density") +
          geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p5 <- ggplot(as.data.frame(itns[,17,as.numeric(input$algo)]), aes(x=1:10000)) +
          geom_line(aes(y=as.vector(itns[,17,as.numeric(input$algo)]))) + 
          xlab("Iteration") +
          ggtitle("Trace Plot") +
          ylab("2nd Entry - 2nd Row") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        q5 <- ggAcf(as.data.frame(itns[,17,as.numeric(input$algo)])) +
          ggtitle("ACF Plot") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p105 <- ggplot(as.data.frame(itns[,17,as.numeric(input$algo)]), 
                      aes(x=itns[,17,as.numeric(input$algo)]))+
          ggtitle("Density Plot") +
          xlab("2nd Entry - 2nd Row") +
          ylab("Density") +
          geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p6 <- ggplot(as.data.frame(itns[,18,as.numeric(input$algo)]), aes(x=1:10000)) +
          geom_line(aes(y=as.vector(itns[,18,as.numeric(input$algo)]))) + 
          xlab("Iteration") +
          ggtitle("Trace Plot") +
          ylab("2nd Entry - 3rd Row") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        q6 <- ggAcf(as.data.frame(itns[,18,as.numeric(input$algo)])) +
          ggtitle("ACF Plot") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p125 <- ggplot(as.data.frame(itns[,18,as.numeric(input$algo)]), 
                      aes(x=itns[,18,as.numeric(input$algo)]))+
          ggtitle("Density Plot") +
          xlab("2nd Entry - 3rd Row") +
          ylab("Density") +
          geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p7 <- ggplot(as.data.frame(itns[,19,as.numeric(input$algo)]), aes(x=1:10000)) +
          geom_line(aes(y=as.vector(itns[,19,as.numeric(input$algo)]))) + 
          xlab("Iteration") +
          ggtitle("Trace Plot") +
          ylab("3rd Entry - 1st Row") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        q7 <- ggAcf(as.data.frame(itns[,19,as.numeric(input$algo)])) +
          ggtitle("ACF Plot") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p145 <- ggplot(as.data.frame(itns[,19,as.numeric(input$algo)]), 
                      aes(x=itns[,19,as.numeric(input$algo)]))+
          ggtitle("Density Plot") +
          xlab("3rd Entry - 1st Row") +
          ylab("Density") +
          geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p8 <- ggplot(as.data.frame(itns[,20,as.numeric(input$algo)]), aes(x=1:10000)) +
          geom_line(aes(y=as.vector(itns[,20,as.numeric(input$algo)]))) + 
          xlab("Iteration") +
          ggtitle("Trace Plot") +
          ylab("3rd Entry - 2nd Row") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        q8 <- ggAcf(as.data.frame(itns[,20,as.numeric(input$algo)])) +
          ggtitle("ACF Plot") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p165 <- ggplot(as.data.frame(itns[,20,as.numeric(input$algo)]), 
                      aes(x=itns[,20,as.numeric(input$algo)]))+
          ggtitle("Density Plot") +
          xlab("3rd Entry - 2nd Row") +
          ylab("Density") +
          geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        p9 <- ggplot(as.data.frame(itns[,21,as.numeric(input$algo)]), aes(x=1:10000)) +
          geom_line(aes(y=as.vector(itns[,21,as.numeric(input$algo)]))) + 
          xlab("Iteration") +
          ggtitle("Trace Plot") +
          ylab("3rd Entry - 3rd Row") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
        q9 <- ggAcf(as.data.frame(itns[,21,as.numeric(input$algo)])) +
          ggtitle("ACF Plot") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
      
        p185 <- ggplot(as.data.frame(itns[,21,as.numeric(input$algo)]), 
                      aes(x=itns[,21,as.numeric(input$algo)]))+
          ggtitle("Density Plot") +
          xlab("3rd Entry - 3rd Row") +
          ylab("Density") +
          geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
          theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
                axis.title=element_text(size=8))
        
      p <- grid.arrange(p1, q1, p25, p2, q2, p45, p3,q3, p65, p4, q4, p85, p5, q5, p105,
                        p6, q6, p125, p7,q7, p145, p8, q8, p165, p9,q9, p185, ncol=3)
      p 
      
    }
    
    if (input$lev==10 & input$type==2){
      
      p1 <- ggplot(as.data.frame(itns[,22,as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,22,as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("1st Entry - 1st Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      q1 <- ggAcf(as.data.frame(itns[,22,as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p25 <- ggplot(as.data.frame(itns[,22,as.numeric(input$algo)]), 
                    aes(x=itns[,22,as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("1st Entry - 1st Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p2 <- ggplot(as.data.frame(itns[,23,as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,23,as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("1st Entry - 2nd Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      q2 <- ggAcf(as.data.frame(itns[,23,as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p45 <- ggplot(as.data.frame(itns[,23,as.numeric(input$algo)]), 
                    aes(x=itns[,23,as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("1st Entry - 2nd Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p3 <- ggplot(as.data.frame(itns[,24,as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,24,as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("1st Entry - 3rd Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      q3 <- ggAcf(as.data.frame(itns[,24,as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p65 <- ggplot(as.data.frame(itns[,24,as.numeric(input$algo)]), 
                    aes(x=itns[,24,as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("1st Entry - 3rd Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p4 <- ggplot(as.data.frame(itns[,25,as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,25,as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("2nd Entry - 1st Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      q4 <- ggAcf(as.data.frame(itns[,25,as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p85 <- ggplot(as.data.frame(itns[,25,as.numeric(input$algo)]), 
                    aes(x=itns[,1,as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("2nd Entry - 1st Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p5 <- ggplot(as.data.frame(itns[,26,as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,26,as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("2nd Entry - 2nd Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      q5 <- ggAcf(as.data.frame(itns[,26,as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p105 <- ggplot(as.data.frame(itns[,26,as.numeric(input$algo)]), 
                     aes(x=itns[,26,as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("2nd Entry - 2nd Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p6 <- ggplot(as.data.frame(itns[,27,as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,27,as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("2nd Entry - 3rd Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      q6 <- ggAcf(as.data.frame(itns[,18,as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p125 <- ggplot(as.data.frame(itns[,27,as.numeric(input$algo)]), 
                     aes(x=itns[,27,as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("2nd Entry - 3rd Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p7 <- ggplot(as.data.frame(itns[,28,as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,28,as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("3rd Entry - 1st Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      q7 <- ggAcf(as.data.frame(itns[,28,as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p145 <- ggplot(as.data.frame(itns[,28,as.numeric(input$algo)]), 
                     aes(x=itns[,28,as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("3rd Entry - 1st Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p8 <- ggplot(as.data.frame(itns[,29,as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,29,as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("3rd Entry - 2nd Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      q8 <- ggAcf(as.data.frame(itns[,29,as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p165 <- ggplot(as.data.frame(itns[,29,as.numeric(input$algo)]), 
                     aes(x=itns[,29,as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("3rd Entry - 2nd Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p9 <- ggplot(as.data.frame(itns[,30,as.numeric(input$algo)]), aes(x=1:10000)) +
        geom_line(aes(y=as.vector(itns[,30,as.numeric(input$algo)]))) + 
        xlab("Iteration") +
        ggtitle("Trace Plot") +
        ylab("3rd Entry - 3rd Row") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      q9 <- ggAcf(as.data.frame(itns[,30,as.numeric(input$algo)])) +
        ggtitle("ACF Plot") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p185 <- ggplot(as.data.frame(itns[,30,as.numeric(input$algo)]), 
                     aes(x=itns[,30,as.numeric(input$algo)]))+
        ggtitle("Density Plot") +
        xlab("3rd Entry - 3rd Row") +
        ylab("Density") +
        geom_density(fill="dodgerblue", alpha=0.5, color="cornflowerblue") +
        theme(plot.title = element_text(size=12, color="deepskyblue3", hjust = 0.5),
              axis.title=element_text(size=8))
      
      p <- grid.arrange(p1, q1, p25, p2, q2, p45, p3,q3, p65, p4, q4, p85, p5, q5, p105,
                        p6, q6, p125, p7,q7, p145, p8, q8, p165, p9,q9, p185, ncol=3)
      p 
      
    }
    
  }, height=1800)
  
}

shinyApp(ui = ui, server = server)

