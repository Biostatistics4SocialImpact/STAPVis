#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(purrr)
# Define UI for application that draws a histogram
ui <- navbarPage("STAPVis",
                 tabPanel("Step 1",
                          sidebarLayout(
                              sidebarPanel(
                                  sliderInput("d_max",
                                              "Maximum Distance",
                                              min = 1,
                                              max = 20,
                                              value = 5),
                                  sliderInput("theta",
                                              "Spatial parameter",
                                              min = 0.001,
                                              max = 20,
                                              value = 1),
                                  sliderInput("theta_2",
                                              "Spatial shape parameter (Weibull only)",
                                              min = 0.001,
                                              max = 20,
                                              value = 1),
                                  radioButtons("K","Exposure Function",
                                               c("Complementary Error Function" = "erfc",
                                                 "Exponential" = "exp",
                                                 "Weibull" = "wei")),
                                  sliderInput("p",
                                              "Exposure Precision",
                                              min = 0.001,
                                              max = 1,
                                              value = 0.05),
                              ),
                              mainPanel(
                                  plotOutput("expPlot")
                              )
                          )
                 ),
                 tabPanel("Step 2",
                          sidebarLayout(
                              sidebarPanel(
                                  sliderInput("d_max2",
                                              "Maximum Distance",
                                              min = 1,
                                              max = 20,
                                              value = 5),
                                  radioButtons("prior_dist_1","Spatial Parameter Prior Distribution",
                                               c("|Normal|" = "rnorm",
                                                 "Log Normal" = "rlnorm")),
                                  sliderInput("mean_1",
                                              "Mean",
                                              min = -10,
                                              max = 10,value = 1),
                                  sliderInput("scale_1",
                                              "Scale",
                                              min = 0.01,
                                              max = 5,
                                              value = 1),
                                  radioButtons("prior_dist_2","Spatial Shape Prior Distribution (Weibull only)",
                                               c("|Normal|" = "rnorm",
                                                 "Log Normal" = "rlnorm")),
                                  sliderInput("mean_2",
                                              "Mean",
                                              min = -10,
                                              max = 10,
                                              value = 1),
                                  sliderInput("scale_2",
                                              "Scale",
                                              min = 0.01,
                                              max = 5, value = 1),
                                  radioButtons("K_2","Exposure Function",
                                               c("Complementary Error Function" = "erfc",
                                                 "Exponential" = "exp",
                                                 "Weibull" = "wei")),
                                  sliderInput("p_2",
                                              "Exposure Precision",
                                              min = 0.001,
                                              max = 1,
                                              value = 0.05),
                              ),
                              mainPanel(
                                  plotOutput("expPlot2"),
                                  tableOutput("sumTable")
                              )
                          ), class = "rightAlign")
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$expPlot <- renderPlot({
        expf <- switch(input$K,
                       erfc = function(x,y) pracma::erfc(x/y),
                       exp = function(x,y) pexp(x,rate = 1/y, lower.tail = FALSE),
                       wei = function(x,y) pweibull(x,shape = input$theta_2, scale = y, lower.tail = FALSE))
        
        df <- tibble(Distance = seq(from = 0.01, to = input$d_max,by=0.01),
                     Exposure = expf(Distance,input$theta))
        d_star <- uniroot(function(x) expf(x,input$theta) - input$p,interval = c(0,input$d_max))$root
        
        p <- df %>% ggplot(aes(x=Distance,y=Exposure)) + geom_line() + ggthemes::theme_hc()
        if(!is.null(d_star))
            p <- p + geom_vline(aes(xintercept = d_star),
                                linetype = 2,
                                color='red') + annotate("text",label = paste0("Spatial Scale: ",round(d_star,2)), x = d_star-0.35, y = 0.77,
                                                        size = 5)
        return(p)
    })
    
    output$expPlot2 <- renderPlot({
        expf <- switch(input$K_2,
                       erfc = function(x,y,z) pracma::erfc(x/y),
                       exp = function(x,y,z) pexp(x,rate = 1/y, lower.tail = FALSE),
                       wei = function(x,y,z) pweibull(x,shape = z, scale = y, lower.tail = FALSE))
        
        dist1 <- switch(input$prior_dist_1,
                        rnorm = function(x,y,z) abs(rnorm(x,y,z)),
                        rlnorm = rlnorm)
        dist2 <- switch(input$prior_dist_2,
                        rnorm = function(x,y,z) abs(rnorm(x,y,z)),
                        rlnorm = rlnorm)
        num_samples <- 1E3
        theta_samples <- dist1(num_samples,input$mean_1,input$scale_1)
        shape_samples <- dist2(num_samples,input$mean_2,input$scale_2)
        d <- seq(from = 0.01, to = input$d_max2,by=0.01)
        
        df <- purrr::map_dfr(d,function(x){tibble(Distance = x,
                                                  Spatial_par = theta_samples,
                                                  Spatial_shape = shape_samples,
                                                  Exposure = expf(x,theta_samples,shape_samples)
                                                  )})
        df <- df %>% group_by(Distance) %>% 
            summarise(lower = quantile(Exposure,0.025),
                      median = median(Exposure),
                      upper = quantile(Exposure,0.975))
        
        
        p <- df %>% ggplot(aes(x=Distance,y=median)) + geom_line() + ggthemes::theme_hc() + 
            geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3) + ylab("Exposure") + 
            labs(title = "Prior Spatial Exposure Function Visualization",
                 subtitle = "Shaded Area is 95% Credible Interval")
            
        
        return(p)
        
    })
    output$sumTable <- renderTable({
        expf <- switch(input$K_2,
                       erfc = function(x,y,z) pracma::erfc(x/y),
                       exp = function(x,y,z) pexp(x,rate = 1/y, lower.tail = FALSE),
                       wei = function(x,y,z) pweibull(x,shape = z, scale = y, lower.tail = FALSE))
        
        dist1 <- switch(input$prior_dist_1,
                        rnorm = function(x,y,z) abs(rnorm(x,y,z)),
                        rlnorm = rlnorm)
        dist2 <- switch(input$prior_dist_2,
                        rnorm = function(x,y,z) abs(rnorm(x,y,z)),
                        rlnorm = rlnorm)
        num_samples <- 1E3
        theta_samples <- dist1(num_samples,input$mean_1,input$scale_1)
        shape_samples <- dist2(num_samples,input$mean_2,input$scale_2)
        d <- seq(from = 0.01, to = input$d_max2,by=0.01)
        
        df <- purrr::map_dfr(d,function(x){tibble(Distance = x,
                                                  Spatial_par = theta_samples,
                                                  Spatial_shape = shape_samples,
                                                  Exposure = expf(x,theta_samples,shape_samples)
        )})
        df <- df %>% group_by(Distance) %>% 
            summarise(lower = quantile(Exposure,0.025),
                      median = median(Exposure),
                      upper = quantile(Exposure,0.975))
        df %>% ungroup() %>% mutate(ix = 1:n(),
                      `2.5%` = ((lower - input$p_2) <= 0.01),
                      `50%` = ((median - input$p_2) <= 0.01),
                      `97.5%` = ((upper - input$p_2) <= 0.01)) %>% 
            tidyr::gather(`2.5%`,`50%`,`97.5%`,key="Quantile",value="Indicator") %>% 
            filter(Indicator==1) %>% 
            group_by(Quantile) %>% filter(Distance==min(Distance)) %>% 
            select(Distance,Quantile) %>% tidyr::spread(Quantile,Distance) -> df
        df <- df %>% mutate(Quantile= "Spatial Scale") %>% select(Quantile,`2.5%`,`50%`,`97.5%`)
        return(df)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
