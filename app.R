library(shiny)
library(ggplot2)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      sliderInput("r0", "R0", min = 0, max = 10, value = 4, step = 0.1),
      sliderInput("i0", "Initial infected", min = 0, max = 20, value = 20, step = 1),
    ),
    mainPanel(
      plotOutput("sirPlot"),
      textOutput("ssq")
    )
  )
)

server <- function(input, output) {

library(tibble)

  uncontrolled_period <- tribble(
    ~Date,        ~Cases, ~Cumulative_cases,
    "2020-03-07",     2,   50,
    "2020-03-08",     0,   50,
    "2020-03-09",     0,   50,
    "2020-03-10",     3,   53,
    "2020-03-11",     6,   59,
    "2020-03-12",    11,   70,
    "2020-03-13",     5,   75,
    "2020-03-14",     7,   82,
    "2020-03-15",    32,  114,
    "2020-03-16",    33,  147,
    "2020-03-17",    30,  177,
    "2020-03-18",    35,  212,
    "2020-03-19",    60,  272,
    "2020-03-20",    50,  322,
    "2020-03-21",    89,  411,
    "2020-03-22",   188,  599,
    "2020-03-23",   122,  721,
    "2020-03-24",   106,  827,
    "2020-03-25",   107,  934,
    "2020-03-26",   111, 1045
  )
  uncontrolled_period$Date <- as.Date(uncontrolled_period$Date)

  incidence_frame <- reactive({
    N <- 6.6e7
    gamma <- 1/6
    beta <- input$r0 * gamma
    omega <- 1/5
    I0 <- input$i0
    S0 <- N - I0
    R0 <- 0
    days <- nrow(uncontrolled_period)
    dt <- 0.1

    n <- ceiling(days / dt)
    time <- seq(0, by = dt, length.out = n)

    S <- numeric(n)
    E <- numeric(n)
    I <- numeric(n)
    R <- numeric(n)

    S[1] <- S0
    E[1] <- 0
    I[1] <- I0
    R[1] <- R0

    for (t in 2:n) {
      dS <- -beta * S[t-1] * I[t-1] / N
      dE <- beta * S[t-1] * I[t-1] / N - omega * E[t-1]
      dI <- omega * E[t-1] - gamma * I[t-1]
      dR <- gamma * I[t-1]

      S[t] <- S[t-1] + dt * dS
      E[t] <- E[t-1] + dt * dE
      I[t] <- I[t-1] + dt * dI
      R[t] <- R[t-1] + dt * dR
    }

    incidence <- E * omega

    incidence_frame <- data.frame(time = uncontrolled_period$Date[1] + time, incidence = incidence, t = time)

    incidence_frame
  })

  incidence_frame_whole <- reactive({
    incidence_frame_whole <- incidence_frame()
    indexes <- which(incidence_frame_whole$t == floor(incidence_frame_whole$t))
    incidence_frame_whole <- incidence_frame_whole[indexes,]
    incidence_frame_whole$actual <- uncontrolled_period$Cases
    incidence_frame_whole$difference <- incidence_frame_whole$actual - incidence_frame_whole$incidence

    incidence_frame_whole
  })

  ssq <- reactive({
    sum(incidence_frame_whole()$difference ^ 2)
  })
  
  output$sirPlot <- renderPlot({
    
    # matplot(time, cbind(S, I, E, R),
    #         type = "l", lty = 1, col = c("blue", "red", "green"),
    #         xlab = "Time (days)", ylab = "Proportion",
    #         main = "SIR Model (Euler's Method)")
    # legend("right", legend = c("Susceptible", "Infected", "Recovered"),
    #        col = c("blue", "red", "green", "orange"), lty = 1)

    ggplot(incidence_frame(), aes(x=time, y=incidence)) +
    geom_col(aes(x=Date, y=Cases), data = uncontrolled_period, fill = "red") +
    geom_line() +
    geom_point(data = incidence_frame_whole(), aes(x=time, y=actual)) +
    geom_point(data = incidence_frame_whole()) +
    geom_segment(aes(x = time, xend = time, y = actual, yend = incidence), data = incidence_frame_whole()) +
    labs(x = "Date", y = "Incidence")
  })

  output$ssq <- renderText({
    paste0("Sum of Squared Error (SSQ): ", prettyNum(round(ssq(), 2), big.mark=","))
  })
}

shinyApp(ui, server)