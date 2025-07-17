library(shiny)
library(ggplot2)
if (FALSE) {
  library(munsell)
}

ui <- fluidPage(
  titlePanel(
    tags$div(
  style = "
    display: flex;
    justify-content: center;
    align-items: center;
    gap: 15px;
    margin-bottom: 20px;
  ",
  tags$img(src = "logo.png", height = "60px"),
  tags$h2("Fitting an Epidemic Model with Least Squares", style = "margin: 0;")
  ),
  "Fitting an Epidemic Model with Least Squares"),


  tags$style(HTML("
  
    @media (min-width: 992px) {
        body::before {
            content: '';
            position: fixed;
            top: -100vh;
            left: 20;
            width: 30px;
            height: 200vh;
            background: #DB4436;
            transform: rotate(45deg);
            z-index: -1;
            pointer-events: none;

            /* Simulate 3 lines using box-shadow */
            box-shadow:
            35px 0 #4BA5AB,
            70px 0 #74C2CE;
        }
    }"
  )),

    tags$div(
    style = "
      background-color: white;
      padding: 30px;
      border-radius: 8px;
      box-shadow: 0 0 10px rgba(0,0,0,0.05);
    ",
  
  
  tags$p("This app simulates a basic SEIR model and compares the model's incidence to observed COVID-19 case data from March 2020."),
  tags$p("Use the sliders to adjust the basic reproduction number (R₀) and the initial number of infected individuals. The plot shows observed vs. simulated incidence, and the Sum of Squared Errors (SSQ) quantifies the fit."),
  tags$p("Try and find the best solution by moving around the sliders."),

  br(),
  
  fluidRow(
  tags$style(HTML("
    .flex-container {
      display: flex;
      align-items: center;
    }
    .slider-column {
      flex: 0 0 300px;
      padding-right: 30px;
      padding-left: 30px;
    }
    .plot-column {
      flex: 1;
    }
  ")),
  
  div(class = "flex-container",
      div(class = "slider-column",
          sliderInput("r0", "Basic Reproduction Number (R₀)", min = 0, max = 8, value = 4, step = 0.1),
          sliderInput("i0", "Initial Infected Individuals", min = 0, max = 40, value = 20, step = 1)
      ),
      div(class = "plot-column",
          h4(textOutput("ssq"), style = "color:#2C3E50; font-weight:bold; margin-bottom: 20px;"),
          plotOutput("sirPlot", height = "500px")
      )
  )
)
    )
)

server <- function(input, output) {
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
    gamma <- 1 / 6
    beta <- input$r0 * gamma
    omega <- 1 / 5
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
      dS <- -beta * S[t - 1] * I[t - 1] / N
      dE <- beta * S[t - 1] * I[t - 1] / N - omega * E[t - 1]
      dI <- omega * E[t - 1] - gamma * I[t - 1]
      dR <- gamma * I[t - 1]

      S[t] <- S[t - 1] + dt * dS
      E[t] <- E[t - 1] + dt * dE
      I[t] <- I[t - 1] + dt * dI
      R[t] <- R[t - 1] + dt * dR
    }

    incidence <- E * omega
    data.frame(time = uncontrolled_period$Date[1] + time, incidence = incidence, t = time)
  })

  incidence_frame_whole <- reactive({
    df <- incidence_frame()
    indexes <- which(df$t == floor(df$t))
    df <- df[indexes, ]
    df$actual <- uncontrolled_period$Cases
    df$difference <- df$actual - df$incidence
    df
  })

  ssq <- reactive({
    sum(incidence_frame_whole()$difference ^ 2)
  })

  output$sirPlot <- renderPlot({
    ggplot() +
      geom_col(data = uncontrolled_period, aes(x = Date, y = Cases), fill = "#E74C3C", alpha = 0.6, width = 0.8) +
      geom_line(data = incidence_frame(), aes(x = time, y = incidence), color = "#2C3E50", size = 1.2) +
      geom_point(data = incidence_frame_whole(), aes(x = time, y = actual), color = "#E74C3C", size = 2) +
      geom_point(data = incidence_frame_whole(), aes(x = time, y = incidence), color = "#2C3E50", size = 2) +
      geom_segment(data = incidence_frame_whole(),
                   aes(x = time, xend = time, y = actual, yend = incidence),
                   color = "gray40", linetype = "dashed") +
      labs(x = "Date", y = "Daily Incidence") +
      theme_minimal(base_size = 14)
  })

  output$ssq <- renderText({
    paste0("Sum of Squared Error (SSQ): ", prettyNum(round(ssq(), 2), big.mark = ","))
  })
}

shinyApp(ui, server)