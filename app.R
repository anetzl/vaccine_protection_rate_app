library(shiny)
library(DT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readxl)
library(rstatix)

# calculate total protection rates
total_pr <- function(df) {
  
  # calculate protection rates per strain and vaccine
  df %>%
    separate(expected_circulation, into = c("circulation_lower", "circulation_upper"), sep = "-", convert = TRUE, fill= "right") %>%
    mutate(circulation_upper = ifelse(is.na(circulation_upper), circulation_lower, circulation_upper)) %>%
    pivot_longer(
      cols = starts_with("vaccine"),
      names_to = "vaccine",
      values_to = "protection_rate"
    ) %>%
    separate(protection_rate, into = c("pr_lower", "pr_upper"), sep = "-", convert = TRUE, fill= "right") %>%
    rowwise() %>%
    mutate(pr_upper = ifelse(is.na(pr_upper), pr_lower, pr_upper),
           circulation_mean = mean(c(circulation_upper, circulation_lower)),
           pr_mean = mean(c(pr_lower, pr_upper))) %>%
    ungroup() %>%
    pivot_longer(
      cols = starts_with("circulation"),
      names_to = "circulation_estimate",
      values_to = "circulation"
    ) %>%
    pivot_longer(
      cols = starts_with("pr"),
      names_to = "pr_estimate",
      values_to = "protection_rate"
    ) %>%
    mutate(strain_circ_pr = protection_rate * circulation) %>%
    group_by(vaccine, circulation_estimate, pr_estimate) %>%
    mutate(vaccine_circ_pr = sum(strain_circ_pr),
           vaccine_circ_pr = ifelse(vaccine_circ_pr > 1, 1, vaccine_circ_pr),
           strain_circ_pr = ifelse(strain_circ_pr > 1, 1, strain_circ_pr),
           total_circulation = sum(circulation),
           vaccine = gsub("_pr", "", vaccine )) %>%
    ungroup() %>%
    mutate(max_circulation = max(total_circulation)) -> df
  
  if(unique(df$max_circulation) > 1){
    warning("The estimated strain circulation values you provided result in some cases where the total strain circulation exceeds 1, and 
            hence total protection rates (PR) exceed 1. Please provide estimated circulation 
            rates where the upper limits of all strains combined do not exceed 1.")
  }
  
  
  return(df)
}

pr_protection_plot <- function(df){
  
  df_strain <- df %>%
    select(strain_name, vaccine, strain_circ_pr, circulation_estimate, pr_estimate) %>%
    group_by(strain_name, vaccine) %>%
    mutate(min_pr = min(strain_circ_pr),
           max_pr = max(strain_circ_pr),
           mean_pr = strain_circ_pr[circulation_estimate == "circulation_mean" & pr_estimate == "pr_mean"]) %>%
    ungroup() %>%
    select(strain_name, vaccine, min_pr, max_pr, mean_pr) %>%
    mutate(Strain = strain_name) %>%
    unique()
  
  
  df_vacc <- df %>%
    select(vaccine, vaccine_circ_pr, circulation_estimate, pr_estimate) %>%
    unique() %>%
    group_by(vaccine) %>%
    mutate(min_pr = min(vaccine_circ_pr),
           max_pr = max(vaccine_circ_pr),
           mean_pr = vaccine_circ_pr[circulation_estimate == "circulation_mean" & pr_estimate == "pr_mean"],
           strain_name = "Total") %>%
    ungroup() %>%
    select(strain_name, vaccine, min_pr, max_pr, mean_pr) %>%
    unique()
  
  df_vacc %>%
    ggplot(aes(x = vaccine)) +
    geom_pointrange(aes(y = mean_pr, ymin = min_pr, ymax = max_pr)) +
    geom_pointrange(data = df_strain, aes(y = mean_pr, ymin = min_pr, ymax = max_pr, x = vaccine, color = Strain, group = Strain), position = position_dodge(width = 0.2), alpha = 0.6) + 
    labs(
      x = "Vaccine",
      y = "Protection Rate"
    ) +
    ylim(c(0,1)) + 
    theme_bw() -> p
  
  return(p)
}


# functions for sampling from distribution
## simulate from distribution
simulate_posterior <- function(lower_circulation, upper_circulation,
                               mean_pr = 0.85, sd_pr = 0.10,
                               n_sim = 10000) {
  # Prior distributions
  circulation_samples <- runif(n_sim, min = lower_circulation, max = upper_circulation)
  protection_samples <- rnorm(n_sim, mean = mean_pr, sd = sd_pr)
  
  # as protection rates can only be between 0 and 1, censor them like this
  protection_samples[protection_samples > 1] <- 1
  protection_samples[protection_samples < 0] <- 0
  
  # Simulated posterior: product of sampled circulation and protection
  total_protection_samples <- circulation_samples * protection_samples
  
  return(data.frame(
    circulation_samples = circulation_samples,
    protection_samples = protection_samples,
    protection_rate = total_protection_samples,
    sample_nr = 1:n_sim
  ))
}

format_dist_input <- function(df_dist){
  
  df_dist %>%
    separate(expected_circulation, into = c("circulation_lower", "circulation_upper"), sep = "-", convert = TRUE, fill= "right") %>%
    mutate(circulation_upper = ifelse(is.na(circulation_upper), circulation_lower, circulation_upper)) %>%
    pivot_longer(cols = starts_with("vaccine"), names_to = "vaccine", values_to = "value") %>%
    separate(value, into = c("mean", "sd"), sep = ";", convert = TRUE) %>%
    separate(vaccine, into = c("vaccine", "extra"), sep = "_", convert = TRUE) %>%
    select(!extra) -> df_dist_long
  
  return(df_dist_long)
  
}

sample_vaccine_posteriors <- function(df, n_sim = 1000) {
  
  df_dist_long <- format_dist_input(df)
  
  n_vacc <- unique(df_dist_long$vaccine)
  n_strain <- unique(df_dist_long$strain_name)
  full_pr <- list()
  for(v in n_vacc){
    
    for(s in n_strain){
      
      temp_df <- df_dist_long %>%
        filter(strain_name == s) %>%
        filter(vaccine == v) %>%
        unique()
      
      full_pr[[paste0(v, "_", s)]] <- simulate_posterior(temp_df$circulation_lower, 
                                                         temp_df$circulation_upper, 
                                                         temp_df$mean,
                                                         temp_df$sd,
                                                         n_sim = n_sim) %>%
        mutate(strain = s,
               vaccine = v)
      
    }
    
  }
  
  full_pr <- do.call(rbind, full_pr) 
  
  full_pr %>%
    group_by(vaccine, sample_nr) %>%
    mutate(total_circulation = sum(circulation_samples)) %>%
    ungroup() %>%
    filter(total_circulation > 1) %>%
    pull(sample_nr) %>%
    unique() -> repeat_samples
  
  full_pr %>%
    group_by(vaccine, sample_nr) %>%
    mutate(total_pr = sum(protection_rate),
           total_circulation = sum(circulation_samples)) %>%
    ungroup() %>%
    filter(!sample_nr %in% repeat_samples) -> full_pr
  
  
  return(list("full_pr" = full_pr, "repeat_samples" = repeat_samples))
  
}

write_sample_discard <- function(df, n_sim){
  
  repeat_samples <- sample_vaccine_posteriors(df, n_sim)$repeat_samples
  n_repeat_samples <- length(repeat_samples)
  
  if(n_repeat_samples > 0){
    output_text <- paste0('The estimated strain circulation values you provided result in ', n_repeat_samples, ' cases where the total strain circulation exceeds 1. 
            These samples will be discarded, and hence the effective sample size for calculating posterior protection rates is ', n_sim - n_repeat_samples, ', not your input ', n_sim,'. 
            Please provide estimated circulation 
            rates where the upper limits of all strains combined do not exceed one to avoid this warning.')
  } else {
    output_text <- ""
  }
  
  return(output_text)
  
}

# plot posterior distributions
plot_posterior_pr <- function(df, n_sim = 1000){
  
  full_pr <- sample_vaccine_posteriors(df, n_sim)$full_pr
  
  # calculate mean and sd pr
  full_pr %>%
    group_by(strain, vaccine) %>%
    reframe(mean = Rmisc::CI(protection_rate)["mean"],
            lower = Rmisc::CI(protection_rate)["lower"],
            upper = Rmisc::CI(protection_rate)["upper"]) -> means_strain
  
  full_pr %>%
    group_by(vaccine) %>%
    reframe(mean = Rmisc::CI(total_pr)["mean"],
            lower = Rmisc::CI(total_pr)["lower"],
            upper = Rmisc::CI(total_pr)["upper"]) -> means_vacc
  
  full_pr %>%
    select(total_pr, vaccine, sample_nr) %>%
    unique() %>%
    ggplot(aes(x = total_pr, fill = vaccine)) + 
    geom_histogram(alpha=0.6, position = 'identity') + 
    geom_vline(data = means_vacc, aes(xintercept = mean, color = vaccine)) +
    geom_vline(data = means_vacc, aes(xintercept = lower, color = vaccine), linetype = "dashed") +
    geom_vline(data = means_vacc, aes(xintercept = upper, color = vaccine), linetype = "dashed") +
    theme_bw() + 
    xlim(c(0, 1)) + 
    xlab(paste0("Total PR (", "\u2211", "Strain PR)")) -> p_vacc
  
  
  full_pr %>%
    ggplot(aes(x = protection_rate, fill = vaccine)) + 
    geom_histogram(alpha=0.6, position = 'identity') + 
    geom_vline(data = means_strain, aes(xintercept = mean, color = vaccine)) +
    geom_vline(data = means_strain, aes(xintercept = lower, color = vaccine), linetype = "dashed") +
    geom_vline(data = means_strain, aes(xintercept = upper, color = vaccine), linetype = "dashed") +
    facet_wrap(~strain) + 
    theme_bw() + 
    xlab("Strain PR (circulation x Vaccine PR)") + 
    xlim(c(0, 1)) + 
    theme(strip.background.x = element_blank()) -> p_strain
  
  p_vacc / p_strain + plot_layout(widths = c(1, 1/length(unique(full_pr$strain)))) -> p_comb
  
  return(p_comb)
  
}

calculate_dist_t_test <- function(df, n_sim){
  
  round_three_digits <- function(x){
    round(x, 3)
  }
  
  full_pr <- sample_vaccine_posteriors(df, n_sim)$full_pr
  
  full_pr %>%
    select(sample_nr, total_pr, vaccine) %>%
    unique() %>%
    mutate(strain = "Total") %>%
    group_by(strain) %>%
    t_test(total_pr~vaccine, detailed = TRUE) %>%
    rbind(.,
          full_pr %>%
            group_by(strain) %>%
            t_test(protection_rate~vaccine, detailed = TRUE)) %>%
    add_significance() %>%
    mutate(across(where(is.double), round_three_digits)) %>%
    mutate("Mean PR group 1 - group 2" = estimate,
           "Mean PR group 1" = estimate1,
           "Mean PR group 2" = estimate2) %>%
    select(!starts_with("estimate")) %>%
    select(!`.y.`) -> tab_out
  
  return(tab_out)
}


ui <- fluidPage(
  titlePanel("Interactive Total Protection Rate calculator for vaccine strains"),
  
  HTML(paste('<font size="2">This app calculates total Protection Rates (PR) elicited by various vaccine strains <i>V</i> from supplied estimates 
             using the following formula: <br><br>')),
             
  withMathJax(HTML("\\( PR(v) = \\sum_{s}^{S} PR_{s}(v) \\cdot PE_s \\)<br><br><br>")),
  withMathJax(HTML(paste("\\( PR_{s}(v) = \\)", "Protection rate by vaccine strain <i>v</i> against virus strain <i>s</i> <br><br>",
                         "\\( PE_s = \\)", "Probability of exposure to strain <i>s</i> (~ expected circulation probability) <br><br><br>")
                         )),      
  
  
  HTML(paste('Two evaluation modes are possible: A frequentist approach, where the user supplies discrete values or ranges for 
             which the PRs should be calcualated, and a Baeysian approach, where strain circulation values are sampled from 
             a uniform disitribution and strain protection rates from a normal distribution, the distribution parameters 
             are user specified input. The Bayesian approach results in a distribution of total protection rates per vaccine.<br><br>')),
  
  sidebarLayout(
    sidebarPanel(
      HTML(paste('<B><font size="2">Frequentist data file (.csv, .xlsx):</B><br>')),
      HTML(paste('<n><font size="1">
                  File of format S+1 rows, V+2 columns. <br>
                  File requires these columns in this order: strain_name, expected_circulation, vaccine1_pr, ..., vaccineV_pr <br>
                  strain_name: Any format, S rows for S strains to be included in the protection rate calculation <br>
                  expected_circulation: value between 0-1 per strain s, can be specified as "lower-upper" range (e.g: 0.2-0.4). <i>A warning will be issued
                   if the sum of expected_circulation of all strains is larger than 1</i> <br>
                  vaccineV_pr: estimated protection rate from vaccine V for strain s, between 0-1, can be specified as "lower-upper" range (e.g: 0.6-0.8)
                  ')),
      fileInput("file_upload", "",
                accept = c(".csv", ".xlsx")),
      
      HTML(paste('<B><font size="2">Bayesian data file (.csv, .xlsx):</B><br>')),
      HTML(paste('<n><font size="1">
                  File of format S+1 rows, V+2 columns. <br>
                  File requires these columns in this order: strain_name, expected_circulation, vaccine1_pr (mean;sd), ..., vaccineV_pr (mean;sd)<br>
                  strain_name: Any format, S rows for S strains to be included in the protection rate calculation <br>
                  expected_circulation: value between 0-1 per strain s, can be specified as "lower-upper" range (e.g: 0.2-0.4). Samples will be 
                  drawn randomly from a uniform distribution between "lower" and "upper" range. <i>A warning will be issued
                   if the sum of expected_circulation of all strains is larger than 1. Samples with total circulation > 1 will be discarded.</i> <br>
                  vaccineV_pr: estimated mean protection rate from vaccine V for strain s with standard deviation, values between 0 and 1 are required. Samples will be drawn from 
                  a normal distribution with the specified mean and sd (e.g: 0.8;0.1). Any samples > 1 will be set to 1, any samples < 0 will be set to 0.
                  ')),
      fileInput("dist_file_upload", "",
                accept = c(".csv", ".xlsx")),
      
      HTML(paste('<B><font size="2">Number of samples for the Bayesian summary:</B><br>')),
      HTML(paste('<n><font size="1">
                  Integer between 10 and 10000. Samples with total strain circulation > 1 will be removed, leading to an effective sample size lower than the 
                  input value in such cases.
                  <B><font size="2">
                  ')),
      numericInput(
        inputId = "sim_n",         # ID to reference in server
        label = "",  # Label in the UI
        value = 1000,              # Default value
        min = 10,                 # Minimum allowed
        max = 10000,              # Maximum allowed
        step = 1                 # Step size
      ),
      numericInput(
        inputId = "sim_seed",         # ID to reference in server
        label = "Random seed for reproducible Bayesian outcome",  # Label in the UI
        value = 666,              # Default value
        min = 1,                 # Minimum allowed
        max = 10000,              # Maximum allowed
        step = 1                 # Step size
      )
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Frequentist Summary",
                 DTOutput("csv_table"),
                 plotOutput("protection_plot"),
                 textOutput("error_msg")),
        tabPanel("Bayesian Summary",
                 DTOutput("csv_dist_table"),
                 textOutput("dist_sample_discard_text"),
                 plotOutput("protection_dist_plot"),
                 DTOutput("stat_dist_table"),
                 textOutput("error_msg"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  csv_data <- reactive({
    req(input$file_upload)  # Ensure a file is selected
    tryCatch({
      if(grepl(".csv", input$file_upload$datapath)){
        read.csv(input$file_upload$datapath, stringsAsFactors = FALSE)
      } else {
        read_excel(input$file_upload$datapath)
      }
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      NULL
    })
  })
  
  output$csv_table <- renderDT({
    dat <- csv_data()
    req(dat)  # Only render if data is available
    datatable(dat, options = list(pageLength = 10))
  })
  
  
  output$protection_plot <- renderPlot({
    pr_test <- total_pr(csv_data())
    validate(
      need((unique(pr_test$max_circulation) <= 1), paste0("The maximum circulation probability is ", unique(pr_test$max_circulation), ". 
      Please lower the expected circulation values such that the sum of all strains is <=1.")),
      need((max(pr_test$protection_rate) <= 1), paste0("The maximum protection rate you provided is ", max(pr_test$protection_rate), ". 
      Please lower the maximum protection rate such that it is <=1."))
    )
    pr_protection_plot(pr_test)
  })
  
  # here for distribution data
  csv_dist_data <- reactive({
    req(input$dist_file_upload)  # Ensure a file is selected
    tryCatch({
      if(grepl(".csv", input$dist_file_upload$datapath)){
        read.csv(input$dist_file_upload$datapath, stringsAsFactors = FALSE)
      } else {
        read_excel(input$dist_file_upload$datapath)
      }
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      NULL
    })
  })
  
  observe({set.seed(input$sim_seed)})
  
  output$csv_dist_table <- renderDT({
    dat <- csv_dist_data()
    req(dat)  # Only render if data is available
    datatable(dat, options = list(pageLength = 10))
  })
  
  output$protection_dist_plot <- renderPlot({
    dat <- csv_dist_data()
    req(dat)  # Only render if data is available
    plot_posterior_pr(dat, input$sim_n)
  })
  
  output$dist_sample_discard_text <- renderText({
    dat <- csv_dist_data()
    req(dat)  # Only render if data is available
    write_sample_discard(dat, input$sim_n)
  })
    
  output$stat_dist_table <- renderDT({
    tab_out <- calculate_dist_t_test(csv_dist_data(), input$sim_n)
    req(tab_out)  # Only render if data is available
    datatable(tab_out, options = list(pageLength = 10))
    
  })
  
}

shinyApp(ui, server)


# Feed in serological data, antibody titers
# look more and more at serology, comparative tests
# convert titer into protection