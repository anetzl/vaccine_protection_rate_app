# This app takes a google sheet as input
library(shiny)
library(DT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readxl)
library(rstatix)
library(googlesheets4)
library(RColorBrewer)
library(scales)
library(tibble)
library(colorspace)

theme_set(theme_bw() + 
            theme(strip.background = element_blank(),
                #  axis.text.x = element_text(angle = 45, hjust = 1)),
                  panel.spacing.x = unit(1.3, "lines")))


# sampling posteriors
sample_posterior <- function(mean_norm = 0.1, sd_norm = 0.1, n_sim = 100){
  res <- rnorm(n_sim, mean = mean_norm, sd = sd_norm)
  res[is.na(res)] <- 0
  res[res < 0] <- 0
  res
}

get_color_pallette <- function(variable, pallette = NULL){
  if(is.null(pallette)){
    target_cols <- hue_pal()(length(variable))
  } else {
   
    base_pallette <- brewer.pal(n = 8, name = pallette)
    target_cols <- colorRampPalette(base_pallette)(12)
  }
  names(target_cols) <- variable
  
  return(target_cols)
}

sample_normalised_circulation_scenario <- function(data, target_scenario, sd_norm, n_sim){
  
  lapply(data[[target_scenario]], function(x){
    sample_posterior(x, sd_norm, n_sim)
  }) -> test_circ
  names(test_circ) <- data$variant
  do.call(rbind, test_circ) -> test_circ
  test_circ / colSums(test_circ)[col(test_circ)] -> circ_norm
  as.data.frame(circ_norm) %>%
    rownames_to_column(var = "variant") %>%
    pivot_longer(cols = colnames(.)[2:ncol(.)], names_to = "sample", values_to = "circulation_rate") %>%
    mutate(scenario = target_scenario) -> circ_norm
  
  return(circ_norm)
}


sample_pr_scenario <- function(data, target_scenario, sd_norm, n_sim){
  
  lapply(data[[target_scenario]], function(x){
    res <- sample_posterior(x, sd_norm, n_sim)
    res[res>1] <- 1
    res
  }) -> test_circ
  names(test_circ) <- data$variant
  do.call(rbind, test_circ) -> test_circ
  as.data.frame(test_circ) %>%
    rownames_to_column(var = "variant") %>%
    pivot_longer(cols = colnames(.)[2:ncol(.)], names_to = "sample", values_to = "protection_rate") %>%
    mutate(vaccine = target_scenario) -> circ_norm
  
  return(circ_norm)
}

simulate_total_pr_posterior <- function(data_circulation, data_pr, sd_circ, sd_pr, sd_data, n_sim){
  
  # rename reactive containers
  data_circulation <- data_circulation %>%
    rename(., variant = colnames(.)[1]) %>%
    filter(!is.na(variant))
  
  data_pr <- data_pr %>%
    rename(., variant = colnames(.)[1]) %>%
    filter(!is.na(variant)) %>%
    filter(variant %in% data_circulation$variant)
  
  data_pr <- data_pr %>%
    filter(variant %in% data_circulation$variant)
  
  all_norm_circulation <- lapply(colnames(data_circulation)[2:ncol(data_circulation)], function(x){
    sample_normalised_circulation_scenario(data_circulation, x, sd_circ, n_sim)
  })
  
  all_pr <- lapply(colnames(data_pr)[2:ncol(data_pr)], function(x){
    sample_pr_scenario(data_pr, x, sd_pr, n_sim)
  })
  
  lapply(all_norm_circulation, function(x){
    temp <- lapply(all_pr, function(pr){
      x %>%
        left_join(., pr, by = c("variant", "sample")) %>%
        mutate(strain_pr = circulation_rate * protection_rate,
               full_scenario = paste0(scenario,", v: ", vaccine))
    })
    do.call(rbind, temp)
  }) -> combined_pr
  
  combined_pr <- do.call(rbind, combined_pr)
  
  total_pr <- combined_pr %>%
    group_by(full_scenario, sample) %>%
    mutate(total_pr = sum(strain_pr)) %>%
    ungroup()
  
  return(total_pr)
}


calc_mean_vacc_pr <- function(total_pr){
  # summarise prs
  total_pr %>%
    filter(!is.na(total_pr)) %>%
    group_by(vaccine, scenario, full_scenario) %>%
    reframe(mean = Rmisc::CI(total_pr)["mean"],
            lower = Rmisc::CI(total_pr)["lower"],
            upper = Rmisc::CI(total_pr)["upper"],
            ymin = 0,
            ymax = Inf,
            total_pr = mean) -> means_vacc
  
  return(means_vacc)
}

calc_mean_strain_pr <- function(total_pr){
  # summarise prs
  total_pr %>%
    filter(!is.na(total_pr)) %>%
    group_by(vaccine, scenario, full_scenario, variant) %>%
    reframe(mean = Rmisc::CI(strain_pr)["mean"],
            lower = Rmisc::CI(strain_pr)["lower"],
            upper = Rmisc::CI(strain_pr)["upper"],
            ymin = 0,
            ymax = Inf) -> means_strain
  
  return(means_strain)
}

plot_total_pr <- function(total_pr, means_vacc, plot_colors){
  
  # plot total pr
  total_pr %>%
    filter(!is.na(total_pr)) %>%
    select(scenario, vaccine, full_scenario, total_pr) %>%
    unique() %>%
    ggplot(aes(x = total_pr, fill = vaccine)) + 
    geom_rect(data = means_vacc, aes(ymin = ymin, ymax = ymax, xmin = lower, xmax = upper, fill = vaccine), color = NA, alpha = 0.3) +
    geom_vline(data = means_vacc, aes(xintercept = mean, color = vaccine)) +
    geom_histogram(alpha=0.6, position = 'identity', color = "grey40") + 
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    facet_wrap(~scenario) + 
    # theme_bw() + 
    theme(strip.background.x = element_blank()) +
    xlim(c(-0.02, 1.02)) + 
    coord_cartesian(expand = FALSE) +
    xlab(paste0("Total PR (", "\u2211", "Strain PR)")) -> p_vacc
  
  return(p_vacc)
  
}

plot_specific_scenario <- function(pr_data, target_scenario, means_vacc, means_strain, plot_colors){
  
  temp_data <- pr_data %>%
    filter(scenario == target_scenario)
  
  temp_means <- means_vacc %>%
    filter(scenario == target_scenario)
  
  temp_strains <- means_strain %>%
    filter(scenario == target_scenario) %>%
    mutate(strain_pr = mean)
  
  full_plot <- plot_total_pr(temp_data, temp_means, plot_colors)
  
  strain_plot <- temp_data %>%
    select(scenario, vaccine, variant, full_scenario, strain_pr) %>%
    unique() %>%
    ggplot(aes(x = strain_pr, fill = vaccine)) + 
    geom_rect(data = temp_strains, aes(ymin = ymin, ymax = ymax, xmin = lower, xmax = upper, fill = vaccine), color = NA, alpha = 0.3) +
    geom_vline(data = temp_strains, aes(xintercept = mean, color = vaccine)) +
    geom_histogram(alpha=0.6, position = 'identity', color = "grey40") + 
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    facet_wrap(~variant) + 
    # theme_bw() + 
    theme(strip.background.x = element_blank()) +
    xlim(c(-0.02, 1.02)) + 
    coord_cartesian(expand = FALSE) +
    xlab(paste0("Strain PR")) 
  
  comb_plot <- full_plot / strain_plot + plot_layout(guides = 'collect')
  
  return(comb_plot)
}

plot_all_scenarios_per_strain <- function(total_pr, means_vacc, means_strain, plot_colors){
  
  all_scenarios <- total_pr %>%
    filter(!is.na(total_pr)) %>%
    pull(scenario) %>%
    unique()
  
  all_scenarios_per_strain <- lapply(all_scenarios,
                                     function(x){plot_specific_scenario(total_pr, x, means_vacc, means_strain, plot_colors)})
  
  comb_plot <- wrap_plots(all_scenarios_per_strain, ncol = 1)
  
  return(comb_plot)
}

plot_pr_per_vaccine <- function(total_pr, means_vacc, plot_colors){
  
  # plot total pr
  total_pr %>%
    filter(!is.na(total_pr)) %>%
    select(scenario, vaccine, full_scenario, total_pr) %>%
    unique() %>%
    ggplot(aes(x = total_pr, fill = scenario)) + 
    geom_rect(data = means_vacc, aes(ymin = ymin, ymax = ymax, xmin = lower, xmax = upper, fill = vaccine), color = NA, alpha = 0.3) +
    geom_vline(data = means_vacc, aes(xintercept = mean, color = scenario)) +
    geom_histogram(alpha=0.6, position = 'identity', color = "grey40") + 
    facet_wrap(~vaccine) + 
    # theme_bw() + 
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    theme(strip.background.x = element_blank()) +
    scale_x_continuous(name = paste0("Total PR (", "\u2211", "Strain PR)"),
                       limits = c(-0.02, 1.02)) + 
    coord_cartesian(expand = FALSE)-> p_vacc
  
  return(p_vacc)
  
}

plot_base_circ_distributions <- function(total_pr, scenario_order, variant_order, vacc_order){
 
  total_pr %>%
    ungroup() %>%
    filter(circulation_rate != "NaN") %>%
    select(sample, vaccine, scenario, variant, circulation_rate, protection_rate) %>%
    mutate(vaccine = factor(vaccine, levels = vacc_order),
           scenario = factor(scenario, levels = scenario_order),
           variant = factor(variant, levels = variant_order[!is.na(variant_order)])) -> sub_sim
  
  sub_sim %>%
    filter(!is.na(circulation_rate)) %>%
    ggplot(aes(x = circulation_rate)) + 
    geom_histogram(alpha=0.6, position = 'identity', color = "grey40") + 
    facet_grid(variant~scenario,
               labeller = label_wrap_gen(width=10)) + 
    # theme_bw() + 
    theme(strip.background = element_blank()) +
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       limits = c(-0.02, 1.02)) +
    coord_cartesian(expand = FALSE) +
    xlab(paste0("Simulated Circulation Rates")) -> plot_sim_circ
  
  
  return(plot_sim_circ)
  
}

plot_base_pr_distributions <- function(total_pr, scenario_order, variant_order, vacc_order,
                                       plot_colors){
  
  total_pr %>%
    ungroup() %>%
    filter(circulation_rate != "NaN") %>%
    select(sample, vaccine, scenario, variant, circulation_rate, protection_rate) %>%
    mutate(vaccine = factor(vaccine, levels = vacc_order),
           scenario = factor(scenario, levels = scenario_order),
           variant = factor(variant, levels = variant_order[!is.na(variant_order)])) -> sub_sim
  
  
  sub_sim %>%
    filter(!is.na(protection_rate)) %>%
    ggplot(aes(x = protection_rate, fill = vaccine)) + 
    geom_histogram(alpha=0.6, position = 'identity', color = "grey40") + 
    facet_grid(variant ~ vaccine,
               labeller = label_wrap_gen(width=10)) + 
    # theme_bw() + 
    theme(strip.background = element_blank(),
          legend.position = "none") +
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       limits = c(-0.02, 1.02)) +
    scale_fill_manual(values = plot_colors) +
    coord_cartesian(expand = FALSE) +
    xlab(paste0("Simulated Protection Rates")) -> plot_sim_pr
  
  return(plot_sim_pr)
  
}

subset_to_target_vaccine <- function(total_pr, target_vaccines){
  
  target_vacc <- strsplit(target_vaccines, ",")[[1]]
  
  total_pr %>%
    filter(vaccine %in% target_vacc) -> sub_pr
  
  return(sub_pr)
  
}

calculate_dist_t_test <- function(total_pr){
  
  round_three_digits <- function(x){
    round(x, 3)
  }
  
  total_pr %>%
    ungroup() %>%
    filter(!is.na(total_pr)) %>%
    select(scenario, total_pr, vaccine) %>%
    unique() %>%
    group_by(scenario) %>%
    t_test(total_pr~vaccine, detailed = TRUE) %>%
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
  
  
  HTML(paste('Strain circulation values and protection rates are sampled from 
             a normal distribution, the distribution parameters 
             are user specified input. The Bayesian approach results in a distribution of total protection rates per vaccine.<br><br>')),
  
  
  sidebarLayout(
    sidebarPanel(
      textInput("selected_sheet", "Google sheet name",
                value = "Overview"),
      textInput("circ_range", "Range for circulation values, including strain names and circulation scenarios",
                value = "B6:N20"),
      textInput("pr_range", "Range for protection rate values, including strain names and vaccine names",
                value = "P7:AA21"),
      HTML(paste('<B><font size="2">URL to Google sheet:</B><br>')),
      HTML(paste('<n><font size="1">
                  Estimated mean protection rate and strain circulation, values between 0 and 1 are required. Samples will be drawn from 
                  a normal distribution with the specified mean from the google sheet and sd. Any samples > 1 will be set to 1, any samples < 0 will be set to 0.
                  ')),
      textInput("file_upload", "",
                value = ""),
      HTML(paste('<br><B><font size="2">Number of samples for the Bayesian summary:</B><br>')),
      HTML(paste('<n><font size="1">
                  Integer between 10 and 10000.
                  <B><font size="2">
                  ')),
      numericInput(
        inputId = "sim_n",         # ID to reference in server
        label = "",  # Label in the UI
        value = 50,              # Default value
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
      ),
      numericInput(
        inputId = "pr_sd",         # ID to reference in server
        label = "Standard deviation for protection rate (normal distribution)",  # Label in the UI
        value = 0.1,              # Default value
        min = 0.01,                 # Minimum allowed
        max = 1,              # Maximum allowed
        step = 0.01                 # Step size
      ),
      numericInput(
        inputId = "circ_sd",         # ID to reference in server
        label = "Standard deviation for expected circulation (normal distribution)",  # Label in the UI
        value = 0.1,              # Default value
        min = 0.01,                 # Minimum allowed
        max = 1,              # Maximum allowed
        step = 0.01                 # Step size
      ),
      HTML(paste('<B><font size="2">Target vaccine strain:</B><br>')),
      HTML(paste('<n><font size="1">
                 Put in your target vaccine strain as named in the google sheet. Separate multiple strains with a "," and no spaces. 
                  <B><font size="2">')),
      textInput("target_vaccine", "",
                value = "")
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Total protection rates by circulation scenario",
                 plotOutput("protection_dist_plot", height = 500),
             #    HTML(paste('<br><B><font size="2">Input mean circulation:</B><br>')),
                textOutput("stat_table"),
                 DTOutput("stat_dist_table"),
                 textOutput("error_msg")),
        tabPanel("Total protection rates by vaccine",
                 plotOutput("protection_dist_vacc_plot", height = 500),
                 textOutput("error_msg")),
        tabPanel("Strain specific protection rates by vaccine and circulation scenario",
                # DTOutput("csv_dist_table"),
                 plotOutput("protection_strain_plot"),
             #    DTOutput("stat_dist_table"),
                 textOutput("error_msg")),
        tabPanel("Target vaccine",
                 plotOutput("target_vaccine_plot", height = 6000),
                 textOutput("stat_table"),
                 DTOutput("stat_dist_table_target"),
                 textOutput("error_msg")),
        tabPanel("Input data",
                 textOutput("circulation_table"),
                 DTOutput("csv_table_circ"),
                 plotOutput("base_circ_distribution_plot", height = 800, width = "150%"),
                 textOutput("pr_table"),
                 DTOutput("csv_table_pr"),
                 plotOutput("base_pr_distribution_plot", height = 800, width ="180%"),
                 textOutput("error_msg"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  csv_data_circulation <- reactive({
    req(input$file_upload)
    req(input$selected_sheet)
    req(input$circ_range)
  #  observeEvent(input$submit_button, {
    tryCatch({
      googlesheets4::range_read(input$file_upload, range = paste0(input$selected_sheet,"!", input$circ_range))
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      NULL
  #  })
    })
  })
  
  csv_data_pr <- reactive({
    req(input$file_upload)
    req(input$selected_sheet)
    req(input$pr_range)
   # observeEvent(input$submit_button, {
    tryCatch({
      googlesheets4::range_read(input$file_upload, range = paste0(input$selected_sheet,"!", input$pr_range))
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      NULL
  #  })
    })
  })
  
  observe({set.seed(input$sim_seed)})
  
  observe(input$sim_n)
  
  observe(input$circ_sd)
  
  observe(input$pr_sd)
  
  total_pr <- reactive({
    req(csv_data_circulation())  # Ensure a file is selected
    tryCatch({
      simulate_total_pr_posterior(csv_data_circulation(),
                                  csv_data_pr(),
                                  input$circ_sd, input$pr_sd,
                                  n_sim = input$sim_n)
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      NULL
    })
  })
  
  # get factor levels
  levels_vacc <- reactive({
    req(csv_data_pr())
    colnames(csv_data_pr())[2:ncol(csv_data_pr())]
  })
  
  levels_scenario <- reactive({
    req(csv_data_circulation())
    colnames(csv_data_circulation())[2:ncol(csv_data_circulation())]
  })
  
  levels_variant <- reactive({
    req(csv_data_circulation())
    csv_data_circulation()[,1]
  })
  
  # get colors
  vacc_colors <- reactive({
    req(levels_vacc())
    get_color_pallette(levels_vacc(), NULL)
  })
  circ_colors <- reactive({
    req(levels_scenario())
    get_color_pallette(levels_scenario(), "Dark2")
  })
  
  output$protection_dist_plot <- renderPlot({
    req(total_pr())  # Ensure a file is selected
    req(vacc_colors())
    tryCatch({
      means_vacc <- calc_mean_vacc_pr(total_pr())
      plot_total_pr(total_pr(), means_vacc, vacc_colors())
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      NULL
    })
  },
  height = 500)
  
  output$protection_dist_vacc_plot <- renderPlot({
    req(total_pr())  # Ensure a file is selected
    req(circ_colors())
    tryCatch({
      means_vacc <- calc_mean_vacc_pr(total_pr())
      plot_pr_per_vaccine(total_pr(), means_vacc, circ_colors())
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      NULL
    })
  },
  height = 500)
  
  output$protection_strain_plot <- renderPlot({
    req(total_pr())  # Ensure a file is selected
    req(vacc_colors())
    tryCatch({
      means_vacc <- calc_mean_vacc_pr(total_pr())
      means_strain <- calc_mean_strain_pr(total_pr())
      plot_all_scenarios_per_strain(total_pr(), means_vacc, means_strain, vacc_colors())
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      NULL
    })
  },
  height = 6000)
  
  output$target_vaccine_plot <- renderPlot({
    req(total_pr())  # Ensure a file is selected
    req(input$target_vaccine)
    tryCatch({
      sub_pr <- subset_to_target_vaccine(total_pr(), input$target_vaccine)
      means_vacc <- calc_mean_vacc_pr(sub_pr)
      means_strain <- calc_mean_strain_pr(sub_pr)
      plot_all_scenarios_per_strain(sub_pr, means_vacc, means_strain, vacc_colors())
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      NULL
    })
  },
  height = 6000)
  
  # create plot for simulated expected circulation and protection rates
  output$base_circ_distribution_plot <- renderPlot({
    req(total_pr())
  #  req(levels_variant())
    tryCatch({
      plot_base_circ_distributions(total_pr(), 
                              variant_order = csv_data_circulation()[,1], 
                              vacc_order = colnames(csv_data_pr())[2:ncol(csv_data_pr())],
                              scenario_order = colnames(csv_data_circulation())[2:ncol(csv_data_circulation())])
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      NULL
    })
  },
  height = 800)
  
  output$base_pr_distribution_plot <- renderPlot({
    req(total_pr())
    #req(levels_variant())
    req(vacc_colors())
    tryCatch({
      plot_base_pr_distributions(total_pr(), 
                                 variant_order = csv_data_circulation()[,1], 
                                 vacc_order = colnames(csv_data_pr())[2:ncol(csv_data_pr())],
                                 scenario_order = colnames(csv_data_circulation())[2:ncol(csv_data_circulation())],
                                 plot_colors = vacc_colors())
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      NULL
    })
  },
  height = 800)
  
  output$circulation_table <- renderText("Input mean circulation:")
  output$csv_table_circ <- renderDT({
    dat <- csv_data_circulation()
    req(dat)  # Only render if data is available
    dat <- dat %>%
      rename(., variant = colnames(.)[1]) %>%
      filter(!is.na(variant))
    datatable(dat, options = list(pageLength = 15))
  })
  
  output$pr_table <- renderText("Input mean vaccine protection rates:")
  output$csv_table_pr <- renderDT({
    dat <- csv_data_pr()
    req(dat)  # Only render if data is available
    dat <- dat %>%
    rename(., variant = colnames(.)[1]) %>%
      filter(!is.na(variant))
    datatable(dat, options = list(pageLength = 15))
  })
  
  output$stat_table <- renderText("\n \n \nT test statistical comparison:")
  output$stat_dist_table <- renderDT({
    tab_out <- calculate_dist_t_test(total_pr())
    req(tab_out)  # Only render if data is available
    datatable(tab_out, options = list(pageLength = 20))
  })
  
  output$stat_dist_table_target <- renderDT({
    req(input$target_vaccine)
    sub_pr <- subset_to_target_vaccine(total_pr(), input$target_vaccine)
    tab_out <- calculate_dist_t_test(sub_pr)
    req(tab_out)  # Only render if data is available
    datatable(tab_out, options = list(pageLength = 20))
  })
  
}

# lets google sheets know that you don't need authentication
gs4_deauth()
shinyApp(ui, server)
