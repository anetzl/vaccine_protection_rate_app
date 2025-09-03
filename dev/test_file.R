library(shiny)
library(DT)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(readxl)
library(rstatix)
library(googlesheets4)

test_path <- "https://docs.google.com/spreadsheets/d/1HN37wI3sHbwiIgfdDKIGT_IUMwspR1kmv4DODpEvY7Y/edit?gid=1659497864#gid=1659497864"

# these work for reading
data_circulation <- googlesheets4::range_read(test_path, range = "Overview!B6:N19") %>%
  rename(., variant = colnames(.)[1]) %>%
  filter(!is.na(variant))

data_pr <- googlesheets4::range_read(test_path, range = "Overview!P7:AA21") %>%
  rename(., variant = colnames(.)[1]) %>%
  filter(!is.na(variant)) %>%
  filter(variant %in% data_circulation$variant)

  
# input data by person
target_scenario <- "Approx current"
sd_pr <- 0.1
sd_circ <- 0.1
n_sim <- 20

# sampling posteriors
sample_posterior <- function(mean_norm = 0.1, sd_norm = 0.1, n_sim = 100){
  res <- rnorm(n_sim, mean = mean_norm, sd = sd_norm)
  res[is.na(res)] <- 0
  res[res < 0] <- 0
  res
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
    mutate(total_pr = sum(strain_pr))
  
  return(total_pr)
}


calc_mean_vacc_pr <- function(total_pr){
  # summarise prs
  total_pr %>%
    filter(!is.na(total_pr)) %>%
    group_by(vaccine, scenario, full_scenario) %>%
    reframe(mean = Rmisc::CI(total_pr)["mean"],
            lower = Rmisc::CI(total_pr)["lower"],
            upper = Rmisc::CI(total_pr)["upper"]) -> means_vacc
  
  return(means_vacc)
}

calc_mean_strain_pr <- function(total_pr){
  # summarise prs
  total_pr %>%
    filter(!is.na(total_pr)) %>%
    group_by(vaccine, scenario, variant, full_scenario) %>%
    reframe(mean = Rmisc::CI(strain_pr)["mean"],
            lower = Rmisc::CI(strain_pr)["lower"],
            upper = Rmisc::CI(strain_pr)["upper"]) -> means_strain
  
  return(means_strain)
}

plot_total_pr <- function(total_pr, means_vacc, means_strain){
  
  # plot total pr
  total_pr %>%
    filter(!is.na(total_pr)) %>%
    select(scenario, vaccine, full_scenario, total_pr) %>%
    unique() %>%
    ggplot(aes(x = total_pr, fill = vaccine)) + 
    geom_vline(data = means_vacc, aes(xintercept = mean, color = vaccine)) +
    geom_vline(data = means_vacc, aes(xintercept = lower, color = vaccine), linetype = "dashed") +
    geom_vline(data = means_vacc, aes(xintercept = upper, color = vaccine), linetype = "dashed") +
    geom_histogram(alpha=0.6, position = 'identity', color = "grey40") + 
    facet_wrap(~scenario) + 
    theme_bw() + 
    theme(strip.background.x = element_blank()) +
    xlim(c(0, 1)) + 
    xlab(paste0("Total PR (", "\u2211", "Strain PR)")) -> p_vacc
  
  return(p_vacc)
  
}

plot_specific_scenario <- function(pr_data, target_scenario, means_vacc, means_strain){
  
  temp_data <- pr_data %>%
    filter(scenario == target_scenario)
  
  temp_means <- means_vacc %>%
    filter(scenario == target_scenario)

  temp_strains <- means_strain %>%
    filter(scenario == target_scenario)  
  
  full_plot <- plot_total_pr(temp_data, temp_means, temp_strains)
  
  strain_plot <- temp_data %>%
    select(scenario, vaccine, variant, full_scenario, strain_pr) %>%
    unique() %>%
    ggplot(aes(x = strain_pr, fill = vaccine)) + 
    geom_vline(data = temp_strains, aes(xintercept = mean, color = vaccine)) +
    geom_vline(data = temp_strains, aes(xintercept = lower, color = vaccine), linetype = "dashed") +
    geom_vline(data = temp_strains, aes(xintercept = upper, color = vaccine), linetype = "dashed") +
    geom_histogram(alpha=0.6, position = 'identity', color = "grey40") + 
    facet_wrap(~variant) + 
    theme_bw() + 
    theme(strip.background.x = element_blank()) +
    xlim(c(0, 1)) + 
    xlab(paste0("Strain PR")) 
  
  comb_plot <- full_plot / strain_plot + plot_layout(guides = 'collect')
  
  return(comb_plot)
}

plot_all_scenarios_per_strain <- function(total_pr, means_vacc, means_strain){
  
  all_scenarios <- total_pr %>%
    filter(!is.na(total_pr)) %>%
    pull(scenario) %>%
    unique()
  
  all_scenarios_per_strain <- lapply(all_scenarios,
                                     function(x){plot_specific_scenario(total_pr, x, means_vacc, means_strain)})
  
  comb_plot <- wrap_plots(all_scenarios_per_strain[1:2], ncol = 1)
  
  return(comb_plot)
}

total_pr <- simulate_total_pr_posterior(data_circulation, data_pr,
                                        sd_circ, sd_pr, 
                                        n_sim = n_sim)

means_vacc <- calc_mean_vacc_pr(total_pr)

means_strain <- calc_mean_strain_pr(total_pr)

final_full_plot <- plot_total_pr(total_pr, means_vacc, means_strain)

final_strain_plot <- plot_all_scenarios_per_strain(total_pr, means_vacc, means_strain)

plot_specific_scenario(total_pr, "Approx current", means_vacc, means_strain)

total_pr %>%
  filter(variant == "145S") %>%
  View()

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

  
# plot simulated distributions


