# Here, I just want to test simulation of the new pgs corr methodology
set.seed(320)
'%>%' <- magrittr::'%>%'
library(ggplot2)

Sigma <- list()
#Sigma[[1]] <- (.25 ^ 2) * matrix(data = 1, nrow = 2, ncol = 2)
#Sigma[[2]] <- (.25 ^ 2) * matrix(data = c(1, 0, 0, 0), nrow = 2, ncol = 2)
#Sigma[[3]] <- (.25 ^ 2) * matrix(data = c(0, 0, 0, 1), ncol = 2)

Sigma[[1]] <- matrix(data = 1, nrow = 2, ncol = 2)
Sigma[[2]] <- matrix(data = c(3, sqrt(3), sqrt(3), 1), ncol = 2)

maf_simulator <- function(n) {

  pmax(pmin(rbeta(n, .5, 10), .5 - 1e-05), .005)

}

metric_vec <- c()
model_type_vec <- c()
selection_method_vec <- c()
value_vec <- c()
amp_vec <- c()

for (amp in seq(0, 1, length.out = 8)) {

  print(amp)

  sim_out <- gamp::simulate_pgs_corr_fast_v2(
    n_e0 = ceiling(10000 / .75),
    n_e1 = ceiling(10000 / .75),
    n_snps = 1000,
    maf_simulator,
    e0_h2 = .5,
    e1_h2 = .5,
    Sigma = Sigma,
    pi = c(1 - amp, amp),
    num_sims = 25,
    s = .5,
    pval_thresh_additive = .1,
    pval_thresh_GxE = .1,
    selection_methods = c("additive", "GxE")
  )

  for (model_type in c("additive", "GxE")) {

    for (metric in c("mse", "power", "type_1_error", "accuracy", "corr")) {

      for (selection_method in c("pval", "same_number")) {

        metric_vec <- c(metric_vec, metric)
        model_type_vec <- c(model_type_vec, model_type)
        value_vec <- c(value_vec, sim_out[[selection_method]][[model_type]][[metric]])
        selection_method_vec <- c(selection_method_vec, selection_method)
        amp_vec <- c(amp_vec, amp)

      }

    }

  }

  readr::write_rds(sim_out, glue::glue("rds_data/polymetric_{amp}_zhu2_test.rds"))

}

# out_df <- data.frame(
#   metric = metric_vec,
#   amp = amp_vec,
#   model_type = model_type_vec,
#   selection_method = selection_method_vec,
#   value = value_vec
# )

# Now, want to make plots of the power using:
# (1) Threshold - Additive
# (2) Threshold - GxE
# (3) Same Number

# out_df <- out_df %>%
#   dplyr::filter(!(selection_method == "same_number" & model_type == "additive"))
#
# out_df <- out_df %>%
#   dplyr::mutate(model_method = paste(model_type, selection_method))
#
# ggplot(data = out_df %>% dplyr::filter(metric == "power")) +
#   geom_point(aes(x = amp, y = value, color = model_method)) +
#   geom_line(aes(x = amp, y = value, color = model_method)) +
#   xlab("Percent Independent") +
#   ylab("Power")
#
# ggplot(data = out_df %>% dplyr::filter(metric == "mse")) +
#   geom_point(aes(x = amp, y = value, color = model_method)) +
#   geom_line(aes(x = amp, y = value, color = model_method)) +
#   xlab("Percent Independent") +
#   ylab("MSE")
#
# ggplot(data = out_df %>% dplyr::filter(metric == "corr")) +
#   geom_point(aes(x = amp, y = value, color = model_method)) +
#   geom_line(aes(x = amp, y = value, color = model_method)) +
#   xlab("Percent Independent") +
#   ylab("Corr")

# Here, I want to extract the different dataframes so that I can plot them
# I don't think that this should be too hard

metric_vec <- c()
value_vec <- c()
model_type_vec <- c()
amp_vec <- c()

for (amp in seq(0, 1, .25)) {

  sim_data <- readr::read_rds(
    glue::glue("rds_data/polymetric_{amp}_zhu_test.rds")
  )
  for (metric in c("mse", "power", "type_1_error")) {

    for (model_type in c("GxE", "additive", "mash")) {

      amp_vec <- c(amp_vec, amp)
      model_type_vec <- c(model_type_vec, model_type)
      metric_vec <- c(metric_vec, metric)
      value_vec <- c(value_vec, sim_data[[model_type]][[metric]])

    }

  }

}

out_df <- data.frame(
  model_type = model_type_vec,
  amp = amp_vec,
  metric = metric_vec,
  value = value_vec
)

ggplot(data = out_df %>% dplyr::filter(metric == "power")) +
  geom_point(aes(x = amp, y = value, color = model_type)) +
  geom_line(aes(x = amp, y = value, color = model_type)) +
  xlab("Percent Independent") +
  ylab("Power")

ggplot(data = out_df %>% dplyr::filter(metric == "mse")) +
  geom_point(aes(x = amp, y = value, color = model_type)) +
  geom_line(aes(x = amp, y = value, color = model_type)) +
  xlab("Percent Independent") +
  ylab("MSE")

ggplot(data = out_df %>% dplyr::filter(metric == "type_1_error")) +
  geom_point(aes(x = amp, y = value, color = model_type)) +
  geom_line(aes(x = amp, y = value, color = model_type)) +
  xlab("Percent Independent") +
  ylab("Type 1 Error")

# Now, want to make a correlation dataframe
# I don't think that this should be too hard as long as I pay attention
# to each of the different variables going in

corr_vec <- c()
amp_vec <- c()
selection_method_vec <- c()
beta_method_vec <- c()

for (selection_method in c("additive", "GxE")) {

  for (beta_method in c("additive", "GxE")) {

    for (amp in seq(0, 1, length.out = 8)) {

      sim_data <- readr::read_rds(
        glue::glue("rds_data/polymetric_{amp}_zhu2_test.rds")
      )

      amp_vec <- c(amp_vec, amp)
      selection_method_vec <- c(selection_method_vec, selection_method)
      beta_method_vec <- c(beta_method_vec, beta_method)
      corr_vec <- c(corr_vec, sim_data[[selection_method]][[beta_method]][['corr']])

    }

  }

}

corr_df <- data.frame(
  selection_method = selection_method_vec,
  beta_method = beta_method_vec,
  amp = amp_vec,
  corr = corr_vec
)

corr_df <- corr_df %>%
  dplyr::filter(beta_method == selection_method)

corr_df <- corr_df %>%
  dplyr::mutate(
    selection_method = dplyr::if_else(selection_method == "mash", "all", selection_method)
  )

ggplot(data =corr_df) +
  geom_point(aes(x = amp, y = corr, color = selection_method), size = 1) +
  geom_line(aes(x = amp, y = corr, color = selection_method)) +
  xlab("Percent of SNPs 3x Amplified in Environment 1") +
  ylab("PGS Test Set Correlation in Both Environments") +
  labs(color="Model")

ggplot(data =corr_df) +
  geom_point(aes(x = amp, y = corr, color = selection_method), size = 1) +
  geom_line(aes(x = amp, y = corr, color = selection_method, linetype = beta_method)) +
  xlab("Percent Amplified") +
  ylab("Corr") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"))
