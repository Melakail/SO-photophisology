
# ============================================================
# 1. Load packages and data
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(tibble)

fluorescence <- read_excel(
  "md_R_analysis_ready.xlsx",
  sheet = "fluorescence_replicates"
)

summaries <- read_excel(
  "md_R_analysis_ready.xlsx",
  sheet = "treatment_summaries"
)

deltas <- read_excel(
  "md_R_analysis_ready.xlsx",
  sheet = "delta_comparisons"
)

parameters <- read_excel(
  "md_R_analysis_ready.xlsx",
  sheet = "station_parameters"
)


# ============================================================
# 2. Clean delta comparison data
# ============================================================

deltas <- deltas %>%
  filter(comparison != "Calculation") %>%
  mutate(
    across(
      c(delta_Fv_Fm, mean_delta_Fv_Fm, sd_delta_Fv_Fm),
      ~ parse_number(as.character(.x))
    )
  )


# ============================================================
# 3. Summarise raw fluorescence measurements
# ============================================================

treatment_summary <- fluorescence %>%
  filter(treatment %in% c("Control", "Fe")) %>%
  group_by(station, treatment) %>%
  summarise(
    mean_Fv_Fm = mean(Fv_Fm, na.rm = TRUE),
    sd_Fv_Fm = sd(Fv_Fm, na.rm = TRUE),
    mean_sigmaPSII = mean(sigmaPSII, na.rm = TRUE),
    sd_sigmaPSII = sd(sigmaPSII, na.rm = TRUE),
    n = sum(!is.na(Fv_Fm)),
    .groups = "drop"
  )

treatment_summary


# ============================================================
# 4. Identify incubation replicates
# ============================================================

# Check the number of observations per station and treatment

fluorescence %>%
  count(station, treatment)

summaries %>%
  count(station, group)


# Extract the three Control and Fe incubation replicates

incubation_replicates <- summaries %>%
  filter(group != "Tzero") %>%
  mutate(
    treatment = if_else(
      grepl("^Control", group),
      "Control",
      "Fe"
    ),
    replicate = as.integer(
      sub(".* ", "", group)
    )
  )

incubation_replicates %>%
  count(station, treatment)


# ============================================================
# 5. Calculate station-level treatment means and SD
# ============================================================

treatment_comparison <- incubation_replicates %>%
  group_by(station, treatment) %>%
  summarise(
    Fv_Fm_mean = mean(mean_Fv_Fm, na.rm = TRUE),
    Fv_Fm_sd = sd(mean_Fv_Fm, na.rm = TRUE),
    sigmaPSII_mean = mean(mean_sigmaPSII, na.rm = TRUE),
    sigmaPSII_sd = sd(mean_sigmaPSII, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

treatment_comparison


# ============================================================
# 6. Calculate Fe treatment differences
# ============================================================

treatment_differences <- treatment_comparison %>%
  select(station, treatment, Fv_Fm_mean, sigmaPSII_mean) %>%
  pivot_wider(
    names_from = treatment,
    values_from = c(Fv_Fm_mean, sigmaPSII_mean)
  ) %>%
  mutate(
    delta_Fv_Fm = Fv_Fm_mean_Fe - Fv_Fm_mean_Control,
    delta_sigmaPSII = sigmaPSII_mean_Fe - sigmaPSII_mean_Control
  )

treatment_differences


# ============================================================
# 7. Plot change in Fv/Fm following Fe addition
# ============================================================

ggplot(
  treatment_differences,
  aes(x = station, y = delta_Fv_Fm)
) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Station",
    y = expression(Delta * F[v] / F[m]),
    title = "Response of Fv/Fm to iron addition"
  ) +
  theme_classic()


# ============================================================
# 8. Plot changes in Fv/Fm and sigmaPSII
# ============================================================

# Scale sigmaPSII changes to the Fv/Fm plotting range

scale_factor <- max(abs(treatment_differences$delta_Fv_Fm)) /
  max(abs(treatment_differences$delta_sigmaPSII))

plot_data <- treatment_differences %>%
  select(station, delta_Fv_Fm, delta_sigmaPSII) %>%
  pivot_longer(
    cols = c(delta_Fv_Fm, delta_sigmaPSII),
    names_to = "parameter",
    values_to = "delta"
  ) %>%
  mutate(
    scaled_delta = if_else(
      parameter == "delta_sigmaPSII",
      delta * scale_factor,
      delta
    )
  )

ggplot(
  plot_data,
  aes(x = station, y = scaled_delta, fill = parameter)
) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c(
      "delta_Fv_Fm" = "#0072B2",
      "delta_sigmaPSII" = "#D55E00"
    ),
    labels = c(
      "delta_Fv_Fm" = expression(Delta * F[v] / F[m]),
      "delta_sigmaPSII" = expression(Delta * sigma[PSII])
    )
  ) +
  scale_y_continuous(
    name = expression(Delta * F[v] / F[m]),
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = expression(Delta * sigma[PSII])
    )
  ) +
  labs(
    x = "Station",
    fill = "Parameter",
    title = "Photophysiological response to iron addition"
  ) +
  theme_classic() +
  theme(
    legend.position = "top"
  )


# ============================================================
# 9. Clean environmental and phytoplankton data
# ============================================================

parameters_clean <- parameters %>%
  mutate(
    value_mean = if_else(
      parameter == "Z",
      NA_real_,
      value_original %>%
        str_replace("±.*$", "") %>%
        str_replace_all("×10", "e") %>%
        str_replace_all("\\s+", "") %>%
        parse_number()
    )
  )


# Check nutrient concentrations

parameters_clean %>%
  filter(
    parameter %in% c("Nitrate", "Phosphate", "Silicate", "dFe")
  ) %>%
  select(station, parameter, value_mean) %>%
  pivot_wider(
    names_from = parameter,
    values_from = value_mean
  ) %>%
  arrange(station) %>%
  print(n = Inf)


# ============================================================
# 10. Prepare environmental data
# ============================================================

environmental_variables <- c(
  "Nitrate",
  "Phosphate",
  "Silicate",
  "dFe",
  "SST",
  "MLD",
  "zeu"
)

environmental_data <- parameters_clean %>%
  filter(parameter %in% environmental_variables) %>%
  select(station, parameter, value_mean) %>%
  pivot_wider(
    names_from = parameter,
    values_from = value_mean
  ) %>%
  left_join(
    treatment_differences %>%
      select(station, delta_Fv_Fm, delta_sigmaPSII),
    by = "station"
  )


# ============================================================
# 11. Environmental correlation analysis
# ============================================================

environmental_cor <- environmental_data %>%
  select(-station) %>%
  cor(
    method = "spearman",
    use = "pairwise.complete.obs"
  )

round(environmental_cor, 3)


# Calculate the number of observations used for each
# pairwise correlation

environmental_matrix <- as.matrix(
  environmental_data %>%
    select(-station)
)

environmental_n <- crossprod(
  !is.na(environmental_matrix)
)

environmental_n


# Check the number of observations, unique values and SD

environmental_data %>%
  summarise(
    across(
      -station,
      list(
        n = ~ sum(!is.na(.x)),
        unique = ~ n_distinct(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      )
    )
  )


# Convert correlation matrix to long format

environmental_cor_long <- as.data.frame(environmental_cor) %>%
  rownames_to_column("Variable1") %>%
  pivot_longer(
    cols = -Variable1,
    names_to = "Variable2",
    values_to = "Correlation"
  )


# Environmental correlation heatmap

ggplot(
  environmental_cor_long,
  aes(x = Variable1, y = Variable2, fill = Correlation)
) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = ifelse(
        is.na(Correlation),
        "NA",
        round(Correlation, 2)
      )
    ),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey80"
  ) +
  coord_fixed() +
  labs(
    title = "Environmental and photophysiological correlations",
    x = NULL,
    y = NULL,
    fill = "Spearman's rho"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


# ============================================================
# 12. Prepare phytoplankton–environment data
# ============================================================

phytoplankton_groups <- c(
  "Diatoms",
  "Haptophytes",
  "Prasinophytes",
  "Dinoflagellates",
  "Chlorophytes",
  "Cyanobacteria",
  "Cryptophytes",
  "Pelagophytes"
)

taxa_environment_data <- parameters_clean %>%
  filter(
    parameter %in% c(
      phytoplankton_groups,
      environmental_variables
    )
  ) %>%
  select(station, parameter, value_mean) %>%
  pivot_wider(
    names_from = parameter,
    values_from = value_mean
  ) %>%
  left_join(
    treatment_differences %>%
      select(station, delta_Fv_Fm, delta_sigmaPSII),
    by = "station"
  )


# Check phytoplankton proportions

taxa_environment_data %>%
  mutate(
    phytoplankton_total = rowSums(
      across(all_of(phytoplankton_groups)),
      na.rm = TRUE
    )
  ) %>%
  select(station, phytoplankton_total)


# Inspect phytoplankton composition by station

taxa_environment_data %>%
  select(station, all_of(phytoplankton_groups)) %>%
  print(n = Inf)


# ============================================================
# 13. Phytoplankton–environment correlation analysis
# ============================================================

taxa_environment_cor <- taxa_environment_data %>%
  select(-station) %>%
  cor(
    method = "spearman",
    use = "pairwise.complete.obs"
  )


# Extract phytoplankton correlations with environmental
# variables and photophysiological responses

taxa_environment_cor_long <- as.data.frame(
  taxa_environment_cor
) %>%
  rownames_to_column("Variable1") %>%
  pivot_longer(
    cols = -Variable1,
    names_to = "Variable2",
    values_to = "Correlation"
  ) %>%
  filter(
    Variable1 %in% phytoplankton_groups,
    Variable2 %in% c(
      environmental_variables,
      "delta_Fv_Fm",
      "delta_sigmaPSII"
    )
  )


# Print correlations, ordered by absolute strength

taxa_environment_cor_long %>%
  arrange(Variable1, desc(abs(Correlation))) %>%
  print(n = Inf)


# Phytoplankton–environment correlation heatmap

ggplot(
  taxa_environment_cor_long,
  aes(
    x = Variable2,
    y = Variable1,
    fill = Correlation
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = ifelse(
        is.na(Correlation),
        "NA",
        round(Correlation, 2)
      )
    ),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey80"
  ) +
  labs(
    title = "Phytoplankton–environment correlations",
    x = NULL,
    y = NULL,
    fill = "Spearman's rho"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

