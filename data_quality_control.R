library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# merged_data = reference(AF2) structures and ranked_0 sample(AF3) structures merged in a table
# merged_data_avg = reference(AF2) structures and averaged samples(AF3) structures merged in a table
# sample = AF3 ref = AF2

#---LOADING THE FILES AND PREPARING THE MERGED TABLE WITH AF2 AND AF2---

# Load both files. samples are AF3 structures while references are AF2 structures.
samples <- read_excel("/Users/imb/Desktop/AF3_hydrogens_metrics.xlsx")
reference <- read_excel("/Users/imb/Desktop/reference_metrics.xlsx", sheet = "AF-MMv2.2 result")

# Filter only ranked_0 from samples for the best ranked samples for all known extension samples
ranked_0_table <- samples %>%
  filter(model_id == "ranked_0")

# Add prediction group column to the ranked_0 table. Prediction groups are minimal names 
# for the prediction_names and they are for easy access
ranked_0_table <- ranked_0_table %>%
  mutate(prediction_group = str_extract(prediction_name, "^[^M]+"))


# Merge regarding both the samples (best ranks) and reference
merged_data <- merge(ranked_0_table, reference,
                     by = c("prediction_name"),
                     suffixes = c("_sample", "_ref"))

#Calculate differences
merged_data <- merged_data %>%
  mutate(
    diff_RMSD_backbone_peptide = RMSD_backbone_peptide_sample - RMSD_backbone_peptide_ref,
    perc_diff_RMSD_backbone_peptide = 100 * (RMSD_backbone_peptide_sample - RMSD_backbone_peptide_ref) / RMSD_backbone_peptide_ref,
    
    diff_RMSD_all_atom_peptide = RMSD_all_atom_peptide_sample - RMSD_all_atom_peptide_ref,
    perc_diff_RMSD_all_atom_peptide = 100 * (RMSD_all_atom_peptide_sample - RMSD_all_atom_peptide_ref) / RMSD_all_atom_peptide_ref,
    
    diff_RMSD_domain = RMSD_domain_sample - RMSD_domain_ref,
    perc_diff_RMSD_domain = 100 * (RMSD_domain_sample - RMSD_domain_ref) / RMSD_domain_ref,
    
    diff_DockQ = DockQ_sample - DockQ_ref,
    perc_diff_DockQ = 100 * (DockQ_sample - DockQ_ref) / DockQ_ref,
    
    diff_iPAE = iPAE_sample - iPAE_ref,
    perc_diff_iPAE = 100 * (iPAE_sample - iPAE_ref) / iPAE_ref,
    
    diff_pDockQ = pDockQ_sample - pDockQ_ref,
    perc_diff_pDockQ = 100 * (pDockQ_sample - pDockQ_ref) / pDockQ_ref
    
  )


# Check the intersection of columns to be sure
intersect(names(samples), names(reference))

# Check the names in the merged data to see if it has both ref and sample
names(merged_data)

# How many matched samples you now have
nrow(merged_data)  

# Add numbers as prediction names
merged_data$PredictionID <- seq_len(nrow(merged_data))

# Exporting the ID mapping if needed
write.csv(merged_data[, c("PredictionID", "prediction_name")], "prediction_id_mapping.csv", row.names = FALSE)




#---PLOTS---

# BAR PLOTS of percentage differences for all metrics

# 1. RMSD Backbone Peptide
ggplot(merged_data, aes(x = PredictionID, y = perc_diff_RMSD_backbone_peptide)) +
  geom_col(fill = "#5F9EA0") +
  labs(title = "Percent Difference in RMSD Backbone Peptide",
       y = "% Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. RMSD All Atom Peptide
ggplot(merged_data, aes(x = PredictionID, y = perc_diff_RMSD_all_atom_peptide)) +
  geom_col(fill = "#9ACD32") +
  labs(title = "Percent Difference in RMSD All Atom Peptide",
       y = "% Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. RMSD Domain
ggplot(merged_data, aes(x = PredictionID, y = perc_diff_RMSD_domain)) +
  geom_col(fill = "#D2691E") +
  labs(title = "Percent Difference in RMSD Domain",
       y = "% Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 4. DockQ
ggplot(merged_data, aes(x = PredictionID, y = perc_diff_DockQ)) +
  geom_col(fill = "#7B68EE") +
  labs(title = "Percent Difference in DockQ",
       y = "% Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 5. iPAE
ggplot(merged_data, aes(x = PredictionID, y = perc_diff_iPAE)) +
  geom_col(fill = "#F08080") +
  labs(title = "Percent Difference in iPAE",
       y = "% Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 6. pDockQ
ggplot(merged_data, aes(x = PredictionID, y = perc_diff_pDockQ)) +
  geom_col(fill = "#00CED1") +
  labs(title = "Percent Difference in pDockQ",
       y = "% Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# BAR PLOTS for absolute differences

# RMSD Backbone Peptide
ggplot(merged_data, aes(x = PredictionID, y = diff_RMSD_backbone_peptide)) +
  geom_col(fill = "#4682B4") +
  labs(title = "Absolute Difference in RMSD Backbone Peptide",
       y = "Absolute Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# RMSD All Atom Peptide
ggplot(merged_data, aes(x = PredictionID, y = diff_RMSD_all_atom_peptide)) +
  geom_col(fill = "#6B8E23") +
  labs(title = "Absolute Difference in RMSD All Atom Peptide",
       y = "Absolute Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# RMSD Domain
ggplot(merged_data, aes(x = PredictionID, y = diff_RMSD_domain)) +
  geom_col(fill = "#D2691E") +
  labs(title = "Absolute Difference in RMSD Domain",
       y = "Absolute Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# DockQ
ggplot(merged_data, aes(x = PredictionID, y = diff_DockQ)) +
  geom_col(fill = "#8A2BE2") +
  labs(title = "Absolute Difference in DockQ",
       y = "Absolute Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# iPAE
ggplot(merged_data, aes(x = PredictionID, y = diff_iPAE)) +
  geom_col(fill = "#CD5C5C") +
  labs(title = "Absolute Difference in iPAE",
       y = "Absolute Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# pDockQ
ggplot(merged_data, aes(x = PredictionID, y = diff_pDockQ)) +
  geom_col(fill = "#20B2AA") +
  labs(title = "Absolute Difference in pDockQ",
       y = "Absolute Difference", x = "Prediction #") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# SCATTER PLOTS

# RMSD Backbone Peptide
ggplot(merged_data, aes(x = RMSD_backbone_peptide_ref, y = RMSD_backbone_peptide_sample)) +
  geom_point(color = "#4682B4", size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "RMSD Backbone Peptide: Sample vs Reference",
       x = "Reference", y = "Sample") +
  theme_minimal()

# RMSD All Atom Peptide
ggplot(merged_data, aes(x = RMSD_all_atom_peptide_ref, y = RMSD_all_atom_peptide_sample)) +
  geom_point(color = "#6B8E23", size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "RMSD All Atom Peptide: Sample vs Reference",
       x = "Reference", y = "Sample") +
  theme_minimal()

# RMSD Domain
ggplot(merged_data, aes(x = RMSD_domain_ref, y = RMSD_domain_sample)) +
  geom_point(color = "#D2691E", size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "RMSD Domain: Sample vs Reference",
       x = "Reference", y = "Sample") +
  theme_minimal()

# DockQ
ggplot(merged_data, aes(x = DockQ_ref, y = DockQ_sample)) +
  geom_point(color = "#8A2BE2", size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "DockQ: Sample vs Reference",
       x = "Reference", y = "Sample") +
  theme_minimal()

# iPAE
ggplot(merged_data, aes(x = iPAE_ref, y = iPAE_sample)) +
  geom_point(color = "#CD5C5C", size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "iPAE: Sample vs Reference",
       x = "Reference", y = "Sample") +
  theme_minimal()

# pDockQ
ggplot(merged_data, aes(x = pDockQ_ref, y = pDockQ_sample)) +
  geom_point(color = "#20B2AA", size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "pDockQ: Sample vs Reference",
       x = "Reference", y = "Sample") +
  theme_minimal()

# Extreme RMSD for AF2
extreme_rmsd <- merged_data[merged_data$diff_RMSD_backbone_peptide <= -100,
                            c("prediction_name", "RMSD_backbone_peptide_sample", "RMSD_backbone_peptide_ref")]
print(extreme_rmsd)
View(extreme_rmsd)

# Extreme cases  for AF3
extreme_rmsd_failures <- merged_data %>%
  filter(RMSD_backbone_peptide_sample > 50,
         RMSD_backbone_peptide_ref < 5) %>%
  select(prediction_name,
         RMSD_backbone_peptide_sample,
         RMSD_backbone_peptide_ref,
         diff_RMSD_backbone_peptide)

# View it
View(extreme_rmsd_failures)


# See unmatched references to samples, if needed
unmatched_refs <- anti_join(reference, samples, by = c("prediction_name", "model_id"))
unmatched_prediction_names <- unmatched_refs$prediction_name
print(unmatched_prediction_names)





# DIFFERENT PLOT TYPES FOR EACH CLASS

# Scatter plot
ggplot(merged_data, aes(x = domain_ext_index, y = DockQ_sample)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, color = "#2c7fb8") +  
  geom_smooth(method = "loess", color = "#de2d26", se = TRUE) +       
  labs(title = "DockQ Accuracy vs Domain Extension Length",
       x = "Domain Extension Index",
       y = "DockQ (Sample)") +
  theme_minimal()

# Boxplot
ggplot(merged_data, aes(x = factor(domain_ext_index), y = DockQ_sample)) +
  geom_boxplot(fill = "#41b6c4") +
  labs(title = "DockQ Accuracy Distribution by Domain Extension Index",
       x = "Domain Extension Index",
       y = "DockQ (Sample)") +
  theme_minimal()

# Faceted scatter plot 

ggplot(
  merged_data %>% filter(motif_ext_index == 0),
  aes(x = domain_ext_index, y = DockQ_sample)
) +
  geom_point(color = "#2c7fb8") +
  geom_smooth(method = "loess", se = FALSE, color = "#de2d26") +
  facet_wrap(~ prediction_group) +
  labs(
    title = "DockQ Accuracy vs Domain Extension (Motif Extension Fixed at 0)",
    x = "Domain Extension Index",
    y = "DockQ (Sample)"
  ) +
  theme_minimal()

ggplot(
  merged_data %>% filter(domain_ext_index == 0),
  aes(x = motif_ext_index, y = DockQ_sample)
) +
  geom_point(color = "#2c7fb8", alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, color = "#de2d26") +
  facet_wrap(~ prediction_group) +
  labs(
    title = "DockQ Accuracy vs Motif Extension (Domain Extension Fixed at 0)",
    x = "Motif Extension Index",
    y = "DockQ (Sample)"
  ) +
  theme_minimal()



# PRACTICAL FUNCTION TO CALL FOR FACETED SCATTER PLOT

plot_extension_vs_metric <- function(data, metric, extension_col, title_prefix = "") {
  # Automatically identify the other extension column
  other_extension_col <- ifelse(extension_col == "domain_ext_index", "motif_ext_index", "domain_ext_index")
  
  # Filter data: fix the other extension to 0 and exclude RMSD values > 15
  filtered_data <- data %>%
    filter(!!sym(other_extension_col) == 0) %>%
    filter(!!sym(metric) <= 15)
  
  ggplot(filtered_data, aes_string(x = extension_col, y = metric)) +
    geom_point(color = "#2c7fb8", alpha = 0.7) +
    geom_smooth(method = "loess", se = FALSE, color = "#de2d26") +
    facet_wrap(~ prediction_group) +
    labs(
      title = paste0(title_prefix, metric, " vs ", extension_col),
      x = extension_col,
      y = metric
    ) +
    theme_minimal()
}

# RMSD backbone peptide
plot_extension_vs_metric(merged_data, "RMSD_backbone_peptide_sample", "domain_ext_index", "Domain Extension: ")
plot_extension_vs_metric(merged_data, "RMSD_backbone_peptide_sample", "motif_ext_index", "Motif Extension: ")

# RMSD all atom peptide
plot_extension_vs_metric(merged_data, "RMSD_all_atom_peptide_sample", "domain_ext_index", "Domain Extension: ")
plot_extension_vs_metric(merged_data, "RMSD_all_atom_peptide_sample", "motif_ext_index", "Motif Extension: ")

# RMSD domain
plot_extension_vs_metric(merged_data, "RMSD_domain_sample", "domain_ext_index", "Domain Extension: ")
plot_extension_vs_metric(merged_data, "RMSD_domain_sample", "motif_ext_index", "Motif Extension: ")


# It should be noted that the names for the prediction groups are not prediction_names but prediction_groups,
# and the reason they look shorter (DEG_ for example is DEG_MDM2_SWIB_1) is because the code classifies them
# as they first encounter the letter 'M'. This worked for these 30 protein groups but it is not a
# very generalizable classification approach. So you might need to change it for different groups.




# TABLE FOR AVG(SAMPLE)-REF DIFFERENCE

# First lets create the average sample table 

# Columns to average
cols_to_avg <- c(
  "RMSD_domain",
  "RMSD_backbone_peptide",
  "RMSD_all_atom_peptide",
  "RMSD_all_atom",
  "iPAE",
  "pDockQ",
  "DockQ"
)

average_table <- samples %>%
  group_by(prediction_name) %>%
  summarise(across(all_of(cols_to_avg), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  left_join(
    samples %>% distinct(prediction_name, .keep_all = TRUE) %>% select(-all_of(cols_to_avg)),
    by = "prediction_name"
  )

# Add prediction group column to the average table. Prediction groups are minimal names 
# for the prediction_names and they are for easy access
average_table <- average_table %>%
  mutate(prediction_group = str_extract(prediction_name, "^[^M]+"))


# Merge regarding both the samples (averaged) and reference
merged_data_avg <- merge(average_table, reference,
                     by = c("prediction_name"),
                     suffixes = c("_sample", "_ref"))

#Calculate differences
merged_data_avg <- merged_data_avg %>%
  mutate(
    diff_RMSD_backbone_peptide = RMSD_backbone_peptide_sample - RMSD_backbone_peptide_ref,
    perc_diff_RMSD_backbone_peptide = 100 * (RMSD_backbone_peptide_sample - RMSD_backbone_peptide_ref) / RMSD_backbone_peptide_ref,
    
    diff_RMSD_all_atom_peptide = RMSD_all_atom_peptide_sample - RMSD_all_atom_peptide_ref,
    perc_diff_RMSD_all_atom_peptide = 100 * (RMSD_all_atom_peptide_sample - RMSD_all_atom_peptide_ref) / RMSD_all_atom_peptide_ref,
    
    diff_RMSD_domain = RMSD_domain_sample - RMSD_domain_ref,
    perc_diff_RMSD_domain = 100 * (RMSD_domain_sample - RMSD_domain_ref) / RMSD_domain_ref,
    
    diff_DockQ = DockQ_sample - DockQ_ref,
    perc_diff_DockQ = 100 * (DockQ_sample - DockQ_ref) / DockQ_ref,
    
    diff_iPAE = iPAE_sample - iPAE_ref,
    perc_diff_iPAE = 100 * (iPAE_sample - iPAE_ref) / iPAE_ref,
    
    diff_pDockQ = pDockQ_sample - pDockQ_ref,
    perc_diff_pDockQ = 100 * (pDockQ_sample - pDockQ_ref) / pDockQ_ref
    
  )



# --- Define manual selection list ---
# I had to do this for the groups where MFL_DFL didnt exist, so that I could take the longest motif
# or domain extension possible. 

manual_dfl_list <- c(
  'LIG_TRFH_1_M10_M451_D1_D437',
  'LIG_WW_1_M865_M895_D2925_D3362',
  'LIG_NRBOX_M453_M826_DFL',
  'TRG_AP2beta_CARGO_1_MFL_D4_D937',
  'LIG_Vh1_VBS_1_M494_M740_DFL',
  'LIG_CAP-Gly_1_M1406_M1438_DFL',
  'DOC_MAPK_JIP1_4_MFL_D1_D464',
  'DOC_USP7_UBL2_3_M922_M1603_DFL',
  'LIG_PAM2_2_M1248_M1522_DFL',
  'DOC_SPAK_OSR1_1_MFL_D7_D527',
  'LIG_OCRL_FandH_1_M12_M249_DFL',
  'LIG_NBox_RRM_1_MFL_D1_D559',
  'LIG_REV1ctd_RIR_1_M401_M870_DFL',
  'LIG_RPA_C_Vert_M1_M309_D1_D270',
  'LIG_PDZ_Class_1_MFL_D1185_D1492',
  'LIG_PAM2_1_M1_M127_DFL',
  'LIG_CNOT1_NIM_1_M1_M276_D1377_D2376'
)


# --- Select representatives per group following MFL_DFL if exists > or manual list ---
average_dfl_reps <- merged_data_avg %>%
  group_by(prediction_group) %>%
  filter(
    if (any(str_detect(prediction_name, "MFL_DFL"))) {
      str_detect(prediction_name, "MFL_DFL")
    } else {
      prediction_name %in% manual_dfl_list
    }
  ) %>%
  slice(1) %>%
  ungroup()

# --- Prepare for plotting ---

# After you create average_dfl_reps
# I discarded these 2 structures because they had poor accurracy AF3 predictions, and 
# I didnt want to include them in the comparison since they wouldnt be meaningful

average_dfl_reps <- average_dfl_reps %>%
  filter(prediction_name != "DEG_Kelch_Keap1_1_MFL_DFL",
         prediction_name != "LIG_LIR_Gen_1_MFL_DFL" )


plot_data <- average_dfl_reps %>%
  select(prediction_group,
         prediction_name,
         RMSD_backbone_peptide_ref, RMSD_backbone_peptide_sample,
         RMSD_all_atom_peptide_ref, RMSD_all_atom_peptide_sample,
         RMSD_domain_ref, RMSD_domain_sample,
         diff_RMSD_backbone_peptide,
         diff_RMSD_all_atom_peptide,
         diff_RMSD_domain,
         diff_DockQ,
         diff_iPAE,
         diff_pDockQ) %>%
  pivot_longer(
    cols = starts_with("diff_"),
    names_to = "metric",
    values_to = "difference"
  )

# Add columns to identify reference and sample RMSD values for zeroing
plot_data <- plot_data %>%
  mutate(
    ref_value = case_when(
      metric == "diff_RMSD_backbone_peptide" ~ RMSD_backbone_peptide_ref,
      metric == "diff_RMSD_all_atom_peptide" ~ RMSD_all_atom_peptide_ref,
      metric == "diff_RMSD_domain" ~ RMSD_domain_ref,
      TRUE ~ NA_real_
    ),
    sample_value = case_when(
      metric == "diff_RMSD_backbone_peptide" ~ RMSD_backbone_peptide_sample,
      metric == "diff_RMSD_all_atom_peptide" ~ RMSD_all_atom_peptide_sample,
      metric == "diff_RMSD_domain" ~ RMSD_domain_sample,
      TRUE ~ NA_real_
    ),
    zeroed = !is.na(ref_value) & !is.na(sample_value) &
      ref_value > 15 & sample_value > 15,
    difference = ifelse(zeroed, 0, difference)
  )

# Factor levels and labels for metrics
plot_data$metric <- factor(plot_data$metric,
                           levels = c("diff_RMSD_backbone_peptide", "diff_RMSD_all_atom_peptide",
                                      "diff_RMSD_domain", "diff_DockQ", "diff_iPAE", "diff_pDockQ"),
                           labels = c("RMSD Backbone Peptide", "RMSD All Atom Peptide",
                                      "RMSD Domain", "DockQ", "iPAE", "pDockQ"))



# --- Plot ---
ggplot(plot_data, aes(x = prediction_group, y = difference,
                      fill = ifelse(zeroed, "Zeroed", "Normal"))) +
  geom_col(show.legend = TRUE) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c("Zeroed" = "gray70", "Normal" = "steelblue")) +
  labs(title = "Average Table Sample - Reference Differences (DFL/MFL Representatives)\nGray = Both AF2 & AF3 > 15 Å, difference set to 0",
       x = "Prediction Group", y = "Difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# TABLE FOR RANKED_0 SAMPLE - REF DIFFERENCE

# --- Define manual selection list ---
manual_dfl_list <- c(
  'LIG_TRFH_1_M10_M451_D1_D437',
  'LIG_PCNA_PIPBox_1_M11_M164_D1_D261',
  'LIG_WW_1_M865_M895_D2925_D3362',
  'LIG_NRBOX_M453_M826_DFL',
  'TRG_AP2beta_CARGO_1_MFL_D4_D937',
  'LIG_Vh1_VBS_1_M494_M740_DFL',
  'LIG_CAP-Gly_1_M1406_M1438_DFL',
  'DOC_MAPK_JIP1_4_MFL_D1_D464',
  'DOC_USP7_UBL2_3_M922_M1603_DFL',
  'LIG_PAM2_2_M1248_M1522_DFL',
  'DOC_SPAK_OSR1_1_MFL_D7_D527',
  'LIG_OCRL_FandH_1_M12_M249_DFL',
  'LIG_NBox_RRM_1_MFL_D1_D559',
  'LIG_REV1ctd_RIR_1_M401_M870_DFL',
  'LIG_RPA_C_Vert_M1_M309_D1_D270',
  'LIG_PDZ_Class_1_MFL_D1185_D1492',
  'LIG_PAM2_1_M1_M127_DFL',
  'LIG_CNOT1_NIM_1_M1_M276_D1377_D2376'
)

# --- Select representatives per group ---
dfl_representatives <- merged_data %>%
  group_by(prediction_group) %>%
  filter(
    # Step 1: Prefer MFL_DFL in prediction_name
    if (any(str_detect(prediction_name, "MFL_DFL"))) {
      str_detect(prediction_name, "MFL_DFL")
    } else {
      # Step 3: Multiple DFLs → keep only those in manual list
      prediction_name %in% manual_dfl_list
    }
  ) %>%
  slice(1) %>%
  ungroup()



# --- Prepare for plotting ---

# After you create average_dfl_reps
dfl_representatives <- dfl_representatives %>%
  filter(prediction_name != "DEG_Kelch_Keap1_1_MFL_DFL",
         prediction_name != "LIG_LIR_Gen_1_MFL_DFL" )


plot_data <- dfl_representatives %>%
  select(prediction_group,
         RMSD_backbone_peptide_ref, RMSD_backbone_peptide_sample,
         RMSD_all_atom_peptide_ref, RMSD_all_atom_peptide_sample,
         RMSD_domain_ref, RMSD_domain_sample,
         diff_RMSD_backbone_peptide,
         diff_RMSD_all_atom_peptide,
         diff_RMSD_domain,
         diff_DockQ,
         diff_iPAE,
         diff_pDockQ) %>%
  pivot_longer(
    cols = starts_with("diff_"),
    names_to = "metric",
    values_to = "difference"
  ) %>%
  mutate(
    # Add reference and sample values for RMSD metrics for zeroing check
    ref_value = case_when(
      metric == "diff_RMSD_backbone_peptide" ~ RMSD_backbone_peptide_ref,
      metric == "diff_RMSD_all_atom_peptide" ~ RMSD_all_atom_peptide_ref,
      metric == "diff_RMSD_domain" ~ RMSD_domain_ref,
      TRUE ~ NA_real_
    ),
    sample_value = case_when(
      metric == "diff_RMSD_backbone_peptide" ~ RMSD_backbone_peptide_sample,
      metric == "diff_RMSD_all_atom_peptide" ~ RMSD_all_atom_peptide_sample,
      metric == "diff_RMSD_domain" ~ RMSD_domain_sample,
      TRUE ~ NA_real_
    ),
    zeroed = !is.na(ref_value) & !is.na(sample_value) &
      ref_value > 15 & sample_value > 15,
    difference = ifelse(zeroed, 0, difference)
  )

plot_data$metric <- factor(plot_data$metric,
                           levels = c("diff_RMSD_backbone_peptide", "diff_RMSD_all_atom_peptide",
                                      "diff_RMSD_domain", "diff_DockQ", "diff_iPAE", "diff_pDockQ"),
                           labels = c("RMSD Backbone Peptide", "RMSD All Atom Peptide",
                                      "RMSD Domain", "DockQ", "iPAE", "pDockQ"))


# --- Plot ---
ggplot(plot_data, aes(x = prediction_group, y = difference,
                      fill = ifelse(zeroed, "Zeroed", "Normal"))) +
  geom_col(show.legend = TRUE) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c("Zeroed" = "gray70", "Normal" = "steelblue")) +
  labs(title = "Ranked_0 Sample - Reference Differences (MFL_DFL Representatives)\nGray = Both AF2 & AF3 > 15 Å, difference set to 0",
       x = "Prediction Group", y = "Difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))









# HEATMAP FOR AF3 RESULTS WITH AVERAGE

# Get the miminal structures first
af3_data <- read.delim("/Users/imb/Desktop/AF3_metrics.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# If you need both conditions together
minimal_data <- af3_data[af3_data$benchmark_set == "known_DMI" & af3_data$model_id == "ranked_0", ]

# Keep only the specified columns
minimal_data <- minimal_data[, c("prediction_name", "benchmark_set", "RMSD_domain", "RMSD_backbone_peptide", "RMSD_all_atom_peptide", "RMSD_all_atom", "ELM_instance")]

# Add prediction group 
minimal_data <- minimal_data %>%
  mutate(prediction_group = str_extract(prediction_name, "^[^M]+"))

# Change the column name of minimal_data to prepare it to merge with merged_data_avg
cols_to_rename <- c(
  "RMSD_domain",
  "RMSD_backbone_peptide",
  "RMSD_all_atom_peptide",
  "RMSD_all_atom"
)

# Append "_minimal" to their names
minimal_data <- minimal_data %>%
  rename_with(~ paste0(.x, "_minimal"), all_of(cols_to_rename))

# Join these minimal columns to the merged_data_avg table
merged_data_avg <- merged_data_avg %>%
  left_join(minimal_data, by = "ELM_instance")


# Calculate log2 fold change (min_frag / sample_avg)
merged_data_avg <- merged_data_avg %>%
  mutate(log2_fold_change = log2(RMSD_all_atom_peptide_minimal / RMSD_all_atom_peptide_sample))

# Aggregate mean log2 fold change by extension coordinates
heatmap_data <- merged_data_avg %>%
  group_by(domain_ext_index, motif_ext_index) %>%
  summarise(mean_log2_fold_change = mean(log2_fold_change, na.rm = TRUE)) %>%
  ungroup()

# Optional: Complete the grid if some coordinates are missing
heatmap_data_complete <- heatmap_data %>%
  complete(domain_ext_index = 0:3,
           motif_ext_index = 0:6,
           fill = list(mean_log2_fold_change = NA))

# Plot heatmap with red-blue diverging palette
ggplot(heatmap_data_complete, aes(x = as.factor(domain_ext_index), 
                                  y = as.factor(motif_ext_index), 
                                  fill = mean_log2_fold_change)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_log2_fold_change, 2)), size = 3, na.rm = TRUE) +
  scale_fill_distiller(palette = "RdBu", direction = -1,  # red to blue, direction 1
                       limits = c(-2.5, 2.5),
                       oob = scales::squish,
                       name = expression(log[2](RMSD[min] / RMSD[ext]))) +
  labs(
    title = "AF3 Average — log2 Fold Change (RMSD)",
    x = "Domain Extension Step",
    y = "Motif Extension Step"
  ) +
  scale_y_discrete(limits = rev(as.character(0:6))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )



# HEATMAP FOR AF2 RESULTS (We only had best ranked ones, so there is no plot for the average)

# Group and calculate mean log2 fold change
ref_heatmap_data <- reference %>%
  group_by(domain_ext_index, motif_ext_index) %>%
  summarise(mean_rmsd_fc = mean(RMSD_all_atom_peptide_fold_change, na.rm = TRUE)) %>%
  ungroup()

# (Optional) Complete grid for missing combinations
ref_heatmap_data_complete <- ref_heatmap_data %>%
  complete(domain_ext_index = 0:3,
           motif_ext_index = 0:6,
           fill = list(mean_rmsd_fc = NA))

# Plot
ggplot(ref_heatmap_data_complete, aes(x = as.factor(domain_ext_index), 
                                      y = as.factor(motif_ext_index), 
                                      fill = mean_rmsd_fc)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_rmsd_fc, 2)), size = 3, na.rm = TRUE) +
  scale_fill_distiller(palette = "RdBu", direction = -1,  # red to blue, direction 1
                       limits = c(-2.5, 2.5),
                       oob = scales::squish,
                       name = expression(log[2](RMSD[min] / RMSD[ext]))) +
  labs(
    title = "AF2 Best Ranked — log2 Fold Change (RMSD)",
    x = "Domain Extension Step",
    y = "Motif Extension Step"
  ) +
  scale_y_discrete(limits = rev(as.character(0:6))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

                                  
                                  
                                  
                                  
# HEATMAP FOR AF3 RESULTS WITH BEST RANKED 

# Join the minimal columns calculated above to the merged_data table
merged_data <- merged_data %>%
  left_join(minimal_data, by = "ELM_instance")

# Calculate log2 fold change (min_frag / sample_avg)
merged_data <- merged_data %>%
  mutate(log2_fold_change = log2(RMSD_all_atom_peptide_minimal / RMSD_all_atom_peptide_sample))

# Aggregate mean log2 fold change by extension coordinates
heatmap_data <- merged_data %>%
  group_by(domain_ext_index, motif_ext_index) %>%
  summarise(mean_log2_fold_change = mean(log2_fold_change, na.rm = TRUE)) %>%
  ungroup()


# Complete grid to show all combinations (fill missing with NA)
heatmap_data_complete <- heatmap_data %>%
  complete(domain_ext_index = 0:3,
           motif_ext_index = 0:6,
           fill = list(mean_log2_fold_change = NA))

# Plot heatmap
ggplot(heatmap_data_complete, aes(x = as.factor(domain_ext_index),
                                  y = as.factor(motif_ext_index),
                                  fill = mean_log2_fold_change)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_log2_fold_change, 2)), size = 3, na.rm = TRUE) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       limits = c(-2.5, 2.5),
                       oob = scales::squish,
                       name = expression(log[2](RMSD[min] / RMSD[ext]))) +
  labs(
    title = "AF3 Top-Ranked Samples — log2 Fold Change (RMSD)",
    x = "Domain Extension Step",
    y = "Motif Extension Step"
  ) +
  scale_y_discrete(limits = rev(as.character(0:6))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

