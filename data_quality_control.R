library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

merged_data <- merged_data %>%
  mutate(prediction_group = str_extract(prediction_name, "^[^M]+"))

# Load both files
samples <- read_excel("/Users/imb/Desktop/AF3_hydrogens_metrics.xlsx")
reference <- read_excel("/Users/imb/Desktop/reference_metrics.xlsx", sheet = "AF-MMv2.2 result")

# Add suffixes to the same column names
# Merge regarding both the prection name and model_id
merged_data <- merge(samples, reference,
                     by = c("prediction_name", "model_id"),
                     suffixes = c("_sample", "_ref"))

# Check the intersection of columns to be sure
intersect(names(samples), names(reference))

# Check the names in the merged data to see if it has both ref and sample
names(merged_data)

nrow(merged_data)  # how many matched samples you now have

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

# Add numbers as prediction names
merged_data$PredictionID <- seq_len(nrow(merged_data))

# Exporting the ID mapping if needed
write.csv(merged_data[, c("PredictionID", "prediction_name")], "prediction_id_mapping.csv", row.names = FALSE)


# Bar plots of percentage differences for all metrics

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

# Bar plots of absolute differences

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

# Scatter plots 

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


extreme_rmsd <- merged_data[merged_data$diff_RMSD_backbone_peptide <= -100,
                            c("prediction_name", "RMSD_backbone_peptide_sample", "RMSD_backbone_peptide_ref")]
print(extreme_rmsd)
View(extreme_rmsd)


unmatched_refs <- anti_join(reference, samples, by = c("prediction_name", "model_id"))
unmatched_prediction_names <- unmatched_refs$prediction_name
print(unmatched_prediction_names)


#Now calculate the difference between extension groups

# Earliest and latest domain extension versions per group
domain_earliest <- merged_data %>%
  group_by(prediction_group) %>%
  filter(domain_ext_index == min(domain_ext_index, na.rm = TRUE)) %>%
  slice(1) %>%   # if multiple, pick first
  ungroup()

domain_latest <- merged_data %>%
  group_by(prediction_group) %>%
  filter(domain_ext_index == max(domain_ext_index, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()

# Join domain extension comparison
domain_extension_comparison <- domain_earliest %>%
  select(prediction_group,
         prediction_name_earliest = prediction_name,
         model_id_earliest = model_id,
         domain_ext_earliest = domain_ext_index,
         DockQ_earliest = DockQ_sample,
         RMSD_earliest = RMSD_backbone_peptide_sample,
         pDockQ_earliest = pDockQ_sample,
         iPAE_earliest = iPAE_sample) %>%
  left_join(
    domain_latest %>%
      select(prediction_group,
             prediction_name_latest = prediction_name,
             model_id_latest = model_id,
             domain_ext_latest = domain_ext_index,
             DockQ_latest = DockQ_sample,
             RMSD_latest = RMSD_backbone_peptide_sample,
             pDockQ_latest = pDockQ_sample,
             iPAE_latest = iPAE_sample),
    by = "prediction_group"
  ) %>%
  mutate(
    DockQ_change = DockQ_latest - DockQ_earliest,
    RMSD_change = RMSD_latest - RMSD_earliest,
    pDockQ_change = pDockQ_latest - pDockQ_earliest,
    iPAE_change = iPAE_latest - iPAE_earliest
  )

# Similarly for motif extension

motif_earliest <- merged_data %>%
  group_by(prediction_group) %>%
  filter(motif_ext_index == min(motif_ext_index, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()

motif_latest <- merged_data %>%
  group_by(prediction_group) %>%
  filter(motif_ext_index == max(motif_ext_index, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()

motif_extension_comparison <- motif_earliest %>%
  select(prediction_group,
         prediction_name_earliest = prediction_name,
         model_id_earliest = model_id,
         motif_ext_earliest = motif_ext_index,
         DockQ_earliest = DockQ_sample,
         RMSD_earliest = RMSD_backbone_peptide_sample,
         pDockQ_earliest = pDockQ_sample,
         iPAE_earliest = iPAE_sample) %>%
  left_join(
    motif_latest %>%
      select(prediction_group,
             prediction_name_latest = prediction_name,
             model_id_latest = model_id,
             motif_ext_latest = motif_ext_index,
             DockQ_latest = DockQ_sample,
             RMSD_latest = RMSD_backbone_peptide_sample,
             pDockQ_latest = pDockQ_sample,
             iPAE_latest = iPAE_sample),
    by = "prediction_group"
  ) %>%
  mutate(
    DockQ_change = DockQ_latest - DockQ_earliest,
    RMSD_change = RMSD_latest - RMSD_earliest,
    pDockQ_change = pDockQ_latest - pDockQ_earliest,
    iPAE_change = iPAE_latest - iPAE_earliest
  )

View(domain_extension_comparison)
View(motif_extension_comparison)









# Scatter plot
ggplot(merged_data, aes(x = domain_ext_index, y = DockQ_sample)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, color = "#2c7fb8") +  # jitter to avoid overplotting
  geom_smooth(method = "loess", color = "#de2d26", se = TRUE) +          # smooth trend line with confidence interval
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

ggplot(merged_data, aes(x = domain_ext_index, y = DockQ_sample)) +
  geom_point(color = "#2c7fb8") +
  geom_smooth(method = "loess", se = FALSE, color = "#de2d26") +
  facet_wrap(~ prediction_group) +
  labs(title = "DockQ Accuracy vs Domain Extension by Prediction Group",
       x = "Domain Extension Index",
       y = "DockQ (Sample)") +
  theme_minimal()

ggplot(merged_data, aes(x = motif_ext_index, y = DockQ_sample)) +
  geom_point(color = "#2c7fb8", alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, color = "#de2d26") +
  facet_wrap(~ prediction_group) +
  labs(title = "DockQ Accuracy vs Motif Extension by Prediction Group",
       x = "Motif Extension Index",
       y = "DockQ (Sample)") +
  theme_minimal()



# PRACTICAL FUNCTION TO CALL FOR FACETED SCATTER PLOT
plot_extension_vs_metric <- function(data, metric, extension_col, title_prefix = "") {
  ggplot(data, aes_string(x = extension_col, y = metric)) +
    geom_point(color = "#2c7fb8", alpha = 0.7) +
    geom_smooth(method = "loess", se = FALSE, color = "#de2d26") +
    facet_wrap(~ prediction_group) +
    labs(title = paste(title_prefix, metric, "vs", extension_col),
         x = extension_col,
         y = metric) +
    theme_minimal()
}


# RMSD backbone peptide vs domain extension
plot_extension_vs_metric(merged_data, "RMSD_backbone_peptide_sample", "domain_ext_index", "Domain Extension: ")

# RMSD backbone peptide vs motif extension
plot_extension_vs_metric(merged_data, "RMSD_backbone_peptide_sample", "motif_ext_index", "Motif Extension: ")

# RMSD all atom peptide vs domain extension
plot_extension_vs_metric(merged_data, "RMSD_all_atom_peptide_sample", "domain_ext_index", "Domain Extension: ")

# RMSD all atom peptide vs motif extension
plot_extension_vs_metric(merged_data, "RMSD_all_atom_peptide_sample", "motif_ext_index", "Motif Extension: ")

# RMSD domain vs motif extension
plot_extension_vs_metric(merged_data, "RMSD_domain_sample", "motif_ext_index", "Motif Extension: ")

# RMSD domain vs domain extension
plot_extension_vs_metric(merged_data, "RMSD_domain_sample", "domain_ext_index", "Domain Extension: ")


# extreme cases vol2
extreme_rmsd_failures <- merged_data %>%
  filter(RMSD_backbone_peptide_sample > 50,
         RMSD_backbone_peptide_ref < 5) %>%
  select(prediction_name,
         RMSD_backbone_peptide_sample,
         RMSD_backbone_peptide_ref,
         diff_RMSD_backbone_peptide,
         domain_ext_index,
         motif_ext_index)

# View it
View(extreme_rmsd_failures)



# HEATMAP
# Install necessary packages
install.packages(c("ggplot2", "dplyr", "tidyr", "RColorBrewer"))

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# STEP 1: Summarise from merged_data
heatmap_data <- merged_data %>%
  group_by(domain_ext_index, motif_ext_index) %>%
  summarise(mean_rmsd_fc = mean(RMSD_all_atom_peptide_fold_change, na.rm = TRUE)) %>%
  ungroup()

# STEP 2: Ensure full grid is represented
heatmap_data_complete <- heatmap_data %>%
  complete(domain_ext_index = 0:3,
           motif_ext_index = 0:6,
           fill = list(mean_rmsd_fc = NA))

# STEP 3: Plot the heatmap
ggplot(heatmap_data_complete, aes(x = as.factor(domain_ext_index), 
                                  y = as.factor(motif_ext_index), 
                                  fill = mean_rmsd_fc)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_rmsd_fc, 2)), size = 3, na.rm = TRUE) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       limits = c(-2.5, 0.5),
                       oob = scales::squish,
                       name = expression(log[2](RMSD[min]/RMSD[ext]))) +
  labs(
    title = "AF3 known extension",
    x = "Domain extension step",
    y = "Motif extension step"
  ) +
  scale_y_discrete(limits = rev(as.character(0:6))) +  # Y axis from 6 to 0
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )


#HEATMAP 2 WITH HIGHEST RANK SAMPLES


# STEP 1: Filter only top-ranked model per prediction
top_models <- merged_data %>%
  group_by(domain_ext_index, motif_ext_index) %>%
  slice_max(order_by = ranking_score, n = 1, with_ties = FALSE) %>%
  ungroup()

# STEP 2: Ensure full grid is represented
heatmap_data_complete <- top_models %>%
  complete(domain_ext_index = 0:3,
           motif_ext_index = 0:6,
           fill = list(RMSD_all_atom_peptide_fold_change = NA))

# STEP 3: Plot the heatmap directly from top model data
ggplot(heatmap_data_complete, aes(x = as.factor(domain_ext_index), 
                                  y = as.factor(motif_ext_index), 
                                  fill = RMSD_all_atom_peptide_fold_change)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(RMSD_all_atom_peptide_fold_change, 2)), size = 3, na.rm = TRUE) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       limits = c(-2.5, 0.5),
                       oob = scales::squish,
                       name = expression(log[2](RMSD[min]/RMSD[ext]))) +
  labs(
    title = "AF3 known extension",
    x = "Domain extension step",
    y = "Motif extension step"
  ) +
  scale_y_discrete(limits = rev(as.character(0:6))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

# MIRRORED BAR PLOT FOR AF2 AND AF3
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# 0. Add `base_name` column by splitting before 'M'
merged_data <- merged_data %>%
  mutate(base_name = str_extract(prediction_name, "^[^M]+"))

# 1. Filter only _DFL rows
dfl_data <- merged_data %>%
  filter(str_detect(prediction_name, "_DFL$"))

# 2. Take one sample per extension group
representatives <- dfl_data %>%
  group_by(base_name, domain_ext_index, motif_ext_index) %>%
  slice(1) %>%  # take the first one
  ungroup()

# 3. Create a long format table with AF2 and AF3 values
long_data <- representatives %>%
  mutate(PredictionID = paste(base_name, domain_ext_index, motif_ext_index, sep = "_")) %>%
  select(PredictionID, base_name, domain_ext_index, motif_ext_index,
         AF2 = RMSD_all_atom_peptide_ref,
         AF3 = RMSD_all_atom_peptide_sample) %>%
  pivot_longer(cols = c("AF2", "AF3"), names_to = "Method", values_to = "RMSD") %>%
  mutate(RMSD = ifelse(Method == "AF3", -RMSD, RMSD))  # Mirror AF3 bars downward

# 4. Plot mirrored bar plot
ggplot(long_data, aes(x = PredictionID, y = RMSD, fill = Method)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = c("AF2" = "#1f78b4", "AF3" = "#e31a1c")) +
  labs(title = "Mirrored RMSD All Atom Peptide: AF2 vs AF3",
       x = "Prediction Extension Group",
       y = "RMSD (Å)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    plot.title = element_text(size = 14, face = "bold")
  )

representatives %>%
  select(prediction_name, base_name, domain_ext_index, motif_ext_index)

# RANDOM SEEDS BARPLOT MIRROR
library(dplyr)
library(ggplot2)

# Step 1: Extract base_name from prediction_name
merged_data <- merged_data %>%
  mutate(base_name = sub("M.*", "", prediction_name))

# Step 2: Randomly select one prediction_name per base_name
set.seed(42)  # for reproducibility
representative_names <- merged_data %>%
  group_by(base_name) %>%
  slice_sample(n = 1) %>%
  pull(prediction_name)

# Step 3: Filter merged_data to keep only selected representative predictions
representatives <- merged_data %>%
  filter(prediction_name %in% representative_names)

# Step 4: Prepare tidy data for plotting
representatives_long <- representatives %>%
  mutate(
    RMSD_AF2 = RMSD_all_atom_peptide_ref,
    RMSD_AF3 = RMSD_all_atom_peptide_sample
  ) %>%
  select(base_name, RMSD_AF2, RMSD_AF3) %>%
  pivot_longer(cols = c(RMSD_AF2, RMSD_AF3), names_to = "method", values_to = "RMSD") %>%
  mutate(
    method = ifelse(method == "RMSD_AF2", "AF2", "AF3"),
    RMSD_signed = ifelse(method == "AF3", -RMSD, RMSD)
  )

# Step 5: Plot mirror bar chart
ggplot(representatives_long, aes(x = reorder(base_name, RMSD_signed), y = RMSD_signed, fill = method)) +
  geom_col(position = "identity") +
  scale_y_continuous(labels = abs) +
  labs(title = "AF2 vs AF3 - RMSD All Atom Peptide (Mirror Plot)",
       x = "Prediction (Base Name)",
       y = "RMSD (Å)",
       fill = "Method") +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6),
    plot.title = element_text(size = 14, face = "bold")
  )
# Show the selected representative prediction names and their base names
representatives %>%
  select(base_name, prediction_name) %>%
  arrange(base_name)


