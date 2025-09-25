
library(tidyverse)
library(MSnbase)
library(PNNL.DMS.utils)
library(MSnSet.utils)
library(ggpubr)
library(ggrepel)

#this script is only for plotting and does not perform any analysis-------------

load("res.RData") #load DEA data made in 2_differential_abundance.R

#Volcano plot of all quantifiable histone---------------------------------------
panela <- res %>%
  dplyr::rename(feature_name = PF) %>%
  left_join(fData(m)) %>%
  filter(!is.na(P.Value)) %>%
  ggplot(aes(x = logFC, y = -log10(P.Value))) +
  geom_point(size = 5, alpha=0.5) +
  # xlim(-2.1, 2.1)+
  geom_hline(yintercept = -log10(0.05), alpha=0.2) +
  geom_vline(xintercept = c(1, -1), alpha=0.2) +
  # labs(size = "Completeness") +
  geom_text_repel(
    aes(x = logFC, y = -log10(P.Value), label = ifelse(P.Value < 0.05, feature_name, NA)),
    inherit.aes = FALSE,
    force_pull = 0.01, #Lower values make labels less "attached" to their points.
    nudge_y = 0.5,
    # nudge_x = -1,
    na.rm = TRUE,
    segment.alpha = 0.2
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key.size = unit(0.5, "cm"),  # Shrinking the size of the legend keys
    legend.text = element_text(size = 10),  # Shrinking the size of the legend text
    legend.title = element_text(size = 12),   # Shrink the size of the legend title
    plot.margin = unit(c(1, 1, 3.5, 1), "lines") # Add margin space at the bottom
  ) +
  labs(x = expression(Log[2] * " fold change"), y = expression(-Log[10] * " (P-Value)"))+
  annotate("text", x = 2, y = 0.1, label = "Up in Virus", color = "black", size = 5, hjust = 1) +
  annotate("text", x = -3.5, y = 0.1, label = "Down in Virus", color = "black", size = 5, hjust = 0)
panela

#Volcano plot of by mod type for PTM histones-----------------------------------
panelptm <- res %>%
  filter(!grepl("RP", PF)) %>%
  dplyr::rename(feature_name = PF) %>%
  left_join(fData(m)) %>%
  filter(!is.na(P.Value)) %>%
  mutate(modtype = str_extract(MIScore, "^[^\\[]+")) %>%
  mutate(modtype = ifelse(modtype %in% c("Acetyl", "Trimethyl"), "Acetyl_Trimethyl", modtype)) %>%
  filter(modtype %in% c("Methyl", "Dimethyl", "Acetyl_Trimethyl", "Phosph")) %>%
  mutate(modtype = factor(modtype, levels = c("Methyl", "Dimethyl", "Acetyl_Trimethyl", "Phosph"))) %>%  # Convert `modtype` to a factor with specified order
  mutate(Gene = str_extract(feature_name, "^[^_]+")) %>%
  ggplot(aes(x = logFC, y = -log10(P.Value), color = modtype)) +
  geom_point(size = 5, alpha=0.5) +
  geom_vline(xintercept = c(-1,1), alpha=0.2) +
  # geom_vline(xintercept = c(0), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), alpha=0.2) +
  labs(size = "Completeness") +
  # geom_text_repel(
  #   aes(x = logFC, y = -log10(P.Value), label = feature_name),
  #   inherit.aes = FALSE,  # Ignore default aesthetics for this layer
  #   force_pull = 0.01, #Lower values make labels less "attached" to their points.
  #   seed = 123,  # Use a random seed for reproducibility
  #   max.overlaps = Inf,  # Allow labels to overlap more freely
  #   segment.color = "black",  # Color of the connecting line
  #   segment.size = 0.5,  # Thickness of the connecting line
  #   nudge_y = 0.03,  # Adjust vertical distance from the point
  #   nudge_x = 0.05,  # Adjust horizontal distance from the point
  #   na.rm = TRUE
# ) +
geom_text_repel(
  aes(x = logFC, y = -log10(P.Value), label = ifelse(P.Value < 0.05, feature_name, NA)),
  inherit.aes = FALSE,
  force_pull = 0.01, #Lower values make labels less "attached" to their points.
  nudge_y = 0.5,
  # nudge_x = -1,
  na.rm = TRUE,
  segment.alpha = 0.2
) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key.size = unit(0.5, "cm"),  # Shrinking the size of the legend keys
    legend.text = element_text(size = 10),  # Shrinking the size of the legend text
    legend.title = element_text(size = 12)  # Shrink the size of the legend title
    # plot.margin = unit(c(1, 1, 3.5, 1), "lines") # Add margin space at the bottom
  ) +
  labs(
    x = expression(Log[2] * " fold change"),
    y = expression(-Log[10] * " (P-Value)"),
    color = "Modification"
  ) +
  annotate("text", x = 0.5, y = 0.03, label = "Up in Virus", color = "black", size = 5, hjust = 1) +
  annotate("text", x = -0.75, y = 0.03, label = "Down in Virus", color = "black", size = 5, hjust = 0) +
  scale_color_manual(values = c("black", "#E69F00", "#56b4e9", "#009e73"))
panelptm

#Volcano plot of H2As colored by C-term truncation------------------------------
panelbcont <- res %>%
  dplyr::rename(feature_name = PF) %>%
  left_join(fData(m)) %>%
  filter(grepl("H2A", feature_name)) %>%
  mutate(Cdelta = protLength - lastAA) %>%
  mutate(deltal = if_else(Cdelta > 15, "Cleaved (+15 a.a.)", "Full length")) %>%
  ggplot(aes(x = logFC, y = -log10(P.Value), color = Cdelta)) + # Using Cdelta as a continuous variable
  geom_point(size = 5) +
  # xlim(-2.1, 2.1) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.2) +
  geom_vline(xintercept = c(1, -1), alpha = 0.2) +
  geom_text_repel(
    aes(x = logFC, y = -log10(P.Value), label = ifelse(P.Value < 0.05, feature_name, NA)),
    inherit.aes = FALSE,
    force_pull = 0.01, #Lower values make labels less "attached" to their points.
    nudge_y = 0.5,
    # nudge_x = -1,
    na.rm = TRUE,
    segment.alpha = 0.2
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
    # plot.margin = unit(c(1, 1, 3.5, 1), "lines")
  ) +
  labs(
    x = expression(Log[2] * " fold change"),
    y = expression(-Log[10] * " (P-Value)")
  ) +
  annotate("text", x = 2, y = 0.1, label = "Up in Virus", color = "black", size = 5, hjust = 1) +
  annotate("text", x = -2, y = 0.1, label = "Down in Virus", color = "black", size = 5, hjust = 0) +
  guides(
    color = guide_colorbar(title = "How many residues \n cleaved from C-term?")  # Add a continuous gradient legend for Cdelta
  ) +
  scale_color_gradient2(
    low = "gray",       # Color for low Cdelta
    mid = "blue",       # Midpoint color
    high = "red",       # Color for high Cdelta
    midpoint = 25       # Set the midpoint of Cdelta gradient
  ) # Continuous gradient with gray, blue (midpoint), and red
panelbcont

#Volcano plot of H3s colored by N-term truncation------------------------------
h3si <- res %>%
  dplyr::rename(feature_name = PF) %>%
  left_join(fData(m)) %>%
  filter(grepl("H3", feature_name)) %>%
  filter(firstAA > 2) %>% #remove canonical forms 
  ggplot(aes(x = logFC, y = -log10(P.Value), color = firstAA)) + # Using Cdelta as a continuous variable
  geom_point(size = 5) +
  scale_color_gradient2(
    low = "gray",       # Color for low Cdelta
    mid = "blue",       # Midpoint color
    high = "red",       # Color for high Cdelta
    midpoint = 35       # Set the midpoint of Cdelta gradient
  ) +# Continuous gradient with gray, blue (midpoint), and red
  geom_hline(yintercept = -log10(0.05), alpha = 0.2) +
  geom_vline(xintercept = c(1, -1), alpha = 0.2) +
  geom_text_repel(
    aes(x = logFC, y = -log10(P.Value), label = firstAA),
    inherit.aes = FALSE,
    force_pull = 0.1, #Lower values make labels less "attached" to their points.
    nudge_y = 0.05,
    nudge_x = 0.05,
    na.rm = TRUE,
    segment.alpha = 0.2
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
    # plot.margin = unit(c(1, 1, 3.5, 1), "lines")
  ) +
  labs(
    x = expression(Log[2] * " fold change"),
    y = expression(-Log[10] * " (P-Value)")
    # title = "Differential abundance of N-terminally truncated H3 proteoforms"
  ) +
  annotate("text", x = 0.9, y = 0.1, label = "Up in Virus", color = "black", size = 5, hjust = 1) +
  annotate("text", x = -0.9, y = 0.1, label = "Down in Virus", color = "black", size = 5, hjust = 0) +
  guides(
    color = guide_colorbar(title = "First amino \nacid position")  # Add a continuous gradient legend for Cdelta
  ) 
h3si



