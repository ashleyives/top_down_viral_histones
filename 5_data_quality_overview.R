
library(MSnSet.utils)
library(tidyverse)
library(proteomicsCV) 

load("normalized_td_msnset.RData")

dim(m) #number of proteoform clusters

masses <- fData(m) %>%
  filter(grepl("Histone", fData(m)$'Protein description')) %>%
  group_by(Gene) %>%
  mutate(mass_round = round(`Proteoform mass`)) %>% # Round 'Proteoform mass' to the nearest whole integer
  distinct(mass_round) %>%
  tally() %>%
  dplyr::rename(masses = n) 

pfrs <- fData(m) %>%
  filter(grepl("Histone", fData(m)$'Protein description')) %>%
  group_by(Gene) %>%
  tally() %>% #it is counting the number of rows per gene 
  dplyr::rename(pfrs = n) 

mf <- filter_by_occurrence(m, 0.5) 

masses50 <- fData(mf) %>%
  filter(grepl("Histone", fData(mf)$'Protein description')) %>%
  group_by(Gene) %>%
  mutate(mass_round = round(`Proteoform mass`)) %>% # Round 'Proteoform mass' to the nearest whole integer
  distinct(mass_round) %>%
  tally() %>%
  dplyr::rename(masses = n) 

pfrs50 <- fData(mf) %>%
  filter(grepl("Histone", fData(mf)$'Protein description')) %>%
  group_by(Gene) %>%
  tally() %>% #it is counting the number of rows per gene 
  dplyr::rename(pfrs = n) 

#number of masses and proteoforms per gene plot---------------------------------
panela <- masses %>%
  left_join(pfrs) %>%
  pivot_longer(
    cols = -Gene,
    names_to = "data",
    values_to = "value"
  ) %>%
  ggplot(aes(x = reorder(Gene, -value), y = value, fill = data)) + 
  geom_bar(
    stat = "identity", 
    position = position_dodge(width = 0.8) 
  ) +
  labs(
    x = "Gene", 
    y = "Total Observations (n)", 
    fill = "Observation" 
  ) +
  scale_fill_manual(
    values = c("gray", "black"), 
    labels = c("Mass", "Proteoform") 
  ) +
  theme_bw(base_size = 16) +
  scale_y_continuous(
    expand = c(0, 0), 
    limits = c(0, 200)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.8,0.7)
  ) +
  geom_point(data = masses50 %>%
               left_join(pfrs50) %>%
               pivot_longer(
                 cols = -Gene,
                 names_to = "data",
                 values_to = "value"
               ),
             aes(x = reorder(Gene, -value), y = value, color = data), 
             position = position_dodge(width = 0.8), # Align with bars
             size = 2,
             show.legend = FALSE) + # Adjust size of points
  scale_color_manual(
    values = c("red", "blue")
  )
panela

#data completeness plot---------------------------------------------------------
panelb <- fData(m) %>%
  mutate(Completeness = (non_na_229E + non_na_Mock) / (ncol(m))) %>%
  mutate(Completeness_bin = case_when(
    Completeness <= 1 & Completeness > 0.80 ~ "100-80%",
    Completeness <= 0.80 & Completeness > 0.6 ~ "80-60%",
    Completeness <= 0.6 & Completeness > 0.4 ~ "60-40%",
    Completeness <= 0.4 & Completeness > 0.2 ~ "40-20%",
    Completeness <= 0.2 & Completeness >= 0 ~ "20-0%"
  )) %>%
  mutate(Completeness_bin = factor(Completeness_bin, levels = c("100-80%", "80-60%", "60-40%", "40-20%", "20-0%"))) %>% # Set factor levels explicitly
  group_by(Completeness_bin) %>%
  tally(name = "n") %>% 
  ggplot(aes(x = Completeness_bin, y = n)) +
  geom_bar(stat = "identity", fill = "black") +
  labs(x = "Data Completeness", y = "Proteoforms (n)") +
  theme_bw(base_size = 16)+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 250)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
    # legend.key.size = unit(0.5, "cm"),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 12),
    # plot.margin = unit(c(1, 1, 3.5, 1), "lines")
  )
panelb



#CV plots-----------------------------------------------------------------------

mvirus <- m[,grepl("229E",sampleNames(m))] 

mmock <-  m[,grepl("Mock",sampleNames(m))] 

raw_cvs1<-data.frame(cv = protLogCV(data.frame(exprs(mvirus)), 2)) %>%
  mutate(Condition = "Virus")

raw_cvs2<-data.frame(cv = protLogCV(data.frame(exprs(mmock)), 2)) %>%
  mutate(Condition = "Mock")

raw_cvs <- rbind(raw_cvs1, raw_cvs2)

cv <- raw_cvs %>%
  filter(!is.na(cv)) %>%
  ggplot(aes(x = Condition, y = cv, fill = Condition)) +
  geom_violin(alpha = 0.5) +
scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) +
  scale_fill_manual(values = c("gray", "#56b4e9")) +  # Specify fill colors
  # Add mean (red dot)
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 21, 
    size = 3, 
    color = "red", 
    fill = "red"
  ) +
  # Add median (blue diamond)
  stat_summary(
    fun = median, 
    geom = "point", 
    shape = 23, 
    size = 3, 
    color = "blue", 
    fill = "blue"
  ) + 
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  labs(y = "RSD (%)")
cv
