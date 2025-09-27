
library(tidyverse)

source("parse_ miscore.R")

load("x_grp_td.RData")

x <- parse_miscore(x_grp)

x <- x %>%
  filter(mi_level == 1) %>%
  filter(mi_location != "-") %>%
  mutate(ID = paste0(Gene, "_", mi_mod, "_", mi_adj_location_list)) %>%
  distinct(ID, .keep_all = TRUE) %>%
  dplyr::select(Gene, 'Protein accession', mi_mod, mi_location, mi_location_list, mi_adj_location_list, mi_level) %>%
  separate(mi_adj_location_list, into = c("mi_residue", "mi_adj_position"), 
           sep = "(?<=\\D)(?=\\d)", convert = TRUE) %>% # Splits at the boundary between letters and numbers
  dplyr::rename(Accession = 'Protein accession') %>%
  mutate(site_id = paste0(Accession, "_", mi_adj_position)) %>%
  dplyr::select(-Gene)

#load in clustal alignments made in uniprot for each core histone 
h2a <- Biostrings::readAAStringSet(
  file.path("h2a_align.txt"))

h2b <- Biostrings::readAAStringSet(
  file.path("h2b_align.txt"))

h3 <- Biostrings::readAAStringSet(
  file.path("h3_align.txt"))

h4 <- Biostrings::readAAStringSet(
  file.path("h4_align.txt"))

wrangle_fasta <- function(x){
  
  y <- data.frame(
    Accession = names(x),                          # Extracts protein accession
    FullSequence = as.character(x),                              # Converts sequence data to character
    row.names = NULL,
    check.names = FALSE
  ) %>%
    mutate(Accession = str_extract(Accession, "^\\S+"))  # Extract first word if needed
  
  return(y)
}

h2a_df <- wrangle_fasta(h2a)
h2b_df <- wrangle_fasta(h2b)
h3_df <- wrangle_fasta(h3)
h4_df <- wrangle_fasta(h4)

all_histones <- bind_rows(h2a_df, h2b_df, h3_df, h4_df)

long_fasta <- all_histones %>%
  mutate(letter = strsplit(FullSequence, "")) %>%  # Split sequences into individual characters
  unnest(letter) %>%                              # Expand list-column into rows
  group_by(Accession) %>%               # Group by sequence ID
  mutate(alignment_position = row_number(),       # Add position numbers (1 to n for each sequence)
         sequence_position = cumsum(ifelse(letter != "-", 1, 0)) %>% # Increment position for non "-" letters
           ifelse(letter == "-", NA, .)) %>%     # Assign NA if letter is "-"                  
  ungroup()  %>%
  mutate(core_histone = case_when(
    grepl("H2A", Accession) ~ "H2A",
    grepl("H2B", Accession) ~ "H2B",
    grepl("H3", Accession) ~ "H3",
    grepl("H4", Accession) ~ "H4",
    TRUE ~ NA_character_ # Default value if no match is found
  ))%>%
  mutate(site_id = paste0(Accession, "_", sequence_position)) %>%
  dplyr::select(-Accession)

x_grp$Accession <- x_grp$`Protein accession`

genenames <- x_grp %>%
  dplyr::select(Accession, Gene) %>%
  separate(Accession, into = c("Accession"), sep="_", remove = TRUE) %>%
  distinct()

forplot <- long_fasta %>%
  left_join(x, by = join_by(site_id)) %>%
  separate(site_id, into = c("Accession"), sep="_", remove = FALSE) %>%
  left_join(genenames) %>%
  mutate(symbol = case_when(
    mi_mod == "CTrmAmid" ~ "CTrmAmid",       # Assign for Carbamyl
    mi_mod == "Deamide" ~ "Deamide",          # Assign for Deamide
    mi_mod == "Acetyl" ~ "Acetyl/Trimethyl",        # Assign for Acetyl
    mi_mod == "Carbamyl" ~ "Carbamyl",        # Assign for Carbamyl
    mi_mod == "Phosph" ~ "Phosph",           # Assign for Phospho
    mi_mod == "Trimethyl" ~ "Acetyl/Trimethyl",     # Assign for Trimethyl
    mi_mod == "IronAdduct" ~ "IronAdduct",      # Assign for IronAdduct
    mi_mod == "Dimethyl" ~ "Dimethyl",        # Assign for Dimethyl
    mi_mod == "Methyl" ~ "Methyl",           # Assign for Methyl
    mi_mod == "Plus1Oxy" ~ "Ox",         # Assign for Plus1Oxy
    mi_mod == "Plus2Oxy" ~ "Ox2",        # Assign for Plus2Oxy
    mi_mod == "Glutahion" ~ "Ox2",        # Assign for Plus2Oxy
    TRUE ~ NA_character_                # Default: NA if no match
  )) %>%
  mutate(letter = if_else(is.na(mi_mod), "-", letter)) %>%
  mutate(plot_id = paste(alignment_position, "_", core_histone, "_", symbol)) %>% 
  distinct(plot_id, .keep_all = TRUE) %>%
  group_by(core_histone, alignment_position) %>%
  # Create `symbol_stagger`
  mutate(
    symbol_stagger = cumsum(!duplicated(symbol) & !is.na(symbol))
  ) %>% 
  ungroup() %>%
  mutate(facet = case_when(
    alignment_position >= 1 & alignment_position <= 70 ~ 1,
    alignment_position > 70 ~ 2
  )) %>%
  group_by(core_histone, alignment_position) %>%
  # Remove rows where `letter` is "-" in case of multiple rows per grouping
  filter(!(n() > 1 & letter == "-")) %>%
  ungroup()

forplot <- forplot %>%
  mutate(y_stagger = as.numeric(factor(core_histone)) + symbol_stagger * 0.2)

ptms <- forplot %>%
  ggplot(aes(x = alignment_position)) +
  # Plot letters with red/black coloring based on `mi_mod`
  geom_text(
    aes(y = as.numeric(factor(core_histone)), label = letter, color = ifelse(!is.na(mi_mod), "red", "black")),
    size = 5
  ) +
  scale_y_continuous(
    breaks = c(1, 2, 3, 4), # Define positions on the y-axis
    labels = c("H2A", "H2B", "H3", "H4") # Manually set labels
  ) +
  # Map symbol to specific shapes and control their staggered position using `y_stagger`
  geom_point(
    aes(y = y_stagger, shape = symbol),
    size = 4,
    color = "blue"
  ) +
  theme_minimal() +
  # Ensure facet labels are removed and specify x-axis breaks every 20 residues
  labs(x = "Alignment Position", y = "Core Histone") +
  theme_classic(base_size = 24) +
  theme(
    strip.text.x = element_blank(), # Remove facet labels (strip.text)
    panel.spacing = unit(1, "lines"), # Ensure consistent spacing between facets
    axis.text.y = element_text(color = "black"), # Ensure readable y-axis text
    legend.position = "top" # Change legend position to bottom
  ) +
  scale_x_continuous(breaks = seq(0, max(forplot$alignment_position, na.rm = TRUE), by = 20)) +
  scale_color_identity() +
  # Map each symbol to a specific shape
  scale_shape_manual(
    values = c(
      "CTrmAmid" = 18,         # Diamond for CTrmAmid
      "Deamide" = 0,           # Square for Deamide
      "Acetyl/Trimethyl" = 8,         # Star for Acetyl/Trimethyl
      "Carbamyl" = 4,          # Cross for Carbamyl
      "Phosph" = 1,           # Circle for Phospho
      "IronAdduct" = 2,          # Triangle for IronAdduct
      "Dimethyl" = 6,          # Star shape for Dimethyl
      "Methyl" = 3,           # Plus for Methyl
      "Ox" = 16,          # Filled circle for Plus1Oxy
      "Ox2" = 17          # Filled triangle for Plus2Oxy
    ),
    name = "PTM", # Change legend title to "PTM"
    na.translate = FALSE # Remove NA values from the legend
  ) +
  # Divide into facets with free scales (separate x/y scales for facets)
  facet_wrap(~ facet, ncol = 1, scales = "free")
ptms
