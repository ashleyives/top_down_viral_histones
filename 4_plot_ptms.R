
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
           sep = "(?<=\\D)(?=\\d)", convert = TRUE) %>% 
  dplyr::rename(Accession = 'Protein accession') %>%
  mutate(site_id = paste0(Accession, "_", mi_adj_position)) %>%
  dplyr::select(-Gene)%>%
  filter(mi_mod != "Carbamyl")

h2a <- Biostrings::readAAStringSet(
  file.path("C:/Users/ives435/OneDrive - PNNL/histones/tong_mrc5_24hr/td/h2a_notzv_align.txt"))

h2azv <- Biostrings::readAAStringSet(
  file.path("C:/Users/ives435/OneDrive - PNNL/histones/tong_mrc5_24hr/td/h2azv_align.txt"))

h2b <- Biostrings::readAAStringSet(
  file.path("C:/Users/ives435/OneDrive - PNNL/histones/tong_mrc5_24hr/td/h2b_align.txt"))

h3 <- Biostrings::readAAStringSet(
  file.path("C:/Users/ives435/OneDrive - PNNL/histones/tong_mrc5_24hr/td/h3_align.txt"))

h4 <- Biostrings::readAAStringSet(
  file.path("C:/Users/ives435/OneDrive - PNNL/histones/tong_mrc5_24hr/td/h4_align.txt"))


wrangle_fasta <- function(x){
  
  y <- data.frame(
    Accession = names(x),                          
    FullSequence = as.character(x),                              
    row.names = NULL,
    check.names = FALSE
  ) %>%
    mutate(Accession = str_extract(Accession, "^\\S+"))  
  
  return(y)
}

h2a_df <- wrangle_fasta(h2a)
h2azv_df <- wrangle_fasta(h2azv)
h2b_df <- wrangle_fasta(h2b)
h3_df <- wrangle_fasta(h3)
h4_df <- wrangle_fasta(h4)

all_histones <- bind_rows(h2a_df, h2azv_df, h2b_df, h3_df, h4_df)

long_fasta <- all_histones %>%
  mutate(letter = strsplit(FullSequence, "")) %>%  
  unnest(letter) %>%                              
  group_by(Accession) %>%               
  mutate(alignment_position = row_number(),       
         sequence_position = cumsum(ifelse(letter != "-", 1, 0)) %>% 
           ifelse(letter == "-", NA, .)) %>%                     
  ungroup()  %>%
  mutate(core_histone = case_when(
    grepl("H2AV", Accession) ~ "H2Azv",
    grepl("H2AZ", Accession) ~ "H2Azv",
    grepl("H2AX", Accession) ~ "H2A",
    grepl("H2A2B", Accession) ~ "H2A",
    grepl("H2A2A", Accession) ~ "H2A",
    grepl("H2A1C", Accession) ~ "H2A",
    grepl("H2A1J", Accession) ~ "H2A",
    grepl("H2A1", Accession) ~ "H2A",
    grepl("H2A1D", Accession) ~ "H2A",
    grepl("H2A1B", Accession) ~ "H2A",
    grepl("H2AJ", Accession) ~ "H2A",
    grepl("H2A2C", Accession) ~ "H2A",
    grepl("H2B", Accession) ~ "H2B",
    grepl("H3", Accession) ~ "H3",
    grepl("H4", Accession) ~ "H4",
    TRUE ~ NA_character_ 
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
  group_by(core_histone, alignment_position) %>%
  summarize(
    mi_mod = paste(unique(na.omit(mi_mod)), collapse = "_"), 
    letter = paste(unique(na.omit(letter[letter != "-"])), collapse = "/"), 
    .groups = "drop"
  ) %>%
  ungroup() %>% 
  group_by(core_histone) %>%
  mutate(
    x = row_number() %% 10,    
    x = ifelse(x == 0, 10, x), 
    y = (row_number() - 1) %/% 10 + 1 
  ) %>%
  ungroup() 


forplot <- forplot %>%
  tidyr::separate_rows(mi_mod, sep = "_") %>%
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
  group_by(core_histone, alignment_position) %>%   #removes redundant rows for acetyl/trimethyl 
  distinct(symbol, .keep_all = TRUE) %>%          
  ungroup()  %>%                                     
  group_by(core_histone, alignment_position) %>%
  mutate(
    symbol_stagger = case_when(
      n() <= 1 ~ 0,  # 0 or 1 symbol → 0
      n() == 2 ~ ifelse(row_number() == 1, -0.5, 0.5), 
      n() == 3 ~ case_when(                 
        row_number() == 1 ~ -1,
        row_number() == 2 ~ 0,
        row_number() == 3 ~ 1
      ),
      TRUE ~ NA_real_   
    )
  ) %>% 
  mutate(letter = case_when(nchar(letter) > 1 ~ "X",
                            TRUE ~ letter)) 

ptms <- forplot %>%
  filter(core_histone != "H2Azv") %>%
  ggplot(aes(x=x, y=y, size=7))+
  geom_text(aes(label = letter, 
                color = ifelse(!is.na(symbol), "red", "black")), 
            fontface = "bold",  
            show.legend = FALSE)+
  scale_color_identity() +
  theme_minimal()+
  theme(
    legend.text = element_text(size = 16),  # Font size for the legend labels
    legend.title = element_text(size = 16), 
    panel.grid = element_blank(),           # Remove panel grid
    axis.ticks = element_blank(),          # Remove axis ticks
    axis.text.x = element_blank(),         # Remove x-axis text
    axis.text.y = element_blank(),         # Remove y-axis text
    axis.title.x = element_blank(),        # Remove x-axis label
    axis.title.y = element_blank(),        # Remove y-axis label
    strip.text.x = element_text(size = 24, hjust = 0, fontface="bold"),  # Left-align labels
    legend.position = "top", 
    panel.border = element_rect(
      colour = "black",                    # Black outer border
      fill = NA, 
      size = 1                             # Thickness of border
    )
  ) +
  geom_point(
    aes(y = y-0.5, x= x-symbol_stagger/2,shape = symbol),
    size = 4,
    color = "blue"
  ) +
  scale_y_reverse()  +
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
  guides(
    shape = guide_legend(
      nrow = 1,                # Force a single row in the legend
      byrow = TRUE,            # Arrange items by row (useful for multiple rows)
      title.position = "left",  # Position legend title on top of the keys
      title.hjust = 0.5       # Center-align the title
    )
  ) +
  facet_wrap(~core_histone, scales = "free",
             labeller = as_labeller(c("H2A" = "(A) H2A", "H2B" = "(B) H2B", "H3" = "(C) H3", "H4" = "(D) H4")))
ptms
