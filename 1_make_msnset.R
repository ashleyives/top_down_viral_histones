
library(tidyverse)
library(TopPICR) 
library(MSnbase)
library(PNNL.DMS.utils)
library(MSnSet.utils)

load("toppic_output_td_5573.RData") #PrSM.txt and ms1.feature files 

#MS2 PSMs
ids <- toppic_output$ms2identifications%>%
  filter(!Dataset %in% badruns) %>%
  filter(`Retention time` < 3300) #remove anything eluting after 55 minutes

#MS1 Feature Data 
feat <- toppic_output$ms1features %>%
  filter(!Dataset %in% badruns) %>%
  filter(Time_apex < 55) %>% #remove anything eluting after 55 minutes
  mutate(Time_apex = Time_apex*60) #convert time to seconds


#A given feature can occur in multiple genes. This function selects a single gene for each feature. The gene with the minimum E-value is kept and all other genes are removed.
rm_false_gene_nocv <- function (x) 
{
  y <- x %>% dplyr::select(Dataset, `Feature apex`, `Feature intensity`, 
                           Gene, `E-value`) %>% dplyr::distinct() %>% dplyr::group_by(Dataset, 
                                                                                      `Feature apex`, `Feature intensity`) %>% dplyr::slice_min(`E-value`) %>% 
    dplyr::ungroup() %>% dplyr::select(Dataset, `Feature apex`, 
                                       `Feature intensity`, Gene)
  return(dplyr::semi_join(x, y))
}

x_err <- rm_false_gene_nocv(ids) 

#add length of gene to MS2 assignments------------------------------------------
fst <- Biostrings::readAAStringSet(
  file.path("ID_008569_DD63CF99.fasta"))

add_protein_length <- function(x, fst_obj){
  fst_df <- data.frame(`Protein accession` = names(fst_obj),
                       protLength = BiocGenerics::width(fst_obj), 
                       row.names = NULL, 
                       check.names = FALSE) %>%
    mutate(`Protein accession` = word(`Protein accession`))
  x <- inner_join(x, fst_df)
  return(x)
}

x_err <- add_protein_length(x_err, fst)

#check FDR 
compute_fdr(x_err)

#Filter, proteoform must have 2 PSMs--------------------------------------------
x_fbc <- filter_by_count(x = x_err,
                         count_within = c("Dataset", "Scan(s)"),
                         count = "cleanSeq",
                         threshold = 2) 

#check FDR 
compute_fdr(x_fbc)

# FDR control ------------------------------------------------------------------

# First find the E-value threshold for each of the three annotation types.
# FDR is calculated at Gene level.
the_cutoff <- find_evalue_cutoff(x = x_fbc,
                                 fdr_threshold = 0.01)

#Apply E-Value cutoff 
x_fdr <- x_fbc %>%
  filter(`E-value` < the_cutoff)

compute_fdr(x_fdr)

#The infer_prot function finds the smallest set of UniProt accessions (UniProtAcc) that maps to the largest number of amino acid sequences (cleanSeq).
x_ipf <- infer_prot(x_fdr)

compute_fdr(x_ipf)

#define proteoforms as species with same mass, gene, first and lastaa, and rounded masses for later clustering 
x_ipf <- x_ipf %>%
  mutate(simpProteoform = paste(Gene, firstAA, lastAA, round_any(`Adjusted precursor mass`,2, f = round), sep="_"))

#LOESS model for RT alignment across files--------------------------------------
form_model_ani <- function (x, ref_ds, ...) 
{
  x_ref <- x %>% dplyr::filter(Dataset == ref_ds) %>% dplyr::distinct(`Feature intensity`, 
                                                                      `Feature apex`, simpProteoform) %>% dplyr::group_by(simpProteoform) %>% 
    dplyr::slice_max(`Feature intensity`) %>% dplyr::slice_min(`Feature apex`)
  ds <- unique(x$Dataset)
  x_model <- vector(mode = "list", length = length(ds))
  for (e in 1:length(ds)) {
    x_cur <- x %>% dplyr::filter(Dataset == ds[[e]]) %>% 
      dplyr::distinct(`Feature intensity`, `Feature apex`, 
                      simpProteoform) %>% dplyr::group_by(simpProteoform) %>% 
      dplyr::slice_max(`Feature intensity`) %>% dplyr::slice_min(`Feature apex`)
    both <- dplyr::intersect(x_ref$simpProteoform, x_cur$simpProteoform)
    idx_ref <- which(x_ref$simpProteoform %in% both)
    idx_cur <- which(x_cur$simpProteoform %in% both)
    model_loess <- stats::loess(x_ref$`Feature apex`[idx_ref] ~ 
                                  x_cur$`Feature apex`[idx_cur], ...)
    x_model[[e]] <- model_loess
  }
  names(x_model) <- ds
  return(x_model)
}

the_model <- form_model_ani(
  x = x_ipf,
  ref_ds = find_ref_ds(x = x_ipf),
  control = loess.control(surface = "direct"), # Use direct to avoid NAs.
  span = 0.5,
  family = "symmetric")

# Align retention times according to the model created previously---------------
x_art <- align_rt(
  x = x_ipf,
  model = the_model,
  var_name = "Feature apex")

#sanity check alignment for two files, line should be x=y ideally-----------------------------
x_art %>%
  filter(Dataset == "TD_MRC5_Mock_24h_1_29Mar24_Arwen_WBEH-23-02-03") %>%
  ggplot(aes(x=RTalign, y=`Feature apex`))+
  geom_point()+
  geom_abline(slope=1, intercept = 0)

x_art %>%
  filter(Dataset == "TD_MRC5_Mock_24h_5_29Mar24_Arwen_WBEH-23-02-03") %>%
  ggplot(aes(x=RTalign, y=`Feature apex`))+
  geom_point()+
  geom_abline(slope=1, intercept = 0)

# Recalibrate the masses -------------------------------------------------------

# Calculate the error between the `Precursor mass` and the `Adjusted precursor
# mass`. This acts as the model for recalibrating the mass. A reference data set
# is not used when calculating the mass error model.
the_error <- calc_error(x = x_art,
                        ref_ds = find_ref_ds(x = x_ipf))

# Recalibrate the mass according the the errors computed previously.
x_rcm <- recalibrate_mass(x = x_art,
                          errors = the_error,
                          var_name = "Precursor mass")

#Cluster the data using hierarchical clustering with the hclust function. The data are clustered with the normalized recalibrated mass and the normalized aligned retention times.
x_clu <- TopPICR::cluster(x = x_rcm,
                          errors = the_error,
                          method = "single",
                          height = 1.5, 
                          min_size = 3) 

# Group clusters ---------------------------------------------------------------
#originally n_rt_sd = 1.5 
x_grp <- create_pcg(x = x_clu,
                    errors = the_error,
                    n_mme_sd = 4/the_error$ppm_sd, 
                    n_Da = 2, 
                    n_rt_sd = 5) 

#used for PTM plot
save(x_grp, file = "x_grp_td.RData")

x_meta0 <- create_mdata(x = x_grp,
                        errors = the_error,
                        n_mme_sd = 4/the_error$ppm_sd,
                        n_rt_sd = 150/the_error$rt_sd
                        # ppm_cutoff = 4,
                        # rt_cutoff = 150
)

x_meta <- x_meta0 %>%
  mutate(feature_name = paste0(Gene, "_", pcGroup))  %>%
  dplyr::select(collision, feature_name) %>%
  filter(collision != "-")

# Align the unidentified feature retention times with the model created by the
# identified feature retention times.
feat_art <- align_rt(
  x = feat,
  model = the_model,
  var_name = "Time_apex")

# Recalibrate unidentified feature mass -----------------------------------------

# Recalibrate the unidentified feature masses with the model created by the identified feature masses.
feat_rcm <- recalibrate_mass(
  x = feat_art,
  errors = the_error,
  var_name = "Mass")

# Perform match-between-runs----------------------------------------------------

#add a fake CV to group by for sake of function
feat_rcm <- feat_rcm %>%
  mutate(CV = 1)

x_grp <- x_grp%>%
  mutate(CV = 1)

feat_rtv <- match_features(ms2 = x_grp,
                           ms1 = feat_rcm,
                           errors = the_error,
                           # n_mme_sd = 4/the_error$ppm_sd,
                           # n_rt_sd = 150/the_error$rt_sd
                           n_mme_sd = 3,
                           n_rt_sd = 5
                           # ppm_cutoff = 4,
                           # rt_cutoff = 150
)

#remove fake CV
feat_rtv <- feat_rtv%>%
  dplyr::select(-CV)

#Select lowest E-value per proteoform to assign to meta data, could potentially be isomeric species for post-translationally modified proteoforms. 
meta0 <- x_grp %>%
  mutate(feature_name = paste(Gene, pcGroup, sep = "_"))%>%
  group_by(feature_name) %>%
  slice_min(`E-value`, n = 1) %>%
  ungroup() %>%
  dplyr::select(-Dataset, -Gene, -pcGroup)

#Prepare MSnSet object----------------------------------------------------------

x <- feat_rtv %>%
  mutate(feature_name = paste(Gene, pcGroup, sep = "_"),
         sample_name = Dataset) %>%
  dplyr::select(-Dataset)

x_expr <- x %>%
  pivot_wider(id_cols = "feature_name",
              names_from = "sample_name",
              values_from = "Intensity") %>%
  as.data.frame() %>%
  {rownames(.) <- .$feature_name;.} %>%
  dplyr::select(-feature_name) %>%
  as.matrix()

#count non-NA values for virus versus mock, must be in N=3 per condition for downstream quantification 

columns_229E <- grep("229E", colnames(x_expr)) 
columns_Mock <- grep("Mock", colnames(x_expr)) 

non_na_229E <- rowSums(!is.na(x_expr[, columns_229E]))
non_na_Mock <- rowSums(!is.na(x_expr[, columns_Mock]))

result <- data.frame(
  non_na_229E = non_na_229E,
  non_na_Mock = non_na_Mock
) %>%
  rownames_to_column(var="feature_name")

x_feat <- x %>%
  group_by(Gene,pcGroup, feature_name) %>%
  dplyr::summarize(median_intensity = median(Intensity),
                   count = n(),
                   .groups = "keep") %>%
  ungroup() %>%
  left_join(meta0, by=c("feature_name")) %>%
  left_join(result, by=c("feature_name")) %>%
  left_join(x_meta, by=c("feature_name")) %>%
  mutate(proteoform_id = paste(Gene, pcGroup, sep="_")) %>%
  as.data.frame() %>%
  {rownames(.) <- .$feature_name;.}

pd <- data.frame(SampleID =  colnames(x_expr)) %>%
  separate(SampleID, c("TD", "Cell", "Virus", "Time", "Biorep"), sep = "_", remove = FALSE) 

rownames(pd) <- pd$SampleID

m <- MSnSet(x_expr, x_feat[rownames(x_expr),], pd[colnames(x_expr),])

#Median normalize log-2 transformed intensities---------------------------------

#Log2 transform applied to expression crosstab of ExpressionSet or MSnSet object followed by centering the median around zero. 
m <- log2_zero_center(m)

#Check values before normalization
boxplot(exprs(m))

## Normalize data ----
m <- MSnSet.utils::normalizeByGlob(m, method = "once")

#Check values after normalization
boxplot(exprs(m))


#Annotate unassigned mass shifts using UniMod repository------------------------

annotate_Nterm_acetyls <- function(x, nterm_tol = 3, acetyl_id = "Acetyl"){
  
  temp_mod_names <- "mods"
  if("mod_names" %in% names(x$mods[[1]]))
    temp_mod_names <- "mod_names"
  
  for(i in 1:nrow(x)){
    if(x[i,"firstAA"] <= nterm_tol){
      y <- x[i,"mods"][[1]]
      if(length(y$mods) > 0){
        if(y$mods[1] == acetyl_id & y$mods_left_border[1] == 1){
          y$mods[1] <- paste("N-", acetyl_id, sep="")
          y[[temp_mod_names]][1] <- y$mods[1]
          posi_str <- map2_chr(y$mods_left_border, y$mods_right_border, paste, sep="-")
          y$mods_str <- paste(map2_chr(y[[temp_mod_names]], posi_str, paste, sep="@"), collapse=", ")
          x[i,"mods"][[1]] <- list(y)
        }
      }
    }
  }
  return(x)
}

unimods <- TopPICR::create_mod_data(
  mod_file = "TopPIC_Dynamic_Mods.txt", 
  mod_path = ".",
  use_unimod = TRUE
)

x <- fData(m)
x <- x %>%
  mutate(mods = map(Proteoform, TopPICR:::extract_mods))

mass_annotation_table <- TopPICR:::get_mass_annotation_table(x, unimods, 0.3, 0.1)
x <- x %>%
  mutate(mods = map(mods, TopPICR:::annotate_masses, mass_annotation_table, matching_tol = .Machine$double.eps))

# x <- annotate_Nterm_acetyls(as.data.frame(x), nterm_tol = 3, acetyl_id = "Acetyl")

rownames(x) <- x$proteoform_id

#reassign fData for MsnSet
fData(m) <- x

#save MsnSet object
save(m, file = "normalized_td_msnset.RData")

