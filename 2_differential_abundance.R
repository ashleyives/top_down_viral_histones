library(tidyverse)
library(MSnbase)
library(PNNL.DMS.utils)
library(MSnSet.utils)

load("normalized_td_msnset.RData")

#select for histone proteoforms only
m <- m[grepl("Histone", fData(m)$'Protein description'),] 

#only quantify features that have N=3 in each condition (mock and virus aka 229E)
m <- m[fData(m)$non_na_229E > 2,] 
m <- m[fData(m)$non_na_Mock > 2,] 

#needed to specify direction of comparison for limma 
m <- m[, with(pData(m), order(Virus))]
pData(m)$Virus <- factor(pData(m)$Virus, levels = c("Mock", "229E")) #you explicitly set "Mock" as the reference level. The first level in the levels argument is used as the reference.

res <- limma_gen(m, model.str = "~ Virus", coef.str = "Virus229E") %>%
  rownames_to_column(var="PF")

#add data completeness by counting number of NA obsverations  
x_final <- as.data.frame(m) %>%
  t()

rowwise_na_summary <- apply(x_final, 1, function(x) sum(!is.na(x))) %>%
  as.data.frame() %>%
  rownames_to_column(var = "PF")

rowwise_na_summary$notna <- rowwise_na_summary$.

rowwise_na_summary <- rowwise_na_summary %>%
  mutate(MissingPercent = (notna/ncol(x_final)))

res <- res %>%
  left_join(rowwise_na_summary, by="PF")  %>%
  mutate(MissingPercent = (notna/ncol(x_final))) %>%
  filter(MissingPercent > 0.5) 

save(res, file = "res.RData")
