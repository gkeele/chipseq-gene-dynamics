library(lme4)
library(tidyverse)
library(brms)
library(ggplot2)

setwd("/proj/strahllb/users/Austin/Lerner_LANS_set2d/chipseq-gene-dynamics")
gene_chip_dat <- read.table("data/scaled_matrix_for_LMM_model_nonOverlappingGenes.txt") # path to data
names(gene_chip_dat) <- c("gene", "assay", "replicate", "timepoint", "value")

gene_chip_dat <- gene_chip_dat %>%
  mutate(assay = factor(assay),
         assay = recode(assay, "1" = "writer_loss", "2" = "writer_eraser_loss", "3" = "writer_add"),
         sample_unit = paste(gene, assay, replicate, sep = "_"),
         value = ifelse(value == 0, value + 0.0001, value),
         value = ifelse(value == 1, value - 0.0001, value),
         timepoint_cat = factor(timepoint))
gene_chip_dat$timepoint[gene_chip_dat$assay %in% c("writer_loss", "writer_eraser_loss") & gene_chip_dat$timepoint_cat == 0] <- 0
gene_chip_dat$timepoint[gene_chip_dat$assay %in% c("writer_loss", "writer_eraser_loss") & gene_chip_dat$timepoint_cat == 1] <- 30/60
gene_chip_dat$timepoint[gene_chip_dat$assay %in% c("writer_loss", "writer_eraser_loss") & gene_chip_dat$timepoint_cat == 2] <- 60/60
gene_chip_dat$timepoint[gene_chip_dat$assay %in% c("writer_loss", "writer_eraser_loss") & gene_chip_dat$timepoint_cat == 3] <- 90/60
gene_chip_dat$timepoint[gene_chip_dat$assay == "writer_add" & gene_chip_dat$timepoint_cat == 0] <- 0
gene_chip_dat$timepoint[gene_chip_dat$assay == "writer_add" & gene_chip_dat$timepoint_cat == 1] <- 20/60
gene_chip_dat$timepoint[gene_chip_dat$assay == "writer_add" & gene_chip_dat$timepoint_cat == 2] <- 40/60
gene_chip_dat$timepoint[gene_chip_dat$assay == "writer_add" & gene_chip_dat$timepoint_cat == 3] <- 60/60
 
## Multi-gene model
multigene_fit <- brms::brm(value ~ 1 + timepoint_cat + assay + timepoint_cat:assay + (1 + timepoint_cat + assay + timepoint_cat:assay | sample_unit + gene), 
                           data = gene_chip_dat %>% filter(assay %in% c("writer_loss", "writer_eraser_loss")), 
                           family = "beta",
                           iter = 2000,
                           control = list(adapt_delta = 0.8),
						   cores=288)
						   
saveRDS(multigene_fit, "test_multigene.rds")
