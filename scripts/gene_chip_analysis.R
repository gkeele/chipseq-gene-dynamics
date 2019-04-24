library(lme4)
library(tidyverse)
library(brms)
library(ggplot2)

setwd("~/projects/chipseq-gene-dynamics/")
#setwd("~/Documents/git_repositories/chipseq-gene-dynamics/")


## Read in the data
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


## ggplot theme
gg_theme <- theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  axis.line = element_line(colour = "black"),
                  plot.title = element_text(hjust = 0.5), 
                  axis.text = element_text(size = 12, face = "bold"),
                  axis.title = element_text(size = 12, face = "bold"),
                  axis.text.x = element_text(hjust = 1, face = "bold"))

###### Individual gene plots
practice <- gene_chip_dat %>%
  filter(gene %in% unique(gene_chip_dat$gene)[1]) %>%
  group_by(sample_unit)

## Linear
g <- ggplot(data = practice, aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme
g

## LOESS (curves)
g <- ggplot(data = practice, aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "loess")
g <- g + gg_theme + guides(color = FALSE)
g



######## Example Stan analysis
## Try YAL012W and YAL066W
## Function to run brms
run_brms_on_chipseq <- function(chipseq_dat, 
                                this_gene, 
                                adapt_delta = 0.8,
                                iter = 2000) {
  use_dat <- chipseq_dat %>% filter(gene == this_gene)#,
                                    #!(value %in% c(0, 1)))
  # writer loss
  wl_fit <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                      data = use_dat %>% filter(assay == "writer_loss"), 
                      family = "beta",
                      iter = iter,
                      control = list(adapt_delta = adapt_delta))
  wl_fixed <- summary(wl_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
  wl_random <- summary(wl_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
  # writer-eraser loss
  wel_fit <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                       data = use_dat %>% filter(assay == "writer_eraser_loss"), 
                       family = "beta",
                       iter = iter,
                       control = list(adapt_delta = adapt_delta))
  wel_fixed <- summary(wel_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
  wel_random <- summary(wel_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
  # writer add
  wa_fit <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                      data = use_dat %>% filter(assay == "writer_add"), 
                      family = "beta",
                      iter = iter,
                      control = list(adapt_delta = adapt_delta))
  wa_fixed <- summary(wa_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
  wa_random <- summary(wa_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")

  names(wl_fixed) <- names(wl_random) <- names(wel_fixed) <- names(wel_random) <- names(wa_fixed) <- names(wa_random) <- c("parameter", "estimate", "est_error", "lower_95", "upper_95", "eff_sample", "rhat")
  results <- list(writer_loss = bind_rows(wl_fixed, wl_random),
                  writer_eraser_loss = bind_rows(wel_fixed, wel_random),
                  writer_add = bind_rows(wa_fixed, wa_random))
  results
}

yal012w_models <- run_brms_on_chipseq(chipseq_dat = gene_chip_dat, this_gene = "YAL012W")
yal066w_models <- run_brms_on_chipseq(chipseq_dat = gene_chip_dat, this_gene = "YAL066W")


## Gene with extreme low effect
yal041w_models <- run_brms_on_chipseq(chipseq_dat = gene_chip_dat, this_gene = "YAL041W", adapt_delta = 0.99, iter = 5000)


## Parsing full set of Stan results
#stan_fit_paths <- list.files("individual_gene_Stan_models/", full.names = TRUE)
#stan_fit_paths <- list.files("individual_gene_Stan_models_alpha_0.99_iter_5000//", full.names = TRUE)
stan_fit_paths <- list.files("individual_gene_Stan_models_alpha_0.8_iter_10000/", full.names = TRUE)
stan_fit_paths <- stan_fit_paths[!(grepl(x = stan_fit_paths, pattern = "WL_vs_WEL"))]
#file_path <- "individual_gene_Stan_models//sacCer3_nonOverlappingGenes_noChrMGenes_"
#file_path <- "individual_gene_Stan_models_alpha_0.99_iter_5000///sacCer3_nonOverlappingGenes_noChrMGenes_"
file_path <- "individual_gene_Stan_models_alpha_0.8_iter_10000//sacCer3_nonOverlappingGenes_noChrMGenes_"

for (i in 1:length(stan_fit_paths)) {
  this_fit <- readRDS(stan_fit_paths[i])

  this_gene <- gsub(x = stan_fit_paths[i], pattern = file_path, replacement = "", fixed = TRUE) %>%
    gsub(x = ., replacement = "", pattern = "_Stan.rds", fixed = TRUE)
  
  
  if (i == 1) {
    timepoint_dat <- bind_rows(data.frame(gene = this_gene, assay = "writer_loss", this_fit$writer_loss %>% filter(parameter == "timepoint")),
                               data.frame(gene = this_gene, assay = "writer_eraser_loss", this_fit$writer_eraser_loss %>% filter(parameter == "timepoint")),
                               data.frame(gene = this_gene, assay = "writer_add", this_fit$writer_add %>% filter(parameter == "timepoint")))
  }
  else {
    timepoint_dat <- bind_rows(timepoint_dat,
                               data.frame(gene = this_gene, assay = "writer_loss", this_fit$writer_loss %>% filter(parameter == "timepoint")),
                               data.frame(gene = this_gene, assay = "writer_eraser_loss", this_fit$writer_eraser_loss %>% filter(parameter == "timepoint")),
                               data.frame(gene = this_gene, assay = "writer_add", this_fit$writer_add %>% filter(parameter == "timepoint")))
  }
}
timepoint_dat <- timepoint_dat %>% 
  select(-parameter) %>%
  mutate(is_pos = ifelse(lower_95 > 0, TRUE, FALSE),
         is_neg = ifelse(upper_95 < 0, TRUE, FALSE))
timepoint_dat$category[timepoint_dat$lower_95 > 0] <- "positive"
timepoint_dat$category[timepoint_dat$upper_95 < 0] <- "negative"
timepoint_dat$category[timepoint_dat$lower_95 < 0 & timepoint_dat$upper_95 > 0] <- "zero"

  
ggplot(data = timepoint_dat, aes(x = estimate)) + geom_histogram() + facet_grid(assay~.)

## Grabbing extreme negative genes
low_genes <- timepoint_dat %>% 
  filter(estimate < -10,
         category == "negative") %>%
  pull(gene) %>%
  unique

## Grabbing high genes
high_genes <- timepoint_dat %>% 
  filter(category == "positive") %>%
  pull(gene) %>%
  unique

## Grabbing low genes
high_genes <- timepoint_dat %>% 
  filter(category == "positive") %>%
  pull(gene) %>%
  unique

g <- ggplot(data = gene_chip_dat %>% filter(gene %in% high_genes[1:6]), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + guides(color = FALSE)
g
g <- ggplot(data = gene_chip_dat %>% filter(gene %in% high_genes[1:6]), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + guides(color = FALSE)
g
g <- ggplot(data = gene_chip_dat %>% filter(gene %in% genes_all_zero[1:6]), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
#g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + guides(color = FALSE)
g
g <- ggplot(data = gene_chip_dat %>% filter(gene %in% genes_all_nonzero[1:6]), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
#g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + guides(color = FALSE)
g


genes_all <- unique(gene_chip_dat$gene)
genes_pos_wa <- timepoint_dat %>% 
  filter(assay == "writer_add",
         category == "positive") %>%
  pull(gene)
genes_neg_wl <- timepoint_dat %>% 
  filter(assay == "writer_loss",
         category == "negative") %>%
  pull(gene)
genes_neg_wel <- timepoint_dat %>% 
  filter(assay == "writer_eraser_loss",
         category == "negative") %>%
  pull(gene)
genes_neg_loss <- intersect(genes_neg_wl, genes_neg_wel)
genes_overlap <- intersect(genes_pos_wa, genes_neg_wl) %>% intersect(genes_neg_wel)

genes_neg_loss_only <- genes_neg_loss[!(genes_neg_loss %in% genes_pos_wa)]

library(UpSetR)
upset_dat <- data.frame(genes = genes_all,
                        WA_positive = sapply(1:length(genes_all), function(i) as.numeric(genes_all[i] %in% genes_pos_wa)),
                        WL_negative = sapply(1:length(genes_all), function(i) as.numeric(genes_all[i] %in% genes_neg_wl)),
                        WEL_negative = sapply(1:length(genes_all), function(i) as.numeric(genes_all[i] %in% genes_neg_wel)))
upset(upset_dat, order.by = "freq", 
      text.scale = 1.2,
      sets.bar.color = c("coral", "magenta", "seagreen1"))

## Time trend estimate by error
plot(timepoint_dat$est_error, timepoint_dat$estimate, pch = ifelse(timepoint_dat$category == "zero", 1, 19), col = c("coral", "magenta", "seagreen1")[as.factor(timepoint_dat$assay)], 
     las = 1, ylab = "Trend with time", xlab = "Error on trend")
abline(h = 0, lty = 2)

time_alt_fit <- lm(estimate ~ 1 + assay, data = timepoint_dat, weights = 1/timepoint_dat$est_error)
time_null_fit <- lm(estimate ~ 1, data = timepoint_dat, weights = 1/timepoint_dat$est_error)
anova(time_alt_fit, time_null_fit)

## Time trend for loss categories only
weights <- 1/timepoint_dat %>% filter(assay %in% c("writer_loss", "writer_eraser_loss")) %>% pull(est_error)
time_alt_fit <- lm(estimate ~ 1 + assay, data = timepoint_dat %>% filter(assay %in% c("writer_loss", "writer_eraser_loss")), weights = weights)
time_null_fit <- lm(estimate ~ 1, data = timepoint_dat %>% filter(assay %in% c("writer_loss", "writer_eraser_loss")), weights = weights)
anova(time_alt_fit, time_null_fit)


## Gene level
gene_timepoint_dat <- timepoint_dat %>% 
  select(gene, assay, estimate) %>%
  spread(key = assay, value = estimate) 

genes_all_zero <- timepoint_dat %>% 
  select(gene, assay, category) %>%
  spread(key = assay, value = category) %>%
  filter(writer_add == "zero", 
         writer_eraser_loss == "zero",
         writer_loss == "zero") %>%
  pull(gene)
genes_all_nonzero <- timepoint_dat %>% 
  select(gene, assay, category) %>%
  spread(key = assay, value = category) %>%
  filter(writer_add == "positive", 
         writer_eraser_loss == "negative",
         writer_loss == "negative") %>%
  pull(gene)
genes_all_negative <- timepoint_dat %>% 
  select(gene, assay, category) %>%
  spread(key = assay, value = category) %>%
  filter(writer_eraser_loss == "negative",
         writer_loss == "negative") %>%
  pull(gene)

plot(gene_timepoint_dat$writer_eraser_loss, gene_timepoint_dat$writer_loss)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)


## Looking at differences between wl and wel
loss_compare_dat <- gene_chip_dat %>%
  filter(assay %in% c("writer_loss", "writer_eraser_loss"))

practice_dat <- loss_compare_dat %>% 
  filter(gene == gene_chip_dat$gene[10])
compare_fit <- brms::brm(value ~ 1 + timepoint_cat + assay + timepoint_cat:assay + (1 + timepoint_cat + assay + timepoint_cat:assay | sample_unit), 
                    data = practice_dat, 
                    family = "beta",
                    iter = 2000,
                    control = list(adapt_delta = 0.8))

compare_wl_to_wel_brms <- function(chipseq_dat, 
                                   this_gene, 
                                   adapt_delta = 0.8,
                                   iter = 2000) {
  use_dat <- chipseq_dat %>% filter(gene == this_gene,
                                    assay %in% c("writer_loss", "writer_eraser_loss"))
  
  compare_fit <- brms::brm(value ~ 1 + timepoint_cat + assay + timepoint_cat:assay + (1 + timepoint_cat + assay + timepoint_cat:assay | sample_unit), 
                           data = use_dat, 
                           family = "beta",
                           iter = iter,
                           control = list(adapt_delta = adapt_delta))
  compare_fixed <- summary(compare_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
  compare_fixed
}

practice <- compare_wl_to_wel_brms(chipseq_dat = gene_chip_dat,
                                   this_gene = gene_chip_dat$gene[100])

loss_compare_path <- "individual_gene_Stan_models_alpha_0.8_iter_10000//"
loss_compare_files <- list.files(loss_compare_path, full.names = TRUE)
loss_compare_files <- grep(x = loss_compare_files, pattern = "WL_vs_WEL", value = TRUE)
loss_compare_remove <- "individual_gene_Stan_models_alpha_0.8_iter_10000///sacCer3_nonOverlappingGenes_noChrMGenes_"
for (i in 1:length(loss_compare_files)) {
  this_fit <- readRDS(loss_compare_files[i])
  
  this_gene <- gsub(x = loss_compare_files[i], pattern = loss_compare_remove, replacement = "", fixed = TRUE) %>%
    gsub(x = ., replacement = "", pattern = "_WL_vs_WEL_Stan.rds", fixed = TRUE)
  
  if (i == 1) {
    loss_compare_dat <- data.frame(gene = rep(this_gene, 8), this_fit)
  }
  else {
    loss_compare_dat <- bind_rows(loss_compare_dat,
                                  data.frame(gene = rep(this_gene, 8), this_fit))
  }
}

## Attempt to grab genes with "significant" lag
total_lag <- loss_compare_dat %>%
  filter(grepl(x = parameter, pattern = "assay")) %>%
  mutate(parameter = recode(parameter, 
                            'assaywriter_eraser_loss' = "0",
                            'timepoint_cat1:assaywriter_eraser_loss' = "1",
                            'timepoint_cat2:assaywriter_eraser_loss' = "2",
                            'timepoint_cat3:assaywriter_eraser_loss' = "3")) %>%
  rename(timepoint = parameter) %>%
  mutate(lag = l.95..CI > 0)

total_lag <- lag_dat %>%
  group_by(gene) %>%
  summarize(lag_count = sum(lag)) %>%
  ungroup

lagged_genes <- total_lag %>%
  filter(lag_count > 0) %>%
  pull(gene)

g <- ggplot(data = gene_chip_dat %>% filter(assay != "writer_add", gene %in% lagged_genes), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
#g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + guides(color = FALSE)
g
# Mostly detect ugly things, don't have power at the single gene level

## Overall picture of lag
plot_loss_compare_dat <- loss_compare_dat %>%
  filter(grepl(x = parameter, pattern = "assay")) %>%
  mutate(parameter = recode(parameter, 
                            'assaywriter_eraser_loss' = "0",
                            'timepoint_cat1:assaywriter_eraser_loss' = "1",
                            'timepoint_cat2:assaywriter_eraser_loss' = "2",
                            'timepoint_cat3:assaywriter_eraser_loss' = "3")) %>%
  rename(timepoint = parameter)
g <- ggplot(data = plot_loss_compare_dat, aes(x = timepoint, y = Estimate)) + geom_point(col = "gray") + geom_line(aes(group = gene), col = "gray")
g <- g + scale_color_grey()
g <- g + geom_boxplot(aes(x = timepoint, y = Estimate), color = "tomato")
g <- g + geom_hline(yintercept = 0, linetype = "dashed", color = "red") + gg_theme
g

g <- ggplot(data = plot_loss_compare_dat %>% filter(gene %in% genes_all_negative), aes(x = timepoint, y = Estimate)) + geom_point(col = "gray") + geom_line(aes(group = gene), col = "gray")
g <- g + scale_color_grey()
g <- g + geom_boxplot(aes(x = timepoint, y = Estimate), color = "cadetblue") + stat_summary(fun.y=mean, geom="line", aes(group=1), col = "black", size = 1.5) + ylab("WEL - WL timepoint effects")
g <- g + geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.5) + gg_theme
g

## Multi-gene model
multigene_fit <- brms::brm(value ~ 1 + timepoint_cat + assay + timepoint_cat:assay + (1 + timepoint_cat + assay + timepoint_cat:assay | sample_unit + gene), 
                           data = gene_chip_dat %>% filter(assay %in% c("writer_loss", "writer_eraser_loss")), 
                           family = "beta",
                           iter = 2000,
                           control = list(adapt_delta = 0.8))


##############################################
##
## Additional gene information
##
##############################################
wa_covar_raw_dat <- data.table::fread("~/projects/chipseq-gene-dynamics/data/DtL_setd2_RNA_combined_sf.txt", data.table = FALSE)
wa_covar_dat <- wa_covar_raw_dat %>%
  select(GENE, DtL_0min_set2d_RNA_Rep1, DtL_0min_set2d_RNA_Rep2, DtL_0min_set2d_RNA_Rep3) %>%
  rename(rep1 = DtL_0min_set2d_RNA_Rep1,
         rep2 = DtL_0min_set2d_RNA_Rep2,
         rep3 = DtL_0min_set2d_RNA_Rep3) %>%
  gather(key = rep, value = transcript, -GENE) %>%
  group_by(GENE) %>%
  summarize(transcript = mean(log(transcript))) %>%
  ungroup %>%
  rename(gene = GENE)

wl_covar_raw_dat <- data.table::fread("~/projects/chipseq-gene-dynamics/data/LtD_setd2_RNA_combined_sf.txt", data.table = FALSE)
wl_covar_dat <- wl_covar_raw_dat %>%
  select(GENE, LtD_0min_set2d_RNA_Rep1, LtD_0min_set2d_RNA_Rep2, LtD_0min_set2d_RNA_Rep3) %>%
  rename(rep1 = LtD_0min_set2d_RNA_Rep1,
         rep2 = LtD_0min_set2d_RNA_Rep2,
         rep3 = LtD_0min_set2d_RNA_Rep3) %>%
  gather(key = rep, value = transcript, -GENE) %>%
  group_by(GENE) %>%
  summarize(transcript = mean(log(transcript))) %>%
  ungroup %>%
  rename(gene = GENE)


wa_full_dat <- timepoint_dat %>%
  filter(assay == "writer_add") %>%
  left_join(wa_covar_dat)

wl_full_dat <- timepoint_dat %>%
  filter(assay == "writer_loss") %>%
  left_join(wl_covar_dat)

plot(wa_full_dat$transcript, wa_full_dat$est_error)
plot(wa_full_dat$estimate, wa_full_dat$est_error)

plot(wl_full_dat$transcript, wl_full_dat$est_error)
plot(wl_full_dat$estimate, wl_full_dat$est_error)

par(mfrow = c(1, 2))
plot(wa_full_dat$transcript, wa_full_dat$estimate)
plot(wl_full_dat$transcript, wl_full_dat$estimate)


##############################################
##
## Assay specific raw data plots
##
##############################################
mean_gene_chip_dat <- gene_chip_dat %>%
  group_by(gene, assay, timepoint) %>%
  summarise(value = mean(value)) %>%
  ungroup
# Writer loss
writer_loss <- ggplot(data = mean_gene_chip_dat %>% filter(assay == "writer_loss"), aes(x = timepoint, y = value, color = gene)) + geom_point() + geom_line()
writer_loss <- writer_loss + scale_color_grey()
writer_loss <- writer_loss + geom_smooth(aes(y = value, x = timepoint), method = "lm", col = "magenta") + ggtitle("Writer loss")
writer_loss <- writer_loss + theme(panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), 
                                   axis.line = element_line(colour = "black"),
                                   plot.title = element_text(hjust = 0.5), 
                                   axis.text = element_text(size = 12, face = "bold"),
                                   axis.title = element_text(size = 12, face = "bold"),
                                   axis.text.x = element_text(hjust = 1, face = "bold")) + guides(color = FALSE)
writer_loss
# Writer and eraser loss
writer_eraser_loss <- ggplot(data = mean_gene_chip_dat %>% filter(assay == "writer_eraser_loss"), aes(x = timepoint, y = value, color = sample_unit)) + geom_point() + geom_line()
writer_eraser_loss <- writer_eraser_loss + scale_color_grey()
writer_eraser_loss <- writer_eraser_loss + geom_smooth(aes(y = value, x = timepoint), method = "lm", col = "seagreen1") + ggtitle("Writer and eraser loss")
writer_eraser_loss <- writer_eraser_loss + theme(panel.grid.major = element_blank(), 
                                                 panel.grid.minor = element_blank(),
                                                 panel.background = element_blank(), 
                                                 axis.line = element_line(colour = "black"),
                                                 plot.title = element_text(hjust = 0.5), 
                                                 axis.text = element_text(size = 12, face = "bold"),
                                                 axis.title = element_text(size = 12, face = "bold"),
                                                 axis.text.x = element_text(hjust = 1, face = "bold")) + guides(color = FALSE)
writer_eraser_loss
# Writer add
writer_add <- ggplot(data = mean_gene_chip_dat %>% filter(assay == "writer_add"), aes(x = timepoint, y = value, color = sample_unit)) + geom_point() + geom_line()
writer_add <- writer_add + scale_color_grey()
writer_add <- writer_add + geom_smooth(aes(y = value, x = timepoint), method = "lm", col = "coral") + ggtitle("Writer add")
writer_add <- writer_add + theme(panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), 
                                 axis.line = element_line(colour = "black"),
                                 plot.title = element_text(hjust = 0.5), 
                                 axis.text = element_text(size = 12, face = "bold"),
                                 axis.title = element_text(size = 12, face = "bold"),
                                 axis.text.x = element_text(hjust = 1, face = "bold")) + guides(color = FALSE)
writer_add

gridExtra::grid.arrange(writer_loss, writer_eraser_loss, writer_add, nrow = 1)


