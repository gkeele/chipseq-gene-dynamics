#######################################################################
##
##            Histone tri-methyl dynamics project
##
##            Authors: Greg Keele & Austin Hepperla
##
#######################################################################

## R packages
library(lme4)
library(tidyverse)
library(brms)
library(ggplot2)
library(gridExtra)
library(UpSetR)
library(scales)

## Colors
wa_col <- "purple"
wl_col <- "black"
wel_col <- "red"
# wa_col <- "coral"
# wl_col <- "magenta"
# wel_col <- "seagreen1"

## Set directory for local environment to the git repos
setwd("~/projects/chipseq-gene-dynamics/")
#setwd("~/Documents/git_repositories/chipseq-gene-dynamics/")

## ggplot theme
gg_theme <- theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  axis.line = element_line(colour = "black"),
                  plot.title = element_text(hjust = 0.5), 
                  axis.text = element_text(size = 12, face = "bold"),
                  axis.title = element_text(size = 12, face = "bold"),
                  axis.text.x = element_text(hjust = 1, face = "bold"))

###################################################
##
##  Read in and process the data
##
##  ChIP-seq proprotions
##
###################################################
raw_gene_chip_dat <- read.table("data/scaled_matrix_for_LMM_model_nonOverlappingGenes.txt") # path to data
names(raw_gene_chip_dat) <- c("gene", "assay", "replicate", "timepoint", "value")

gene_chip_dat <- raw_gene_chip_dat %>%
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

##############################################
##
##  Individual gene plots
##
##############################################
single_gene_dat <- gene_chip_dat %>%
  filter(gene %in% unique(gene_chip_dat$gene)[1]) %>%
  group_by(sample_unit)

## Linear
g <- ggplot(data = single_gene_dat, aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme
g

##############################################
##
##  Assay for all genes raw data plots
##
##############################################
mean_gene_chip_dat <- gene_chip_dat %>%
  group_by(gene, assay, timepoint) %>%
  summarise(value = mean(value)) %>%
  ungroup
# Writer loss
writer_loss <- ggplot(data = mean_gene_chip_dat %>% filter(assay == "writer_loss"), aes(x = timepoint, y = value, color = gene)) + geom_point() + geom_line()
writer_loss <- writer_loss + scale_color_grey()
writer_loss <- writer_loss + geom_smooth(aes(y = value, x = timepoint), method = "lm", col = wl_col) + ggtitle("Writer loss")
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
writer_eraser_loss <- ggplot(data = mean_gene_chip_dat %>% filter(assay == "writer_eraser_loss"), aes(x = timepoint, y = value, color = gene)) + geom_point() + geom_line()
writer_eraser_loss <- writer_eraser_loss + scale_color_grey()
writer_eraser_loss <- writer_eraser_loss + geom_smooth(aes(y = value, x = timepoint), method = "lm", col = wel_col) + ggtitle("Writer and eraser loss")
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
writer_add <- ggplot(data = mean_gene_chip_dat %>% filter(assay == "writer_add"), aes(x = timepoint, y = value, color = gene)) + geom_point() + geom_line()
writer_add <- writer_add + scale_color_grey()
writer_add <- writer_add + geom_smooth(aes(y = value, x = timepoint), method = "lm", col = wa_col) + ggtitle("Writer add")
writer_add <- writer_add + theme(panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), 
                                 axis.line = element_line(colour = "black"),
                                 plot.title = element_text(hjust = 0.5), 
                                 axis.text = element_text(size = 12, face = "bold"),
                                 axis.title = element_text(size = 12, face = "bold"),
                                 axis.text.x = element_text(hjust = 1, face = "bold")) + guides(color = FALSE)
writer_add

grid.arrange(writer_loss, writer_eraser_loss, writer_add, nrow = 1)

##############################################
##
##  Beta regression: Stan model
##
##############################################
## Function to run brms for beta regression
#### ChIP-seq normalized to proportions within a replicate
#### Per gene per assay
run_brms_on_chipseq <- function(chipseq_dat, 
                                this_gene, 
                                adapt_delta = 0.8,
                                iter = 2000) {
  use_dat <- chipseq_dat %>% filter(gene == this_gene)#,
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
## Try YAL012W and YAL066W
yal012w_models <- run_brms_on_chipseq(chipseq_dat = gene_chip_dat, this_gene = "YAL012W")
yal066w_models <- run_brms_on_chipseq(chipseq_dat = gene_chip_dat, this_gene = "YAL066W")
## Gene with extreme low effect
yal041w_models <- run_brms_on_chipseq(chipseq_dat = gene_chip_dat, this_gene = "YAL041W", adapt_delta = 0.99, iter = 5000)

##############################################
##
##  Beta regression: aggregate results
##
##############################################
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
## Histogram
ggplot(data = timepoint_dat, aes(x = estimate)) + geom_histogram() + facet_grid(assay~.)
## All genes
all_genes <- unique(gene_chip_dat$gene)
## Grabbing high genes
positive_genes <- timepoint_dat %>% 
  filter(category == "positive") %>%
  pull(gene) %>%
  unique
## Grabbing low genes
negative_genes <- timepoint_dat %>% 
  filter(category == "negative") %>%
  pull(gene) %>%
  unique
## Grabbing more categories of genes
genes_neg_wl <- timepoint_dat %>% 
  filter(assay == "writer_loss",
         category == "negative") %>%
  pull(gene)
genes_neg_wel <- timepoint_dat %>% 
  filter(assay == "writer_eraser_loss",
         category == "negative") %>%
  pull(gene)
trend_genes <- intersect(positive_genes, genes_neg_wl) %>% intersect(genes_neg_wel)
genes_neg_loss_only <- negative_genes[!(negative_genes %in% genes_pos_wa)]
## Grabbing genes with no trends
genes_all_zero <- timepoint_dat %>% 
  select(gene, assay, category) %>%
  spread(key = assay, value = category) %>%
  filter(writer_add == "zero", 
         writer_eraser_loss == "zero",
         writer_loss == "zero") %>%
  pull(gene)

##############################################
##
##  Beta regression: timecourse plots
##
##############################################
## Genes with positive trends (writer add)
g <- ggplot(data = gene_chip_dat %>% filter(gene %in% positive_genes[1:6]), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + guides(color = FALSE)
g
## Genes with no trends
g <- ggplot(data = gene_chip_dat %>% filter(gene %in% genes_all_zero[1:6]), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
#g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + guides(color = FALSE)
g
## Genes with all trends
g <- ggplot(data = gene_chip_dat %>% filter(gene %in% trend_genes[1:6]), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
#g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + guides(color = FALSE)
g
##############################################
##
##  Beta regression: Upset plot
##
##############################################
upset_dat <- data.frame(genes = all_genes,
                        WA_positive = sapply(1:length(all_genes), function(i) as.numeric(all_genes[i] %in% genes_pos_wa)),
                        WL_negative = sapply(1:length(all_genes), function(i) as.numeric(all_genes[i] %in% genes_neg_wl)),
                        WEL_negative = sapply(1:length(all_genes), function(i) as.numeric(all_genes[i] %in% genes_neg_wel)))
upset(upset_dat, order.by = "freq", 
      text.scale = 1.2,
      sets.bar.color = c(wa_col, wl_col, wel_col))
##############################################
##
##  Beta regression: Trend estimate by error
##
##############################################
plot(timepoint_dat$est_error, timepoint_dat$estimate, pch = ifelse(timepoint_dat$category == "zero", 1, 19), col = c(wa_col, wl_col, wel_col)[as.factor(timepoint_dat$assay)], 
     las = 1, ylab = "Histone trend with time", xlab = "Error on trend", frame.plot = FALSE)
abline(h = 0, lty = 2)
legend("bottomright", 
       legend = c("writer add", "writer loss", "writer eraser loss"),
       col = c(wa_col, wl_col, wel_col),
       fill = c(wa_col, wl_col, wel_col),
       bty = "n")
##############################################
##
##  Beta regression: Formal test of assay 
##                   on trend
##
##############################################
time_alt_fit <- lmer(estimate ~ (1 | gene) + assay, 
                     data = timepoint_dat, 
                     weights = 1/timepoint_dat$est_error,
                     REML = FALSE)
time_null_fit <- lmer(estimate ~ (1 | gene), 
                      data = timepoint_dat, 
                      weights = 1/timepoint_dat$est_error,
                      REML = FALSE)
anova(time_alt_fit, time_null_fit)
# p-value < 2.2e-16

## Time trend for loss categories only
loss_dat <- timepoint_dat %>%
  filter(assay %in% c("writer_loss", "writer_eraser_loss"))
loss_weights <- 1/loss_dat$est_error
time_loss_alt_fit <- lmer(estimate ~ (1 | gene) + assay, 
                     data = loss_dat, 
                     weights = loss_weights,
                     REML = FALSE)
time_loss_null_fit <- lmer(estimate ~ (1 | gene), 
                      data = loss_dat, 
                      weights = loss_weights,
                      REML = FALSE)
anova(time_loss_alt_fit, time_loss_null_fit)
# p-value < 2.2e-16
# Writer loss is found to have greater loss than writer eraser loss

## Gene level data
gene_timepoint_dat <- timepoint_dat %>% 
  select(gene, assay, estimate) %>%
  spread(key = assay, value = estimate) 

##################################################################
##
##  Joint loss beta regression: Stan model
##
##################################################################
## Function to run brms for beta regression of loss assays
#### ChIP-seq normalized to proportions within a replicate
#### Per gene modeling of loss categories
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
# Try YAL067W-A
yal067w_a_compare_model <- compare_wl_to_wel_brms(chipseq_dat = gene_chip_dat,
                                                  this_gene = "YAL067W-A",
                                                  iter = 10000)

#####################################################
##
##  Joint loss beta regression: aggregate results
##
#####################################################
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
#############################################################
##
## Joint loss beta regression: characterizing lag in trend
##
#############################################################
#### At timepoints
## timepoint 0
loss_compare_dat %>% 
  filter(parameter == "assaywriter_eraser_loss") %>%
  mutate(check = Estimate > 0) %>%
  pull(check) %>%
  mean # 0.6255844
## timepoint 1
loss_compare_dat %>% 
  filter(parameter == "timepoint_cat1:assaywriter_eraser_loss") %>%
  mutate(check = Estimate > 0) %>%
  pull(check) %>%
  mean # 0.9130353
## timepoint 2
loss_compare_dat %>% 
  filter(parameter == "timepoint_cat2:assaywriter_eraser_loss") %>%
  mutate(check = Estimate > 0) %>%
  pull(check) %>%
  mean # 0.9762484
## timepoint 3
loss_compare_dat %>% 
  filter(parameter == "timepoint_cat3:assaywriter_eraser_loss") %>%
  mutate(check = Estimate > 0) %>%
  pull(check) %>%
  mean # 0.9064896
#### Transitions
## timepoint 0 to timepoint 1
loss_compare_dat %>%
  select(gene, parameter, Estimate) %>%
  filter(parameter %in% c("assaywriter_eraser_loss", "timepoint_cat1:assaywriter_eraser_loss")) %>%
  spread(key = parameter, value = Estimate) %>%
  mutate(check = `timepoint_cat1:assaywriter_eraser_loss` > assaywriter_eraser_loss) %>%
  pull(check) %>%
  mean # 0.8913409
## timepoint 1 to timepoint 2
loss_compare_dat %>%
  select(gene, parameter, Estimate) %>%
  filter(parameter %in% c("timepoint_cat2:assaywriter_eraser_loss", "timepoint_cat1:assaywriter_eraser_loss")) %>%
  spread(key = parameter, value = Estimate) %>%
  mutate(check = `timepoint_cat2:assaywriter_eraser_loss` > `timepoint_cat1:assaywriter_eraser_loss`) %>%
  pull(check) %>%
  mean # 0.7637928
## timepoint 2 to timepoint 3
loss_compare_dat %>%
  select(gene, parameter, Estimate) %>%
  filter(parameter %in% c("timepoint_cat3:assaywriter_eraser_loss", "timepoint_cat2:assaywriter_eraser_loss")) %>%
  spread(key = parameter, value = Estimate) %>%
  mutate(check = `timepoint_cat2:assaywriter_eraser_loss` > `timepoint_cat3:assaywriter_eraser_loss`) %>%
  pull(check) %>%
  mean # 0.9452029
## timepoint 0 to timepoint 3
loss_compare_dat %>%
  select(gene, parameter, Estimate) %>%
  filter(parameter %in% c("timepoint_cat3:assaywriter_eraser_loss", "assaywriter_eraser_loss")) %>%
  spread(key = parameter, value = Estimate) %>%
  mutate(check = assaywriter_eraser_loss < `timepoint_cat3:assaywriter_eraser_loss`) %>%
  pull(check) %>%
  mean # 0.8782495
## Grab genes with "significant" lag
lag_dat <- loss_compare_dat %>%
  filter(grepl(x = parameter, pattern = "assay")) %>%
  mutate(parameter = recode(parameter, 
                            'assaywriter_eraser_loss' = "0",
                            'timepoint_cat1:assaywriter_eraser_loss' = "1",
                            'timepoint_cat2:assaywriter_eraser_loss' = "2",
                            'timepoint_cat3:assaywriter_eraser_loss' = "3")) %>%
  rename(timepoint = parameter) %>%
  mutate(lag = l.95..CI > 0)
## Looking at genes with consecutive lags
total_lag_dat <- lag_dat %>%
  group_by(gene) %>%
  summarize(lag_count = sum(lag)) %>%
  ungroup
## Grab genes with some lag between WL and WEL
lagged_genes <- total_lag_dat %>%
  filter(lag_count > 0) %>%
  pull(gene)
## Plot genes with significant lags between WL and WEL
g <- ggplot(data = gene_chip_dat %>% filter(assay != "writer_add", gene %in% lagged_genes), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
#g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + guides(color = FALSE)
g
# Mostly detect ugly things, don't have power at the single gene level
#############################################################
##
## Joint loss beta regression: overall picture of lag
##
#############################################################
## Processing data
plot_loss_compare_dat <- loss_compare_dat %>%
  filter(grepl(x = parameter, pattern = "assay")) %>%
  mutate(parameter = recode(parameter, 
                            'assaywriter_eraser_loss' = "0",
                            'timepoint_cat1:assaywriter_eraser_loss' = "1",
                            'timepoint_cat2:assaywriter_eraser_loss' = "2",
                            'timepoint_cat3:assaywriter_eraser_loss' = "3")) %>%
  rename(timepoint = parameter)
## Plot
g <- ggplot(data = plot_loss_compare_dat %>% filter(gene %in% negative_genes), aes(x = timepoint, y = Estimate)) + geom_line(aes(group = gene, col = gene), alpha = 0.2) + scale_color_grey(guide = FALSE)
g <- g + geom_boxplot(aes(x = timepoint, y = Estimate), outlier.alpha = 0, color = wel_col, size = 0.7, alpha = 0.3) + ylab("WEL - WL timepoint effects")
g <- g + stat_summary(fun.y = mean, geom="line", aes(group=1), col = wel_col, size = 1)
g <- g + gg_theme# + guides(fill=FALSE)
g <- g + geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)
g



##############################################
##
## Comparison of histone mark loss
## with additional factors (covariates)
##
##############################################
## Transcript levels and gene length
wa_covar_raw_dat <- data.table::fread("data/DtL_setd2_RNA_combined_sf.txt", data.table = FALSE) 
wa_covar_dat <- wa_covar_raw_dat %>%
  rename(gene_length = "Gene Length",
         gene = GENE) %>%
  select(gene, gene_length, DtL_0min_set2d_RNA_Rep1, DtL_0min_set2d_RNA_Rep2, DtL_0min_set2d_RNA_Rep3) %>%
  rename(rep1 = DtL_0min_set2d_RNA_Rep1,
         rep2 = DtL_0min_set2d_RNA_Rep2,
         rep3 = DtL_0min_set2d_RNA_Rep3) %>%
  gather(key = rep, value = transcript, -c(gene, gene_length)) %>%
  group_by(gene, gene_length) %>%
  summarize(transcript = mean(log(transcript + 1))) %>%
  ungroup

wl_covar_raw_dat <- data.table::fread("data/LtD_setd2_RNA_combined_sf.txt", data.table = FALSE)
wl_covar_dat <- wl_covar_raw_dat %>%
  rename(gene_length = "Gene Length",
         gene = GENE) %>%
  select(gene, gene_length, LtD_0min_set2d_RNA_Rep1, LtD_0min_set2d_RNA_Rep2, LtD_0min_set2d_RNA_Rep3) %>%
  rename(rep1 = LtD_0min_set2d_RNA_Rep1,
         rep2 = LtD_0min_set2d_RNA_Rep2,
         rep3 = LtD_0min_set2d_RNA_Rep3) %>%
  gather(key = rep, value = transcript, -c(gene, gene_length)) %>%
  group_by(gene, gene_length) %>%
  summarize(transcript = mean(log(transcript + 1))) %>%
  ungroup

## Absolute H3K36me3 levels
wa_histone_covar <- data.table::fread("data/absolute_matrix_for_LMM_model_nonOverlappingGenes_noChrMGenes.txt", data.table = FALSE) %>%
  filter(V4 == 0,
         V2 == 3) %>%
  select(V1, V5) %>%
  group_by(V1) %>%
  summarize(V5 = mean(V5)) %>%
  ungroup %>%
  rename(gene = V1,
         histone = V5)
wl_histone_covar <- data.table::fread("data/absolute_matrix_for_LMM_model_nonOverlappingGenes_noChrMGenes.txt", data.table = FALSE) %>%
  filter(V4 == 0,
         V2 == 1) %>%
  select(V1, V5) %>%
  group_by(V1) %>%
  summarize(V5 = mean(V5)) %>%
  ungroup %>%
  rename(gene = V1,
         histone = V5)
wel_histone_covar <- data.table::fread("data/absolute_matrix_for_LMM_model_nonOverlappingGenes_noChrMGenes.txt", data.table = FALSE) %>%
  filter(V4 == 0,
         V2 == 1) %>%
  select(V1, V5) %>%
  group_by(V1) %>%
  summarize(V5 = mean(V5)) %>%
  ungroup %>%
  rename(gene = V1,
         histone = V5)
## Merging transcript and histone mark data together
### Do not have transcripts for WEL
wa_full_dat <- timepoint_dat %>%
  filter(assay == "writer_add") %>%
  left_join(wa_covar_dat) %>%
  left_join(wa_histone_covar)
wl_full_dat <- timepoint_dat %>%
  filter(assay == "writer_loss") %>%
  left_join(wl_covar_dat) %>%
  left_join(wl_histone_covar)
## Reduce to only genes with trends in all assays
wa_reduced_dat <- wa_full_dat %>% filter(gene %in% trend_genes)
wl_reduced_dat <- wl_full_dat %>% filter(gene %in% trend_genes)

## Time 0 transcript by time 0 histone marks
### Full data
par(mfrow = c(1, 2))
plot(x = wl_full_dat$hist, y = wl_full_dat$transcript, 
     xlab = "Time 0 mean histone marks", ylab = "Time 0 mean transcription", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.4), las = 1)
plot(wa_full_dat$hist, wa_full_dat$transcript, 
     xlab = "Time 0 mean histone marks", ylab = "Time 0 mean transcription", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.4), las = 1)
### Trend genes only
par(mfrow = c(1, 2))
plot(x = wl_reduced_dat$hist, y = wl_reduced_dat$transcript, 
     xlab = "Time 0 mean histone marks", ylab = "Time 0 mean transcription", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.5), las = 1)
plot(wa_reduced_dat$hist, wa_reduced_dat$transcript, 
     xlab = "Time 0 mean histone marks", ylab = "Time 0 mean transcription", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.5), las = 1)

## Histone mark trends by time 0 marks
### Full data
par(mfrow = c(2, 2))
plot(x = wa_full_dat$histone, y = wa_full_dat$estimate, 
     xlab = "Mean transcription at time 0", ylab = "Histone gain estimate", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.5), las = 1)
plot(x = wl_full_dat$histone, y = abs(wl_full_dat$estimate), 
     xlab = "Mean transcription at time 0", ylab = "Histone loss estimate", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.5), las = 1)
plot(x = wa_full_dat$histone, y = wa_full_dat$est_error, 
     xlab = "Mean transcription at time 0", ylab = "Error on histone gain estimate", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.5), las = 1)
plot(x = wl_full_dat$histone, y = abs(wl_full_dat$est_error), 
     xlab = "Mean transcription at time 0", ylab = "Error on histone loss estimate", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.5), las = 1)
### Trend genes only
par(mfrow = c(2, 2))
plot(x = wa_reduced_dat$histone, y = wa_reduced_dat$estimate, 
     xlab = "Mean transcription at time 0", ylab = "Histone gain estimate", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.5), las = 1)
plot(x = wl_reduced_dat$histone, y = abs(wl_reduced_dat$estimate), 
     xlab = "Mean transcription at time 0", ylab = "Histone loss estimate", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.5), las = 1)
plot(x = wa_reduced_dat$histone, y = wa_reduced_dat$est_error, 
     xlab = "Mean transcription at time 0", ylab = "Error on histone gain estimate", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.5), las = 1)
plot(x = wl_reduced_dat$histone, y = abs(wl_reduced_dat$est_error), 
     xlab = "Mean transcription at time 0", ylab = "Error on histone loss estimate", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.5), las = 1)

## Histone mark trends by time 0 transcripts
### Full data
par(mfrow = c(2, 2))
plot(x = wa_full_dat$transcript, y = wa_full_dat$estimate, 
     xlab = "Mean time 0 transcript", ylab = "Histone gain estimate", 
     main = "Writer add", las = 1, 
     pch = 20, col = alpha(wa_col, 0.3))
abline(lm(estimate ~ transcript, data = wa_full_dat), col = "gray")
text(x = 4, y = 4, paste("r = ", 
                         round(cor(wa_full_dat$transcript, wa_full_dat$estimate), 3)),
     col = "goldenrod", font = 2, cex = 2)
plot(x = wl_full_dat$transcript, y = abs(wl_full_dat$estimate), 
     xlab = "Mean time 0 transcript", ylab = "Histone loss estimate", 
     main = "Writer loss", las = 1,
     pch = 20, col = alpha(wl_col, 0.3))
abline(lm(abs(estimate) ~ transcript, data = wl_full_dat), col = "gray")
text(x = 4, y = 3, paste("r = ", 
                         round(cor(wl_full_dat$transcript, abs(wl_full_dat$estimate)), 3)),
     col = "goldenrod", font = 2, cex =2)
plot(x = wa_full_dat$histone, y = wa_full_dat$est_error, 
     xlab = "Mean time 0 transcript", ylab = "Error on histone loss estimate", 
     main = "Writer add", las = 1,
     pch = 20, col = alpha(wa_col, 0.3))
plot(x = wl_full_dat$histone, y = abs(wl_full_dat$est_error), 
     xlab = "Mean time0 H3K36me3", ylab = "Error on histone loss estimate", 
     main = "Writer loss", las = 1,
     pch = 20, col = alpha(wl_col, 0.3))
cor.test(x = wa_full_dat$transcript, y = wa_full_dat$estimate)
cor.test(x = wl_full_dat$transcript, y = wl_full_dat$estimate)
### Trend genes only
par(mfrow = c(2, 2))
plot(x = wa_reduced_dat$transcript, y = wa_reduced_dat$estimate, 
     xlab = "Mean time 0 transcript", ylab = "Histone gain estimate", 
     main = "Writer add", las = 1, 
     pch = 20, col = alpha(wa_col, 0.3))
abline(lm(estimate ~ transcript, data = wa_reduced_dat), col = "gray")
text(x = 4, y = 4, paste("r = ", 
                         round(cor(wa_reduced_dat$transcript, wa_reduced_dat$estimate), 3)),
     col = "goldenrod", font = 2, cex = 2)
plot(x = wl_reduced_dat$transcript, y = abs(wl_reduced_dat$estimate), 
     xlab = "Mean time 0 transcript", ylab = "Histone loss estimate", 
     main = "Writer loss", las = 1,
     pch = 20, col = alpha(wl_col, 0.3))
abline(lm(abs(estimate) ~ transcript, data = wl_reduced_dat), col = "gray")
text(x = 4, y = 3, paste("r = ", 
                         round(cor(wl_reduced_dat$transcript, abs(wl_reduced_dat$estimate)), 3)),
     col = "goldenrod", font = 2, cex =2)
plot(x = wa_reduced_dat$histone, y = wa_reduced_dat$est_error, 
     xlab = "Mean time 0 transcript", ylab = "Error on histone loss estimate", 
     main = "Writer add", las = 1,
     pch = 20, col = alpha(wa_col, 0.3))
plot(x = wl_reduced_dat$histone, y = abs(wl_reduced_dat$est_error), 
     xlab = "Mean time0 H3K36me3", ylab = "Error on histone loss estimate", 
     main = "Writer loss", las = 1,
     pch = 20, col = alpha(wl_col, 0.3))
cor.test(x = wa_reduced_dat$transcript, y = wa_reduced_dat$estimate)
cor.test(x = wl_reduced_dat$transcript, y = wl_reduced_dat$estimate)

## Histone mark trends by gene length
### Full data
par(mfrow = c(2, 2))
plot(wa_full_dat$gene_length, wa_full_dat$estimate, 
     xlab = "Gene length", ylab = "Trend with time", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.3), las = 1)
plot(wl_full_dat$gene_length, wl_full_dat$estimate, 
     xlab = "Gene length", ylab = "Trend with time", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.3), las = 1)
plot(wa_full_dat$gene_length, wa_full_dat$est_error,
     xlab = "Gene length", ylab = "Error on trend with time", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.3), las = 1)
plot(wl_full_dat$gene_length, wl_full_dat$est_error, 
     xlab = "Gene length", ylab = "Error on trend with time", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.3), las = 1)
### Trend genes only
par(mfrow = c(2, 2))
plot(wa_reduced_dat$gene_length, wa_reduced_dat$estimate, 
     xlab = "Gene length", ylab = "Trend with time", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.3), las = 1)
plot(wl_reduced_dat$gene_length, wl_reduced_dat$estimate, 
     xlab = "Gene length", ylab = "Trend with time", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.3), las = 1)
plot(wa_reduced_dat$gene_length, wa_reduced_dat$est_error,
     xlab = "Gene length", ylab = "Error on trend with time", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.3), las = 1)
plot(wl_reduced_dat$gene_length, wl_reduced_dat$est_error, 
     xlab = "Gene length", ylab = "Error on trend with time", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.3), las = 1)

## Gene length by Time 0 transcript
par(mfrow = c(2, 2))
plot(wa_full_dat$histone, wa_full_dat$gene_length, 
     xlab = "Mean H3K36me3 at time 0", ylab = "Gene length", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.3), las = 1)
plot(wl_full_dat$histone, wl_full_dat$gene_length, 
     xlab = "Mean H3K36me3 at time 0", ylab = "Gene length", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.3), las = 1)
plot(wa_reduced_dat$histone, wa_reduced_dat$gene_length, 
     xlab = "Mean H3K36me3 at time 0", ylab = "Gene length", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.3), las = 1)
plot(wl_reduced_dat$histone, wl_reduced_dat$gene_length, 
     xlab = "Mean H3K36me3 at time 0", ylab = "Gene length", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.3), las = 1)

## Gene length, histone mark trend, and time 0 transcript
g1 <- ggplot(data = wa_full_dat, aes(x = gene_length, y = estimate)) + geom_point(aes(color = transcript)) + ggtitle("Writer add") + gg_theme
g2 <- ggplot(data = wl_full_dat, aes(x = gene_length, y = estimate)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss") + gg_theme
grid.arrange(g1, g2)
g1 <- ggplot(data = wa_full_dat, aes(x = gene_length, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer add") + gg_theme
g2 <- ggplot(data = wl_full_dat, aes(x = gene_length, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss") + gg_theme
grid.arrange(g1, g2)

## Time 0 histone, histone mark trend, and time 0 transcript
g1 <- ggplot(data = wa_full_dat, aes(x = histone, y = estimate)) + geom_point(aes(color = transcript)) + ggtitle("Writer add") + gg_theme
g2 <- ggplot(data = wl_full_dat, aes(x = histone, y = estimate)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss") + gg_theme
grid.arrange(g1, g2)
g1 <- ggplot(data = wa_full_dat, aes(x = histone, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer add") + gg_theme
g2 <- ggplot(data = wl_full_dat, aes(x = histone, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss") + gg_theme
grid.arrange(g1, g2)

###################################################
##
##  Read in and process the data
##
##  Transcripts
##
###################################################
# WA
wa_tx_raw_dat <- data.table::fread("data/DtL_setd2_RNA_combined_sf.txt", data.table = FALSE) %>%
  rename(gene_length = "Gene Length",
         gene = GENE) %>%
  gather(key = "category", value = "transcript", -c(gene, gene_length)) %>%
  separate(category, c("V1", "timepoint", "V2", "V3", "replicate")) %>%
  select(-c(V1, V2, V3)) %>%
  mutate(timepoint = as.numeric(gsub(x = timepoint, pattern = "min", replacement = ""))/60,
         replicate = as.numeric(gsub(x = replicate, pattern = "Rep", replacement = "")),
         assay = 2) %>%
  arrange(gene, timepoint, replicate)
# WL
wl_tx_raw_dat <- data.table::fread("data/LtD_setd2_RNA_combined_sf.txt", data.table = FALSE) %>%
  rename(gene_length = "Gene Length",
         gene = GENE) %>%
  gather(key = "category", value = "transcript", -c(gene, gene_length)) %>%
  separate(category, c("V1", "timepoint", "V2", "V3", "replicate")) %>%
  select(-c(V1, V2, V3)) %>%
  mutate(timepoint = as.numeric(gsub(x = timepoint, pattern = "min", replacement = ""))/60,
         replicate = as.numeric(gsub(x = replicate, pattern = "Rep", replacement = "")),
         assay = 1) %>%
  arrange(gene, timepoint, replicate)
# No transcription data for WEL
## Bind together
tx_raw_dat <- bind_rows(wa_tx_raw_dat,
                        wl_tx_raw_dat) %>%
  mutate(assay = factor(assay),
         assay = recode(assay, "1" = "writer_loss", "2" = "writer_add"),
         sample_unit = paste(gene, assay, replicate, sep = "_"))

##############################################
##
##  Raw data plots
##  Compare transcripts with histone marks
##
##############################################
overlap_genes <- intersect(unique(gene_chip_dat$gene), unique(tx_raw_dat$gene))

de_genes <- c("YAL037W", "YAL012W", "YBL106C", "YBL098W", "YBL078C")
# Transcripts
g_transcript <- ggplot(data = tx_raw_dat %>% filter(gene %in% de_genes), aes(x = timepoint, y = log(transcript + 1), col = assay)) + scale_color_manual(values = c(wl_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g_transcript <- g_transcript + geom_smooth(aes(y = log(transcript + 1), x = timepoint), method = "lm", size = 2)
g_transcript <- g_transcript + gg_theme + facet_wrap(~gene)
# Histone marks
g_histone <- ggplot(data = gene_chip_dat %>% filter(gene %in% de_genes), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g_histone <- g_histone + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g_histone <- g_histone + gg_theme + facet_wrap(~gene)
grid.arrange(g_transcript, g_histone)

# YFR044C
g <- ggplot(data = tx_raw_dat %>% filter(gene %in% "YFR044C"), aes(x = timepoint, y = log(transcript + 1), col = assay)) + scale_color_manual(values = c(wl_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = log(transcript + 1), x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + facet_wrap(~gene)
g
# YAL037W
g <- ggplot(data = tx_raw_dat %>% filter(gene %in% "YAL037W"), aes(x = timepoint, y = log(transcript + 1), col = assay)) + scale_color_manual(values = c(wl_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = log(transcript + 1), x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + facet_wrap(~gene)
g
# YPR191W
g <- ggplot(data = tx_raw_dat %>% filter(gene %in% "YPR191W"), aes(x = timepoint, y = log(transcript + 1), col = assay)) + scale_color_manual(values = c(wl_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = log(transcript + 1), x = timepoint), method = "lm", size = 2)
g <- g + gg_theme + facet_wrap(~gene)
g

###################################################
##
##  Read in and process the data
##
##  Raw histones data
##
###################################################
histone_raw_dat <- data.table::fread("data/absolute_matrix_for_LMM_model_nonOverlappingGenes_noChrMGenes.txt", data.table = FALSE) %>%
  rename(gene = V1,
         assay = V2,
         replicate = V3,
         timepoint = V4,
         histone = V5) %>%
  mutate(assay = factor(assay),
         assay = recode(assay, "1" = "writer_loss", "2" = "writer_eraser_loss", "3" = "writer_add"),
         timepoint_cat = factor(timepoint),
         timepoint = 0,
         sample_unit = paste(gene, assay, replicate, sep = "_")) 
histone_raw_dat$timepoint[histone_raw_dat$assay %in% c("writer_loss", "writer_eraser_loss") & histone_raw_dat$timepoint_cat == 0] <- 0
histone_raw_dat$timepoint[histone_raw_dat$assay %in% c("writer_loss", "writer_eraser_loss") & histone_raw_dat$timepoint_cat == 1] <- 30/60
histone_raw_dat$timepoint[histone_raw_dat$assay %in% c("writer_loss", "writer_eraser_loss") & histone_raw_dat$timepoint_cat == 2] <- 60/60
histone_raw_dat$timepoint[histone_raw_dat$assay %in% c("writer_loss", "writer_eraser_loss") & histone_raw_dat$timepoint_cat == 3] <- 90/60
histone_raw_dat$timepoint[histone_raw_dat$assay == "writer_add" & histone_raw_dat$timepoint_cat == 0] <- 0
histone_raw_dat$timepoint[histone_raw_dat$assay == "writer_add" & histone_raw_dat$timepoint_cat == 1] <- 20/60
histone_raw_dat$timepoint[histone_raw_dat$assay == "writer_add" & histone_raw_dat$timepoint_cat == 2] <- 40/60
histone_raw_dat$timepoint[histone_raw_dat$assay == "writer_add" & histone_raw_dat$timepoint_cat == 3] <- 60/60

## Merge transcript and histone data
tx_histone_merge_dat <- tx_raw_dat %>% 
  select(gene, assay, replicate, timepoint, sample_unit, transcript) %>%
  inner_join(histone_raw_dat %>%
               select(sample_unit, timepoint, histone))

##############################################
##
##  Log-normal regression: Stan model
##
##############################################
## Function to run brms for log-normal regression
#### Per gene per assay
run_brms_on_transcript <- function(transcript_dat, 
                                   this_gene, 
                                   adapt_delta = 0.8,
                                   iter = 2000) {
  use_dat <- transcript_dat %>% filter(gene == this_gene)#,
  # writer loss
  wl_fit <- brms::brm(transcript ~ 1 + timepoint + histone + (1 + timepoint + histone | sample_unit), 
                      data = use_dat %>% filter(assay == "writer_loss"), 
                      family = "hurdle_lognormal",
                      iter = iter,
                      control = list(adapt_delta = adapt_delta))
  wl_fixed <- summary(wl_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
  wl_random <- summary(wl_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
  # writer add
  wa_fit <- brms::brm(transcript ~ 1 + timepoint + histone + (1 + timepoint + histone | sample_unit), 
                      data = use_dat %>% filter(assay == "writer_add"), 
                      family = "hurdle_lognormal",
                      iter = iter,
                      control = list(adapt_delta = adapt_delta))
  wa_fixed <- summary(wa_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
  wa_random <- summary(wa_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
  
  names(wl_fixed) <- names(wl_random) <- names(wa_fixed) <- names(wa_random) <- c("parameter", "estimate", "est_error", "lower_95", "upper_95", "eff_sample", "rhat")
  results <- list(writer_loss = bind_rows(wl_fixed, wl_random),
                  writer_add = bind_rows(wa_fixed, wa_random))
  results
}

ypr191w_models = run_brms_on_transcript(transcript_dat = tx_histone_merge_dat,
                                        this_gene = "YPR191W",
                                        iter = 10000)

################################################
##
##  Transcription regression: aggregate results
##
################################################
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






## Un-scaled histone
raw_unscaled_gene_chip_dat <- read.table("data/absolute_matrix_for_LMM_model_nonOverlappingGenes.txt") # path to data

unscaled_gene_chip_dat <- raw_unscaled_gene_chip_dat %>%
  rename(gene = V1,
         assay = V2,
         replicate = V3,
         timepoint = V4,
         value = V5) %>%
  mutate(assay = factor(assay),
         assay = recode(assay, "1" = "writer_loss", "2" = "writer_eraser_loss", "3" = "writer_add"),
         sample_unit = paste(gene, assay, replicate, sep = "_"),
         timepoint_cat = factor(timepoint),
         timepoint = 0,
         transformed_value = round(value * 1000))
unscaled_gene_chip_dat$timepoint[unscaled_gene_chip_dat$assay %in% c("writer_loss", "writer_eraser_loss") & unscaled_gene_chip_dat$timepoint_cat == 0] <- 0
unscaled_gene_chip_dat$timepoint[unscaled_gene_chip_dat$assay %in% c("writer_loss", "writer_eraser_loss") & unscaled_gene_chip_dat$timepoint_cat == 1] <- 30/60
unscaled_gene_chip_dat$timepoint[unscaled_gene_chip_dat$assay %in% c("writer_loss", "writer_eraser_loss") & unscaled_gene_chip_dat$timepoint_cat == 2] <- 60/60
unscaled_gene_chip_dat$timepoint[unscaled_gene_chip_dat$assay %in% c("writer_loss", "writer_eraser_loss") & unscaled_gene_chip_dat$timepoint_cat == 3] <- 90/60
unscaled_gene_chip_dat$timepoint[unscaled_gene_chip_dat$assay == "writer_add" & unscaled_gene_chip_dat$timepoint_cat == 0] <- 0
unscaled_gene_chip_dat$timepoint[unscaled_gene_chip_dat$assay == "writer_add" & unscaled_gene_chip_dat$timepoint_cat == 1] <- 20/60
unscaled_gene_chip_dat$timepoint[unscaled_gene_chip_dat$assay == "writer_add" & unscaled_gene_chip_dat$timepoint_cat == 2] <- 40/60
unscaled_gene_chip_dat$timepoint[unscaled_gene_chip_dat$assay == "writer_add" & unscaled_gene_chip_dat$timepoint_cat == 3] <- 60/60

unscaled_gene_chip_dat %>% filter(assay == "writer_add",
                                  gene == unscaled_gene_chip_dat$gene[1])

g <- ggplot(data = unscaled_gene_chip_dat %>% filter(gene == unscaled_gene_chip_dat$gene[1]), aes(x = timepoint, y = value * 1000, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme
g

g <- ggplot(data = gene_chip_dat %>% filter(gene == unscaled_gene_chip_dat$gene[1]), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme
g

hist(unscaled_gene_chip_dat$value * 1000)
hist(unscaled_gene_chip_dat$transformed_value)

## Arbitrary, but use transformed_value = round(value * 1000)
run_brms_on_raw_chipseq <- function(chipseq_dat, 
                                    this_gene, 
                                    adapt_delta = 0.8,
                                    iter = 2000) {
  use_dat <- chipseq_dat %>% filter(gene == this_gene)#,
  #!(value %in% c(0, 1)))
  # writer loss
  wl_fit <- brms::brm(transformed_value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                      data = use_dat %>% filter(assay == "writer_loss"), 
                      family = "zero_inflated_negbinomial",
                      iter = iter,
                      control = list(adapt_delta = adapt_delta))
  wl_fixed <- summary(wl_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
  wl_random <- summary(wl_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
  # writer-eraser loss
  wel_fit <- brms::brm(transformed_value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                       data = use_dat %>% filter(assay == "writer_eraser_loss"), 
                       family = "zero_inflated_negbinomial",
                       iter = iter,
                       control = list(adapt_delta = adapt_delta))
  wel_fixed <- summary(wel_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
  wel_random <- summary(wel_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
  # writer add
  wa_fit <- brms::brm(transformed_value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                      data = use_dat %>% filter(assay == "writer_add"), 
                      family = "zero_inflated_negbinomial",
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

g <- ggplot(data = unscaled_gene_chip_dat %>% filter(gene == "YAL012W"), aes(x = timepoint, y = transformed_value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = transformed_value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme
g

g <- ggplot(data = gene_chip_dat %>% filter(gene == "YAL012W"), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme
g

raw_yal012w_models <- run_brms_on_chipseq(chipseq_dat = unscaled_gene_chip_dat, this_gene = "YAL012W", iter = 10000)


g <- ggplot(data = unscaled_gene_chip_dat %>% filter(gene == "YAL066W"), aes(x = timepoint, y = transformed_value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = transformed_value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme
g

g <- ggplot(data = gene_chip_dat %>% filter(gene == "YAL066W"), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + gg_theme
g

g1 <- ggplot(data = unscaled_gene_chip_dat %>% filter(gene == "YAL041W"), aes(x = timepoint, y = transformed_value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g1 <- g1 + geom_smooth(aes(y = transformed_value, x = timepoint), method = "lm", size = 2)
g1 <- g1 + gg_theme
g1

g2 <- ggplot(data = gene_chip_dat %>% filter(gene == "YAL041W"), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g2 <- g2 + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g2 <- g2 + gg_theme
g2

gene_chip_dat %>% filter(gene == "YAL041W", 
                         assay == "writer_eraser_loss")
unscaled_gene_chip_dat %>% filter(gene == "YAL041W", 
                         assay == "writer_eraser_loss")

grid.arrange(g1, g2)










####################################
##
##  Etc
##
####################################

####################################
##
## Permutations
##
####################################
all_perms <- gtools::permutations(n = 4, r = 4) - 1
perm_wl_effects <- perm_wel_effects <- perm_wa_effects <- rep(NA, nrow(all_perms) - 1)
for (i in 2:nrow(all_perms)) {
  perm_dat <- gene_chip_dat
  perm_dat$timepoint = all_perms[i, gene_chip_dat$timepoint + 1]
  
  perm_fit <- run_brms_on_chipseq(chipseq_dat = perm_dat, this_gene = positive_genes[1])
  perm_wl_effects[i - 1] <- perm_fit$writer_loss %>% 
    filter(parameter == "timepoint") %>%
    pull(estimate)
  perm_wel_effects[i - 1] <- perm_fit$writer_eraser_loss %>% 
    filter(parameter == "timepoint") %>%
    pull(estimate)
  perm_wa_effects[i - 1] <- perm_fit$writer_add %>% 
    filter(parameter == "timepoint") %>%
    pull(estimate)
}
actual_fit <- run_brms_on_chipseq(chipseq_dat = gene_chip_dat, this_gene = positive_genes[1])
mean(actual_fit$writer_add %>% filter(parameter == "timepoint") %>% pull(estimate) <= c(perm_wa_effects,
                                                                                        actual_fit$writer_add %>% filter(parameter == "timepoint") %>% pull(estimate)))
mean(actual_fit$writer_loss %>% filter(parameter == "timepoint") %>% pull(estimate) >= c(perm_wl_effects,
                                                                                         actual_fit$writer_add %>% filter(parameter == "timepoint") %>% pull(estimate)))
mean(actual_fit$writer_eraser_loss %>% filter(parameter == "timepoint") %>% pull(estimate) >= c(perm_wel_effects,
                                                                                                actual_fit$writer_add %>% filter(parameter == "timepoint") %>% pull(estimate)))
####################################
##
## Multi-gene model
##
####################################
multigene_fit <- brms::brm(value ~ 1 + timepoint_cat + assay + timepoint_cat:assay + (1 + timepoint_cat + assay + timepoint_cat:assay | sample_unit + gene), 
                           data = gene_chip_dat %>% filter(assay %in% c("writer_loss", "writer_eraser_loss")), 
                           family = "beta",
                           iter = 2000,
                           control = list(adapt_delta = 0.8))
# Seems intractible for running


