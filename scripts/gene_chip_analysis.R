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
wa_col <- "plum"
wl_col <- "darkseagreen2"
wel_col <- "tomato"
# wa_col <- "coral"
# wl_col <- "magenta"
# wel_col <- "seagreen1"

## Set directory for local environment to the git repos
#setwd("~/projects/chipseq-gene-dynamics/")
setwd("~/Documents/git_repositories/chipseq-gene-dynamics/")

## ggplot theme
gg_theme <- theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  axis.line = element_line(colour = "black"),
                  plot.title = element_text(hjust = 0.5, face = "bold"), 
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
         raw_value = value,
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
  filter(gene %in% gene_chip_dat$gene[1]) %>%
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
run_beta_brms_on_chipseq <- function(chipseq_dat, 
                                     this_gene, 
                                     adapt_delta = 0.8,
                                     iter = 2000,
                                     zo_inflation = TRUE,
                                     seed = 123,
                                     rhat_cutoff = 1.1,
                                     run_limit = 10) {
  use_dat <- chipseq_dat %>% filter(gene == this_gene)#,
  
  stop_now <- FALSE
  num_runs <- 0
  while(!stop_now) {
    # writer loss
    wl_fit <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                        data = use_dat %>% filter(assay == "writer_loss"), 
                        family = ifelse(zo_inflation, "zero_one_inflated_beta", "beta"),
                        iter = iter,
                        control = list(adapt_delta = adapt_delta),
                        seed = seed)
    wl_fixed <- summary(wl_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
    wl_random <- summary(wl_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
    # writer-eraser loss
    wel_fit <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                         data = use_dat %>% filter(assay == "writer_eraser_loss"), 
                         family = ifelse(zo_inflation, "zero_one_inflated_beta", "beta"),
                         iter = iter,
                         control = list(adapt_delta = adapt_delta),
                         seed = seed)
    wel_fixed <- summary(wel_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
    wel_random <- summary(wel_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
    # writer add
    wa_fit <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                        data = use_dat %>% filter(assay == "writer_add"), 
                        family = ifelse(zo_inflation, "zero_one_inflated_beta", "beta"),
                        iter = iter,
                        control = list(adapt_delta = adapt_delta), 
                        seed = seed)
    wa_fixed <- summary(wa_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
    wa_random <- summary(wa_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
    
    num_runs <- num_runs + 1
    if ((all(wl_fixed$Rhat < rhat_cutoff) & all(wel_fixed$Rhat < rhat_cutoff) & all(wa_fixed$Rhat < rhat_cutoff)) | num_runs == run_limit) {
      stop_now <- TRUE
    }
    else {
      seed <- seed + 1
    }
  }
  
  names(wl_fixed) <- names(wl_random) <- names(wel_fixed) <- names(wel_random) <- names(wa_fixed) <- names(wa_random) <- c("parameter", "estimate", "est_error", "lower_95", "upper_95", "eff_sample", "rhat")
  results <- list(writer_loss = bind_rows(wl_fixed, wl_random),
                  writer_eraser_loss = bind_rows(wel_fixed, wel_random),
                  writer_add = bind_rows(wa_fixed, wa_random),
                  seed = seed,
                  run_limit_stop = ifelse(num_runs == run_limit, TRUE, FALSE),
                  zo_inflation = zo_inflation)
  results
}
## Try YAL012W and YAL066W
yal012w_beta_models <- run_beta_brms_on_chipseq(chipseq_dat = gene_chip_dat, 
                                                this_gene = "YAL012W", 
                                                iter = 10000,
                                                zo_inflation = FALSE)
yal066w_beta_models <- run_beta_brms_on_chipseq(chipseq_dat = gene_chip_dat, 
                                                this_gene = "YAL066W",
                                                iter = 10000,
                                                zo_inflation = FALSE)
## Gene with extreme low effect
yal041w_beta_models <- run_beta_brms_on_chipseq(chipseq_dat = gene_chip_dat, 
                                                this_gene = "YAL041W", 
                                                iter = 5000,
                                                zo_inflation = FALSE)

zoi_gene_chip_dat <- gene_chip_dat %>%
  mutate(value = raw_value)
## Try YAL012W
yal012w_zoibeta_models <- run_beta_brms_on_chipseq(chipseq_dat = zoi_gene_chip_dat, 
                                                   this_gene = "YAL012W", 
                                                   iter = 10000,
                                                   zo_inflation = TRUE)


##############################################
##
##  Beta regression: aggregate results
##
##############################################
## Parsing full set of Stan results
stan_fit_paths <- list.files("individual_gene_Stan_models_alpha_0.8_iter_10000/", full.names = TRUE)
stan_fit_paths <- stan_fit_paths[!(grepl(x = stan_fit_paths, pattern = "WL_vs_WEL"))]
file_path <- "individual_gene_Stan_models_alpha_0.8_iter_10000//sacCer3_nonOverlappingGenes_noChrMGenes_"
for (i in 1:length(stan_fit_paths)) {
  this_fit <- readRDS(stan_fit_paths[i])

  this_gene <- gsub(x = stan_fit_paths[i], pattern = file_path, replacement = "", fixed = TRUE) %>%
    gsub(x = ., replacement = "", pattern = "_Stan.rds", fixed = TRUE)
  
  
  if (i == 1) {
    timepoint_dat <- bind_rows(data.frame(gene = this_gene, assay = "writer_loss", this_fit$writer_loss %>% filter(parameter == "timepoint")),
                               data.frame(gene = this_gene, assay = "writer_eraser_loss", this_fit$writer_eraser_loss %>% filter(parameter == "timepoint")),
                               data.frame(gene = this_gene, assay = "writer_add", this_fit$writer_add %>% filter(parameter == "timepoint")))
    timepoint_dat$seed <- this_fit$seed
    timepoint_dat$run_limit_stop <- this_fit$run_limit_stop
  }
  else {
    holder_dat <- bind_rows(data.frame(gene = this_gene, assay = "writer_loss", this_fit$writer_loss %>% filter(parameter == "timepoint")),
                            data.frame(gene = this_gene, assay = "writer_eraser_loss", this_fit$writer_eraser_loss %>% filter(parameter == "timepoint")),
                            data.frame(gene = this_gene, assay = "writer_add", this_fit$writer_add %>% filter(parameter == "timepoint")))
    holder_dat$seed <- this_fit$seed
    holder_dat$run_limit_stop <- this_fit$run_limit_stop
    timepoint_dat <- bind_rows(timepoint_dat,
                               holder_dat)
  }
}
timepoint_dat <- timepoint_dat %>% 
  select(-parameter) %>%
  mutate(is_pos = ifelse(lower_95 > 0, TRUE, FALSE),
         is_neg = ifelse(upper_95 < 0, TRUE, FALSE))
timepoint_dat$category[timepoint_dat$lower_95 > 0] <- "positive"
timepoint_dat$category[timepoint_dat$upper_95 < 0] <- "negative"
timepoint_dat$category[timepoint_dat$lower_95 < 0 & timepoint_dat$upper_95 > 0] <- "zero"
timepoint_dat$attempts <- (timepoint_dat$seed - 123) + 1

## Histogram
ggplot(data = timepoint_dat, aes(x = estimate)) + geom_histogram() + gg_theme + facet_grid(assay~.)
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
genes_pos_wa <- timepoint_dat %>% 
  filter(assay == "writer_add",
         category == "positive") %>%
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

g <- ggplot(data = gene_chip_dat %>% filter(gene %in% "Q0032"), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
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
#color circle if non-zero estimate within each assay
plot(timepoint_dat$est_error, timepoint_dat$estimate, pch = ifelse(timepoint_dat$category == "zero", 1, 19), col = c(wa_col, wl_col, wel_col)[as.factor(timepoint_dat$assay)], 
     las = 1, ylab = "Histone trend with time", xlab = "Error on trend", frame.plot = FALSE)
abline(h = 0, lty = 2)
legend("bottomright", 
       legend = c("writer add", "writer loss", "writer eraser loss", "confident within assay", "not confident within assay"),
       col = c(wa_col, wl_col, wel_col, "gray", "gray"),
	   pch=c(15, 15, 15, 19, 1),
       bty = "n")

#color circle if non-zero across all assays
plot(timepoint_dat$est_error, timepoint_dat$estimate, pch = ifelse(timepoint_dat$gene %in% trend_genes, 19, 1), col = c(wa_col, wl_col, wel_col)[as.factor(timepoint_dat$assay)], 
     las = 1, ylab = "Histone trend with time", xlab = "Error on trend", frame.plot = FALSE)
abline(h = 0, lty = 2)
legend("bottomright", 
       legend = c("writer add", "writer loss", "writer eraser loss", "confident within assay", "not confident within assay"),
       col = c(wa_col, wl_col, wel_col, "gray", "gray"),
	   pch=c(15, 15, 15, 19, 1),
       bty = "n")

##############################################
##
##  Beta regression: Formal test of assay 
##                   on trend
##
##############################################
time_alt_fit <- lmer(estimate ~ -1 + (1 | gene) + assay, 
                     data = timepoint_dat, 
                     weights = 1/timepoint_dat$est_error,
                     REML = FALSE)
time_null_fit <- lmer(estimate ~ (1 | gene), 
                      data = timepoint_dat, 
                      weights = 1/timepoint_dat$est_error,
                      REML = FALSE)
anova(time_alt_fit, time_null_fit)
# p-value < 2.2e-16
emmeans::emmeans(time_alt_fit, list(pairwise ~ assay), adjust = "tukey")

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
gene_error_data <- timepoint_dat %>%
  select(gene, assay, est_error) %>%
  spread(key = assay, value = est_error) %>%
  rename(writer_loss_error = writer_loss,
         writer_eraser_loss_error = writer_eraser_loss,
         writer_add_error = writer_add)
gene_timepoint_dat <- bind_cols(gene_timepoint_dat, gene_error_data)

##################################################################
##
##  Joint loss beta regression: Stan model
##
##################################################################
## Function to run brms for beta regression of loss assays
#### ChIP-seq normalized to proportions within a replicate
#### Per gene modeling of loss categories
compare_wl_to_wel_beta_brms <- function(chipseq_dat, 
                                         this_gene, 
                                         adapt_delta = 0.8,
                                         iter = 2000,
                                        zo_inflation = TRUE,
                                         seed = 123,
                                         rhat_cutoff = 1.1,
                                         run_limit = 10) {
  use_dat <- chipseq_dat %>% filter(gene == this_gene,
                                    assay %in% c("writer_loss", "writer_eraser_loss"))
  
  stop_now <- FALSE
  num_runs <- 0
  while(!stop_now) {
    compare_fit <- brms::brm(value ~ 1 + timepoint_cat + assay + timepoint_cat:assay + (1 + timepoint_cat + assay + timepoint_cat:assay | sample_unit), 
                             data = use_dat, 
                             family = ifelse(zo_inflation, "zero_one_inflated_beta", "beta"),
                             iter = iter,
                             control = list(adapt_delta = adapt_delta))
    compare_fixed <- summary(compare_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
    
    num_runs <- num_runs + 1
    if (all(compare_fixed$Rhat < rhat_cutoff) | num_runs == run_limit) {
      stop_now <- TRUE
    }
    else {
      seed <- seed + 1
    }
  }
  results <- list(comparison_wel_to_wl = compare_fixed,
                  seed = seed,
                  run_limit_stop = ifelse(num_runs == run_limit, TRUE, FALSE),
                  zo_inflation = zo_inflation)
  results
}
# Try YAL067W-A
yal067w_a_compare_beta_model <- compare_wl_to_wel_beta_brms(chipseq_dat = gene_chip_dat,
                                                       this_gene = "YAL067W-A",
                                                       iter = 10000, 
                                                       seed = seed)
yal067w_a_compare_model_zoibeta_model <- compare_wl_to_wel_beta_brms(chipseq_dat = gene_chip_dat,
                                                                     this_gene = "YAL067W-A",
                                                                     iter = 10000, 
                                                                     zero_inflation = FALSE,
                                                                     seed = seed)
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
    loss_compare_dat$seed <- this_fit$seed
    loss_compare_dat$run_limit_stop <- this_fit$run_limit_stop
  }
  else {
    holder_dat <- data.frame(gene = rep(this_gene, 8), this_fit)
    holder_dat$seed <- this_fit$seed
    holder_dat$run_limit_stop <- this_fit$run_limit_stop
    loss_compare_dat <- bind_rows(loss_compare_dat,
                                  holder_dat)
  }
}
names(loss_compare_dat) <- gsub(x = names(loss_compare_dat), pattern = "comparison_wel_to_wl.", replacement = "", fixed = TRUE)
loss_compare_dat$attempts <- (loss_compare_dat$seed - 123) + 1
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
  mean # 0.6618665
## timepoint 1
loss_compare_dat %>% 
  filter(parameter == "timepoint_cat1:assaywriter_eraser_loss") %>%
  mutate(check = Estimate > 0) %>%
  pull(check) %>%
  mean # 0.9229474
## timepoint 2
loss_compare_dat %>% 
  filter(parameter == "timepoint_cat2:assaywriter_eraser_loss") %>%
  mutate(check = Estimate > 0) %>%
  pull(check) %>%
  mean # 0.9841032
## timepoint 3
loss_compare_dat %>% 
  filter(parameter == "timepoint_cat3:assaywriter_eraser_loss") %>%
  mutate(check = Estimate > 0) %>%
  pull(check) %>%
  mean # 0.9244436
#### Transitions
## timepoint 0 to timepoint 1
loss_compare_dat %>%
  select(gene, parameter, Estimate) %>%
  filter(parameter %in% c("assaywriter_eraser_loss", "timepoint_cat1:assaywriter_eraser_loss")) %>%
  spread(key = parameter, value = Estimate) %>%
  mutate(check = `timepoint_cat1:assaywriter_eraser_loss` > assaywriter_eraser_loss) %>%
  pull(check) %>%
  mean # 0.9137834
## timepoint 1 to timepoint 2
loss_compare_dat %>%
  select(gene, parameter, Estimate) %>%
  filter(parameter %in% c("timepoint_cat2:assaywriter_eraser_loss", "timepoint_cat1:assaywriter_eraser_loss")) %>%
  spread(key = parameter, value = Estimate) %>%
  mutate(check = `timepoint_cat2:assaywriter_eraser_loss` > `timepoint_cat1:assaywriter_eraser_loss`) %>%
  pull(check) %>%
  mean # 0.7766972
## timepoint 2 to timepoint 3
loss_compare_dat %>%
  select(gene, parameter, Estimate) %>%
  filter(parameter %in% c("timepoint_cat3:assaywriter_eraser_loss", "timepoint_cat2:assaywriter_eraser_loss")) %>%
  spread(key = parameter, value = Estimate) %>%
  mutate(check = `timepoint_cat2:assaywriter_eraser_loss` > `timepoint_cat3:assaywriter_eraser_loss`) %>%
  pull(check) %>%
  mean # 0.9672714
## timepoint 0 to timepoint 3
loss_compare_dat %>%
  select(gene, parameter, Estimate) %>%
  filter(parameter %in% c("timepoint_cat3:assaywriter_eraser_loss", "assaywriter_eraser_loss")) %>%
  spread(key = parameter, value = Estimate) %>%
  mutate(check = assaywriter_eraser_loss < `timepoint_cat3:assaywriter_eraser_loss`) %>%
  pull(check) %>%
  mean # 0.9014401
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
g <- ggplot(data = gene_chip_dat %>% filter(assay != "writer_add", gene %in% lagged_genes), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("black", "red")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
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
## Plot all negative trend genes
g <- ggplot(data = plot_loss_compare_dat %>% filter(gene %in% negative_genes), aes(x = timepoint, y = Estimate)) + geom_line(aes(group = gene, col = gene), alpha = 0.2) + scale_color_grey(guide = FALSE)
g <- g + geom_boxplot(aes(x = timepoint, y = Estimate), outlier.alpha = 0, color = wel_col, size = 0.7, alpha = 0.3) + ylab("WEL - WL timepoint effects") + ylim(-7,10)
g <- g + stat_summary(fun.y = mean, geom="line", aes(group=1), col = wel_col, size = 1)
g <- g + gg_theme# + guides(fill=FALSE)
g <- g + geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)
g

## Plot all high confident genes (confident/non-zero across all assays)
g <- ggplot(data = plot_loss_compare_dat %>% filter(gene %in% negative_genes & gene %in% trend_genes), aes(x = timepoint, y = Estimate)) + geom_line(aes(group = gene, col = gene), alpha = 0.2) + scale_color_grey(guide = FALSE)
g <- g + geom_boxplot(aes(x = timepoint, y = Estimate), outlier.alpha = 0, color = wel_col, size = 0.7, alpha = 0.3) + ylab("WEL - WL timepoint effects") +  ylim(-7,10)
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
## Transcript levels at 0min and gene length
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

## Transcript levels at 60min and gene length
wa_covar_raw_dat <- data.table::fread("data/DtL_setd2_RNA_combined_sf.txt", data.table = FALSE) 
wa_covar_dat <- wa_covar_raw_dat %>%
  rename(gene_length = "Gene Length",
         gene = GENE) %>%
  select(gene, gene_length, DtL_60min_set2d_RNA_Rep1, DtL_60min_set2d_RNA_Rep2) %>%
  rename(rep1 = DtL_60min_set2d_RNA_Rep1,
         rep2 = DtL_60min_set2d_RNA_Rep2) %>%
  gather(key = rep, value = transcript, -c(gene, gene_length)) %>%
  group_by(gene, gene_length) %>%
  summarize(transcript = mean(log(transcript + 1))) %>%
  ungroup
  
## Transcript levels at 0min and gene length
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
     xlab = "Mean H3K36me3 at 0min", ylab = "Histone gain estimate", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.5), las = 1)
plot(x = wl_full_dat$histone, y = abs(wl_full_dat$estimate), 
     xlab = "Mean H3K36me3 at 0min", ylab = "Histone loss estimate", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.5), las = 1)
plot(x = wa_full_dat$histone, y = wa_full_dat$est_error, 
     xlab = "Mean H3K36me3 at 0min", ylab = "Error on histone gain estimate", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.5), las = 1)
plot(x = wl_full_dat$histone, y = abs(wl_full_dat$est_error), 
     xlab = "Mean H3K36me3 at 0min", ylab = "Error on histone loss estimate", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.5), las = 1)
### Trend genes only
par(mfrow = c(2, 2))
plot(x = wa_reduced_dat$histone, y = wa_reduced_dat$estimate, 
     xlab = "Mean H3K36me3 at 0min", ylab = "Histone gain estimate", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.5), las = 1)
plot(x = wl_reduced_dat$histone, y = abs(wl_reduced_dat$estimate), 
     xlab = "Mean H3K36me3 at 0min", ylab = "Histone loss estimate", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.5), las = 1)
plot(x = wa_reduced_dat$histone, y = wa_reduced_dat$est_error, 
     xlab = "Mean H3K36me3 at 0min", ylab = "Error on histone gain estimate", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.5), las = 1)
plot(x = wl_reduced_dat$histone, y = abs(wl_reduced_dat$est_error), 
     xlab = "Mean H3K36me3 at 0min", ylab = "Error on histone loss estimate", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.5), las = 1)

## Histone mark trends by time 0 transcripts
### Full data
par(mfrow = c(2, 2))
plot(x = wa_full_dat$transcript, y = wa_full_dat$estimate, 
     xlab = "Mean 0min transcript", ylab = "Histone gain estimate", 
     main = "Writer add", las = 1, 
     pch = 20, col = alpha(wa_col, 0.3))
abline(lm(estimate ~ transcript, data = wa_full_dat), col = "gray")
text(x = 4, y = 4, paste("r = ", 
                         round(cor(wa_full_dat$transcript, wa_full_dat$estimate), 3)),
     col = "goldenrod", font = 2, cex = 2)
plot(x = wl_full_dat$transcript, y = abs(wl_full_dat$estimate), 
     xlab = "Mean 0min transcript", ylab = "Histone loss estimate", 
     main = "Writer loss", las = 1,
     pch = 20, col = alpha(wl_col, 0.3))
abline(lm(abs(estimate) ~ transcript, data = wl_full_dat), col = "gray")
text(x = 4, y = 3, paste("r = ", 
                         round(cor(wl_full_dat$transcript, abs(wl_full_dat$estimate)), 3)),
     col = "goldenrod", font = 2, cex =2)
plot(x = wa_full_dat$transcript, y = wa_full_dat$est_error, 
     xlab = "Mean 0min transcript", ylab = "Error on histone loss estimate", 
     main = "Writer add", las = 1,
     pch = 20, col = alpha(wa_col, 0.3))
plot(x = wl_full_dat$transcript, y = abs(wl_full_dat$est_error), 
     xlab = "Mean 0min transcript", ylab = "Error on histone loss estimate", 
     main = "Writer loss", las = 1,
     pch = 20, col = alpha(wl_col, 0.3))
cor.test(x = wa_full_dat$transcript, y = wa_full_dat$estimate)
cor.test(x = wl_full_dat$transcript, y = wl_full_dat$estimate)
### Trend genes only
par(mfrow = c(2, 2))
plot(x = wa_reduced_dat$transcript, y = wa_reduced_dat$estimate, 
     xlab = "Mean 0min transcript", ylab = "Histone gain estimate", 
     main = "Writer add", las = 1, 
     pch = 20, col = alpha(wa_col, 0.3))
abline(lm(estimate ~ transcript, data = wa_reduced_dat), col = "gray")
text(x = 4, y = 4, paste("r = ", 
                         round(cor(wa_reduced_dat$transcript, wa_reduced_dat$estimate), 3)),
     col = "goldenrod", font = 2, cex = 2)
plot(x = wl_reduced_dat$transcript, y = abs(wl_reduced_dat$estimate), 
     xlab = "Mean 0min transcript", ylab = "Histone loss estimate", 
     main = "Writer loss", las = 1,
     pch = 20, col = alpha(wl_col, 0.3))
abline(lm(abs(estimate) ~ transcript, data = wl_reduced_dat), col = "gray")
text(x = 4, y = 3, paste("r = ", 
                         round(cor(wl_reduced_dat$transcript, abs(wl_reduced_dat$estimate)), 3)),
     col = "goldenrod", font = 2, cex =2)
plot(x = wa_reduced_dat$transcript, y = wa_reduced_dat$est_error, 
     xlab = "Mean 0min transcript", ylab = "Error on histone loss estimate", 
     main = "Writer add", las = 1,
     pch = 20, col = alpha(wa_col, 0.3))
plot(x = wl_reduced_dat$transcript, y = abs(wl_reduced_dat$est_error), 
     xlab = "Mean 0min transcript", ylab = "Error on histone loss estimate", 
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
plot(wl_full_dat$gene_length, wl_full_dat$estimate*-1, 
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
plot(wl_reduced_dat$gene_length, wl_reduced_dat$estimate*-1, 
     xlab = "Gene length", ylab = "Trend with time", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.3), las = 1)
plot(wa_reduced_dat$gene_length, wa_reduced_dat$est_error,
     xlab = "Gene length", ylab = "Error on trend with time", main = "Writer add",
     pch = 20, col = alpha(wa_col, 0.3), las = 1)
plot(wl_reduced_dat$gene_length, wl_reduced_dat$est_error, 
     xlab = "Gene length", ylab = "Error on trend with time", main = "Writer loss",
     pch = 20, col = alpha(wl_col, 0.3), las = 1)

## Histone mark trend comparisons (WA vs WL)
par(mfrow = c(1,2))
plot(x = gene_timepoint_dat$writer_add, y = gene_timepoint_dat$writer_loss*-1, 
     xlab = "Gain in WA", ylab = "Loss in WL", main = "WL loss vs WA gain",
     pch = 20, col = alpha("gray", 0.5), las = 1)
abline(lm(writer_loss*-1 ~ writer_add, data = gene_timepoint_dat), col = "blue")
text(x = 4, y = 4, paste("r = ", 
                        round(cor(gene_timepoint_dat$writer_add, gene_timepoint_dat$writer_loss*-1), 3)),
    col = "goldenrod", font = 2, cex = 2)
reduced_gene_timepoint_dat <- gene_timepoint_dat %>% filter(gene %in% trend_genes)
plot(x = reduced_gene_timepoint_dat$writer_add, y = reduced_gene_timepoint_dat$writer_loss*-1, 
     xlab = "Gain in WA", ylab = "Loss in WL", main = "WL loss vs WA gain in high confidence genes",
     pch = 20, col = alpha("gray", 0.5), las = 1)
abline(lm(writer_loss*-1 ~ writer_add, data = reduced_gene_timepoint_dat), col = "blue")
text(x = 4, y = 4, paste("r = ", 
                        round(cor(reduced_gene_timepoint_dat$writer_add, reduced_gene_timepoint_dat$writer_loss*-1), 3)),
    col = "goldenrod", font = 2, cex = 2)
cor.test(x = gene_timepoint_dat$writer_add, y = gene_timepoint_dat$writer_loss*-1)
#0.387814
cor.test(x = reduced_gene_timepoint_dat$writer_add, y = reduced_gene_timepoint_dat$writer_loss*-1)
#0.284447

# Error coloring
compare_dat <- gene_timepoint_dat 
compare_dat$combined_error <- as.numeric(apply(compare_dat, 1, function(x) max(x["writer_add_error"], x["writer_loss_error"])))
g1 <- ggplot(data = compare_dat, aes(x = writer_add, y = writer_loss*-1, col = combined_error)) + geom_point() + ggtitle("All Genes") + gg_theme + scale_colour_gradient(low="yellow",high="red") + xlab("Gain in WA") + ylab("Loss in WL") + ylim(-0.5,5) + xlim(-1,10)
g2 <- ggplot(data = compare_dat %>% filter(gene %in% trend_genes), aes(x = writer_add, y = writer_loss*-1, col = combined_error)) + geom_point() + ggtitle("High Confidence Genes") + gg_theme + scale_colour_gradient(low="yellow",high="red") + xlab("Gain in WA") + ylab("Loss in WL") + ylim(-0.5,5) + xlim(-1,10)
grid.arrange(g1, g2)

## Histone mark trend error comparisons (WA vs WL)
par(mfrow = c(1,2))
plot(x = gene_timepoint_dat$writer_add_error, y = gene_timepoint_dat$writer_loss_error, 
     xlab = "Error on gain in WA", ylab = "Error on loss in WL", main = "Errors on WL loss vs WA gain",
     pch = 20, col = alpha("gray", 0.5), las = 1)
plot(x = reduced_gene_timepoint_dat$writer_add_error, y = reduced_gene_timepoint_dat$writer_loss_error, 
     xlab = "Error on gain in WA", ylab = "Error on loss in WL", main = "Errors on WL loss vs WA gain in high confidence genes",
     pch = 20, col = alpha("gray", 0.5), las = 1)

## Gene length by initial H3K36me3 levels
par(mfrow = c(2, 2))
plot(wa_full_dat$histone, wa_full_dat$gene_length, 
     xlab = "Mean H3K36me3 at time 0", ylab = "Gene length", main = "Writer add (All Genes)",
     pch = 20, col = alpha(wa_col, 0.3), las = 1)
plot(wl_full_dat$histone, wl_full_dat$gene_length, 
     xlab = "Mean H3K36me3 at time 0", ylab = "Gene length", main = "Writer loss (All Genes)",
     pch = 20, col = alpha(wl_col, 0.3), las = 1)
plot(wa_reduced_dat$histone, wa_reduced_dat$gene_length, 
     xlab = "Mean H3K36me3 at time 0", ylab = "Gene length", main = "Writer add (High Confidence Genes)",
     pch = 20, col = alpha(wa_col, 0.3), las = 1)
plot(wl_reduced_dat$histone, wl_reduced_dat$gene_length, 
     xlab = "Mean H3K36me3 at time 0", ylab = "Gene length", main = "Writer loss (High Confidence Genes)",
     pch = 20, col = alpha(wl_col, 0.3), las = 1)

## Gene length, histone estimate, and time 0 transcript
g1 <- ggplot(data = wa_full_dat, aes(x = gene_length, y = estimate)) + geom_point(aes(color = transcript)) + ggtitle("Writer add (all genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WA Gain Estimate") + ylim(-1, 10) + xlim(0,15000)
g2 <- ggplot(data = wl_full_dat, aes(x = gene_length, y = estimate*-1)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss (all genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WL Loss Estimate") + ylim(-.5, 8) + xlim(0,15000)
g3 <- ggplot(data = wa_reduced_dat, aes(x = gene_length, y = estimate)) + geom_point(aes(color = transcript)) + ggtitle("Writer add (high confidence genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WA Gain Estimate") + ylim(-1, 10) + xlim(0,15000)
g4 <- ggplot(data = wl_reduced_dat, aes(x = gene_length, y = estimate*-1)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss (high confidence genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WL Loss Estimate") + ylim(-.5, 8) + xlim(0,15000)
grid.arrange(g1, g2, g3, g4)

g1 <- ggplot(data = wa_full_dat, aes(x = gene_length, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer add (all genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WA Gain Error") + ylim(0, 6) + xlim(0,15000)
g2 <- ggplot(data = wl_full_dat, aes(x = gene_length, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss (all genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WL Loss Error") + ylim(0, 5) + xlim(0,15000)
g3 <- ggplot(data = wa_reduced_dat, aes(x = gene_length, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer add (high confidence genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WA Gain Error") + ylim(0, 6) + xlim(0,15000)
g4 <- ggplot(data = wl_reduced_dat, aes(x = gene_length, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss (high confidence genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WL Loss Error") + ylim(0, 5) + xlim(0,15000)
grid.arrange(g1, g2, g3, g4)


## Time 0 histone, histone estimate, and time 0 transcript
g1 <- ggplot(data = wa_full_dat, aes(x = histone, y = estimate)) + geom_point(aes(color = transcript)) + ggtitle("Writer add (all genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WA Gain Estimate") + ylim(-1, 10) + xlim(0,0.3)
g2 <- ggplot(data = wl_full_dat, aes(x = histone, y = estimate*-1)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss (all genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WL Loss Estimate") + ylim(-0.5, 5) + xlim(0,1)
g3 <- ggplot(data = wa_full_dat, aes(x = histone, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer add (all genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WA Gain Error") + ylim(0, 6) + xlim(0,0.3)
g4 <- ggplot(data = wl_full_dat, aes(x = histone, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss (all genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WL Loss Error") + ylim(-0.5, 5.5) + xlim(0,1)
grid.arrange(g1, g2, g3, g4)

g1 <- ggplot(data = wa_reduced_dat, aes(x = histone, y = estimate)) + geom_point(aes(color = transcript)) + ggtitle("Writer add (high confidence genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WA Gain Estimate") + ylim(-1, 10) + xlim(0,0.3)
g2 <- ggplot(data = wl_reduced_dat, aes(x = histone, y = estimate*-1)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss (high confidence genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WL Loss Estimate") + ylim(-0.5, 5) + xlim(0,1)
g3 <- ggplot(data = wa_reduced_dat, aes(x = histone, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer add (high confgenes genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WA Gain Error") + ylim(0, 6) + xlim(0,0.3)
g4 <- ggplot(data = wl_reduced_dat, aes(x = histone, y = est_error)) + geom_point(aes(color = transcript)) + ggtitle("Writer loss (high confgenes genes)") + gg_theme + scale_colour_gradient(low="yellow",high="red") + ylab("WL Loss Error") + ylim(-0.5, 5.5) + xlim(0,1)
grid.arrange(g1, g2, g3, g4)


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
##  Raw histones data for time 0 covariate
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
transcription_stan_fit_paths <- list.files("transcription_Stan_models", full.names = TRUE)
file_path <- "transcription_Stan_models/sacCer3_nonOverlappingGenes_noChrMGenes_"
for (i in 1:length(transcription_stan_fit_paths)) {
  this_fit <- readRDS(transcription_stan_fit_paths[i])
  
  this_gene <- gsub(x = transcription_stan_fit_paths[i], pattern = file_path, replacement = "", fixed = TRUE) %>%
    gsub(x = ., replacement = "", pattern = "_transcription_Stan.rds", fixed = TRUE)
  
  
  if (i == 1) {
    transcript_timepoint_dat <- bind_rows(data.frame(gene = this_gene, assay = "writer_loss", this_fit$writer_loss %>% filter(parameter == "timepoint")),
                                          data.frame(gene = this_gene, assay = "writer_add", this_fit$writer_add %>% filter(parameter == "timepoint")))
    histone_transcript_dat <- bind_rows(data.frame(gene = this_gene, assay = "writer_loss", this_fit$writer_loss %>% filter(parameter == "histone")),
                                        data.frame(gene = this_gene, assay = "writer_add", this_fit$writer_add %>% filter(parameter == "histone")))
  }
  else {
    transcript_timepoint_dat <- bind_rows(transcript_timepoint_dat,
                                          data.frame(gene = this_gene, assay = "writer_loss", this_fit$writer_loss %>% filter(parameter == "timepoint")),
                                          data.frame(gene = this_gene, assay = "writer_add", this_fit$writer_add %>% filter(parameter == "timepoint")))
    histone_transcript_dat <- bind_rows(histone_transcript_dat,
                                        data.frame(gene = this_gene, assay = "writer_loss", this_fit$writer_loss %>% filter(parameter == "histone")),
                                        data.frame(gene = this_gene, assay = "writer_add", this_fit$writer_add %>% filter(parameter == "histone")))
  }
}
transcript_timepoint_dat <- transcript_timepoint_dat %>% 
  select(-parameter) %>%
  mutate(is_pos = ifelse(lower_95 > 0, TRUE, FALSE),
         is_neg = ifelse(upper_95 < 0, TRUE, FALSE))
transcript_timepoint_dat$category[transcript_timepoint_dat$lower_95 > 0] <- "positive"
transcript_timepoint_dat$category[transcript_timepoint_dat$upper_95 < 0] <- "negative"
transcript_timepoint_dat$category[transcript_timepoint_dat$lower_95 < 0 & transcript_timepoint_dat$upper_95 > 0] <- "zero"

histone_transcript_dat <- histone_transcript_dat %>% 
  select(-parameter) %>%
  mutate(is_pos = ifelse(lower_95 > 0, TRUE, FALSE),
         is_neg = ifelse(upper_95 < 0, TRUE, FALSE))
histone_transcript_dat$category[histone_transcript_dat$lower_95 > 0] <- "positive"
histone_transcript_dat$category[histone_transcript_dat$upper_95 < 0] <- "negative"
histone_transcript_dat$category[histone_transcript_dat$lower_95 < 0 & histone_transcript_dat$upper_95 > 0] <- "zero"

## Grabbing genes with positive transcription
positive_transcript_genes <- transcript_timepoint_dat %>% 
  filter(category == "positive") %>%
  pull(gene) %>%
  unique
## Grabbing genes with positive histone
positive_histone_genes <- histone_transcript_dat %>% 
  filter(category == "positive") %>%
  pull(gene) %>%
  unique
## Grabbing genes with negative transcription
negative_transcript_genes <- transcript_timepoint_dat %>% 
  filter(category == "negative") %>%
  pull(gene) %>%
  unique
## Grabbing genes with negative histone
negative_histone_genes <- histone_transcript_dat %>% 
  filter(category == "negative") %>%
  pull(gene) %>%
  unique

##############################################
##
##  Transcription regression: Trend estimate by error
##
##############################################
plot(transcript_timepoint_dat$est_error, transcript_timepoint_dat$estimate, pch = ifelse(transcript_timepoint_dat$category == "zero", 1, 19), col = c(wa_col, wl_col)[as.factor(timepoint_dat$assay)], 
     las = 1, ylab = "Transcription estimate with time", xlab = "Error on estimate", frame.plot = FALSE, xlim = c(0, 10), ylim = c(-7, 7))
abline(h = 0, lty = 2)
legend("bottomright", 
       legend = c("writer add", "writer loss"),
       col = c(wa_col, wl_col),
       fill = c(wa_col, wl_col),
       bty = "n")


#YFR044C
g <- ggplot(data = tx_raw_dat %>% filter(gene %in% positive_transcript_genes[8]), aes(x = timepoint, y = log(transcript + 1), col = assay)) + scale_color_manual(values = c(wl_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = log(transcript + 1), x = timepoint), method = "lm", size = 2)
g <- g + gg_theme
g

#YCL002C
g <- ggplot(data = tx_raw_dat %>% filter(gene %in% negative_transcript_genes[6]), aes(x = timepoint, y = log(transcript + 1), col = assay)) + scale_color_manual(values = c(wl_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = log(transcript + 1), x = timepoint), method = "lm", size = 2)
g <- g + gg_theme
g

table(transcript_timepoint_dat$category, transcript_timepoint_dat$assay)
table(histone_transcript_dat$category, histone_transcript_dat$assay)

combined_dat <- inner_join(transcript_timepoint_dat %>%
                             select(gene, assay, category) %>%
                             rename(transcript_category = category),
                           histone_transcript_dat %>%
                             select(gene, assay, category) %>%
                             rename(histone_category = category))

table(combined_dat$transcript_category, combined_dat$histone_category, dnn = c("transcript", "histone"))
chisq.test(combined_dat$transcript_category, combined_dat$histone_category)

shared_positive_genes <- combined_dat %>%
  filter(transcript_category == "positive" & histone_category == "positive") %>%
  pull(gene)
combined_dat %>% filter(gene == shared_positive_genes[1])

#YDL216C
g_tx <- ggplot(data = tx_raw_dat %>% filter(gene %in% shared_positive_genes[1],
                                         assay == "writer_loss"), 
            aes(x = timepoint, y = log(transcript + 1), col = assay)) + scale_color_manual(values = c(wl_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g_tx <- g_tx + geom_smooth(aes(y = log(transcript + 1), x = timepoint), method = "lm", size = 2)
g_tx <- g_tx + gg_theme
#g_tx
g_h <- ggplot(data = histone_raw_dat %>% filter(gene %in% shared_positive_genes[1],
                                         assay == "writer_loss"), 
            aes(x = timepoint, y = histone, col = assay)) + scale_color_manual(values = c(wl_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g_h <- g_h + geom_smooth(aes(y = histone, x = timepoint), method = "lm", size = 2)
g_h <- g_h + gg_theme
#g_h
grid.arrange(g_tx, g_h)



## Wide transcription data
wide_transcript_timepoint_dat <- transcript_timepoint_dat %>%
  select(gene, assay, estimate) %>%
  spread(key = assay, value = estimate)
par(mfrow=c(1,2))
plot(wide_transcript_timepoint_dat$writer_add, wide_transcript_timepoint_dat$writer_loss, xlim = c(-10, 10), ylim = c(-10, 10), main="WA vs WL Transcription GLMM (All genes)", xlab="WA Estimate", ylab="WL Estimate")
wide_reduced_transcript_timepoint_dat  <- wide_transcript_timepoint_dat %>%
  filter(gene %in% unique(c(positive_transcript_genes, negative_transcript_genes)))
plot(wide_reduced_transcript_timepoint_dat$writer_add, wide_reduced_transcript_timepoint_dat$writer_loss, xlim = c(-10, 10), ylim = c(-10, 10), main="WA vs WL Transcription GLMM (High confidence transcription genes)", xlab="WA Estimate", ylab="WL Estimate")


## Modeling unscaled histone marks
raw_unscaled_gene_chip_dat <- read.table("data/absolute_matrix_for_LMM_model_nonOverlappingGenes_noChrMGenes.txt") # path to data

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

#YAL001C
g1 <- ggplot(data = unscaled_gene_chip_dat %>% filter(gene == unscaled_gene_chip_dat$gene[1]), aes(x = timepoint, y = value * 1000, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g1 <- g + geom_smooth(aes(y = value*1000, x = timepoint), method = "lm", size = 2)
g1 <- g + gg_theme
g1

#YAL001C
g2 <- ggplot(data = gene_chip_dat %>% filter(gene == unscaled_gene_chip_dat$gene[1]), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g2 <- g + geom_smooth(aes(y = value*1000, x = timepoint), method = "lm", size = 2)
g2 <- g + gg_theme
g2
grid.arrange(g1,g2)

hist(unscaled_gene_chip_dat$value * 1000)
hist(unscaled_gene_chip_dat$transformed_value)

## Arbitrary, but use transformed_value = round(value * 1000)
run_zinegbin_brms_on_chipseq <- function(chipseq_dat, 
                                         this_gene, 
                                         adapt_delta = 0.8,
                                         iter = 2000,
                                         seed = 123,
                                         rhat_cutoff = 1.1,
                                         run_limit = 10) {
  use_dat <- chipseq_dat %>% filter(gene == this_gene)#,

  num_runs <- 0
  stop_now <- FALSE
  while(!stop_now) {
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
    
    num_runs <- num_runs + 1
    if ((all(wl_fixed$Rhat < rhat_cutoff) & all(wel_fixed$Rhat < rhat_cutoff) & all(wa_fixed$Rhat < rhat_cutoff)) | num_runs == run_limit) {
      stop_now <- TRUE
    }
    else {
      seed <- seed + 1
    }
  }
  names(wl_fixed) <- names(wl_random) <- names(wel_fixed) <- names(wel_random) <- names(wa_fixed) <- names(wa_random) <- c("parameter", "estimate", "est_error", "lower_95", "upper_95", "eff_sample", "rhat")
  results <- list(writer_loss = bind_rows(wl_fixed, wl_random),
                  writer_eraser_loss = bind_rows(wel_fixed, wel_random),
                  writer_add = bind_rows(wa_fixed, wa_random),
                  seed = seed,
                  run_limit_stop = ifelse(num_runs == run_limit, TRUE, FALSE))
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

raw_yal012w_models <- run_zinegbin_brms_on_chipseq(chipseq_dat = unscaled_gene_chip_dat, this_gene = "YAL012W", iter = 10000)


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




########################################################
##
##  ZI Negative Binomial regression: aggregate results
##
########################################################
## Parsing full set of Stan results
stan_negbin_paths <- list.files("individual_gene_Stan_models_absolute_signal_alpha_0.8_iter_10000/", full.names = TRUE)
file_path <- "individual_gene_Stan_models_absolute_signal_alpha_0.8_iter_10000//sacCer3_nonOverlappingGenes_noChrMGenes_"
for (i in 1:length(stan_negbin_paths)) {
  this_fit <- readRDS(stan_negbin_paths[i])
  
  this_gene <- gsub(x = stan_negbin_paths[i], pattern = file_path, replacement = "", fixed = TRUE) %>%
    gsub(x = ., replacement = "", pattern = "_absolute_signal_Stan.rds", fixed = TRUE)
  
  
  if (i == 1) {
    absolute_dat <- bind_rows(data.frame(gene = this_gene, assay = "writer_loss", this_fit$writer_loss %>% filter(parameter == "timepoint")),
                               data.frame(gene = this_gene, assay = "writer_eraser_loss", this_fit$writer_eraser_loss %>% filter(parameter == "timepoint")),
                               data.frame(gene = this_gene, assay = "writer_add", this_fit$writer_add %>% filter(parameter == "timepoint")))
    absolute_dat$seed <- this_fit$seed
    absolute_dat$run_limit_stop <- this_fit$run_limit_stop
  }
  else {
    holder_dat <-bind_rows(data.frame(gene = this_gene, assay = "writer_loss", this_fit$writer_loss %>% filter(parameter == "timepoint")),
                           data.frame(gene = this_gene, assay = "writer_eraser_loss", this_fit$writer_eraser_loss %>% filter(parameter == "timepoint")),
                           data.frame(gene = this_gene, assay = "writer_add", this_fit$writer_add %>% filter(parameter == "timepoint")))
    holder_dat$seed <- this_fit$seed
    holder_dat$run_limit_stop <- this_fit$run_limit_stop
    absolute_dat <- bind_rows(absolute_dat,
                              holder_dat)
  }
}
absolute_dat <- absolute_dat %>% 
  select(-parameter) %>%
  mutate(is_pos = ifelse(lower_95 > 0, TRUE, FALSE),
         is_neg = ifelse(upper_95 < 0, TRUE, FALSE))
absolute_dat$category[absolute_dat$lower_95 > 0] <- "positive"
absolute_dat$category[absolute_dat$upper_95 < 0] <- "negative"
absolute_dat$category[absolute_dat$lower_95 < 0 & absolute_dat$upper_95 > 0] <- "zero"
absolute_dat$attempts <- (absolute_dat$seed - 123) + 1


## Grabbing high genes
absolute_positive_genes <- absolute_dat %>% 
  filter(category == "positive") %>%
  pull(gene) %>%
  unique
## Grabbing low genes
absolute_negative_genes <- absolute_dat %>% 
  filter(category == "negative") %>%
  pull(gene) %>%
  unique
## Grabbing more categories of genes
absolute_genes_neg_wl <- absolute_dat %>% 
  filter(assay == "writer_loss",
         category == "negative") %>%
  pull(gene)
absolute_genes_neg_wel <- absolute_dat %>% 
  filter(assay == "writer_eraser_loss",
         category == "negative") %>%
  pull(gene)
absolute_genes_pos_wa <- absolute_dat %>% 
  filter(assay == "writer_add",
         category == "positive") %>%
  pull(gene)
absolute_trend_genes <- intersect(absolute_positive_genes, absolute_genes_neg_wl) %>% intersect(absolute_genes_neg_wel)
absolute_genes_neg_loss_only <- absolute_negative_genes[!(absolute_negative_genes %in% absolute_genes_pos_wa)]
## Grabbing genes with no trends
absolute_genes_all_zero <- absolute_dat %>% 
  select(gene, assay, category) %>%
  spread(key = assay, value = category) %>%
  filter(writer_add == "zero", 
         writer_eraser_loss == "zero",
         writer_loss == "zero") %>%
  pull(gene)


absolute_upset_dat <- data.frame(genes = all_genes,
                        WA_positive = sapply(1:length(all_genes), function(i) as.numeric(all_genes[i] %in% absolute_genes_pos_wa)),
                        WL_negative = sapply(1:length(all_genes), function(i) as.numeric(all_genes[i] %in% absolute_genes_neg_wl)),
                        WEL_negative = sapply(1:length(all_genes), function(i) as.numeric(all_genes[i] %in% absolute_genes_neg_wel)))
upset(absolute_upset_dat, 
      order.by = "freq", 
      text.scale = 1.2,
      sets.bar.color = c(wa_col, wl_col, wel_col))

plot(absolute_dat$est_error, absolute_dat$estimate, pch = ifelse(absolute_dat$category == "zero", 1, 19), col = c(wa_col, wl_col, wel_col)[as.factor(absolute_dat$assay)], 
     las = 1, ylab = "Histone trend with time", xlab = "Error on trend", frame.plot = FALSE)
abline(h = 0, lty = 2)
legend("bottomright", 
       legend = c("writer add", "writer loss", "writer eraser loss", "confident within assay", "not confident within assay"),
       col = c(wa_col, wl_col, wel_col, "gray", "gray"),
       pch=c(15, 15, 15, 19, 1),
       bty = "n")

confident_absolute_dat <- absolute_dat %>% filter(rhat < 1.1)
plot(confident_absolute_dat$est_error, confident_absolute_dat$estimate, pch = ifelse(confident_absolute_dat$category == "zero", 1, 19), col = c(wa_col, wl_col, wel_col)[as.factor(confident_absolute_dat$assay)], 
     las = 1, ylab = "Histone trend with time", xlab = "Error on trend", frame.plot = FALSE)
abline(h = 0, lty = 2)
legend("bottomright", 
       legend = c("writer add", "writer loss", "writer eraser loss", "confident within assay", "not confident within assay"),
       col = c(wa_col, wl_col, wel_col, "gray", "gray"),
       pch=c(15, 15, 15, 19, 1),
       bty = "n")
sum(absolute_dat$rhat > 1.1)

plot(absolute_dat$est_error, absolute_dat$eff_sample, pch = ifelse(absolute_dat$category == "zero", 1, 19), col = c(wa_col, wl_col, wel_col)[as.factor(absolute_dat$assay)], 
     las = 1, ylab = "Histone trend with time", xlab = "Effective sample size", frame.plot = FALSE)
abline(h = 0, lty = 2)
legend("bottomright", 
       legend = c("writer add", "writer loss", "writer eraser loss", "confident within assay", "not confident within assay"),
       col = c(wa_col, wl_col, wel_col, "gray", "gray"),
       pch=c(15, 15, 15, 19, 1),
       bty = "n")

plot(absolute_dat$est_error, absolute_dat$estimate, pch = ifelse(absolute_dat$category == "zero", 1, 19), col = c(wa_col, wl_col, wel_col)[as.factor(absolute_dat$assay)], 
     las = 1, ylab = "Histone trend with time", xlab = "Error on trend", frame.plot = FALSE, ylim = c(-1.5, 1.5), xlim = c(0, 10))
abline(h = 0, lty = 2)
legend("bottomright", 
       legend = c("writer add", "writer loss", "writer eraser loss", "confident within assay", "not confident within assay"),
       col = c(wa_col, wl_col, wel_col, "gray", "gray"),
       pch=c(15, 15, 15, 19, 1),
       bty = "n")


## Compare negbinom to beta model
compare_rate_dat <- inner_join(timepoint_dat %>%
                                 select(gene, assay, estimate, category) %>%
                                 dplyr::rename(beta_estimate = estimate,
                                               beta_category = category),
                               absolute_dat %>%
                                 select(gene, assay, estimate, category) %>%
                                 dplyr::rename(negbin_estimate = estimate,
                                               negbin_category = category))
par(mfrow = c(1, 2))
plot(compare_rate_dat$negbin_estimate, compare_rate_dat$beta_estimate, pch = ifelse(absolute_dat$category == "zero", 1, 19), col = c(wa_col, wl_col, wel_col)[as.factor(absolute_dat$assay)], 
     las = 1, ylab = "Histone trend with time (proportions)", xlab = "Histone trend with time (absolute)", frame.plot = FALSE)
plot(compare_rate_dat$negbin_estimate, compare_rate_dat$beta_estimate, pch = ifelse(absolute_dat$category == "zero", 1, 19), col = c(wa_col, wl_col, wel_col)[as.factor(absolute_dat$assay)], 
     las = 1, ylab = "Histone trend with time (proportions)", xlab = "Histone trend with time (absolute)", frame.plot = FALSE, xlim = c(-2, 2))

high_wa <- absolute_dat %>%
  filter(assay == "writer_add",
         estimate > 100) %>%
  pull(gene)

g1 <- ggplot(data = unscaled_gene_chip_dat %>% filter(gene == high_wa[1]), aes(x = timepoint, y = transformed_value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g1 <- g1 + geom_smooth(aes(y = transformed_value, x = timepoint), method = "lm", size = 2)
g1 <- g1 + gg_theme
g1

YAL001C_negbin_practice <- run_zinegbin_brms_on_chipseq(chipseq_dat = unscaled_gene_chip_dat,
                                               this_gene = "YAL001C",
                                               iter = 10000)
YER175C_negbin_practice <- run_zinegbin_brms_on_chipseq(chipseq_dat = unscaled_gene_chip_dat,
                                                        this_gene = "YER175C",
                                                        iter = 10000)
g2 <- ggplot(data = unscaled_gene_chip_dat %>% filter(gene == "YAL001C"), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g2 <- g2 + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g2 <- g2 + gg_theme
g2
g2 <- ggplot(data = unscaled_gene_chip_dat %>% filter(gene == "YER175C"), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g2 <- g2 + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g2 <- g2 + gg_theme
g2

high_wel <- absolute_dat %>%
  filter(assay == "writer_eraser_loss") %>%
  filter(estimate == max(estimate)) %>%
  pull(gene)

absolute_dat %>%
  filter(assay == "writer_add",
         estimate > 5) %>%
  head


g2 <- ggplot(data = unscaled_gene_chip_dat %>% filter(gene == "YER175C"), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g2 <- g2 + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g2 <- g2 + gg_theme
g2


####################################
## Genes Austin requested
## 
####################################
YCR008W_prop_g <- ggplot(data = gene_chip_dat %>% filter(gene %in% "YCR008W"), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point(size = 2) + geom_line(aes(group = sample_unit), linetype = "longdash", size = 1.1)
YCR008W_prop_g <- YCR008W_prop_g + gg_theme + guides(color = FALSE) + ggtitle("YCR008W proportional ChIP-seq data")
YCR008W_prop_g <- YCR008W_prop_g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
YCR008W_prop_g

YCR005C_prop_g <- ggplot(data = gene_chip_dat %>% filter(gene %in% "YCR005C"), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point(size = 2) + geom_line(aes(group = sample_unit), linetype = "longdash", size = 1.1)
YCR005C_prop_g <- YCR005C_prop_g + gg_theme + guides(color = FALSE) + ggtitle("YCR005C proportional ChIP-seq data")
YCR005C_prop_g <- YCR005C_prop_g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
YCR005C_prop_g

YCR008W_nonprop_g <- ggplot(data = unscaled_gene_chip_dat %>% filter(gene %in% "YCR008W"), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point(size = 2) + geom_line(aes(group = sample_unit), linetype = "longdash", size = 1.1)
YCR008W_nonprop_g <- YCR008W_nonprop_g + gg_theme + guides(color = FALSE) + ggtitle("YCR008W normalized ChIP-seq data")
YCR008W_nonprop_g <- YCR008W_nonprop_g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
YCR008W_nonprop_g

YCR005C_nonprop_g <- ggplot(data = unscaled_gene_chip_dat %>% filter(gene %in% "YCR005C"), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c(wl_col, wel_col, wa_col)) + geom_point(size = 2) + geom_line(aes(group = sample_unit), linetype = "longdash", size = 1.1)
YCR005C_nonprop_g <- YCR005C_nonprop_g + gg_theme + guides(color = FALSE) + ggtitle("YCR005C normalized ChIP-seq data")
YCR005C_nonprop_g <- YCR005C_nonprop_g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
YCR005C_nonprop_g

pdf("~/Documents/git_repositories/chipseq-gene-dynamics/modeling_figures_gkeele/YCR005C_data.pdf", height = 4, width = 8)
grid.arrange(YCR005C_nonprop_g, YCR005C_prop_g, nrow = 1)
dev.off()

pdf("~/Documents/git_repositories/chipseq-gene-dynamics/modeling_figures_gkeele/YCR008W_data.pdf", height = 4, width = 8)
grid.arrange(YCR008W_nonprop_g, YCR008W_prop_g, nrow = 1)
dev.off()

##################################
##
## Table for paper
##
##################################
table_dat <- inner_join(timepoint_dat %>%
                          mutate(estimate = ifelse(rhat > 1.1, NA, estimate),
                          est_error = ifelse(rhat > 1.1, NA, est_error)) %>%
                          select(gene, assay, estimate, est_error, lower_95, upper_95, category) %>%
                          dplyr::rename(Gene = gene,
                                        ZOIBetaEstimate = estimate,
                                        ZOIBetaError = est_error,
                                        ZOIBetaLower95 = lower_95,
                                        ZOIBetaUpper95 = upper_95,
                                        ZOIBetaTrend = category),
                         absolute_dat %>%
                           mutate(estimate = ifelse(rhat > 1.1, NA, estimate),
                                  est_error = ifelse(rhat > 1.1, NA, est_error)) %>%
                           select(gene, assay, estimate, est_error, lower_95, upper_95, category) %>%
                           dplyr::rename(Gene = gene,
                                         ZINegEstimate = estimate,
                                         ZINegError = est_error,
                                         ZINegLower95 = lower_95,
                                         ZINegUpper95 = upper_95,
                                         ZINegTrend = category))

table_dat %>% 
  filter(assay == "writer_loss") %>%
  select(-assay) %>%
  write.table(., 
              file = "tables/writer_loss_gene_table.txt",
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
table_dat %>% 
  filter(assay == "writer_eraser_loss") %>%
  select(-assay) %>%
  write.table(., 
              file = "tables/writer_eraser_loss_gene_table.txt",
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
table_dat %>% 
  filter(assay == "writer_add") %>%
  select(-assay) %>%
  write.table(., 
              file = "tables/writer_add_gene_table.txt",
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)



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


