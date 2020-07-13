##############################################################################
##############################################################################
###
###       Title: Supplemental_Code.R
###
###       Description: R code to generate findings on the dynamics 
###                    of histone trimethylation in yeast genes
###           
###       Publication: Lerner et al.
###                    An optogenetic switch for the Set2 
###                    methyltransferase provides evidence for 
###                    transcription-dependent and independent 
###                    dynamics of H3K36 methylation
###             
###       Authors: Greg Keele & Austin Hepperla
###
##############################################################################
##############################################################################

### Required R package dependencies
library(tidyverse)
library(brms)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(Vennerable)
library(data.table)

### Plot settings
## Colors
wa_col <- "plum"
wl_col <- "darkseagreen2"
wel_col <- "tomato"

## ggplot theme
plot_theme <- theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), 
                    axis.line = element_line(colour = "black"),
                    plot.title = element_text(hjust = 0.5), 
                    axis.text = element_text(size = 12),
                    axis.title = element_text(size = 12),
                    axis.text.x = element_text(size = 12))

# Note: Stan fits do not perfectly replicate across software versions and computing environments, even when seed is set
#       Exact numbers and plots may vary across different computing environments, but results should be consistent

###########################################
##        
##          Read in and process data
##
###########################################
### Set working directory to one containging
# setwd("~/Documents/git_repositories/chipseq-gene-dynamics/supplemental_data/")
# setwd("~/projects/chipseq-gene-dynamics/supplemental_data/")

###  Relative h3k36 methylation data (quantiles of maximum value)
raw_rel_h3k36_dat <- read.table("Supplemental_Data1.txt")
names(raw_rel_h3k36_dat) <- c("gene", "assay", "replicate", "timepoint", "value")

## Process for further analysis
rel_h3k36_dat <- raw_rel_h3k36_dat %>%
  mutate(assay = as.factor(assay),
         assay = recode(assay, "1" = "writer_loss", "2" = "writer_eraser_loss", "3" = "writer_add"),
         sample_unit = paste(gene, assay, replicate, sep = "_"),
         timepoint_cat = factor(timepoint)) %>%
  mutate(timepoint = case_when(assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 0 ~ 0,
                               assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 1 ~ 30/60,
                               assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 2 ~ 60/60,
                               assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 3 ~ 90/60,
                               assay == "writer_add" & timepoint_cat == 0 ~ 0,
                               assay == "writer_add" & timepoint_cat == 1 ~ 20/60,
                               assay == "writer_add" & timepoint_cat == 2 ~ 40/60,
                               assay == "writer_add" & timepoint_cat == 3 ~ 60/60),
         timepoint_min = case_when(assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 0 ~ 0,
                                   assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 1 ~ 30,
                                   assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 2 ~ 60,
                                   assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 3 ~ 90,
                                   assay == "writer_add" & timepoint_cat == 0 ~ 0,
                                   assay == "writer_add" & timepoint_cat == 1 ~ 20,
                                   assay == "writer_add" & timepoint_cat == 2 ~ 40,
                                   assay == "writer_add" & timepoint_cat == 3 ~ 60))

##  Mean summaries of the relative h3k36 data
mean_rel_h3k36_dat <- rel_h3k36_dat %>%
  group_by(gene, assay, timepoint) %>%
  summarise(value = mean(value)) %>%
  ungroup

###  Absolute h3k36 data
raw_abs_h3k36_dat <- read.table("Supplemental_Data2.txt")
names(raw_abs_h3k36_dat) <- c("gene", "assay", "replicate", "timepoint", "value")

## Process for further analysis
abs_h3k36_dat <- raw_abs_h3k36_dat %>%
  mutate(transformed_value = round(value * 1000),
         assay = as.factor(assay),
         assay = recode(assay, "1" = "writer_loss", "2" = "writer_eraser_loss", "3" = "writer_add"),
         sample_unit = paste(gene, assay, replicate, sep = "_"),
         timepoint_cat = factor(timepoint)) %>%
  mutate(timepoint = case_when(assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 0 ~ 0,
                               assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 1 ~ 30/60,
                               assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 2 ~ 60/60,
                               assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 3 ~ 90/60,
                               assay == "writer_add" & timepoint_cat == 0 ~ 0,
                               assay == "writer_add" & timepoint_cat == 1 ~ 20/60,
                               assay == "writer_add" & timepoint_cat == 2 ~ 40/60,
                               assay == "writer_add" & timepoint_cat == 3 ~ 60/60),
         timepoint_min = case_when(assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 0 ~ 0,
                                   assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 1 ~ 30,
                                   assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 2 ~ 60,
                                   assay %in% c("writer_loss", "writer_eraser_loss") & timepoint_cat == 3 ~ 90,
                                   assay == "writer_add" & timepoint_cat == 0 ~ 0,
                                   assay == "writer_add" & timepoint_cat == 1 ~ 20,
                                   assay == "writer_add" & timepoint_cat == 2 ~ 40,
                                   assay == "writer_add" & timepoint_cat == 3 ~ 60))

##  Mean summaries of the relative h3k36 data
mean_abs_h3k36_dat <- abs_h3k36_dat %>%
  group_by(gene, assay, timepoint) %>%
  summarise(value = mean(value)) %>%
  ungroup

###########################################
##        
##          Functions to fit Bayesian
##          GLMMs with brms 
##          -- R package for running Stan
##
###########################################
### Beta GLMM
## Used with the relative h3k36 data
## Option to use zero-one inflated model
## Fits three models for a given gene
## -- writer add, writer loss, and writer eraser loss
run_beta_brms <- function(rel_dat, 
                          this_gene, 
                          adapt_delta = 0.8,
                          iter = 2000,
                          zo_inflation = TRUE,
                          seed = 123,
                          rhat_cutoff = 1.1,
                          run_limit = 10) {
  
  use_dat <- rel_dat %>% 
    filter(gene == this_gene)
  
  stop_now <- FALSE
  num_runs <- 0
  while(!stop_now) {
    # writer loss
    wl_fit <- brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                  data = use_dat %>% filter(assay == "writer_loss"), 
                  family = ifelse(zo_inflation, "zero_one_inflated_beta", "beta"),
                  iter = iter,
                  control = list(adapt_delta = adapt_delta),
                  seed = seed)
    wl_fixed <- summary(wl_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
    wl_random <- summary(wl_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
    # writer-eraser loss
    wel_fit <- brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                   data = use_dat %>% filter(assay == "writer_eraser_loss"), 
                   family = ifelse(zo_inflation, "zero_one_inflated_beta", "beta"),
                   iter = iter,
                   control = list(adapt_delta = adapt_delta),
                   seed = seed)
    wel_fixed <- summary(wel_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
    wel_random <- summary(wel_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
    # writer add
    wa_fit <- brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
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

# Fit a single beta model and return brms object
run_single_beta_brms <- function(rel_dat, 
                                 this_gene, 
                                 this_assay = c("writer_loss", "writer_eraser_loss", "writer_add"),
                                 adapt_delta = 0.8,
                                 iter = 2000,
                                 zo_inflation = TRUE,
                                 seed = 123,
                                 rhat_cutoff = 1.1,
                                 run_limit = 10) {
  this_assay <- this_assay[1]
  
  use_dat <- rel_dat %>% 
    filter(gene == this_gene,
           assay == this_assay)
  
  stop_now <- FALSE
  num_runs <- 0
  while(!stop_now) {
    brm_fit <- brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                   data = use_dat, 
                   family = ifelse(zo_inflation, "zero_one_inflated_beta", "beta"),
                   iter = iter,
                   control = list(adapt_delta = adapt_delta),
                   seed = seed)
    brm_fixed <- summary(brm_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
    brm_random <- summary(brm_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
    
    num_runs <- num_runs + 1
    if ((all(brm_fixed$Rhat < rhat_cutoff) & all(brm_fixed$Rhat < rhat_cutoff) & all(brm_fixed$Rhat < rhat_cutoff)) | num_runs == run_limit) {
      stop_now <- TRUE
    }
    else {
      seed <- seed + 1
    }
  }
  
  brm_fit
}

## Function to run brms for ZOI beta regression of loss assays
## Allows for comparison of writer loss and writer eraser loss
## Used with the relative h3k36 data
## Option to use zero-one inflated model
compare_wl_to_wel_beta_brms <- function(rel_dat, 
                                        this_gene, 
                                        adapt_delta = 0.8,
                                        iter = 2000,
                                        zo_inflation = TRUE,
                                        seed = 123,
                                        rhat_cutoff = 1.1,
                                        run_limit = 10) {

  use_dat <- rel_dat %>% filter(gene == this_gene,
                                    assay %in% c("writer_loss", "writer_eraser_loss"))
  
  stop_now <- FALSE
  num_runs <- 0
  while(!stop_now) {
    compare_fit <- brm(value ~ 1 + timepoint_cat + assay + timepoint_cat:assay + (1 + timepoint_cat + assay + timepoint_cat:assay | sample_unit), 
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

### Zero-inflated negative binomial GLMM
## Used with the absolute h3k36 data 
## Pseudo-counts: arbitrary, but uses transformed_value = round(value * 1000)
## Fits three models for a given gene
## -- writer add, writer loss, and writer eraser loss
run_zinegbin_brms <- function(abs_dat, 
                              this_gene, 
                              adapt_delta = 0.8,
                              iter = 2000,
                              seed = 123,
                              rhat_cutoff = 1.1,
                              run_limit = 10) {
  
  use_dat <- abs_dat %>% 
    filter(gene == this_gene)
  
  num_runs <- 0
  stop_now <- FALSE
  while(!stop_now) {
    # writer loss
    wl_fit <- brm(transformed_value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                  data = use_dat %>% filter(assay == "writer_loss"), 
                  family = "zero_inflated_negbinomial",
                  iter = iter,
                  control = list(adapt_delta = adapt_delta))
    wl_fixed <- summary(wl_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
    wl_random <- summary(wl_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
    # writer-eraser loss
    wel_fit <- brm(transformed_value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                   data = use_dat %>% filter(assay == "writer_eraser_loss"), 
                   family = "zero_inflated_negbinomial",
                   iter = iter,
                   control = list(adapt_delta = adapt_delta))
    wel_fixed <- summary(wel_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
    wel_random <- summary(wel_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
    # writer add
    wa_fit <- brm(transformed_value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
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

# Fit a single negative binomial model and return brms object
run_single_zinegbin_brms <- function(abs_dat, 
                                     this_gene, 
                                     this_assay,
                                     adapt_delta = 0.8,
                                     iter = 2000,
                                     seed = 123,
                                     rhat_cutoff = 1.1,
                                     run_limit = 10) {
  
  use_dat <- abs_dat %>% 
    filter(gene == this_gene,
           assay == this_assay)
  
  num_runs <- 0
  stop_now <- FALSE
  while(!stop_now) {
    # writer loss
    brm_fit <- brm(transformed_value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                   data = use_dat, 
                   family = "zero_inflated_negbinomial",
                   iter = iter,
                   control = list(adapt_delta = adapt_delta))
    brm_fixed <- summary(brm_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
    brm_random <- summary(brm_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")

    num_runs <- num_runs + 1
    if ((all(brm_fixed$Rhat < rhat_cutoff) & all(brm_fixed$Rhat < rhat_cutoff) & all(brm_fixed$Rhat < rhat_cutoff)) | num_runs == run_limit) {
      stop_now <- TRUE
    }
    else {
      seed <- seed + 1
    }
  }

  brm_fit
}

###########################################
##        
##          Figure 3
##
###########################################
### Plots of unmodeled data summaries across all genes
## DtL for absolute data (Figure 3B)
fig3b_lines <- ggplot(data = mean_abs_h3k36_dat %>%
                  filter(assay == "writer_add"), 
                aes(x = timepoint * 60, y = value)) + 
  geom_line(aes(group = gene), alpha = 0.05, col = wa_col) + 
  ylim(0, 1) +
  xlab("Time (minutes)") + ylab("Mean H3K36me3/H3 (normalized ChIP signal)") +
  plot_theme

outlier_mean_abs_h3k36_dat <- mean_abs_h3k36_dat %>%
  group_by(assay, factor(timepoint * 60)) %>%
  mutate(outlier = value > quantile(value, probs = 0.75, na.rm = TRUE) + 1.5 * IQR(value, na.rm = TRUE)) %>%
  ungroup
fig3b_boxplot <- ggplot(data = mean_abs_h3k36_dat %>%
                        filter(assay == "writer_add"), 
                      aes(x = factor(timepoint * 60), y = value)) + 
  stat_boxplot(geom = "errorbar", width = 0.33, na.rm = TRUE) +
  geom_boxplot(outlier.shape = NA) +  
  geom_point(data = function(x) outlier_mean_abs_h3k36_dat %>%
               filter(assay == "writer_add",
                      outlier), position = position_jitter(w = 0.1, h = 0), size = 0.33) +
  ylim(0, 1) +
  xlab("Time (minutes)") + ylab("Mean H3K36me3/H3\n(normalized ChIP signal)") +
  plot_theme
# fig3b_lines and fig3b_boxplot overlaid in manuscript

## LtD for relative data (Figure 3C)
fig3d_lines <- ggplot(data = mean_abs_h3k36_dat %>%
                        filter(assay == "writer_loss"), 
                      aes(x = timepoint * 60, y = value)) + 
  geom_line(aes(group = gene), alpha = 0.05, col = wl_col) + 
  ylim(0, 1) +
  xlab("Time (minutes)") + ylab("Mean H3K36me3/H3\n(normalized ChIP signal)") +
  plot_theme

fig3c_boxplot <- ggplot(data = mean_abs_h3k36_dat %>%
                          filter(assay == "writer_loss"), 
                        aes(x = factor(timepoint * 60), y = value)) + 
  stat_boxplot(geom = "errorbar", width = 0.33) +
  geom_boxplot(outlier.shape = NA) +  
  geom_point(data = function(x) outlier_mean_abs_h3k36_dat %>%
               filter(assay == "writer_loss",
                      outlier), position = position_jitter(w = 0.1, h = 0), size = 0.33) +
  ylim(0, 1) +
  xlab("Time (minutes)") + ylab("Mean H3K36me3/H3\n(normalized ChIP signal)") +
  plot_theme
# fig3c_lines and fig3c_boxplot overlaid in manuscript

## DtL for absolute data (Figure 3D)
fig3d_lines <- ggplot(data = mean_rel_h3k36_dat %>%
                        filter(assay == "writer_add"), 
                      aes(x = timepoint * 60, y = value)) + 
  geom_line(aes(group = gene), alpha = 0.05, col = wa_col) + 
  ylim(0, 1) +
  xlab("Time (minutes)") + ylab("Proportional H3K36me3/H3\n(as a fraction of each gene max)") +
  plot_theme

outlier_mean_rel_h3k36_dat <- mean_rel_h3k36_dat %>%
  group_by(assay, factor(timepoint * 60)) %>%
  mutate(outlier = value > quantile(value, probs = 0.75, na.rm = TRUE) + 1.5 * IQR(value, na.rm = TRUE) | value < quantile(value, probs = 0.25, na.rm = TRUE) - 1.5 * IQR(value, na.rm = TRUE)) %>%
  ungroup
fig3d_boxplot <- ggplot(data = mean_rel_h3k36_dat %>%
                          filter(assay == "writer_add"), 
                        aes(x = factor(timepoint * 60), y = value)) + 
  stat_boxplot(geom = "errorbar", width = 0.33) +
  geom_boxplot(outlier.shape = NA) +  
  geom_point(data = function(x) outlier_mean_rel_h3k36_dat %>%
               filter(assay == "writer_add",
                      outlier), position = position_jitter(w = 0.1, h = 0), size = 0.33) +
  ylim(0, 1) +
  xlab("Time (minutes)") + ylab("Proportional H3K36me3/H3\n(as a fraction of each gene max)") +
  plot_theme
# fig3d_lines and fig3d_boxplot overlaid in manuscript

## LtD for absolute data (Figure 3E)
fig3e_lines <- ggplot(data = mean_rel_h3k36_dat %>%
                        filter(assay == "writer_loss"), 
                      aes(x = timepoint * 60, y = value)) + 
  geom_line(aes(group = gene), alpha = 0.05, col = wl_col) + 
  ylim(0, 1) +
  xlab("Time (minutes)") + ylab("Mean H3K36me3/H3\n(as a fraction of each gene max)") +
  plot_theme
fig3e_boxplot <- ggplot(data = mean_rel_h3k36_dat %>%
                          filter(assay == "writer_loss"), 
                        aes(x = factor(timepoint * 60), y = value)) + 
  stat_boxplot(geom = "errorbar", width = 0.33) +
  geom_boxplot(outlier.shape = NA) +  
  geom_point(data = function(x) outlier_mean_rel_h3k36_dat %>%
               filter(assay == "writer_loss",
                      outlier), position = position_jitter(w = 0.1, h = 0), size = 0.33) +
  ylim(0, 1) +
  xlab("Time (minutes)") + ylab("Mean H3K36me3/H3\n(as a fraction of each gene max)") +
  plot_theme
# fig3e_lines and fig3e_boxplot overlaid in manuscript

###########################################
##        
##          Figure 4
##
###########################################
## DtL and LtD for non-relative CDS1 (Figure 4A left)
# Run the Stan model fits
cds1_negbin_writer_add_fit <- run_single_zinegbin_brms(abs_dat = abs_h3k36_dat, 
                                                       this_gene = "YBR029C",
                                                       this_assay = "writer_add",
                                                       seed = 124)
cds1_negbin_writer_loss_fit <- run_single_zinegbin_brms(abs_dat = abs_h3k36_dat, 
                                                        this_gene = "YBR029C",
                                                        this_assay = "writer_loss",
                                                        seed = 124)
# Posterior data sets
cds1_negbin_writer_add_post_dat <- conditional_effects(cds1_negbin_writer_add_fit, plot = FALSE)
cds1_negbin_writer_loss_post_dat <- conditional_effects(cds1_negbin_writer_loss_fit, plot = FALSE)

fig4a_left <- ggplot(data = abs_h3k36_dat %>%
                      filter(gene == "YBR029C",
                             assay %in% c("writer_add", "writer_loss")),
                    aes(x = timepoint_min, y = value, col = assay)) +
  geom_line(aes(group = sample_unit), linetype = "dashed") +
  scale_color_manual(values = c(wl_col, wa_col), labels = c("Light to Dark", "Dark to Light")) +
  geom_smooth(data = cds1_negbin_writer_add_post_dat$timepoint, 
              aes(x = `effect1__` * 60, y = `estimate__`/1000, ymin = `lower__`/1000, ymax = `upper__`/1000), 
              stat = "identity", se = TRUE, col = wa_col, method = "loess", size = 2) +
  geom_smooth(data = cds1_negbin_writer_loss_post_dat$timepoint,
              aes(x = `effect1__` * 60, y = `estimate__`/1000, ymin = `lower__`/1000, ymax = `upper__`/1000), 
              stat = "identity", se = TRUE, col = wl_col, method = "loess", size = 2) + 
  xlim(0, 90) +
  ylim(0, 0.5) +
  ylab("Mean H3K36me3/H3\n(normalized ChIP signal)") + xlab("Time (minutes)") +
  annotate(geom = "text", x = 85, y = 0.5, label = "CDS1") +
  plot_theme +
  theme(legend.position = c(0.2, 0.9)) +
  labs(col = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:")))) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

## DtL and LtD for relative CDS1 (Figure 4A right)
## Run the Stan model fits
cds1_beta_writer_add_fit <- run_single_beta_brms(rel_dat = rel_h3k36_dat, 
                                                 this_gene = "YBR029C",
                                                 this_assay = "writer_add",
                                                 seed = 124)
cds1_beta_writer_loss_fit <- run_single_beta_brms(rel_dat = rel_h3k36_dat, 
                                                  this_gene = "YBR029C",
                                                  this_assay = "writer_loss",
                                                  seed = 124)
# Posterior data sets
cds1_beta_writer_add_post_dat <- conditional_effects(cds1_beta_writer_add_fit, plot = FALSE)
cds1_beta_writer_loss_post_dat <- conditional_effects(cds1_beta_writer_loss_fit, plot = FALSE)

fig4a_right <- ggplot(data = rel_h3k36_dat %>%
                      filter(gene == "YBR029C",
                             assay %in% c("writer_add", "writer_loss")),
                    aes(x = timepoint_min, y = value, col = assay)) +
  geom_line(aes(group = sample_unit), linetype = "dashed") +
  scale_color_manual(values = c(wl_col, wa_col), labels = c("Light to Dark", "Dark to Light")) +
  geom_smooth(data = cds1_beta_writer_loss_post_dat$timepoint,
              aes(x = `effect1__` * 60, y = `estimate__`, ymin = `lower__`, ymax = `upper__`), 
              stat = "identity", se = TRUE, col = wl_col, method = "loess", size = 2) + 
  geom_smooth(data = cds1_beta_writer_add_post_dat$timepoint, 
              aes(x = `effect1__` * 60, y = `estimate__`, ymin = `lower__`, ymax = `upper__`), 
              stat = "identity", se = TRUE, col = wa_col, method = "loess", size = 2) +
  ylim(0, 1) + xlim(0, 90) +
  ylab("Proportional H3K36me3/H3\n(as a fraction of each gene max)") + xlab("Time (minutes)") +
  annotate(geom = "text", x = 85, y = 1, label = "CDS1") +
  plot_theme +
  theme(legend.position = c(0.75, 0.15)) +
  labs(col = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:")))) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

## Venn diagram of high confidence genes under DtL and LtD (Figure 4B)
# ZOI Beta GLMM rate data
beta_glmm_dat <- read.csv("Supplemental_Data3.csv", header = TRUE)

# High confidence genes
high_confidence_genes_down <- beta_glmm_dat %>%
  filter(assay == "writer_loss",
         category == "negative") %>%
  pull(gene)
high_confidence_genes_up <- beta_glmm_dat %>%
  filter(assay == "writer_add",
         category == "positive") %>%
  pull(gene)
high_confidence_genes <- intersect(high_confidence_genes_down, high_confidence_genes_up)
venn_list <- list("Light to Dark" = high_confidence_genes_down,
                  "Dark to Light" = high_confidence_genes_up)
fig4b <- compute.Venn(Venn(venn_list))

# Adjusting features of the Venn diagram
gp <- VennThemes(fig4b)
gp$Set$Set1$col <- wl_col
gp$FaceText$`10`$col <- wl_col
gp$SetText$Set1$col <- wl_col
gp$Set$Set2$col <- wa_col
gp$FaceText$`01`$col <- wa_col
gp$SetText$Set2$col <- wa_col

SetLabels <- VennGetSetLabels(fig4b)
SetLabels[SetLabels$Label == "Light to Dark", "x"] <- -50
SetLabels[SetLabels$Label == "Dark to Light", "x"] <- 50
SetLabels[SetLabels$Label == "Light to Dark", "y"] <- 30
SetLabels[SetLabels$Label == "Dark to Light", "y"] <- -40
fig4b <- VennSetSetLabels(fig4b, SetLabels)

plot(fig4b, gp = gp, show = list(Faces = FALSE, Universe = FALSE))

## Comparison of ZIO beta GLMM rates to ZI negbin GLMM rates (Figure 4C)
# ZI Negbin GLMM rate data
negbin_glmm_dat <- read.csv("Supplemental_Data4.csv", header = TRUE)

# Make data.frame to compare rates
compare_rate_dat <- inner_join(beta_glmm_dat %>%
  filter(assay %in% c("writer_add", "writer_loss"),
         !run_limit_stop) %>%
  dplyr::select(gene, assay, estimate, category) %>%
  dplyr::rename(beta_estimate = estimate,
                beta_category = category),
  negbin_glmm_dat %>%
    filter(assay %in% c("writer_add", "writer_loss"),
           !run_limit_stop) %>%
    dplyr::select(gene, assay, estimate, category) %>%
    dplyr::rename(negbin_estimate = estimate,
                  negbin_category = category)) %>%
  mutate(category = ifelse(beta_category != "zero" & negbin_category != "zero", "nonzero", "zero"),
         assay = factor(assay, levels = c("writer_loss", "writer_add")))
  
fig4c <- ggplot(data = compare_rate_dat,
                aes(x = negbin_estimate, y = beta_estimate, col = assay, shape = category)) +
  geom_point() +
  scale_color_manual(values = c(wl_col, wa_col), labels = c("Light to Dark", "Dark to Light")) +
  scale_shape_manual(values = c(20, 21), labels = c("High Confidence", "Low Confidence")) +
  xlim(-2, 3) +
  ylim(-4, 8) + 
  ylab("Proportional H3K36me3/H3\nGLMM Estimate") +
  xlab("Mean H3K36me3/H3\nGLMM Estimate") +
  plot_theme + 
  theme(legend.position = c(0.8, 0.3)) +
  labs(col = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:"))), shape = "")

## Histogram of ZIO beta GLMM rates of overlapping high confidence genes (Figure 4D)
fig4d <- ggplot(data = beta_glmm_dat %>%
                  filter(gene %in% high_confidence_genes,
                         assay %in% c("writer_add", "writer_loss")) %>%
                  mutate(assay = factor(assay, levels = c("writer_loss", "writer_add"))),
                aes(x = estimate, fill = assay)) + 
  geom_histogram(col = "black", bins = 90, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = "dashed", col = "gray") +
  scale_fill_manual(values = c(wl_col, wa_col), labels = c("Light to Dark", "Dark to Light")) +
  ylab("Frequency") + xlab("GLMM Estimate") +
  plot_theme +
  theme(legend.position = c(0.6, 0.5)) +
  labs(fill = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:"))))

## ZIO beta GLMM rates by mean h3k36 at time = 60, colored by RNA abundance for DtL (Figure 4E)
# Grab RNA abundance and average at time = 60
wa_rna_dat <- data.table::fread("Supplemental_Data5.txt", data.table = FALSE) %>%
  dplyr::rename(gene = GENE)
wa_rna_60min_dat <- wa_rna_dat %>%
  dplyr::select(gene, contains("60min")) %>%
  dplyr::rename(rep1 = DtL_60min_set2d_RNA_Rep1,
                rep2 = DtL_60min_set2d_RNA_Rep2) %>%
  gather(key = replicate, value = transcript, -gene) %>%
  group_by(gene) %>%
  summarize(transcript = log(mean(transcript) + 1)) %>%
  ungroup

# Grab h3k36 abundance and average at time = 60
wa_h3k36_60min_dat <- abs_h3k36_dat %>%
  filter(assay == "writer_add",
         timepoint == 1) %>%
  dplyr::select(gene, value) %>%
  group_by(gene) %>%
  summarize(h3k36_value = mean(value)) %>%
  ungroup

wa_glmm_vs_h3k36_dat <- beta_glmm_dat %>%
  filter(assay == "writer_add",
         gene %in% high_confidence_genes) %>%
  dplyr::select(gene, estimate) %>%
  left_join(wa_rna_60min_dat) %>%
  left_join(wa_h3k36_60min_dat)

fig4e <- ggplot(data = wa_glmm_vs_h3k36_dat,
                aes(x = h3k36_value, y = estimate, col = transcript)) +
  geom_point() + 
  ylim(0, 10) +
  xlim(0, 1) +
  xlab("Mean normalized H3K36me3/H3\nChIP signal") +
  ylab("GLMM Estimate") +
  scale_color_gradient(low = "yellow", high = "red") +
  plot_theme +
  theme(legend.position = c(x = 0.75, y = 0.25), legend.direction = "horizontal") +
  guides(color = guide_colorbar(title.position = "top", 
                                title = "RNA Abundance\n(log TPM) at t = 60 min")) +
  annotate(geom = "text", x = 0.25, y = 2, label = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:"), " Dark -> Light")))

## ZIO beta GLMM rates by mean h3k36 at time = 0, colored by RNA abundance for LtD (Figure 4F)
wl_rna_dat <- data.table::fread("Supplemental_Data6.txt", data.table = FALSE) %>%
  dplyr::rename(gene = GENE)

# Grab RNA abundance and average at time = 0
wl_rna_0min_dat <- wl_rna_dat %>%
  dplyr::select(gene, contains("_0min")) %>%
  dplyr::rename(rep1 = LtD_0min_set2d_RNA_Rep1,
                rep2 = LtD_0min_set2d_RNA_Rep2,
                rep3 = LtD_0min_set2d_RNA_Rep3) %>%
  gather(key = replicate, value = transcript, -gene) %>%
  group_by(gene) %>%
  summarize(transcript = log(mean(transcript) + 1)) %>%
  ungroup

# Grab h3k36 abundance and average at time = 0
wl_h3k36_0min_dat <- abs_h3k36_dat %>%
  filter(assay == "writer_loss",
         timepoint == 0) %>%
  dplyr::select(gene, value) %>%
  group_by(gene) %>%
  summarize(h3k36_value = mean(value)) %>%
  ungroup

wl_glmm_vs_h3k36_dat <- beta_glmm_dat %>%
  filter(assay == "writer_loss",
         gene %in% high_confidence_genes) %>%
  dplyr::select(gene, estimate) %>%
  left_join(wl_rna_0min_dat) %>%
  left_join(wl_h3k36_0min_dat)

fig4f <- ggplot(data = wl_glmm_vs_h3k36_dat,
                aes(x = h3k36_value, y = estimate, col = transcript)) +
  geom_point() + 
  xlab("Mean normalized H3K36me3/H3\nChIP signal") +
  ylab("GLMM Estimate") +
  scale_y_reverse(lim = c(0, -5)) +
  scale_color_gradient(low = "yellow", high = "red") +
  plot_theme +
  xlim(0, 1) +
  theme(legend.position = c(x = 0.75, y = 0.25), legend.direction = "horizontal") +
  guides(color = guide_colorbar(title.position = "top", 
                                title = "RNA Abundance\n(log TPM) at t = 0 min")) +
  annotate(geom = "text", x = 0.25, y = -1, label = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:"), " Light to Dark")))

## ZIO beta GLMM rates by RNA abundance for dark to light (Figure 4G)
wa_glmm_vs_rna_dat <- beta_glmm_dat %>%
  filter(assay == "writer_add",
         gene %in% high_confidence_genes) %>%
  dplyr::select(gene, estimate) %>%
  left_join(wa_rna_60min_dat)

fig4g <- ggplot(data = wa_glmm_vs_rna_dat,
                aes(x = transcript, y = estimate)) +
  geom_point(col = wa_col, alpha = 0.4) + 
  geom_smooth(aes(y = estimate, x = transcript), method = "lm", se = FALSE, size = 1.5, linetype = "dashed", col = "gray") + 
  ylim(3, 7) +
  xlab("RNA Abundance (log TPM) at t = 60 min") +
  ylab("GLMM Estimate") +
  annotate(geom = "text", x = 9, y = 7, label = paste("r =", round(cor(wa_glmm_vs_rna_dat$estimate, wa_glmm_vs_rna_dat$transcript), 3))) +
  plot_theme +
  annotate(geom = "text", x = 3, y = 7, label = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:"), " Dark -> Light")))

## ZIO beta GLMM rates by RNA abundance for light to dark (Figure 4H)
wl_glmm_vs_rna_dat <- beta_glmm_dat %>%
  filter(assay == "writer_loss",
         gene %in% high_confidence_genes) %>%
  dplyr::select(gene, estimate) %>%
  left_join(wl_rna_0min_dat)

fig4h <- ggplot(data = wl_glmm_vs_rna_dat,
                aes(x = transcript, y = estimate)) +
  geom_point(col = wl_col, alpha = 0.4) + 
  geom_smooth(aes(y = estimate, x = transcript), method = "lm", se = FALSE, size = 1.5, linetype = "dashed", col = "gray") + 
  scale_y_reverse(lim = c(-1.5, -3.5)) +
  xlab("RNA Abundance (log TPM) at t = 0 min") +
  ylab("GLMM Estimate") +
  annotate(geom = "text", x = 0.5, y = -3.5, label = paste("r =", round(cor(wl_glmm_vs_rna_dat$estimate, wl_glmm_vs_rna_dat$transcript), 3))) +
  plot_theme +
  annotate(geom = "text", x = 5, y = -1.75, label = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:"), " Light to Dark")))


###########################################
##        
##          Figure 6
##
###########################################
## LtD with and without Rph1 for absolute data (Figure 6B)
fig6b_lines <- ggplot(data = mean_abs_h3k36_dat %>%
                        filter(assay %in% c("writer_loss", "writer_eraser_loss")) %>%
                        mutate(unit = paste(gene, assay)), 
                      aes(x = timepoint * 60, y = value, col = assay)) + 
  geom_line(aes(group = unit), alpha = 0.05) + 
  xlab("Time (minutes)") + ylab("Mean H3K36me3/H3\n(normalized ChIP signal)") +
  scale_color_manual(values = c(wl_col, wel_col), 
                     labels = c(expression(paste(italic("set2"), Delta, " + LANS-Set2")),
                                expression(paste(italic("set2"), Delta, italic("rph1"), Delta, " + LANS-Set2")))) +
  theme(legend.position = c(x = 0.5, y = 0.9)) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  labs(color = "") +
  plot_theme

fig6b_boxplot <- ggplot(data = mean_abs_h3k36_dat %>%
                          filter(assay %in% c("writer_loss", "writer_eraser_loss")) %>%
                          mutate(assay = factor(assay, levels = c("writer_loss", "writer_eraser_loss")),
                                 unit = paste(gene, assay)),
                          aes(x = factor(timepoint * 60), y = value, col = assay)) + 
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +  
  geom_point(data = function(x) outlier_mean_abs_h3k36_dat %>%
               filter(assay %in% c("writer_loss", "writer_eraser_loss"),
                      outlier),  
             position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.5), size = 0.33) +
  ylim(0, 1) +
  xlab("Time (minutes)") + ylab("Mean H3K36me3/H3\n(normalized ChIP signal)") +
  scale_color_manual(values = c(wl_col, wel_col), 
                     labels = c(expression(paste(italic("set2"), Delta, " + LANS-Set2")),
                                expression(paste(italic("set2"), Delta, italic("rph1"), Delta, " + LANS-Set2")))) +
  theme(legend.position = c(x = 0.5, y = 0.9)) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  labs(color = "") +
  plot_theme
# fig6b_lines and fig6b_boxplot overlaid in manuscript

## LtD with and without Rph1 for relative data (Figure 6C)
fig6c_lines <- ggplot(data = mean_rel_h3k36_dat %>%
                        filter(assay %in% c("writer_loss", "writer_eraser_loss")) %>%
                        mutate(unit = paste(gene, assay)), 
                      aes(x = timepoint * 60, y = value, col = assay)) + 
  geom_line(aes(group = unit), alpha = 0.05) + 
  xlab("Time (minutes)") + ylab("Proportional H3K36me3/H3\n(as a fraction of each gene max)") +
  scale_color_manual(values = c(wl_col, wel_col), 
                     labels = c(expression(paste(italic("set2"), Delta, " + LANS-Set2")),
                                expression(paste(italic("set2"), Delta, italic("rph1"), Delta, " + LANS-Set2")))) +
  theme(legend.position = c(x = 0.2, y = 0.2)) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  labs(color = "") +
  plot_theme

fig6c_boxplot <- ggplot(data = mean_rel_h3k36_dat %>%
                          filter(assay %in% c("writer_loss", "writer_eraser_loss")) %>%
                          mutate(assay = factor(assay, levels = c("writer_loss", "writer_eraser_loss")),
                                 unit = paste(gene, assay)),
                        aes(x = factor(timepoint * 60), y = value, col = assay)) + 
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +  
  geom_point(data = function(x) outlier_mean_rel_h3k36_dat %>%
               filter(assay %in% c("writer_loss", "writer_eraser_loss"),
                      outlier),  
             position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.5), size = 0.33) +
  ylim(0, 1) +
  xlab("Time (minutes)") + ylab("Proportional H3K36me3/H3\n(as a fraction of each gene max)") +
  scale_color_manual(values = c(wl_col, wel_col), 
                     labels = c(expression(paste(italic("set2"), Delta, " + LANS-Set2")),
                                expression(paste(italic("set2"), Delta, italic("rph1"), Delta, " + LANS-Set2")))) +
  theme(legend.position = c(x = 0.2, y = 0.2)) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  labs(color = "") +
  plot_theme
# fig6c_lines and fig6c_boxplot overlaid in manuscript

## Venn diagram of high confidence genes under LtD with and without Rph1 (Figure 6D)
# High confidence genes
high_confidence_genes_down_eraser <- beta_glmm_dat %>%
  filter(assay == "writer_eraser_loss",
         category == "negative") %>%
  pull(gene)
venn_list <- list("Set2 inactivation" = high_confidence_genes_down,
                  "Set2 inactivation without Rph1" = high_confidence_genes_down_eraser)
fig6d <- compute.Venn(Venn(venn_list))

# Adjusting features of the Venn diagram
gp <- VennThemes(fig6d)
gp$Set$Set1$col <- wl_col
gp$FaceText$`10`$col <- wl_col
gp$SetText$Set1$col <- wl_col
gp$Set$Set2$col <- wel_col
gp$FaceText$`01`$col <- wel_col
gp$SetText$Set2$col <- wel_col

SetLabels <- VennGetSetLabels(fig6d)
SetLabels[SetLabels$Label == "Set2 inactivation", "x"] <- -50
SetLabels[SetLabels$Label == "Set2 inactivation", "y"] <- -45
SetLabels[SetLabels$Label == "Set2 inactivation without Rph1", "x"] <- 50
SetLabels[SetLabels$Label == "Set2 inactivation without Rph1", "y"] <- -45
fig6d <- VennSetSetLabels(fig6d, SetLabels)

plot(fig6d, gp = gp, show = list(Faces = FALSE, Universe = FALSE))

## Histogram of ZIO beta GLMM rates of overlapping high confidence genes (Figure 6E)
high_confidence_genes_down_both <- intersect(high_confidence_genes_down_eraser, high_confidence_genes_down)

fig6e <- ggplot(data = beta_glmm_dat %>%
                  filter(gene %in% high_confidence_genes_down_both,
                         assay %in% c("writer_loss", "writer_eraser_loss")) %>%
                  mutate(assay = factor(assay, levels = c("writer_loss", "writer_eraser_loss"))),
                aes(x = estimate, fill = assay, alpha = 0.6)) + 
  geom_histogram(col = "black", bins = 50, size = 0.5) + 
  scale_x_reverse(lim = c(0, -4)) +
  scale_fill_manual(values = c(wl_col, wel_col)) +
  ylab("Number of Genes") + xlab("GLMM Estimate") +
  plot_theme +
  theme(legend.position = c(0.6, 0.5)) +
  guides(fill = FALSE, alpha = FALSE)

## Difference in modeled timepoints between writer loss and writer eraser loss (Figure 6F)
# Parameter estimates data for model comparing writer loss and writer eraser loss
loss_compare_dat <- read.csv("Supplemental_Data7.csv", header = TRUE)

loss_compare_plot_dat <- loss_compare_dat %>%
  filter(grepl(x = parameter, pattern = "assay")) %>%
  mutate(parameter = recode(parameter, 
                            'assaywriter_eraser_loss' = "0",
                            'timepoint_cat1:assaywriter_eraser_loss' = "1",
                            'timepoint_cat2:assaywriter_eraser_loss' = "2",
                            'timepoint_cat3:assaywriter_eraser_loss' = "3")) %>%
  dplyr::rename(timepoint = parameter) %>%
  mutate(time = recode(timepoint, 
                       '0' = "0",
                       '1' = "30",
                       '2' = "60",
                       '3' = "90"))
  
## Plot all negative trend genes
fig6f <- ggplot(data = loss_compare_plot_dat %>% 
                              filter(gene %in% high_confidence_genes_down_both), 
                            aes(x = time, y = Estimate)) + 
  geom_line(aes(group = gene, col = gene), alpha = 0.1) + 
  geom_boxplot(aes(x = time, y = Estimate), outlier.alpha = 0, color = wel_col, size = 0.7, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  scale_color_grey(guide = FALSE) + 
  ylab(expression(paste(italic("set2"), Delta, italic("rph1"), Delta, " - ", italic("set2"), Delta, " Timepoint Effects"))) +
  xlab("Timepoint (min)") +
  ylim(-5, 5) + 
  stat_summary(fun = mean, geom="line", aes(group=1), col = wel_col, size = 1) + 
  plot_theme


###########################################
##        
##          Supplemental Figure S4
##
###########################################
## DtL and LtD for non-relative YAR009C (Supplemental Figure S4A left)
# Run the actual fits
yar009c_negbin_writer_add_fit <- run_single_zinegbin_brms(abs_dat = abs_h3k36_dat, 
                                                       this_gene = "YAR009C",
                                                       this_assay = "writer_add",
                                                       seed = 125)
yar009c_negbin_writer_loss_fit <- run_single_zinegbin_brms(abs_dat = abs_h3k36_dat, 
                                                        this_gene = "YAR009C",
                                                        this_assay = "writer_loss",
                                                        seed = 125)
# Posterior data sets
yar009c_negbin_writer_add_post_dat <- conditional_effects(yar009c_negbin_writer_add_fit, plot = FALSE)
yar009c_negbin_writer_loss_post_dat <- conditional_effects(yar009c_negbin_writer_loss_fit, plot = FALSE)

figs4a_left <- ggplot(data = abs_h3k36_dat %>%
                       filter(gene == "YAR009C",
                              assay %in% c("writer_add", "writer_loss")),
                     aes(x = timepoint_min, y = value, col = assay)) +
  geom_line(aes(group = sample_unit), linetype = "dashed") +
  scale_color_manual(values = c(wl_col, wa_col), labels = c("Light to Dark", "Dark to Light")) +
  geom_smooth(data = yar009c_negbin_writer_add_post_dat$timepoint, 
              aes(x = `effect1__` * 60, y = `estimate__`/1000, ymin = `lower__`/1000, ymax = `upper__`/1000), 
              stat = "identity", se = TRUE, col = wa_col, method = "loess", size = 2) +
  xlim(0, 90) +
  ylim(0, 1) +
  ylab("Mean H3K36me3/H3\n(normalized ChIP signal)") + xlab("Time (minutes)") +
  annotate(geom = "text", x = 85, y = 1, label = "YAR009C") +
  plot_theme +
  theme(legend.position = c(0.8, 0.4)) +
  labs(col = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:")))) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

## LtD for non-relative YAR009C blown up (Supplemental Figure S4A right)
figs4a_right <- ggplot(data = abs_h3k36_dat %>%
                        filter(gene == "YAR009C",
                               assay %in% c("writer_loss")),
                      aes(x = timepoint_min, y = value, col = assay)) +
  geom_smooth(data = yar009c_negbin_writer_loss_post_dat$timepoint %>%
                mutate(assay = "writer_loss"), 
              aes(x = `effect1__` * 60, y = `estimate__`/1000, ymin = `lower__`/1000, ymax = `upper__`/1000, col = assay), 
            stat = "identity", se = TRUE, method = "loess", size = 2) +
  scale_color_manual(values = wl_col, labels = "Light to Dark") +
  xlim(0, 90) +
  ylim(0, 25) +
  ylab("Mean H3K36me3/H3\n(normalized ChIP signal)") + xlab("Time (minutes)") +
  annotate(geom = "text", x = 85, y = 25, label = "YAR009C") +
  plot_theme +
  theme(legend.position = c(0.25, 0.9)) +
  labs(col = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:")))) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

## DtL and LtD for relative YAR009C (Supplemental Figure S4B)
# Run the actual fits
yar009c_beta_writer_add_fit <- run_single_beta_brms(rel_dat = rel_h3k36_dat, 
                                                    this_gene = "YAR009C",
                                                    this_assay = "writer_add",
                                                    seed = 123)
yar009c_beta_writer_loss_fit <- run_single_beta_brms(rel_dat = rel_h3k36_dat, 
                                                     this_gene = "YAR009C",
                                                     this_assay = "writer_loss",
                                                     seed = 123)
# Posterior data sets
yar009c_beta_writer_add_post_dat <- conditional_effects(yar009c_beta_writer_add_fit, plot = FALSE)
yar009c_beta_writer_loss_post_dat <- conditional_effects(yar009c_beta_writer_loss_fit, plot = FALSE)

figs4b <- ggplot(data = rel_h3k36_dat %>%
                        filter(gene == "YAR009C",
                               assay %in% c("writer_add", "writer_loss")),
                      aes(x = timepoint_min, y = value, col = assay)) +
  geom_line(aes(group = sample_unit), linetype = "dashed") +
  scale_color_manual(values = c(wl_col, wa_col), labels = c("Light to Dark", "Dark to Light")) +
  geom_smooth(data = yar009c_beta_writer_add_post_dat$timepoint, 
              aes(x = `effect1__` * 60, y = `estimate__`, ymin = `lower__`, ymax = `upper__`), 
              stat = "identity", se = TRUE, col = wa_col, method = "loess", size = 2) +
  geom_smooth(data = yar009c_beta_writer_loss_post_dat$timepoint, 
              aes(x = `effect1__` * 60, y = `estimate__`, ymin = `lower__`, ymax = `upper__`), 
              stat = "identity", se = TRUE, col = wl_col, method = "loess", size = 2) +
  xlim(0, 90) +
  ylim(0, 1) +
  ylab("Proportional H3K36me3/H3\n(as a fraction of gen max)") + xlab("Time (minutes)") +
  annotate(geom = "text", x = 85, y = 0.9, label = "YAR009C") +
  plot_theme +
  theme(legend.position = c(0.8, 0.4)) +
  labs(col = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:")))) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

## Comparison of ZIO beta GLMM rates to their error for all genes (Supplemental Figure S4C)
figs4c <- ggplot(data = beta_glmm_dat %>%
                   filter(assay %in% c("writer_loss", "writer_add")) %>%
                   mutate(assay = factor(assay, levels = c("writer_loss", "writer_add")),
                          category = ifelse(category != "zero", "nonzero", "zero")),
                 aes(x = est_error, y = estimate, col = assay, shape = category)) +
  geom_point() +
  scale_color_manual(values = c(wa_col, wl_col), labels = c("Light to Dark", "Dark to Light")) +
  scale_shape_manual(values = c(20, 21), labels = c("High Confidence", "Low Confidence")) +
  ylab("Proportional H3K36me3/H3\nGLMM Estimate") +
  xlab("Proportional H3K36me3/H3\nGLMM Estimate Error") +
  plot_theme +
  theme(legend.position = c(0.8, 0.3)) +
  labs(col = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:"))), shape = "")

## Comparison of ZIO beta GLMM rates to ZI negbin GLMM rates, all genes included (Supplemental Figure S4D)
figs4d <- ggplot(data = compare_rate_dat %>%
                   mutate(assay = factor(assay, levels = c("writer_loss", "writer_add"))),
                 aes(x = negbin_estimate, y = beta_estimate, col = assay, shape = category)) +
  geom_point() +
  scale_color_manual(values = c(wl_col, wa_col), labels = c("Light to Dark", "Dark to Light")) +
  scale_shape_manual(values = c(20, 21), labels = c("High Confidence", "Low Confidence")) +
  ylab("Proportional H3K36me3/H3\nGLMM Estimate") +
  xlab("Mean H3K36me3/H3\nGLMM Estimate") +
  plot_theme +
  theme(legend.position = c(0.8, 0.3)) +
  labs(col = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:"))), shape = "")

## Histogram of ZIO beta GLMM rates of overlapping genes (Supplemental Figure S4E)
figs4e <- ggplot(data = beta_glmm_dat %>%
                   filter(assay %in% c("writer_add", "writer_loss")) %>%
                   mutate(assay = factor(assay, levels = c("writer_loss", "writer_add"))),
                 aes(x = estimate, fill = assay)) + 
  geom_histogram(col = "black", bins = 50, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = "dashed", col = "gray") +
  scale_fill_manual(values = c(wl_col, wa_col), labels = c("Light to Dark", "Dark to Light")) +
  ylab("Frequency") + xlab("GLMM Estimate") +
  plot_theme +
  theme(legend.position = c(0.6, 0.5)) +
  labs(fill = expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:"))))

## Scatterplot of ZIO beta GLMM rates for LtD vs DtL of high confidence genes (Supplemental Figure S4F)
wa_vs_wl_beta_dat <- beta_glmm_dat %>%
  filter(assay %in% c("writer_add", "writer_loss"),
         gene %in% high_confidence_genes) %>%
  dplyr::select(gene, assay, estimate) %>%
  spread(key = assay, value = estimate)

figs4f <- ggplot(data = wa_vs_wl_beta_dat,
                 aes(y = writer_loss, x = writer_add)) + 
  geom_point(col = "black", alpha = 0.5) + 
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", col = "gray") +
  scale_y_reverse() +
  ylab("Proportional H3K36me3/H3\nGLMM Estimate\nLight to Dark") + xlab("Proportional H3K36me3/H3\nGLMM Estimate\nDark to Light") +
  annotate(geom = "text", x = 3.8, y = -3.4, label = paste("r =", round(cor(wa_vs_wl_beta_dat$writer_add, wa_vs_wl_beta_dat$writer_loss), 3))) +
  plot_theme 

## Scatterplot of mean RNA abundance for LtD vs mean RNA abundance for DtL of high confidence genes (Supplemental Figure S4G)
# Grab mean RNA abundances at time = 0 for LtD and time = 60 for DtL
wa_vs_wl_rna_comparison_dat <- inner_join(wl_rna_0min_dat %>%
                                      dplyr::select(gene, wl_min0 = transcript),
                                    wa_rna_60min_dat %>%
                                      dplyr::select(gene, wa_min60 = transcript))

figs4g <- ggplot(data = wa_vs_wl_rna_comparison_dat %>%
                   filter(gene %in% high_confidence_genes),
                 aes(y = wl_min0, x = wa_min60)) + 
  geom_point(col = "black", alpha = 0.5) + 
  ylab("Mean RNA Abundance (t = 0)\n(log TPM)\nLight to Dark") + xlab("Mean RNA Abundance (t = 60)\n(log TPM)\nDark to Light") +
  plot_theme 

## Scatterplot of mean RNA abundance vs mean trimetyl for high confidence genes at DtL (Supplemental Figure S4J left)
wa_rna_vs_h3k36_comparison_dat <- inner_join(wa_h3k36_60min_dat %>%
                                            dplyr::select(gene, h3k36_value),
                                          wa_rna_60min_dat %>%
                                            dplyr::select(gene, transcript)) %>%
  filter(gene %in% high_confidence_genes)

figs4j_left <- ggplot(data = wa_rna_vs_h3k36_comparison_dat,
                      aes(y = transcript, x = h3k36_value)) + 
  geom_point(col = wa_col, alpha = 0.1) + 
  xlim(0, 1) +
  ylab("Mean RNA Abundance (t = 60 min)\n(log TPM)") + xlab("Mean H3K36me3/H3 Signal (t = 60 min)") +
  annotate(geom = "text", x = 0.9, y = 10, label = paste("r =", round(cor(wa_rna_vs_h3k36_comparison_dat$transcript, wa_rna_vs_h3k36_comparison_dat$h3k36_value, use = "pairwise.complete.obs"), 3))) +
  plot_theme 

## Scatterplot of mean RNA abundance vs mean trimetyl for high confidence genes at LtD (Supplemental Figure S4J right)
wl_rna_vs_h3k36_comparison_dat <- inner_join(wl_h3k36_0min_dat %>%
                                               dplyr::select(gene, h3k36_value),
                                            wl_rna_0min_dat %>%
                                              dplyr::select(gene, transcript)) %>%
  filter(gene %in% high_confidence_genes)

figs4j_right <- ggplot(data = wl_rna_vs_h3k36_comparison_dat, 
                       aes(y = transcript, x = h3k36_value)) + 
  geom_point(col = wl_col, alpha = 0.1) + 
  xlim(0, 1) +
  ylab("Mean RNA Abundance (t = 0 min)\n(log TPM)") + xlab("Mean H3K36me3/H3 Signal (t = 0 min)") +
  annotate(geom = "text", x = 0.9, y = 10, label = paste("r =", round(cor(wl_rna_vs_h3k36_comparison_dat$transcript, wl_rna_vs_h3k36_comparison_dat$h3k36_value, use = "pairwise.complete.obs"), 3))) +
  plot_theme 

## Scatterplot of gene length vs mean h3k36 levels for high confidence genes at DtL (Supplemental Figure S4K left)
gene_length_dat <- read.table("Supplemental_Data8.txt", header = FALSE) %>%
  dplyr::rename(gene = V4) %>%
  dplyr::mutate(gene_length = V3 - V2) %>%
  dplyr::select(gene, gene_length)

wa_genelength_vs_h3k36_comparison_dat <- inner_join(wa_h3k36_60min_dat,
                                                    gene_length_dat) %>%
  filter(gene %in% high_confidence_genes)

figs4k_left <- ggplot(data = wa_genelength_vs_h3k36_comparison_dat,
                      aes(y = gene_length, x = h3k36_value)) + 
  geom_point(col = wa_col, alpha = 0.1) + 
  xlim(0, 1) +
  ylim(0, 15000) +
  ylab("Gene Length (bp)") + xlab("Mean H3K36me3/H3 Signal (t = 60 min)") +
  annotate(geom = "text", x = 0.9, y = 15000, label = paste("r =", round(cor(wa_genelength_vs_h3k36_comparison_dat$h3k36_value, wa_genelength_vs_h3k36_comparison_dat$gene_length), 3))) +
  plot_theme 

## Scatterplot of gene length vs mean h3k36 levels for high confidence genes at LtD (Supplemental Figure S4K right)
wl_genelength_vs_h3k36_comparison_dat <- inner_join(wl_h3k36_0min_dat,
                                                    gene_length_dat) %>%
  filter(gene %in% high_confidence_genes)

figs4k_right <- ggplot(data = wa_genelength_vs_h3k36_comparison_dat, 
                       aes(y = gene_length, x = h3k36_value)) + 
  geom_point(col = wl_col, alpha = 0.1) + 
  xlim(0, 1) +
  ylim(0, 15000) +
  ylab("Gene Length (bp)") + xlab("Mean H3K36me3/H3 Signal (t = 0 min)") +
  annotate(geom = "text", x = 0.9, y = 15000, label = paste("r =", round(cor(wl_genelength_vs_h3k36_comparison_dat$h3k36_value, wl_genelength_vs_h3k36_comparison_dat$gene_length), 3))) +
  plot_theme 

## Scatterplot of gene length vs mean RNA abundance for high confidence genes at DtL (Supplemental Figure S4L left)
wa_genelength_vs_rna_comparison_dat <- inner_join(wa_rna_60min_dat,
                                                  gene_length_dat) %>%
  filter(gene %in% high_confidence_genes)

figs4l_left <- ggplot(data = wa_genelength_vs_rna_comparison_dat,
                      aes(y = gene_length, x = transcript)) + 
  geom_point(col = wa_col, alpha = 0.1) + 
  xlim(0, 11) +
  ylim(0, 15000) +
  ylab("Gene Length (bp)") + xlab("Mean RNA Abundance (t = 60 min)\n(log TPM)") +
  annotate(geom = "text", x = 10, y = 15000, label = paste("r =", round(cor(wa_genelength_vs_rna_comparison_dat$transcript, wa_genelength_vs_rna_comparison_dat$gene_length), 3))) +
  plot_theme 

## Scatterplot of gene length vs mean RNA abundance for high confidence genes at LtD (Supplemental Figure S4L right)
wl_genelength_vs_rna_comparison_dat <- inner_join(wl_rna_0min_dat,
                                                  gene_length_dat) %>%
  filter(gene %in% high_confidence_genes)

figs4l_right <- ggplot(data = wl_genelength_vs_rna_comparison_dat, 
                       aes(y = gene_length, x = transcript)) + 
  geom_point(col = wl_col, alpha = 0.1) + 
  xlim(0, 11) +
  ylim(0, 15000) +
  ylab("Gene Length (bp)") + xlab("Mean RNA Abundance (t = 0 min)\n(log TPM)") +
  annotate(geom = "text", x = 10, y = 15000, label = paste("r =", round(cor(wl_genelength_vs_rna_comparison_dat$transcript, wl_genelength_vs_rna_comparison_dat$gene_length, use = "pairwise.complete.obs"), 3))) +
  plot_theme 

## ZIO beta GLMM rates by gene length, colored by RNA abundance for DtL (Supplemental Figure S4M)
# Grab GLMM estimates, errors, average h3k36 and averate RNA at time = 60 for DtL, and gene length
wa_full_60min_dat <- beta_glmm_dat %>%
  filter(assay == "writer_add",
         gene %in% high_confidence_genes) %>%
  dplyr::select(gene, estimate, est_error) %>%
  inner_join(wa_rna_60min_dat) %>%
  inner_join(wa_h3k36_60min_dat) %>%
  inner_join(gene_length_dat)

figs4m <- ggplot(data = wa_full_60min_dat,
                aes(x = gene_length, y = estimate, col = transcript)) +
  geom_point() + 
  ylim(0, 10) +
  xlim(0, 15000) +
  xlab("Gene Length (bp)") +
  ylab("GLMM Estimate\nDark to Light") +
  scale_color_gradient(low = "yellow", high = "red") +
  plot_theme +
  theme(legend.position = c(x = 0.75, y = 0.25), legend.direction = "vertical") +
  guides(color = guide_colorbar(title.position = "top", 
                                title = "RNA Abundance\n(t = 60 min)\n(log TPM)"))

## ZIO beta GLMM rate errors by gene length, colored by RNA abundance for DtL (Supplemental Figure S4N)
figs4n <- ggplot(data = wa_full_60min_dat,
                 aes(x = gene_length, y = est_error, col = transcript)) +
  geom_point() + 
  ylim(0, 3) +
  xlim(0, 15000) +
  xlab("Gene Length (bp)") +
  ylab("GLMM Estimate Error\nDark to Light") +
  scale_color_gradient(low = "yellow", high = "red") +
  plot_theme +
  theme(legend.position = c(x = 0.75, y = 0.25), legend.direction = "vertical") +
  guides(color = guide_colorbar(title.position = "top", 
                                title = "RNA Abundance\n(t = 60 min)\nlog TPM"))

## ZIO beta GLMM rates by mean h3k36 at time = 60, colored by RNA abundance for DtL (Supplemental Figure S4O)
figs4o <- ggplot(data = wa_full_60min_dat,
                 aes(x = h3k36_value, y = estimate, col = transcript)) +
  geom_point() + 
  ylim(0, 10) +
  xlim(0, 1) +
  xlab("Mean H3K36me3/H3 Signal (t = 60 min)") +
  ylab("GLMM Estimate\nDark to Light") +
  scale_color_gradient(low = "yellow", high = "red") +
  plot_theme +
  theme(legend.position = c(x = 0.75, y = 0.25), legend.direction = "vertical") +
  guides(color = guide_colorbar(title.position = "top", 
                                title = "RNA Abundance\n(t = 60 min)\nlog TPM"))

## ZIO beta GLMM rate errors by mean h3k36 at time = 60, colored by RNA abundance for DtL (Supplemental Figure S4P)
figs4p <- ggplot(data = wa_full_60min_dat,
                 aes(x = h3k36_value, y = est_error, col = transcript)) +
  geom_point() + 
  ylim(0, 3) +
  xlim(0, 1) +
  xlab("Mean H3K36me3/H3 Signal (t = 60 min)") +
  ylab("GLMM Estimate Error\nDark to Light") +
  scale_color_gradient(low = "yellow", high = "red") +
  plot_theme +
  theme(legend.position = c(x = 0.75, y = 0.25), legend.direction = "vertical") +
  guides(color = guide_colorbar(title.position = "top", 
                                title = "RNA Abundance\n(t = 60 min)\nlog TPM"))

## ZIO beta GLMM rates by gene length, colored by RNA abundance for LtD (Supplemental Figure S4Q)
# Grab GLMM estimates, errors, average h3k36 and averate RNA at time = 0 for LtD, and gene length
wl_full_0min_dat <- beta_glmm_dat %>%
  filter(assay == "writer_loss",
         gene %in% high_confidence_genes) %>%
  dplyr::select(gene, estimate, est_error) %>%
  inner_join(wl_rna_0min_dat) %>%
  inner_join(wl_h3k36_0min_dat) %>%
  inner_join(gene_length_dat)

figs4q <- ggplot(data = wl_full_0min_dat,
                 aes(x = gene_length, y = estimate, col = transcript)) +
  geom_point() + 
  scale_y_reverse(lim = c(0, -5)) +
  xlim(0, 15000) +
  xlab("Gene Length (bp)") +
  ylab("GLMM Estimate\nLight to Dark") +
  scale_color_gradient(low = "yellow", high = "red") +
  plot_theme +
  theme(legend.position = c(x = 0.75, y = 0.25), legend.direction = "vertical") +
  guides(color = guide_colorbar(title.position = "top", 
                                title = "RNA Abundance\n(t = 0 min)\nlog TPM"))

## ZIO beta GLMM rate errors by gene length, colored by RNA abundance for LtD (Supplemental Figure S4R)
figs4r <- ggplot(data = wl_full_0min_dat,
                 aes(x = gene_length, y = est_error, col = transcript)) +
  geom_point() + 
  ylim(0, 3) +
  xlim(0, 15000) +
  xlab("Gene Length (bp)") +
  ylab("GLMM Estimate Error\nLight to Dark") +
  scale_color_gradient(low = "yellow", high = "red") +
  plot_theme +
  theme(legend.position = c(x = 0.75, y = 0.25), legend.direction = "vertical") +
  guides(color = guide_colorbar(title.position = "top", 
                                title = "RNA Abundance\n(t = 0 min)\nlog TPM"))

## ZIO beta GLMM rates by mean h3k36 at time = 60, colored by RNA abundance for LtD (Supplemental Figure S4S)
figs4s <- ggplot(data = wl_full_0min_dat,
                 aes(x = h3k36_value, y = estimate, col = transcript)) +
  geom_point() + 
  scale_y_reverse(lim = c(0, -5)) +
  xlim(0, 1) +
  xlab("Mean H3K36me3/H3 Signal (t = 0 min)") +
  ylab("GLMM Estimate\nLight to Dark") +
  scale_color_gradient(low = "yellow", high = "red") +
  plot_theme +
  theme(legend.position = c(x = 0.75, y = 0.25), legend.direction = "vertical") +
  guides(color = guide_colorbar(title.position = "top", 
                                title = "RNA Abundance\n(t = 0 min)\nlog TPM"))

## ZIO beta GLMM rate errors by mean h3k36 at time = 60, colored by RNA abundance for LtD (Supplemental Figure S4T)
figs4t <- ggplot(data = wl_full_0min_dat,
                 aes(x = h3k36_value, y = est_error, col = transcript)) +
  geom_point() + 
  ylim(0, 3) +
  xlim(0, 1) +
  xlab("Mean H3K36me3/H3 Signal (t = 0 min)") +
  ylab("GLMM Estimate Error\nLight to Dark") +
  scale_color_gradient(low = "yellow", high = "red") +
  plot_theme +
  theme(legend.position = c(x = 0.75, y = 0.25), legend.direction = "vertical") +
  guides(color = guide_colorbar(title.position = "top", 
                                title = "RNA Abundance\n(t = 0 min)\nlog TPM"))

###########################################
##        
##          Supplemental Figure S6
##
###########################################
## Scatterplot of mean h3k36 for LtD - RPH1 vs mean RNA h3k36 for LtD (Supplemental Figure S6A)
# Grab mean RNA abundances at time = 0 for LtD and time = 60 for DtL
wl_vs_wel_h3k36_comparison_dat <- mean_abs_h3k36_dat %>%
  filter(assay %in% c("writer_loss", "writer_eraser_loss"),
                      timepoint == 0) %>%
  dplyr::select(gene, assay, value) %>%
  spread(key = assay, value = value)

figs6a <- ggplot(data = wl_vs_wel_h3k36_comparison_dat %>%
                   filter(gene %in% high_confidence_genes_down_both),
                 aes(y = writer_eraser_loss, x = writer_loss)) + 
  geom_point(col = "black", size = 0.9) + 
  geom_abline(intercept = 0, slope = 1, col = "gray", linetype = "dashed") +
  xlim(0, 1.2) + ylim(0, 1.2) +
  ylab(expression(paste(italic("set2"), Delta, italic("rph1"), Delta, " + LANS-Set2 Mean H3K36me3/H3"))) + 
  xlab(expression(paste(italic("set2"), Delta, " + LANS-Set2 Mean H3K36me3/H3"))) + 
  plot_theme 

## LtD and LtD + RPH1delta for non-relative PCL6 (Supplemental Figure 6B)
# Run the Stan model fits
pcl6_negbin_writer_loss_fit <- run_single_zinegbin_brms(abs_dat = abs_h3k36_dat, 
                                                        this_gene = "YAL012W",
                                                        this_assay = "writer_loss",
                                                        seed = 125)
pcl6_negbin_writer_eraser_loss_fit <- run_single_zinegbin_brms(abs_dat = abs_h3k36_dat, 
                                                               this_gene = "YAL012W",
                                                               this_assay = "writer_eraser_loss",
                                                               seed = 125)
# Posterior data sets
pcl6_negbin_writer_loss_post_dat <- conditional_effects(pcl6_negbin_writer_loss_fit, plot = FALSE)
pcl6_negbin_writer_eraser_loss_post_dat <- conditional_effects(pcl6_negbin_writer_eraser_loss_fit, plot = FALSE)

figs6b <- ggplot(data = abs_h3k36_dat %>%
                       filter(gene == "YAL012W",
                              assay %in% c("writer_loss", "writer_eraser_loss")),
                     aes(x = timepoint_min, y = value, col = assay)) +
  geom_line(aes(group = sample_unit), linetype = "dashed") +
  scale_color_manual(values = c(wl_col, wel_col), 
                     labels = c(expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:"))), 
                                expression(paste(bolditalic("set2"), Delta, bolditalic("rph1"), Delta, bold(" + LANS-Set2:"))))) +
  geom_smooth(data = pcl6_negbin_writer_eraser_loss_post_dat$timepoint,
              aes(x = `effect1__` * 60, y = `estimate__`/1000, ymin = `lower__`/1000, ymax = `upper__`/1000), 
              stat = "identity", se = TRUE, col = wel_col, method = "loess", size = 2) + 
  geom_smooth(data = pcl6_negbin_writer_loss_post_dat$timepoint, 
              aes(x = `effect1__` * 60, y = `estimate__`/1000, ymin = `lower__`/1000, ymax = `upper__`/1000), 
              stat = "identity", se = TRUE, col = wl_col, method = "loess", size = 2) +
  xlim(0, 90) +
  ylab("Mean H3K36me3/H3\n(normalized ChIP Signal)") + xlab("Time (minutes)") +
  annotate(geom = "text", x = 0, y = 0.6, label = "PCL6") +
  plot_theme +
  theme(legend.position = c(0.8, 0.9)) +
  labs(col = "") +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

## LtD and LtD + RPH1delta for relative PCL6 (Supplemental Figure 6C)
## Run the Stan model fits
pcl6_beta_writer_loss_fit <- run_single_beta_brms(rel_dat = rel_h3k36_dat, 
                                                  this_gene = "YAL012W",
                                                  this_assay = "writer_loss",
                                                  seed = 124)
pcl6_beta_writer_eraser_loss_fit <- run_single_beta_brms(rel_dat = rel_h3k36_dat, 
                                                         this_gene = "YAL012W",
                                                         this_assay = "writer_eraser_loss",
                                                         seed = 124)
# Posterior data sets
pcl6_beta_writer_loss_post_dat <- conditional_effects(pcl6_beta_writer_loss_fit, plot = FALSE)
pcl6_beta_writer_eraser_loss_post_dat <- conditional_effects(pcl6_beta_writer_eraser_loss_fit, plot = FALSE)

figs6c <- ggplot(data = rel_h3k36_dat %>%
                        filter(gene == "YAL012W",
                               assay %in% c("writer_loss", "writer_eraser_loss")),
                      aes(x = timepoint_min, y = value, col = assay)) +
  geom_line(aes(group = sample_unit), linetype = "dashed") +
  scale_color_manual(values = c(wl_col, wel_col),
                     labels = c(expression(paste(bolditalic("set2"), Delta, bold(" + LANS-Set2:"))), 
                                expression(paste(bolditalic("set2"), Delta, bolditalic("rph1"), Delta, bold(" + LANS-Set2:"))))) +
  geom_smooth(data = pcl6_beta_writer_loss_post_dat$timepoint,
              aes(x = `effect1__` * 60, y = `estimate__`, ymin = `lower__`, ymax = `upper__`), 
              stat = "identity", se = TRUE, col = wl_col, method = "loess", size = 2) + 
  geom_smooth(data = pcl6_beta_writer_eraser_loss_post_dat$timepoint, 
              aes(x = `effect1__` * 60, y = `estimate__`, ymin = `lower__`, ymax = `upper__`), 
              stat = "identity", se = TRUE, col = wel_col, method = "loess", size = 2) +
  ylim(0, 1) + xlim(0, 90) +
  ylab("Proportional H3K36me3/H3\n(as a fraction of each gene max)") + xlab("Time (minutes)") +
  annotate(geom = "text", x = 85, y = 1, label = "PCL6") +
  plot_theme +
  theme(legend.position = c(0.15, 0.15)) +
  labs(col = "") +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

## ZIO beta GLMM rates by RNA abundance for dark to light + RPH1detla (Supplemental Figure S6D)
# Grab RNA abundance and average at time = 0
wel_rna_dat <- read.table("Supplemental_Data9.txt", header = TRUE) %>%
  dplyr::rename(gene = Gene) %>%
  gather(key = rep, value = transcript, -gene) %>%
  group_by(gene) %>%
  summarize(transcript = log(mean(transcript) + 1)) %>%
  ungroup

wel_glmm_vs_rna_dat <- beta_glmm_dat %>%
  filter(assay == "writer_eraser_loss") %>%
  dplyr::select(gene, estimate) %>%
  left_join(wel_rna_dat) %>%
  filter(gene %in% high_confidence_genes_down_both)
  
figs6d <- ggplot(data = wel_glmm_vs_rna_dat,
                aes(x = transcript, y = estimate)) +
  geom_point(col = wel_col) + 
  geom_smooth(aes(y = estimate, x = transcript), method = "lm", se = FALSE, size = 1, linetype = "dashed", col = "gray") + 
  xlim(0, 10) +
  scale_y_reverse(lim = c(-1.5, -3.5)) +
  xlab(expression(paste(italic("rph1"), italic(Delta), " RNA Abundance (log TPM)"))) +
  ylab("GLMM Estimate") +
  annotate(geom = "text", x = 9, y = -3.5, label = paste("r =", round(cor(wel_glmm_vs_rna_dat$estimate, wel_glmm_vs_rna_dat$transcript), 3))) +
  plot_theme

