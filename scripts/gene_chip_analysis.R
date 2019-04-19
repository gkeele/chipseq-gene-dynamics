library(lme4)
library(tidyverse)
library(brms)

setwd("~/projects/chipseq-gene-dynamics/")

## Read in the data
gene_chip_dat <- read.table("data/scaled_matrix_for_LMM_model_nonOverlappingGenes.txt") # path to data
names(gene_chip_dat) <- c("gene", "assay", "replicate", "timepoint", "value")

gene_chip_dat <- gene_chip_dat %>%
  mutate(assay = factor(assay),
         assay = recode(assay, "1" = "writer_loss", "2" = "writer_eraser_loss", "3" = "writer_add"),
         sample_unit = paste(gene, assay, replicate, sep = "_"))
gene_chip_dat$timepoint[gene_chip_dat$assay %in% c("writer_loss", "writer_eraser_loss") & gene_chip_dat$timepoint == 1] <- 30
gene_chip_dat$timepoint[gene_chip_dat$assay %in% c("writer_loss", "writer_eraser_loss") & gene_chip_dat$timepoint == 2] <- 60
gene_chip_dat$timepoint[gene_chip_dat$assay %in% c("writer_loss", "writer_eraser_loss") & gene_chip_dat$timepoint == 3] <- 90
gene_chip_dat$timepoint[gene_chip_dat$assay == "writer_add" & gene_chip_dat$timepoint == 1] <- 20
gene_chip_dat$timepoint[gene_chip_dat$assay == "writer_add" & gene_chip_dat$timepoint == 2] <- 40
gene_chip_dat$timepoint[gene_chip_dat$assay == "writer_add" & gene_chip_dat$timepoint == 3] <- 60



###### Individual gene plots
practice <- gene_chip_dat %>%
  filter(gene %in% unique(gene_chip_dat$gene)[1]) %>%
  group_by(sample_unit)

## Linear
g <- ggplot(data = practice, aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash")
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), 
                    axis.line = element_line(colour = "black"),
                    plot.title = element_text(hjust = 0.5), 
                    axis.text = element_text(size = 12, face = "bold"),
                    axis.title = element_text(size = 12, face = "bold"),
                    axis.text.x = element_text(hjust = 1, face = "bold"))
g

## LOESS (curves)
g <- ggplot(data = practice, aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "loess")
g <- g + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
               axis.line = element_line(colour = "black"),
               plot.title = element_text(hjust = 0.5), 
               axis.text = element_text(size = 12, face = "bold"),
               axis.title = element_text(size = 12, face = "bold"),
               axis.text.x = element_text(hjust = 1, face = "bold")) + guides(color = FALSE)
g



######## Example Stan analysis
## Try YAL012W and YAL066W
## Function to run brms
run_brms_on_chipseq <- function(chipseq_dat, this_gene) {
  use_dat <- chipseq_dat %>% filter(gene == this_gene,
                                    !(value %in% c(0, 1)))
  # writer loss
  wl_fit <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                      data = use_dat %>% filter(assay == "writer_loss"), 
                      family = "beta")
  wl_fixed <- summary(wl_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
  wl_random <- summary(wl_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
  # writer-eraser loss
  wel_fit <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                       data = use_dat %>% filter(assay == "writer_eraser_loss"), 
                       family = "beta")
  wel_fixed <- summary(wel_fit)$fixed %>% as.data.frame %>% rownames_to_column("parameter")
  wel_random <- summary(wel_fit)$random %>% as.data.frame %>% rownames_to_column("parameter")
  # writer add
  wa_fit <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), 
                       data = use_dat %>% filter(assay == "writer_add"), 
                       family = "beta")
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


## Parsing full set of Stan results
stan_fit_paths <- list.files("individual_gene_Stan_models/", full.names = TRUE)
for (i in 1:length(stan_fit_paths)) {
  this_fit <- readRDS(stan_fit_paths[i])
  this_gene <- gsub(x = stan_fit_paths[i], pattern = "individual_gene_Stan_models//sacCer3_nonOverlappingGenes_noChrMGenes_", replacement = "", fixed = TRUE) %>%
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
  mutate(estimate = estimate * 60,
         is_pos = ifelse(lower_95 > 0), TRUE, FALSE,
         is_neg = ifelse(upper_95 < 0), TRUE, FALSE)
timepoint_dat$category[timepoint_dat$lower_95 > 0] <- "positive"
timepoint_dat$category[timepoint_dat$upper_95 < 0] <- "negative"
timepoint_dat$category[timepoint_dat$lower_95 < 0 & timepoint_dat$upper_95 > 0] <- "zero"

  
ggplot(data = timepoint_dat, aes(x = estimate)) + geom_histogram() + facet_grid(assay~.)
ggplot(data = timepoint_dat %>% filter(category %in% c("positive", "negative"),
                                       estimate > -10), 
       aes(x = estimate)) + geom_histogram() + facet_grid(assay~.)

## Grabbing extreme negative genes
low_genes <- timepoint_dat %>% 
  filter(estimate < -10,
         category == "negative") %>%
  pull(gene) %>%
  unique

## Grabbing extreme negative genes
zero_genes <- timepoint_dat %>% 
  filter(category == "zero") %>%
  pull(gene) %>%
  unique

g <- ggplot(data = gene_chip_dat %>% filter(gene %in% low_genes[1:25]), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
               axis.line = element_line(colour = "black"),
               plot.title = element_text(hjust = 0.5), 
               axis.text = element_text(size = 12, face = "bold"),
               axis.title = element_text(size = 12, face = "bold"),
               axis.text.x = element_text(hjust = 1, face = "bold")) + guides(color = FALSE)
g
g <- ggplot(data = gene_chip_dat %>% filter(gene %in% zero_genes[1:25]), aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm", size = 2)
g <- g + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
               axis.line = element_line(colour = "black"),
               plot.title = element_text(hjust = 0.5), 
               axis.text = element_text(size = 12, face = "bold"),
               axis.title = element_text(size = 12, face = "bold"),
               axis.text.x = element_text(hjust = 1, face = "bold")) + guides(color = FALSE)
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

library(UpSetR)
upset_dat <- data.frame(genes = genes_all,
                        WA_positive = sapply(1:length(genes_all), function(i) as.numeric(genes_all[i] %in% genes_pos_wa)),
                        WL_negative = sapply(1:length(genes_all), function(i) as.numeric(genes_all[i] %in% genes_neg_wl)),
                        WEL_negative = sapply(1:length(genes_all), function(i) as.numeric(genes_all[i] %in% genes_neg_wel)))
upset(upset_dat, order.by = "freq", text.scale = c(1.5, 1.5, 1.5, 1.2, 1.2), sets.bar.color = c("coral", "magenta", "seagreen1"))

## Estimate by error
plot(timepoint_dat$est_error, timepoint_dat$estimate, pch = ifelse(timepoint_dat$category == "zero", 1, 19), col = c("coral", "magenta", "seagreen1")[as.factor(timepoint_dat$assay)], 
     las = 1, ylab = "Trend with time", xlab = "Error on trend")
plot(timepoint_dat$est_error, timepoint_dat$estimate, ylim = c(-7, 7), pch = ifelse(timepoint_dat$category == "zero", 1, 19), col = c("coral", "magenta", "seagreen1")[as.factor(timepoint_dat$assay)], 
     las = 1, ylab = "Trend with time", xlab = "Error on trend")

time_alt_fit <- lm(estimate ~ 1 + assay, data = timepoint_dat, weights = 1/timepoint_dat$est_error)
time_null_fit <- lm(estimate ~ 1, data = timepoint_dat, weights = 1/timepoint_dat$est_error)
anova(time_alt_fit, time_null_fit)


## Assay specific plots
mean_gene_chip_dat <- gene_chip_dat %>%
  group_by(gene, assay, timepoint, batch_id, rep_id) %>%
  summarise(value = mean(value)) %>%
  ungroup
# Writer loss
writer_loss <- ggplot(data = mean_gene_chip_dat %>% filter(assay == "writer_loss"), aes(x = timepoint, y = value, color = batch_id)) + geom_point() + geom_line()
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
writer_eraser_loss <- ggplot(data = mean_gene_chip_dat %>% filter(assay == "writer_eraser_loss"), aes(x = timepoint, y = value, color = batch_id)) + geom_point() + geom_line()
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
writer_add <- ggplot(data = mean_gene_chip_dat %>% filter(assay == "writer_add"), aes(x = timepoint, y = value, color = batch_id)) + geom_point() + geom_line()
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


