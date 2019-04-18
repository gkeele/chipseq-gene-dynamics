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
  filter(gene %in% unique(gene_chip_dat$gene)[1:16]) %>%
  group_by(sample_unit)

## Linear
g <- ggplot(data = practice, aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line(aes(group = sample_unit), linetype = "longdash") + facet_wrap(~gene)
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



######## Stan analysis
## Try YAL012W and YAL066W

# YAL012W - clean signal
yal012w_wl <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), data = gene_chip_dat %>% filter(gene == "YAL012W",
                                                                                                               assay == "writer_loss", 
                                                                                                               value != 0, 
                                                                                                               value != 1), family = "beta")
yal012w_wel <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), data = gene_chip_dat %>% filter(gene == "YAL012W",
                                                                                                          assay == "writer_eraser_loss", 
                                                                                                          value != 0, 
                                                                                                          value != 1), family = "beta")
yal012w_wa <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), data = gene_chip_dat %>% filter(gene == "YAL012W",
                                                                                                          assay == "writer_add", 
                                                                                                          value != 0, 
                                                                                                          value != 1), family = "beta")

# YAL066W - messy signal
yal066w_wl <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), data = gene_chip_dat %>% filter(gene == "YAL066W",
                                                                                                               assay == "writer_loss", 
                                                                                                               value != 0, 
                                                                                                               value != 1), family = "beta")
yal066w_wel <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), data = gene_chip_dat %>% filter(gene == "YAL066W",
                                                                                                                assay == "writer_eraser_loss", 
                                                                                                                value != 0, 
                                                                                                                value != 1), family = "beta")
yal066w_wa <- brms::brm(value ~ 1 + timepoint + (1 + timepoint | sample_unit), data = gene_chip_dat %>% filter(gene == "YAL066W",
                                                                                                               assay == "writer_add", 
                                                                                                               value != 0, 
                                                                                                               value != 1), family = "beta")

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


## Full model
gene_chip_lmm <- lmer(value ~ (1 + timepoint + assay + timepoint:assay | sample_unit), REML = FALSE, data = gene_chip_dat)

gene_chip_lmm <- lmer(value ~ 1 + timepoint + assay + timepoint:assay + (1 | rep_id) + (assay + timepoint + assay:timepoint - 1 | gene), data = gene_chip_dat, REML = FALSE)
#gene_chip_lmm <- lmer(value ~ 1 + timepoint + assay + (1 | rep_id) + (timepoint | rep_id), data = gene_chip_dat) # Fails to converge

## Test for assay
gene_chip_lmm_assay <- lmer(value ~ 1 + timepoint + assay + (1 | rep_id) + (timepoint - 1 | batch_id), data = gene_chip_dat, REML = FALSE)
gene_chip_lmm_noassay <- lmer(value ~ 1 + timepoint + (1 | rep_id) + (timepoint - 1 | batch_id), data = gene_chip_dat, REML = FALSE)
anova(gene_chip_lmm_assay, gene_chip_lmm_noassay)

## Test for timepoint
gene_chip_lmm_assay <- lmer(value ~ 1 + timepoint + assay + (1 | rep_id) + (timepoint - 1 | batch_id), data = gene_chip_dat, REML = FALSE)
gene_chip_lmm_noassay <- lmer(value ~ 1 + timepoint + (1 | rep_id) + (timepoint - 1 | batch_id), data = gene_chip_dat, REML = FALSE)
anova(gene_chip_lmm_assay, gene_chip_lmm_noassay)

## Test for interaction of assay with time
gene_chip_lmm_full <- lmer(value ~ 1 + timepoint + assay + timepoint:assay + (1 | rep_id), data = gene_chip_dat, REML = FALSE)
gene_chip_lmm_noint <- lmer(value ~ 1 + timepoint + assay + (1 | rep_id), data = gene_chip_dat, REML = FALSE)
anova(gene_chip_lmm_full, gene_chip_lmm_noint)


## Mean model
mean_gene_chip_lmm1 <- lmer(value ~ 1 + assay*timepoint + (1 + assay*timepoint | gene), data = mean_gene_chip_dat, REML = FALSE)
mean_gene_chip_lmm2 <- lmer(value ~ 1 + assay*timepoint + (1 + assay + timepoint | gene), data = mean_gene_chip_dat, REML = FALSE)



lmm1 <- lmer(value ~ 1 + timepoint + assay + (1 + timepoint + assay | gene), data = mean_gene_chip_dat, REML = 0)
lmm2 <- lmer(value ~ 1 + timepoint + assay + (1| gene) + (timepoint - 1 | gene) + (assay - 1 | gene), data = mean_gene_chip_dat, REML = 0)
lmm3 <- lmer(value ~ 1 + timepoint + (1| gene) + (timepoint - 1 | gene), data = mean_gene_chip_dat, REML = 0)
anova(lmm1, lmm2)
anova(lmm3, lmm2)


blup <- ranef(mean_gene_chip_lmm)

