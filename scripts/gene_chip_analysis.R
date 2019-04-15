library(lme4)
library(tidyverse)

setwd("~/projects/chipseq-gene-dynamics/")

## Read in the data
gene_chip_dat <- read.table("data/scaled_matrix_for_LMM_model_nonOverlappingGenes.txt") # path to data
names(gene_chip_dat) <- c("gene", "assay", "replicate", "timepoint", "value")

gene_chip_dat <- gene_chip_dat %>%
  mutate(assay = factor(assay),
         assay = recode(assay, "1" = "writer_loss", "2" = "writer_eraser_loss", "3" = "writer_add"),
         sample_unit = paste(gene, assay, sep = "_"),
         rep_id = paste(gene, assay, timepoint, sep = "_"))

###### Individual gene plots
practice <- gene_chip_dat %>%
  filter(gene %in% unique(gene_chip_dat$gene)[1:3])

## Linear
g <- ggplot(data = practice, aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + geom_line() + facet_wrap(~gene)
g <- g + geom_smooth(aes(y = value, x = timepoint), method = "lm")
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
g <- ggplot(data = practice, aes(x = timepoint, y = value, col = assay)) + scale_color_manual(values = c("magenta", "seagreen1", "coral")) + geom_point() + facet_wrap(~gene)
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

