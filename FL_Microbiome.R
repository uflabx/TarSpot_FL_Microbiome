###Load packages####
pacman::p_load(dada2, ShortRead, Biostrings, phyloseq, tidyverse, igraph, MicEco,
               DECIPHER, phangorn, vegan, edgeR, ggpubr, SpiecEasi, grid, ggthemes,
               pegas, ape, lme4, lmerTest, glmmTMB)


###Set up themes####
theme_Publication <- function(base_size=12, base_family="Arial") {
  (theme_foundation(base_size=base_size, base_family=base_family) + 
      theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.border = element_blank(),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_blank(),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0.5, "cm"),
            legend.justification = "center",
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

###DADA2####
path <- "Reads_for_analysis"

list.files(path)

fnFs <- sort((list.files(path=path,pattern="_1.fastq.gz", recursive = TRUE,full.names = TRUE, include.dirs=TRUE)))
fnRs <- sort((list.files(path=path,pattern="_2.fastq.gz", recursive = TRUE,full.names = TRUE, include.dirs=TRUE)))

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

sample.names <- sapply(substr(basename(fnFs), 1,5), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,minLen = 50,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers) #analagous to an OTU table
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.csv(track, "reads_track.csv")

###Phyloseq####
####Make input files for phyloseq function####
#####1. tax_table####
unite.ref <- "sh_general_release_s_04.04.2024/sh_general_release_dynamic_s_04.04.2024.fasta"  
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE, verbose = TRUE)

#####2. sample_data ####
samples.out <- rownames(seqtab.nochim)
Location <- substr(samples.out, 1,3)
rep <- substr(samples.out, 5,5)
samdf <- data.frame(samples=samples.out,Location=Location, rep=rep)
rownames(samdf) <- samples.out
samdf$Tissue <- ifelse(rep == "C", "Asymptomatic", "Symptomatic")
samdf$Year <-ifelse(samdf$Location == "TSA" | samdf$Location == "TSC", "2022", "2023")
samdf$Group <- paste0(samdf$Year, " ", samdf$Tissue)

#####3. refseq####
ASVs.nochim = DNAStringSet(colnames(seqtab.nochim))
names(ASVs.nochim) = paste0("ASV", 1:ncol(seqtab.nochim))

head(names(ASVs.nochim))
tail(names(ASVs.nochim))

#####4. phy_tree####
alignment = AlignSeqs(ASVs.nochim, anchor=NA, processors=30)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)

#####Format taxa names####
tmp.seqtab = seqtab.nochim
colnames(tmp.seqtab) = names(ASVs.nochim)
tmp.taxa = taxa
rownames(tmp.taxa) = names(ASVs.nochim)
rownames(taxa) = names(ASVs.nochim)

taxa <- gsub("k__", "", taxa)
taxa <- gsub("p__", "", taxa)
taxa <- gsub("c__", "", taxa)
taxa <- gsub("f__", "", taxa)
taxa <- gsub("o__", "", taxa)
taxa <- gsub("g__", "", taxa)
taxa <- gsub("s__", "", taxa)

res <- grep("*_Incertae_sedis", taxa)
GS11 <- grep("GS11", taxa)
taxa <- replace(taxa,res, "Uncharacterized")
taxa <- replace(taxa,GS11, "Uncharacterized")
taxa[is.na(taxa)] <- "Uncharacterized"

###Create phyloseq object####
ps <- phyloseq(otu_table(tmp.seqtab, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa),
               refseq(ASVs.nochim),
               phy_tree(fitGTR$tree))
ps

###Pre-processing and filtering####
#1. Removes singletons and NAs
ps <- prune_taxa(taxa_sums(ps) > 1, ps) 
ps = subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "Uncharacterized"))
together <- psmelt(ps) 
together1 <- together %>% group_by(Sample, Abundance, OTU)%>% filter(OTU != "NA1") 

#2. Remove OTUs with less than 0.1% abundance in at least one sample.
together2 <- together1 %>% 
  group_by(Sample) %>%
  mutate(total_reads=sum(Abundance), 
         percent=Abundance/total_reads*100) %>% 
  filter(percent >=0.1)

keep <- together2$OTU
ps1 = subset_taxa(ps, colnames(otu_table(physeq)) %in% keep)

saveRDS(ps1, "ps1.rds")

###Normalization####
#Normalization with edgeR
edgeRnorm = function(physeq, ...) {
  require("edgeR")
  require("phyloseq")
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  y = edgeR::DGEList(counts = x, remove.zeros = TRUE)
  z = edgeR::calcNormFactors(y, ...)
  if (!all(is.finite(z$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
  }
  return(z)
}
physeq <- ps1

# Extract components from the original phyloseq object
otu_table_original <- otu_table(physeq)
tax_table_original <- tax_table(physeq)
sample_data_original <- sample_data(physeq)
phy_tree_original <- phy_tree(physeq)

write.csv(tax_table_original, "taxtable.csv", row.names = T, quote = F)

ps_norm_counts <- equalizeLibSizes(edgeRnorm(physeq, method = 'RLE'))$pseudo.counts
otu_table_normalized <- otu_table(ps_norm_counts, taxa_are_rows = TRUE)
ps_norm <- phyloseq(otu_table_normalized, tax_table_original, sample_data_original, phy_tree_original)

###Alpha diversity####
otu_matrix <- as.matrix(otu_table(ps_norm))  

shannon_index <- apply(otu_matrix, 2, function(x) diversity(x, index = "shannon"))
simpson_index <- apply(otu_matrix, 2, function(x) diversity(x, index = "simpson"))
observed_richness <- apply(otu_matrix, 2, function(x) sum(x > 0))

richness_diversity_manual <- data.frame(
  SampleID = colnames(otu_matrix),
  Observed = observed_richness,
  Shannon = shannon_index,
  Simpson = simpson_index)

sample_data_df <- data.frame(sample_data(ps_norm))
sample_data_df$SampleID <- rownames(sample_data_df)

richness_diversity_manual <- richness_diversity_manual %>%
  left_join(sample_data_df, by = "SampleID")

write.csv(richness_diversity_manual, "richness_diversity.csv", quote = F, row.names = F)

####Shannon####
model_shannon_lm <- lmerTest::lmer(Shannon ~ Tissue + (1 | Location), data = richness_diversity_manual)
model_shannon_glmm <- glmmTMB(Shannon ~ Tissue + (1 | Location), data = richness_diversity_manual)
AIC(model_shannon_lm, model_shannon_glmm)
AIC(model_shannon_glmm,
    update(model_shannon_glmm, family = Gamma()), 
    update(model_shannon_glmm, family = nbinom2()),
    update(model_shannon_glmm, family = gaussian()))

model_shannon = glmmTMB(Shannon ~ Tissue + (1 | Location),
                         family = gaussian(link = 'identity'),
                         data = richness_diversity_manual)
res_shannon = DHARMa::simulateResiduals(model_shannon) |> 
  resid(quantileFunction = qnorm, outlierValues = c(-5, 5))
fit_shannon <- fitted(model_shannon)
ggResidpanel::resid_auxpanel(res_shannon, fit_shannon)

emmeans::joint_tests(model_shannon)

(emm_shannon = emmeans::emmeans(model_shannon, ~ Tissue))

p_shannon = emm_shannon |> 
  as.data.frame() |> 
  ggplot(aes(x = Tissue, fill=Tissue))+
  geom_col(aes(y = emmean))+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = .2) +
  theme_bw() +
  labs(x = '', y = 'Shannon') +
  scale_fill_manual(values = c("#54278fff","#d95f02")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

####Richness####
model_observed <- glmmTMB(Observed ~ Tissue + (1 | Location), 
                          family = poisson(link = 'log'), 
                          data = richness_diversity_manual)
res_observed = DHARMa::simulateResiduals(model_observed) |> 
  resid(quantileFunction = qnorm, outlierValues = c(-5, 5))
fit_observed <- fitted(model_observed)
ggResidpanel::resid_auxpanel(res_observed, fit_observed)

AIC(model_observed,
    update(model_observed, family = nbinom1()), 
    update(model_observed, family = nbinom2()))

model_observed = glmmTMB(Observed ~ Tissue + (1 | Location),
                         family = nbinom2(link = 'log'),
                         data = richness_diversity_manual)
res_observed = DHARMa::simulateResiduals(model_observed) |> 
  resid(quantileFunction = qnorm, outlierValues = c(-5, 5))
fit_observed <- fitted(model_observed)
ggResidpanel::resid_auxpanel(res_observed, fit_observed)

emmeans::joint_tests(model_observed)

(emm_observed = emmeans::emmeans(model_observed, ~ Tissue, type = 'response'))

p_observed = emm_observed |> 
  as.data.frame() |> 
  ggplot(aes(x = Tissue, fill = Tissue))+
  geom_col(aes(y = response))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .2) +
  theme_bw() +
  labs(x = '', y = 'Richness') +
  scale_fill_manual(values = c("#54278fff","#d95f02")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

###Beta diversity####
####Bray-curtis####
logt  = transform_sample_counts(ps_norm, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "MDS", distance = "bray")

p_pcoa <- plot_ordination(logt, out.pcoa.logt,
                          color = "Tissue", shape = "Tissue") +
  geom_point(size = 2) +
  stat_ellipse( type = "t", linetype = 2, color = "black", alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("#54278fff","#d95f02"))

ps_bray <- phyloseq::distance(logt, method = "bray")
ado_t_bray <- adonis2(ps_bray~Tissue, data=samdf)
ado_l_bray <- adonis2(ps_bray~Location, data=samdf)
ado_tl_bray <- adonis2(ps_bray~Tissue/Location, data=samdf)

beta_disp_bray <- betadisper(ps_bray, sample_data(logt)$Tissue)
anova_disp_bray <- anova(beta_disp_bray)

####Uunifrac####
out.pcoa.prop <- ordinate(ps_norm, method = "MDS", distance = "unifrac", weighted=F)
p_pcoa_prop <- plot_ordination(ps_norm, out.pcoa.prop,
                               color = "Tissue", shape = "Tissue") +
  geom_point(size = 2) + 
  stat_ellipse(type = "t", linetype = 2, color = "black", alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("#54278fff","#d95f02"))

ps_uunifrac <- phyloseq::distance(ps_norm, method = "uunifrac")
ado_t_uu <- adonis2(ps_uunifrac~Tissue, data=samdf)
ado_l_uu <- adonis2(ps_uunifrac~Location, data=samdf)
ado_tl_uu <- adonis2(ps_uunifrac~Tissue/Location, data=samdf)


beta_disp <- betadisper(ps_uunifrac, sample_data(ps_norm)$Tissue)
anova_disp <- anova(beta_disp)

diversity_plots <-ggarrange(p_observed,p_shannon, p_pcoa, p_pcoa_prop, nrow = 2, ncol = 2,
          common.legend = T, legend = "right", align = "hv", 
          labels = c("A", "", "B", ""))

ggsave("diversity_plots.svg", width= 8, height = 5)

###Relative abundance####
####Phylum#####
pie_phylum <- together2 %>% 
  group_by(Tissue, Phylum) %>%
  summarize(total_reads=sum(Abundance))

write_csv(pie_phylum, "phylum_table.csv")

####Class#####
pie_class_df <- psmelt(ps1) %>% 
  group_by(Tissue) %>%
  mutate(total_reads=sum(Abundance), 
         percent=Abundance/total_reads*100)

pie_class_df$Section <- ifelse(pie_class_df$Class == "Sordariomycetes" | pie_class_df$Class == "Dothideomycetes", paste0(pie_class_df$Class), "b")

#####Symptomatic#####
pie_class.ts <- pie_class_df %>% 
  filter(Tissue == "Symptomatic") %>%
  mutate(arc_start = cumsum(lag(Abundance, default = 0)) * 2*pi - 2,
         arc_end   = cumsum(Abundance) * 2*pi - 2,
         x_pos = 0 + cos(arc_start - pi/2),
         y_pos = 1 - sin(arc_start - pi/2))

pie_class.ts2 <- pie_class.ts %>%
  group_by(Section) %>%
  summarize(Percent = sum(percent))

pie_class2 <- pie_class.ts %>% 
  ggplot(aes(x="", y=Abundance))+ 
  geom_bar(aes(fill=Section, color=Section), 
           stat="identity", 
           position="fill", 
           width=0.9) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("#f1a340", "#5e3c99", "#998ec3")) +
  scale_color_manual(values = c("#f1a340", "#5e3c99", "#998ec3")) +
  theme_void()

section_df <- pie_class_df %>%
  filter(Tissue == "Symptomatic" & Section == "b") %>% 
  group_by(Class) %>%
  summarize(Percent = sum(percent))

section_df %>%  
  ggplot(aes(x="", y=Percent, fill=Class, label= Percent))+ 
  geom_col(position="fill", 
           width=0.9) +
  scale_fill_brewer(palette = "BrBG") +
  theme_void() + 
  geom_text(size=3, position = position_stack(vjust = 0.5))

section_plot <- section_df %>%  
  filter(Percent >0) %>%
  ggplot(aes(x="", y=Percent, fill=Class, label= round(Percent, digits = 2)))+ 
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "BrBG") +
  theme_void() +
  geom_text(position = position_stack(vjust = 0.5), color = "black")

pie_legend <- get_legend(pie_class2) 
section_leg <- get_legend(section_plot)
legends <- ggarrange(pie_legend, section_leg, nrow=2)

tar <- ggarrange(pie_class2, section_plot, nrow=1, legend = "none",
          widths = c(2,0.4,1), heights = c(2,0.01,1))

#####Asymptomatic#####
pie_class.c <- pie_class_df %>% 
  filter(Tissue == "Asymptomatic") %>%
  mutate(arc_start = cumsum(lag(Abundance, default = 0)) * 2*pi - 2,
         arc_end   = cumsum(Abundance) * 2*pi - 2,
         x_pos = 0 + cos(arc_start - pi/2),
         y_pos = 1 - sin(arc_start - pi/2))

pie_class.c2 <- pie_class.c %>%
  group_by(Section) %>%
  summarize(Percent = sum(percent))

pie_class2.c <- pie_class.c %>% 
  ggplot(aes(x="", y=Abundance))+ 
  geom_bar(aes(fill=Section, color=Section), 
           stat="identity", 
           position="fill", 
           width=0.9) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("#f1a340", "#5e3c99", "#998ec3")) +
  scale_color_manual(values = c("#f1a340", "#5e3c99", "#998ec3")) +
  theme_void()

section_df.c <- pie_class_df %>%
  filter(Tissue == "Asymptomatic" & Section == "b") %>% 
  group_by(Class) %>%
  summarize(Percent = sum(percent))

section_plot.c <- section_df.c %>%  
  filter(Percent >0) %>%
  ggplot(aes(x="", y=Percent, fill=Class, label= round(Percent, digits = 2)))+ 
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "BrBG") +
  theme_void()+ 
  geom_text(position = position_stack(vjust = 0.5), color = "black")

section_leg.c <- get_legend(section_plot.c)

control <- ggarrange(pie_class2.c, section_plot.c, nrow=1, legend = "none",
          widths = c(2,0.4,1), heights = c(2,0.01,1))

class_plot <- ggarrange(control, tar)

ggarrange(diversity_plots, class_plot, nrow = 2, labels = c("A", "B"))


###Venn diagram####
tiff("venn.tiff",  height=1200, width = 1200, units = "px", 
     compression = "none", pointsize = 12, res = 300)
ps_venn(ps1, group="Tissue", plot = T, weight = F, relative = T) 
dev.off() 

###Phylogenetic tree analysis####
#order by relative abundance
abundance <- together2 %>%
  group_by(Sample) %>%
  arrange(-percent) 

#subet top 50 ASVs
top50 <- unique(abundance$OTU)[1:50]
ps50 = subset_taxa(ps1, colnames(otu_table(ps1)) %in% top50)
tax_table(ps50) <- cbind(tax_table(ps50), Strain=taxa_names(ps50))
myranks = c("Strain", "Genus", "Species")
mylabels = apply(tax_table(ps50)[, myranks], 1, paste, sep="", collapse=" ")
mylabels <- sub("Uncharacterized Uncharacterized", "uncharacterized", mylabels)
mylabels <- sub(" Uncharacterized", " sp.", mylabels)
tax_table(ps50) <- cbind(tax_table(ps50), catglab=mylabels)

plot_tree(ps50, label.tips = "catglab", 
          ladderize="left", 
          color = "Family", 
          shape="Tissue",
          size="abundance",
          plot.margin=0.4,
          base.spacing=0.04)  +
  theme(legend.key=element_blank())

ggsave("fig3.svg",device = "svg", height = 10, width = 10,  dpi=300)

###Proportion data analysis####
tax_dict = ps1@tax_table |> 
  as.data.frame() |> 
  rownames_to_column(var = 'OTU') |> 
  dplyr::select(OTU, Genus) |> 
  mutate(Genus = gsub('Uncharacterized', 'Other', Genus))

ps_prop = apply(ps_norm_counts, 2, \(x){x/sum(x)})
low = apply(ps_prop, 1, \(x){all(x < 0.01)})
other = colSums(ps_prop[low, ])
ps_sub = rbind(ps_prop[!low, ], Other = other)

ps_df = ps_sub |> 
  as.data.frame() |> 
  rownames_to_column(var = 'OTU') |> 
  pivot_longer(cols = 2:16, names_to = 'Sample', values_to = 'Proportion') |> 
  separate('Sample', into = c('Location', 'SampleType'),
           sep = '_', remove = F) |> 
  mutate(SampleType = ifelse(SampleType == 'C', 'Asymptomatic', 'Symptomatic')) |> 
  mutate(across(2:4, factor)) |> 
  left_join(tax_dict, by = join_by(OTU)) |> 
  mutate(Genus = ifelse(is.na(Genus), 'Other', Genus)) |> 
  group_by(Sample, Location, SampleType, Genus) |> 
  summarize(Proportion = sum(Proportion), 
            .groups = 'drop')

ps_wide = ps_df |> 
  mutate(Genus = ifelse(Genus %in% c('Alternaria', 'Epicoccum', 'Exserohilum', 'Phyllachora'), Genus, 'Other'), 
         Genus = factor(Genus, levels = c('Alternaria', 'Epicoccum', 'Exserohilum', 'Phyllachora', 'Other'))) |>  
  group_by(Sample, Location, SampleType, Genus) |> 
  summarise(Proportion = sum(Proportion), .groups = 'drop') |> 
  pivot_wider(names_from = 'Genus', values_from = 'Proportion')

sum(apply(ps_wide[, 4:ncol(ps_wide)], 2, \(x){mean(x == 0)}) == 0) 

# Preparing data for Stan
X = model.matrix(~ 0+SampleType, data = ps_wide)
Z = model.matrix(~ 0+Location, data = ps_wide)
N = nrow(X)
P = ncol(X)
Q = ncol(Z)
Y = ps_wide[, 4:8]
K = ncol(Y)

stan_data = list(N = N, P = P, Q = Q, X = X, Z = Z, K = K, Y = Y)

nchain = 4
nwarm = 1e3
niter = 1e3

stan_mod = stan(file = 'Dirichlet GLMM.stan', data = stan_data, 
                 chains = nchain, cores = nchain, iter = niter+nwarm, warmup = nwarm)

print(stan_mod, pars = c('beta', 'theta'))
samps = rstan::extract(stan_mod)
ctrl_samps = cbind(samps$beta[, 1, ], 1)
ctrl_samps = t(apply(ctrl_samps, 1, \(x){exp(x)/sum(exp(x))}))
tar_samps = cbind(samps$beta[, 2, ], 1)
tar_samps = t(apply(tar_samps, 1, \(x){exp(x)/sum(exp(x))}))

ord = order(rowSums(ctrl_samps[, 2:4]))
ctrl_samps = ctrl_samps[ord, ]
ctrl_df = ctrl_samps |> 
  as.data.frame() |> 
  magrittr::set_colnames(colnames(Y)) |> 
  rownames_to_column(var = 'Iter') |> 
  mutate(Iter = as.numeric(Iter)) |> 
  pivot_longer(cols = 2:6, names_to = 'Genus', values_to = 'Proportion') |> 
  mutate(Genus = factor(Genus, levels = c('Phyllachora', 'Alternaria', 'Epicoccum', 'Exserohilum', 'Other')))

ord = order(tar_samps[, 4])
tar_samps = tar_samps[ord, ]
tar_df = tar_samps |> 
  as.data.frame() |> 
  magrittr::set_colnames(colnames(Y)) |> 
  rownames_to_column(var = 'Iter') |> 
  mutate(Iter = as.numeric(Iter)) |> 
  pivot_longer(cols = 2:6, names_to = 'Genus', values_to = 'Proportion') |> 
  mutate(Genus = factor(Genus, levels = c('Phyllachora', 'Alternaria', 'Epicoccum', 'Exserohilum', 'Other')))

for (s in 1:ncol(tar_samps)){
  d = tar_samps[, s] - ctrl_samps[, s]
  est = round(quantile(d, 0.5), 3)
  ci = round(quantile(d, c(0.025, 0.975)), 3)
  cat(paste0(colnames(Y)[s], ' Tarspot Sample Mean % : ', round(mean(tar_samps[, s]), 3), '\nControl Sample Mean %: ', round(mean(ctrl_samps[, s]), 3)), '\n')
  cat(paste0(colnames(Y)[s], ' Difference: ', est, ' (', ci[1], ' - ', ci[2], ')\n'))
}

###Network analysis####
####Node - abundance and presence####
samples.ts <- subset_samples(ps1, Tissue == "Symptomatic")
samples.c <- subset_samples(ps1, Tissue == "Asymptomatic")

otu_tar_spot <- as.data.frame(otu_table(samples.ts))
otu_control <- as.data.frame(otu_table(samples.c))
otu_full <- as.data.frame(otu_table(ps1))

asv_unique.ts <- colnames(otu_tar_spot)[colSums(otu_tar_spot > 0) > 0 & colSums(otu_control > 0) == 0]
asv_unique.c <- colnames(otu_control)[colSums(otu_control > 0) > 0 & colSums(otu_tar_spot > 0) == 0]
asv_shared <- colnames(otu_tar_spot)[colSums(otu_tar_spot > 0) > 0 & colSums(otu_control > 0) > 0]

rel_abundance.ts <- colSums(otu_tar_spot) / sum(otu_tar_spot)
rel_abundance.c <- colSums(otu_control) / sum(otu_control)

asv_combined <- unique(c(colnames(otu_tar_spot), colnames(otu_control)))

node_df <- data.frame(
  ASV_name = asv_combined,
  Presence = ifelse(asv_combined %in% asv_unique.ts, "Symptomatic",
                    ifelse(asv_combined %in% asv_unique.c, "Asymptomatic", "Shared")),
  rel_abundance_tar_spot = ifelse(asv_combined %in% colnames(otu_tar_spot), 
                                  rel_abundance.ts[asv_combined], 0),
  rel_abundance_asymp = ifelse(asv_combined %in% colnames(otu_control), 
                               rel_abundance.c[asv_combined], 0)
)

write.csv(node_df, "nodes.csv", row.names = F, quote = F)

###Network construction####
otu_table_df <- as.data.frame(otu_table_normalized)

# Subset OTU table by tissue type
otu_c <- otu_table_df[, samdf$Tissue == "Asymptomatic"]
otu_tarspot <- otu_table_df[, samdf$Tissue == "Symptomatic"]

#Filter ASVs unique to each tissue to avoid misleading correlations with 0 abundance ASVs
otu_tarspot_filtered <- otu_tarspot[!(rownames(otu_tarspot) %in% asv_unique.c), ]
otu_control_filtered <- otu_c[!(rownames(otu_c) %in% asv_unique.ts), ]

####Sparcc - tar spot####
sparcc_result.ts <- sparcc(t(otu_tarspot_filtered))

sparcc_cor.ts <- sparcc_result.ts$Cor

ASV_names <- rownames(otu_tarspot_filtered)
rownames(sparcc_cor.ts) <- ASV_names
colnames(sparcc_cor.ts) <- ASV_names

# Create a network from the SparCC correlation matrix
g_sparcc.ts <- graph_from_adjacency_matrix(sparcc_cor.ts, mode = "undirected", weighted = TRUE, diag = FALSE)

####P-values - tar spot####
run_sparcc_permuted <- function(data, n_permutations) {
  perm_cor_list <- list()
  
  for (i in 1:n_permutations) {
    permuted_data <- apply(data, 2, sample) 
    permuted_cor <- sparcc(permuted_data)$Cor
    perm_cor_list[[i]] <- permuted_cor
  }
  return(perm_cor_list)
}

n_permutations <- 1000
perm_correlations.ts <- run_sparcc_permuted(t(otu_tarspot_filtered), n_permutations)
perm_cor_array.ts <- simplify2array(perm_correlations.ts)

p_values.ts <- matrix(NA, nrow = nrow(sparcc_cor.ts), ncol = ncol(sparcc_cor.ts))

for (i in 1:nrow(sparcc_cor.ts)) {
  for (j in 1:ncol(sparcc_cor.ts)) {
    observed_cor <- sparcc_cor.ts[i, j]
    perm_cor <- perm_cor_array.ts[i, j, ]
    # Two-tailed p-value (percentage of permuted correlations greater than or equal to observed)
    p_values.ts[i, j] <- mean(abs(perm_cor) >= abs(observed_cor))
  }
}

p_values_vector.ts <- as.vector(p_values.ts)
p_values_adjusted_vector.ts <- p.adjust(p_values_vector.ts, method = "fdr") 
p_values_adjusted.ts <- matrix(p_values_adjusted_vector.ts, nrow = nrow(p_values.ts), ncol = ncol(p_values.ts))

rownames(p_values_adjusted.ts) <- rownames(sparcc_cor.ts)
colnames(p_values_adjusted.ts) <- colnames(sparcc_cor.ts)

edge_indices.ts <- as.data.frame(as_edgelist(g_sparcc.ts))
row_indices.ts <- edge_indices.ts[, 1]
col_indices.ts <- edge_indices.ts[, 2]

edge_p_values.ts <- mapply(function(row, col) p_values_adjusted.ts[row, col], row_indices.ts, col_indices.ts)

edge_data.ts <- data.frame(
  source = ends(g_sparcc.ts, E(g_sparcc.ts))[, 1],
  target = ends(g_sparcc.ts, E(g_sparcc.ts))[, 2],
  weight = E(g_sparcc.ts)$weight,
  p_value = edge_p_values.ts
)

write.csv(edge_data.ts, "edge_data.ts.csv", quote = F, row.names = F)

edge_data_sig.ts <- subset(edge_data.ts, abs(weight) > 0.6 & p_value < 0.05)

####Taxa source nodes####
tax_data <- as.data.frame(tax_table(ps1))

concatenate_taxa <- function(asv_name) {
  if (asv_name %in% rownames(tax_data)) {
    return(paste(tax_data[asv_name, ], collapse = "; "))
  } else {
    return(NA)
  }
}

edge_data_sig.ts$source_taxa <- sapply(edge_data_sig.ts$source, concatenate_taxa)

write.csv(edge_data_sig.ts, file = "TarSpot_edge_data.csv", row.names = FALSE, quote = FALSE)

###Sparcc - asymptomatic####
sparcc_result.c <- sparcc(t(otu_control_filtered))
sparcc_cor.c <- sparcc_result.c$Cor
ASV_names.c <- rownames(otu_control_filtered)
rownames(sparcc_cor.c) <- ASV_names.c
colnames(sparcc_cor.c) <- ASV_names.c

# Create a network from the SparCC correlation matrix
g_sparcc.c <- graph_from_adjacency_matrix(sparcc_cor.c, mode = "undirected", weighted = TRUE, diag = FALSE)

####P-values - asymptomatic####
perm_correlations.c <- run_sparcc_permuted(t(otu_control_filtered), n_permutations)
perm_cor_array.c <- simplify2array(perm_correlations.c)

p_values.c <- matrix(NA, nrow = nrow(sparcc_cor.c), ncol = ncol(sparcc_cor.c))

for (i in 1:nrow(sparcc_cor.c)) {
  for (j in 1:ncol(sparcc_cor.c)) {
    observed_cor.c <- sparcc_cor.c[i, j]
    perm_cor.c <- perm_cor_array.c[i, j, ]
    # Two-tailed p-value (percentage of permuted correlations greater than or equal to observed)
    p_values.c[i, j] <- mean(abs(perm_cor.c) >= abs(observed_cor.c))
  }
}

p_values_vector.c <- as.vector(p_values.c)
p_values_adjusted_vector.c <- p.adjust(p_values_vector.c, method = "fdr") 
p_values_adjusted.c <- matrix(p_values_adjusted_vector.c, nrow = nrow(p_values.c), ncol = ncol(p_values.c))

rownames(p_values_adjusted.c) <- rownames(sparcc_cor.c)
colnames(p_values_adjusted.c) <- colnames(sparcc_cor.c)

edge_indices.c <- as.data.frame(as_edgelist(g_sparcc.c))
row_indices.c <- edge_indices.c[, 1]
col_indices.c <- edge_indices.c[, 2]

edge_p_values.c <- mapply(function(row, col) p_values_adjusted.c[row, col], row_indices.c, col_indices.c)

edge_data.c <- data.frame(
  source = ends(g_sparcc.c, E(g_sparcc.c))[, 1],
  target = ends(g_sparcc.c, E(g_sparcc.c))[, 2],
  weight = E(g_sparcc.c)$weight,
  p_value = edge_p_values.c
)

write.csv(edge_data.c, "edge_data.c.csv", quote = F, row.names = F)

filtered_edge_data.c <- subset(edge_data.c, abs(weight) > 0.6 & p_value < 0.05)
filtered_edge_data.c$source_taxa <- sapply(filtered_edge_data.c$source, concatenate_taxa)

write.csv(filtered_edge_data.c, "Asymptomatic_edge_data.csv", row.names = FALSE, quote = FALSE)

###P. maydis analysis####
ASV1 <- as.character(ASVs.nochim["ASV1"])
alignmentASV1 <-lapply(ASVs.nochim, function(subject_seq) {
  pairwiseAlignment(pattern = ASV1, subject = subject_seq)
})

pid_values <- sapply(alignmentASV1, function(x) pid(x))
high_identity_ASVs <- which(pid_values >= 97)
high_identity_ASVs

list <- data.frame(name=c("ASV1","ASV4","ASV42","ASV57",
                          "ASV66","ASV84","ASV105","ASV137",
                          "ASV158","ASV192","ASV320","ASV333",
                          "ASV389","ASV396","ASV462","ASV463"))
seqs.Pmaydis <- ASVs.nochim[list$name]

writeXStringSet(seqs.Pmaydis,"pmaydis.fasta",  format="fasta")


selected_ASVs <- tmp.seqtab[, colnames(tmp.seqtab) %in% list$name]
selected_ASVs <-t(selected_ASVs)
tissue_groups <- split(colnames(selected_ASVs), samdf$Tissue)
total_reads <- data.frame(
  ASV = rownames(selected_ASVs),
  Symptomatic = rowSums(selected_ASVs[, tissue_groups$`Tar spot`]),
  Asymptomatic = rowSums(selected_ASVs[, tissue_groups$Asymptomatic])
)

write.csv(total_reads, "pmaydis_asvs_reads.csv", row.names = F, quote = F)

library(msa)

alignment_seqs.Pmaydis <- msa(seqs.Pmaydis)
aligned_seqs <- msaConvert(alignment_seqs.Pmaydis, type = "ape::DNAbin")
nucleotide_diversity <- nuc.div(aligned_seqs)

bootstrap_diversities <- numeric(n_permutations)
for (i in 1:n_permutations) {
  resampled_indices <- sample(1:nrow(aligned_seqs), replace = TRUE)  # Sample row indices
  resampled_seqs <- aligned_seqs[resampled_indices, ]  # Use resampled indices to select rows
  
  # Calculate nucleotide diversity for resampled data
  bootstrap_diversities[i] <- nuc.div(resampled_seqs)
}
p_value <- sum(bootstrap_diversities >= nucleotide_diversity) / n_permutations

###Exserohilum####
ASV17 <- ASVs.nochim["ASV17"]
writeXStringSet(ASV17,"ASV17.fasta",  format="fasta")
