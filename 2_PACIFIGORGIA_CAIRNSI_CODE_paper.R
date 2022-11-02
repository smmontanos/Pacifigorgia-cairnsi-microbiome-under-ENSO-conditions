##Code for the analysis of bacterial community associated with Pacifigorgia cairnsi under ENSO conditions
##By Sandra Montaño Salazar

#Load the libraries

library(phyloseq)
library(devtools)
library(ggplot2)
library(readxl)
library(dplyr)
library(ape)
library(gridExtra)
library(knitr)
library(vegan)
library(tidyverse)
library(scales)
library(grid)
library(reshape2)
library(microbiome)
library(ggpubr)
library(RColorBrewer)
library(microbiomeutilities)
library(viridis)
library(tibble)

#Input files from qiime2 outputs for create a phyloseq object
asv<- read_excel("E:/CAIRNSI/OTU_236.xlsx")
tax<- read_excel("E:/CAIRNSI/TAX_236.xlsx")
metadata <- read_excel("E:/CAIRNSI/metadatosCAIRNSI.xlsx")

asv <- data.frame(asv, row.names = 1)
tax <- data.frame(tax, row.names = 1)
samples <- data.frame(metadata, row.names = 1) 
asv <- as.matrix(asv)
tax <- as.matrix(tax)

#Create a phyloseq object
ASV = otu_table(asv, taxa_are_rows = TRUE)
TAX = tax_table(tax)
samples = sample_data(samples)
phca <- phyloseq(ASV, TAX, samples)
phca

#Modified Protocol from Callahan et al. 2016 (Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses)

# Create table, number of features for each phyla
table(tax_table(phca)[, "Phylum"], exclude = NULL)
ph1 <- subset_taxa(phca, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ph1
# Compute prevalence of each feature, store as data.frame
prevph = apply(X = otu_table(ph1),
               MARGIN = ifelse(taxa_are_rows(ph1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevph = data.frame(Prevalence = prevph,
                    TotalAbundance = taxa_sums(ph1),
                    tax_table(ph1))
plyr::ddply(prevph, "Phylum", function(ph1){cbind(mean(ph1$Prevalence),sum(ph1$Prevalence))})
# Define phyla to filter
filterPhylapp = c("D_1__Epsilonbacteraeota", "D_1__Omnitrophicaeota", "D_1__Fusobacteria")
# Filter entries with unidentified Phylum.
ph2 = subset_taxa(ph1, !Phylum %in% filterPhylapp)
ph2
# Subset to the remaining phyla
prevph1 = subset(prevph, Phylum %in% get_taxa_unique(ph2, "Phylum"))
ggplot(prevph1, aes(TotalAbundance, Prevalence / nsamples(phca),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#For create a Rarefaction curve (supplementary file)
rarecurve(t(otu_table(phca)), step=50, cex=0.5)

#Rarefy the data for calculate alpha and beta diversity

ps.rarefied = rarefy_even_depth(ph2, rngseed=4567, replace=F)

##Alpha diversity and plot

alpha.diversity <- estimate_richness(ps.rarefied, measures = c("Observed","Chao1", "Shannon","Simpson"))
alpha.diversity#table with values

prich <- plot_richness(ps.rarefied, x = "yecat", measures = c("Observed", "Chao1", "Shannon", "Simpson"))
prich <- prich + geom_boxplot(aes(fill = yecat), alpha=0.2)
plot(prich)


##Statistical test to evaluate significant effect of host specie in richness

#for normal data
data <- cbind(sample_data(ps.rarefied), alpha.diversity)
shapiro.test(data$Observed)

phy.anova <- aov(Observed ~ yecat, data)
summary(phy.anova)
TukeyHSD(phy.anova)

#for not normal data

adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(ps.rarefied, measures = "Observed"),
  "Simpson" = phyloseq::estimate_richness(ps.rarefied, measures = "Simpson"),
  "Chao1" = phyloseq::estimate_richness(ps.rarefied, measures = "Chao1"),
  "Shannon" = phyloseq::estimate_richness(ps.rarefied, measures = "Shannon"),
  "yecat" = phyloseq::sample_data(ps.rarefied)$yecat,
  "temp" = phyloseq::sample_data(ps.rarefied)$temperature)

phy.kruskal<- kruskal.test(Simpson ~ yecat, data = adiv)
phy.kruskal
pairwise.wilcox.test(adiv$Chao1.Chao1, adiv$yecat, p.adjust.method="fdr")

##models gml with temperature
glm.temp = glm(Shannon ~ temp, data = adiv)
summary(glm.temp)
plot(Shannon ~ temp, data = adiv)
abline(glm.temp)

##Beta diversity and NMDS ordination 

ps_bray <- phyloseq::distance(ph2, method = "bray")

ps_nmds <- ordinate(
  physeq = ph2, 
  method = "NMDS", 
  distance = "bray"
)
plot_ordination(
  physeq = ph2,
  ordination = ps_nmds,
  color = "yecat",
  title = "NMDS bacterial Communities"
) + 
  scale_color_manual(values = c("blue", "red","green", "#95A900", "#00B81F","#00BF7D","#00C0B8", "#00B8E5", "#00A5FF", "#AC88FF", "#E76BF3", "#FF62BC", "#FF6C90" )
  ) +
  geom_point(aes(color = yecat), alpha = 1, size = 3, shape = 16)

#PERMANOVA for evaluate differences between the four types of samples using Bray Curtis dissimilarities

metbd <- as(sample_data(ps.rarefied), "data.frame")

adonis(distance(ps.rarefied, method="bray") ~ yecat,
       data = metbd)
#for betadisper

dispr <- vegan::betadisper(ps_bray, phyloseq::sample_data(ps.rarefied)$yecat)
dispr
permutest(dispr)

##Taxa barplot exploration

#To graph at phylum level and the relative abundances

ps_phylum <- phyloseq::tax_glom(ph2, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]

phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = yecat, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

#To explore at genus level, relative abundances per year and the sd of the procentages 

top150 <- names(sort(taxa_sums(ph2), decreasing=TRUE)[1:150])
top150 #shows 150 results

dat.aglo = tax_glom(ph2, taxrank = "Genus")
test1 = psmelt(dat.aglo)
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
dat.trans
prune.dat.two = prune_taxa(top150, dat.trans)
dat.dataframe = psmelt(prune.dat.two)
dat.agr = aggregate(Abundance~yecat+Genus, data=dat.dataframe, FUN=mean)

ggplot(dat.agr, aes(x=yecat, y=Abundance, fill=Genus)) + geom_bar(stat="identity", color = "black")
prune.dat.two

by(dat.agr$Abundance,dat.agr$Genus,sd)#sd for genus relative abundance

##For the differential expressed ASV according to their relative abundance using Deseq2 test 

library(DESeq2)
pds <- phyloseq_to_deseq2(ph2, ~yecat)
pds = DESeq(pds)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(pds), 1, gm_mean)

pds = estimateSizeFactors(pds, geoMeans = geoMeans)
pds = DESeq(pds)
res <- results(pds, cooksCutoff = FALSE)
alpha=0.05
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab,"data.frame"), as(tax_table(ph2)[rownames(sigtab), ], "matrix"))
head(sigtab,20)

##For the Heatmap
ps1.c <- format_to_besthit(ph2)

heat.sample <- plot_taxa_heatmap(phca, subset.top = 25,
                                 VariableA = "yecat",
                                 heatcolors = brewer.pal(100, "Blues"),
                                 transformation = "compositional")

brewer.pal(100, "Blues")

## For the Venn diagram

myotutable <- read.csv("~/CAIRNSI/venn_otu_ampvis.csv", check.names = FALSE)
mymetadata <- read.csv("~/CAIRNSI/met_venn_ampvis.csv", check.names = FALSE)

library(ampvis2)
d <- amp_load(otutable = myotutable,
              metadata = mymetadata)
amp_venn(
  d,
  group_by = "Year",
  cut_a = 0.1,
  cut_f = 10,
  text_size = 5,
  normalise = TRUE,
  detailed_output = FALSE
)


