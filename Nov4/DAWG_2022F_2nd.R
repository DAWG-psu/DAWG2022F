### download and load phyloseq object
### run below two lines in terminal
### cd ~/scratch
### wget https://github.com/DAWG-psu/DAWG2022F/raw/main/Nov4/ps_16s.rds
setwd("~/scratch")
phyloseq_16s <- readRDS("ps_16s.rds")

### Remove ASVs from non-bacterial taxa
phyloseq_16s_filt <- subset_taxa(phyloseq_16s, domain == "Bacteria")

### Compositional anaylsis of overall microbiome
### 1. stacked bar plot
### convert phyloseq object into designated level of taxonomy (genus)
library(phyloseq)
rank_names(phyloseq_16s_filt) # check the taxonomy structure
phyloseq_genus <- tax_glom(phyloseq_16s, taxrank = "genus") # combine ASVs into genus level
phyloseq_genus_rel <- transform_sample_counts(phyloseq_genus, function(x) x/sum(x)) # convert data into relative abundance
phyloseq_genus_rel
phyloseq_genus_rel_melt <- psmelt(phyloseq_genus_rel) #melt data
str(phyloseq_genus_rel_melt) # check the data structure
phyloseq_genus_rel_melt[,"genus"] <- as.character(phyloseq_genus_rel_melt[,"genus"])
phyloseq_genus_rel_melt[,"genus"][phyloseq_genus_rel_melt$Abundance < 0.1] <- "Other" #convert low abundance taxa into "Other"

# Check how many genus in your dataset after converting low abundance taxa into "Others"
Count = length(unique(phyloseq_genus_rel_melt[,"genus"]))
if (Count > 25) {
  print("Warning: You have more than 25 taxa to plot, consider using higher cut off")
} else {
  print(Count)
} 

# Sample Annotation - optional (Add orders next to genus)
phyloseq_genus_rel_melt$sample_order <- phyloseq_genus_rel_melt[,"genus"]
phyloseq_genus_rel_melt$sample_order <- as.factor(phyloseq_genus_rel_melt$sample_order)
levels(phyloseq_genus_rel_melt$sample_order) <- seq(1:length(levels(phyloseq_genus_rel_melt$sample_order))) 
phyloseq_genus_rel_melt$text <- paste0(phyloseq_genus_rel_melt[,"genus"], "(", phyloseq_genus_rel_melt$sample_order,")")
#add numbers next to the genus name based on the alphabetic order
head(phyloseq_genus_rel_melt$text)

#getting color code for ploting - randomly assigned color
phyloseq_genus_rel_melt[,"genus"] <- as.factor(phyloseq_genus_rel_melt[,"genus"])
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category =='qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, Count)

#Aggregate all relative abundacneb based on the group
phyloseq_genus_rel_melt_aggr <- aggregate(Abundance ~ Sample+text+sample_order+Type_of_sample+genus, 
                                          data = phyloseq_genus_rel_melt, FUN=sum)

#plot stacked barplot
library(ggplot2)
plot <- ggplot(data=phyloseq_genus_rel_melt_aggr, aes(x=Sample, y=Abundance, fill= text)) + 
  geom_bar(stat="identity", position="stack") + 
  facet_grid(.~Type_of_sample, scales= "free") +
  guides(fill=guide_legend(nrow=5)) + 
  scale_fill_manual(values=col)+
  theme(legend.position = "bottom", axis.text.x = element_text(angle=90), 
        axis.ticks.x = element_blank(), 
        panel.background = element_blank(),
        legend.title = element_blank()) +
  xlab("Samples") + 
  ylab("Relative abundance") + 
  geom_text(aes(label=sample_order), size = 1.5, position = position_stack(vjust = 0.5)) +
  theme(legend.text=element_text(size=5.5))
plot
ggsave("stacked_bar_plot_with_numbers.jpg", plot=plot, width = 12.5, height= 5, units = "in", dpi = 600)

#without numbers
plot_without_numbers <- ggplot(data=phyloseq_genus_rel_melt, aes(x=Sample, y=Abundance, fill= genus)) + 
  geom_bar(stat="identity", position="stack") + 
  facet_grid(.~Type_of_sample, scales= "free") +
  guides(fill=guide_legend(nrow=5)) + 
  scale_fill_manual(values=col)+
  theme(legend.position = "bottom", axis.text.x = element_text(angle=90), 
        axis.ticks.x = element_blank(), 
        panel.background = element_blank(),
        legend.title = element_blank()) +
  xlab("Samples") + 
  ylab("Relative abundance") 
plot_without_numbers
ggsave("stacked_bar_plot_without_numbers.jpg", plot=plot_without_numbers, width = 12.5, height= 5, units = "in", dpi = 600)

## Diversity analysis
## 1. Alpha Diversity
library(ggpubr)
library(reshape2)
alpha_diversity <- estimate_richness(phyloseq_16s_filt, measures= c("Observed", "Shannon", "InvSimpson"))
alpha_comparison <- cbind(alpha_diversity, sample_data(phyloseq_16s_filt))
melt_plot <- melt(alpha_comparison)
plot_alpha <- ggplot(data = melt_plot, aes(y = value, x = Type_of_sample, fill = Type_of_sample)) +
  geom_boxplot() +
  facet_wrap(~variable, scale="free") +
  theme_minimal() + 
  scale_fill_brewer(palette="Dark2") + 
  stat_compare_means(label = "p.format") 
plot_alpha
ggsave("alpha_diversity.jpg", plot = plot_alpha, width = 10, height=7, units= "in", dpi =600)


## 2. Beta diversity
phyloseq_genus
# Data normalization (CLR - Central log ratio transformation)
# https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full
#Following step requires samples on rows and ASVs in columns
otus <- otu_table(phyloseq_genus)
#Replace zero values before clr transformation
#Use CZM method to replace zeros and outputs pseudo-counts (1)
require(zCompositions)
otu.n0 <- t(cmultRepl(otus, label =0, method="CZM", output="p-counts"))
otu.n0 <-ifelse(otu.n0 < 0, otu.n0*(-1), otu.n0)
#Convert data to proportions
otu.n0_prop <- apply(otu.n0, 2, function(x) {x/sum(x)})
#CLR transformation
otu.n0.clr<-t(apply(otu.n0_prop, 2, function(x){log(x)-mean(log(x))}))
phyloesq_genus_CLR <- phyloseq(otu_table(otu.n0.clr, taxa_are_rows=F), 
                               tax_table(phyloseq_genus), 
                               sample_data(phyloseq_genus))
#PCA
pc.clr <- prcomp(otu.n0.clr)
require(compositions)
# Calculate total variance
mvar.clr <- mvar(otu.n0.clr)
row <- rownames(otu.n0.clr)
# extract first two PCs
pc_out <- as.data.frame(pc.clr$x[,1:2])
# combine first two PCs with metadata
pc_out_meta <- as.data.frame(cbind(pc_out, sample_data(phyloesq_genus_CLR)))

# Calculate euclidean distance between sample
dist <- dist(otu.n0.clr, method = "euclidean")

#PERMANOVA (Permutational based ANOVA)
permanova <- vegan::adonis2(dist~Type_of_sample, data=pc_out_meta, perm =999)
permanova
label = paste0("p-value = ",permanova$`Pr(>F)`[1]," (PERMANOVA)") #prepare label for plot
label
#Getting colors based on your comparisons
qual_col_pals = brewer.pal.info[brewer.pal.info$category =='qual',]
col_vector <- unlist(mapply(brewer.pal, 8, "Dark2"))
col=sample(col_vector, length(unique(pc_out_meta[,"Type_of_sample"])))

#Making PCA plot
PCA_plot <- ggplot(pc_out_meta, aes(x=PC1, y=PC2, color = Type_of_sample)) +
  geom_point(size=2)+
  theme(legend.position = 'right')+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr$sdev[1]^2/mvar.clr*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr$sdev[2]^2/mvar.clr*100, digits=1), "%", sep="")) +
  scale_color_manual(values= col) +
  theme(legend.title = element_blank()) +
  annotate("text", x = 15, y =-25, label = label)
PCA_plot
ggsave("pca.jpg", plot = PCA_plot, width= 10, height=7.5, units="in", dpi=600)

