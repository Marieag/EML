#Code written by M. R. Aggerbeck. 
#mrag@envs.au.dk
#Creative commons, 2024. 

###----------------PACKAGES--------------

library("tidyverse"); packageVersion("tidyverse")
library("ggpubr"); packageVersion("ggpubr")
library("data.table"); packageVersion("data.table")
library("igraph"); packageVersion("igraph")
library("plotrix"); packageVersion("plotrix")
library("VennDiagram"); packageVersion("VennDiagram")
library("phyloseq"); packageVersion("phyloseq")
library("devtools"); packageVersion("devtools")
library("biomformat"); packageVersion("biomformat")
library("vegan"); packageVersion("vegan")
library("DESeq2"); packageVersion("DESeq2")
library("mvabund"); packageVersion("mvabund")
library("metacoder"); packageVersion("metacoder")
library("taxa"); packageVersion("taxa")
library("adespatial"); packageVersion("adespatial")
library("microbiome"); packageVersion("microbiome")
library("pairwiseAdonis"); packageVersion("pairwiseAdonis")
library("microbiomeSeq"); packageVersion("microbiomeSeq")

###-----------Functions----------------

#----------ANOVA + pvalues------------------------------
plot_anova_diversity_pval <- function(physeq, method, grouping_column, bonf=TRUE, pValueCutoff=0.05)
{
  
  
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  abund_table <- otu_table(physeq)
  meta_table <- sample_data(physeq)
  
  #get diversity measure using selected methods
  div.df <- alpha_div(physeq,method)
  
  #=add grouping information to alpha diversity measures
  df<-data.frame(div.df,(meta_table[,grouping_column])[as.character(div.df$sample),])
  
  #perform anova of diversity measure between groups
  anova_res <- perform_anova(df,meta_table,grouping_column,pValueCutoff)
  df_pw <- anova_res$df_pw #get pairwise p-values
  
  #Apply Bonferroni correction
  if(bonf){  
    
    numberofsites <- length(unique(c(as.vector(levels(df_pw$from)),as.vector(levels(df_pw$to)))))[1]
    bonf.cor <- as.numeric(as.matrix(df_pw$p))* numberofsites
    temp <- as.factor(bonf.cor)
    df_pw$p <- temp
  }
  
  
  
  
  #Draw the boxplots
  p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
  p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
  p<-p+theme_bw()
  p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Alpha Diversity Measure")+xlab("")
  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("")
  
  
  
  #This loop will generate the lines and signficances
  if(!is.null(df_pw)){ #this only happens when we have significant pairwise anova results
    for(i in 1:dim(df_pw)[1]){
      p<-p+geom_path(inherit.aes=F,aes(x,y),
                     data = data.frame(
                       x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),
                             which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), 
                       y = c(as.numeric(as.character(df_pw[i,"y"])),
                             as.numeric(as.character(df_pw[i,"y"]))), 
                       measure=c(as.character(df_pw[i,"measure"]),
                                 as.character(df_pw[i,"measure"]))), 
                     color="black",lineend = "butt",
                     arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
      p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),
                     data=data.frame(
                       x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+
                            which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,
                       y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),
                       label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                                              breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                                              label=c("***", "**", "*", "")))))
    }
  }
  newlist <- list(p, df_pw)
  return(newlist)
}




###---------START OF CODE---------------

#Set working directory
setwd("O:/Tech_ENVS-EMBI-Afdelingsdrev/Marie/PhD/DATA/Non_targeted_screening/NTs_analysis")
#Set user directory
uzdir <- "O:/Tech_ENVS-EMBI-Afdelingsdrev/Marie/PhD/DATA/Non_targeted_screening/NTs_analysis"


###----------Import cleaned, annotated dataset from .csv---------------------
#Clean headers etc in txt editor/excel - remove regexes, special characters etc. 

LVL1 <- read.csv2(file = "LVL1_ID_metadata_FINAL_.csv", header = TRUE, sep=",")

dim(LVL1)

head(LVL1)
class(rownames(LVL1))
colnames(LVL1)

#Rownames will be unique identifiers for each compound - must be unique value! Use compound ID? 

rowname <- as.data.frame(LVL1[1])
rowname <- as.character(as.matrix(LVL1[1]))
rowname[is.na(rowname)] <- "unidentified"

head(rowname)
rowname <- paste("Compound", rowname, sep = "_")
rowname <- make.unique(rowname)

rawmetadata <- LVL1

#Make compound table (OTU table equivalent - area/counts) 
#Find rows with areas. 
COMP <- as.data.frame(LVL1[32:46])
COMP <- apply(COMP, 2, as.numeric)
head(COMP)
COMP[is.na(COMP)] = 0
colnames(COMP)
rownames(COMP) = rowname

head(COMP)
dim(COMP)

#Make drug class table 
#(Tax table equivalent - "taxonomy" of drug classes - Therapeutics;Illegal drugs;Metabolite etc)
#Remember to pick the proper columns
CLA <- LVL1[5:8]
rownames(CLA) = rowname
head(CLA)

#Create rownames for metadata table

temp <- as.data.frame(colnames(COMP))
temp

### Make changes to sample IDs if too long/weird/etc.:

# temp$`colnames(COMP)` <- gsub('Area..','', temp$`colnames(COMP)`)
# temp$`colnames(COMP)` <- gsub('MRA.','', temp$`colnames(COMP)`)
# temp$`colnames(COMP)` <- gsub('WWTP.','', temp$`colnames(COMP)`)
# temp$`colnames(COMP)` <- gsub('raw..','', temp$`colnames(COMP)`, fixed=TRUE)
# temp #just double checking everything is in order
# temp$`colnames(COMP)` <- gsub('.PB','_PB', temp$`colnames(COMP)`, fixed=TRUE)
# temp$`colnames(COMP)` <- gsub('.QC','_QC', temp$`colnames(COMP)`, fixed=TRUE)
# temp #just double checking again
# temp$`colnames(COMP)` <- gsub('\\..*','', temp$`colnames(COMP)`)
# temp #triple check? 
# 
# #Swap names
# rownames(temp) = temp$`colnames(COMP)`
# colnames(COMP) =rownames(temp)
# # rownames(temp) = colnames(COMP)

#Make small metadata table for initial data exploration. Alternatively, makein excel and read in here. 
#Sample names UST match other matrices

site_abbr = c("A","A","A","A","D","D",
              "D","D","QC", "QC","QC","S","S","S","S")
site = c("Avedøre","Avedøre","Avedøre","Avedøre","Damhusåen","Damhusåen",
         "Damhusåen","Damhusåen","QC_Pool", "QC_Pool","QC_Pool","Skævinge","Skævinge","Skævinge","Skævinge")
type = c("Sample","Sample","Sample","Procedural_Blank","Sample",
         "Sample","Sample","Procedural_Blank","QC","QC","QC", "Sample","Sample","Sample","Procedural_Blank")
PE = c(350000,350000,350000,350000,400000,400000,400000,400000,NA,NA,NA,12000,12000,12000,12000)
temp2 <- cbind(temp,site, site_abbr,type,PE)
colnames(temp2) = c("Sample_ID","Sample_Site","Site","Sample_Type","PE")
temp2
metadata <- sample_data(temp2)
rownames(metadata) = metadata$Sample_ID
head(metadata)

COMP_phy = otu_table(COMP, taxa_are_rows = TRUE)
CLA_phy = tax_table(as.matrix(CLA))

#doublecheck data
any(is.na(COMP))
any(is.na(CLA))
any(is.na(metadata))

#create phyloseq object
physeq = phyloseq(COMP_phy, CLA_phy)
physeq
physeq <- merge_phyloseq(physeq, metadata)
physeq

taxa_names(physeq) = rowname
head(taxa_names(physeq))

sample_names(physeq)

###Code for removal of sample, in case of e.g. contamination:
# physeq = subset_samples(physeq, Sample_ID != "All_samples" & Sample_ID != "All_samples.1")


###backup object in case of accidental erasure/change of primary physeq object
backup <- physeq

###Use this line to restore physeq 
###remember to run everything below this line again, to make sure nothing is messed up.
# physeq <- backup


#Create additional objects for removal of blanks and QC samples for simplified visualisations
physeq_noPB = subset_samples(physeq, Sample_ID != "Area_A_PB" & Sample_ID != "Area_D_PB" & Sample_ID != "Area_S_PB" )
physeq_noPB

physeq_noQC = subset_samples(physeq_noPB, Sample_ID != "Area_Pool_QC_1" &                                Sample_ID != "Area_Pool_QC_2" & Sample_ID != "Area_Pool_QC_3" )
physeq_noQC

drugs = subset_taxa(physeq_noQC, Class=="Therapeutics/Drugs")
drugs

unique(tax_table(physeq)[,2])

#Subset without unknowns
known <- subset_taxa(physeq, !Class=="Unknown structure")
known
unique(tax_table(known)[,2])


#no blanks
known_nopb <- subset_taxa(physeq_noPB, !Class=="Unknown structure")
known_nopb
unique(tax_table(known_nopb)[,2])

#samples only
known_so <- subset_taxa(physeq_noQC, !Class=="Unknown structure")
known_so
unique(tax_table(known_so)[,2])

#subset without unknown metabolites
knownmet <- subset_taxa(known, !Class=="Unidentified Natural Product")
knownmet
knownmet_nopb <- subset_taxa(known_nopb, !Class=="Unidentified Natural Product")
knownmet_nopb
knownmet_so <- subset_taxa(known_so, !Class=="Unidentified Natural Product")
knownmet_so

###Write table with corrected names for the paper table

write.table(as.data.frame(tax_table(physeq)), file='tax_table.tsv', quote=FALSE, sep='\t')

###Extra options: 
###-------log-transform data

#physeq_noQC_log10 <- microbiome::transform(physeq_noQC, "log10p")

###-------Relative abundance

#physeq  = transform_sample_counts(physeq, function(x) x / sum(x) )



###--------Merge samples into single site groups (viz only - no stats on this!)------

temp = prune_taxa(taxa_sums(physeq) > 0, physeq)
temp2 <- as.data.frame(unique(sample_data(temp)[,3]))
groups_phy <- as.vector(as.character(temp2$Site))
sample_data(temp)$group_phy <- get_variable(temp, "Site") %in% groups_phy

merged_phy = merge_samples(temp, "Site"); merged_phy
sample_names(temp)
sample_names(merged_phy)

sample_data(merged_phy)
sample_data(merged_phy)$Site = sample_names(merged_phy)
merged_phy; physeq

tmp <- sample_data(merged_phy)

class(tmp)
write.table(tmp, file='merged_metadata.tsv', quote=FALSE, sep='\t')


# ###------------Core compounds----------------
# NB! Does not work with this type of dataset - maybe set all below a certain threshold to 0? 

# 
# core.taxa.standard <- core_members(known_so, detection = 5/100, prevalence = 90/100)
# core.taxa.standard
# class(core.taxa.standard)
# write.table(core.taxa.standard, file='core_compounds.tsv', quote=FALSE, sep='\t')
# 
# Core_subset <- subset(otu_table(known_so), rownames(otu_table(known_so)) %in% core.taxa.standard)
# Core_subset
# new_physeq <- merge_phyloseq(my_subset, tax_table(physeq), sample_data(physeq), ...)



###----Richness + ANOVA ----
#Compare all experimental samples, 
#calculates diversity indices and tests pairwise differences (bonferroni correction)  

# temp <- physeq_noQC
# head(otu_table(temp))
# otu_table(temp) = round(otu_table(temp))
# 
# plot_richness(temp, x = 'Site', measures=c("Observed", "Shannon", "Simpson"))+ 
#   theme_bw()+
#   theme(legend.title = element_blank())+
#   theme(legend.position = "none")+
#   theme(strip.background = element_rect(fill="white" ))
# 
# p_an <-plot_anova_diversity_pval(physeq_noQC, method = c("richness","shannon", "simpson"),
#                                  grouping_column ="Site", bonf=TRUE, pValueCutoff=0.05)
# 
# print(p_an)
# p <- p_an[[1]]
# p
# ggsave(p, file="LVL1_noQC_richness.pdf",
#        width = 35, height = 20, units = "cm")

#----Barplot----

# Colour palette: 

set_class <- c("#8de985", "#412290","#019d4c","#9a0065","#01b6a8",
               "#ffacfd","#ae7400","#0195fb","#9a5184")
set_class10 <- c("#8de985", "#412290","#019d4c","#9a0065","#01b6a8",
                 "#ffacfd","#ae7400","#0195fb","#9a5184", "black")

# temp <- tax_glom(pseq.rel, taxrank = "Class")


p <- plot_bar(known_so, fill="Class", title="Compounds found in samples, except unknown structures")+
  theme_bw() +
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")+
  scale_fill_manual(values=set_class)+
  scale_color_manual(values=set_class)+
  labs(y= "Abundance", x = "Samples")
print(p)
ggsave(p, file="LVL1_noQC_barplot_Class.pdf", 
       width = 21, height = 14, units = "cm")

pseq.rel <- microbiome::transform(known_so, "compositional")

p <- plot_bar(pseq.rel, fill="Class", title="Compounds found in samples, except unknown structures")+
  theme_bw() +
  labs(y= "Relative Abundance", x = "Samples")+
  scale_fill_manual(values=set_class)+
  scale_color_manual(values=set_class)+
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")
print(p)
ggsave(p, file="LVL1_noQC_barplot_Class_RelAb.pdf", 
       width = 21, height = 14, units = "cm")

drugs.rel = subset_taxa(pseq.rel, Class=="Therapeutics/Drugs")

temp <- tax_glom(drugs, taxrank = "Subclass_I")
p <- plot_bar(drugs, fill="Subclass_I", title="Therapeutics/Drugs found in samples") + 
  theme_bw()+
  labs(y= "Abundance", x = "Samples")
print(p)
ggsave(p, file="LVL1_noQC_barplot_Subclass_Drugs.pdf", 
       width = 35, height = 20, units = "cm")

pseq.rel <- microbiome::transform(drugs, "compositional")
p <- plot_bar(pseq.rel, fill="Subclass_I", title="Compounds found in samples, -unknowns")+
  theme_bw() +
  labs(y= "Relative Abundance", x = "Samples")
# geom_bar(aes(color=Subclass_I, fill=Subclass_I), stat="identity", position="stack")
print(p)
ggsave(p, file="LVL1_noQC_barplot_Subclass_Drugs_RelAb.pdf", 
       width = 30, height = 20, units = "cm")

# p = plot_bar(ent10, "Genus", fill="Genus", facet_grid=SeqTech~Enterotype)
# p <- p + theme(axis.title.x =element_blank())+
#   geom_bar(aes(color="Subclass", fill="Subclass"), stat="identity", position="stack")+
#   theme_bw()+
#   theme(
#     strip.background = element_rect(fill="white" ))
print(p)


#---ordination----

#Choosing bray-curtis, as dataset contains null values and may not be parametric

phy.ord <- ordinate(physeq_noPB, "PCoA", "bray", pvalue.cutoff = 0.05)

p1 = plot_ordination(physeq_noPB, phy.ord, type="sites", color="Site", shape="Sample_Type", 
                     title="PCoA, Bray-Curtis Dissimilarity \n\nA, Samples vs QC")+
  geom_point(size=5)+
  scale_shape_manual(values=c(17, 15))+
  theme_bw()+
  # theme(legend.title = element_blank())+
  theme(#legend.position="bottom",
    strip.background = element_rect(fill="white" ))
print(p1)

p2 = plot_ordination(physeq_noPB, phy.ord, type="species", color="Class", 
                     title="\n\nB, Individual compounds")+
  geom_point(size=2, shape=16)+
  theme_bw()+
  scale_fill_manual(values=set_class10)+
  scale_color_manual(values=set_class10)+
  # theme(legend.title = element_blank())+
  facet_wrap(~Class, 3)+
  theme(#legend.position="bottom",
    strip.background = element_rect(fill="white" ))
print(p2)

###Extra variables for colour, labels etc. 
# p1 <- p1 + scale_color_discrete(labels = c("Compound", "A", "D", "QC", "S")) +
#   scale_size_discrete(labels = c("Sample", "Compound")) +
#   scale_shape_discrete(labels = c("Compound", "QC", "Sample"))
# 
# print(p1)

# # Merging legends
# legend_1 <- get_legend(p1)
# legend_2 <- get_legend(p2)
# 
# legends <- ggarrange(legend_1, legend_2, nrow=2)
# 
# # Combining plots
# rm_legend <- function(p){p + theme(legend.position = "none")}
# plots <- ggarrange(rm_legend(p1), rm_legend(p2), nrow = 1, align = "v")
# 
# # plots + merged legends
# p3 <- ggarrange(plots, legends, widths = c(0.75, 0.25))

p3 <- ggarrange(p1,p2,#legend_1, legend_2,
                widths = 1, heights = 1.1,
                # labels = "PCoA, Bray-Curtis Dissimilarity",
                #label.x = 0.8, label.y = .90, hjust = 0, vjust = 0,
                ncol = 2, nrow = 1
                # common.legend = TRUE, legend = "right"
)
p3

ggsave(p3, file="LVL1_Ordination_SamplesvsQC_comp.pdf", 
       width = 30, height = 15, units = "cm", useDingbats=FALSE)

phy.ord <- ordinate(known_nopb, "PCoA", "bray", pvalue.cutoff = 0.05)
p1 = plot_ordination(known_nopb, phy.ord, type="Biplot", color="Sample_type", shape="Sample_Type", 
                     title="PCoA, bray-curtis dissimilarity, -unknowns")+
  geom_point(size=2)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_rect(fill="white" ))
print(p1)

ggsave(p1, file="LVL1_Ordination_SamplesvsQC_col.pdf", 
       width = 25, height = 20, units = "cm", useDingbats=FALSE)


p2 <- p1 + facet_wrap(~Class, 3)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_rect(fill="white" ))

print(p2)

ggsave(p2, file="LVL1_Ordination_SamplesvsQC_col_wrapped.pdf", 
       width = 25, height = 20, units = "cm", useDingbats=FALSE)


phy.ord <- ordinate(knownmet_so, "PCoA", "bray", pvalue.cutoff = 0.05)
p1 = plot_ordination(knownmet_so, phy.ord, type="Split", color="Class", shape="Site", 
                     title="PCoA, bray-curtis dissimilarity")+
  geom_point(size=2)+
  scale_shape_manual(values=c(16,15,17,8))+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_rect(fill="white" ))
print(p1)

ggsave(p1, file="LVL1_Ordination_PCoA_bray_Class.pdf", 
       width = 35, height = 20, units = "cm", useDingbats=FALSE)


p2 <- p1 + facet_wrap(~Class, 3)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_rect(fill="white" ))

print(p2)

ggsave(p2, file="LVL1_Ordination_PCoA_bray_Class_fw.pdf", 
       width = 35, height = 20, units = "cm", useDingbats=FALSE)

#----ordination, drugs----

phy.ord <- ordinate(drugs, "PCoA", "bray", pvalue.cutoff = 0.05)
p1 = plot_ordination(drugs, phy.ord, type="biplot", color="Subclass_I", shape="Site", 
                     title="Therapeutics/Drugs \nPCoA, bray-curtis dissimilarity")+
  geom_point(size=2)+
  scale_shape_manual(values=c(16,15,17,8))+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_rect(fill="white" ))
print(p1)

ggsave(p1, file="LVL1_drugs_Ordination_PCoA_bray_Subclass.pdf", 
       width = 35, height = 20, units = "cm", useDingbats=FALSE)


p2 <- p1 + facet_wrap(~Subclass_I, 3)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_rect(fill="white" ))

# p2 + geom_polygon(aes(fill=Sample_Type)) + geom_point(size=5)

print(p2)

ggsave(p2, file="LVL1_drugs_Ordination_PCoA_bray_Subclass_fw.pdf", 
       width = 50, height = 20, units = "cm", useDingbats=FALSE)


###----Heatmap-----

tmp <- known_so

tmp <- prune_taxa(names(sort(taxa_sums(tmp),TRUE)[1:50]), tmp)
p <- plot_heatmap(tmp, method="NMDS", distance="bray", sample.label="Site", taxa.label = "Name", 
                  sample.order = "Sample_Type", title="Heatmap, NMDS + Bray-Curtis") +
  theme(axis.title.y=element_blank(), axis.title.x =element_blank()) +
  labs(fill = "Area")
print(p)
ggsave(p, file="LVL1_noQC_Heatmap.pdf", width = 35, height = 20, units = "cm")

tmp <- prune_taxa(names(sort(taxa_sums(tmp),TRUE)[1:50]), tmp)
p <- plot_heatmap(tmp, method="NMDS", distance="bray", sample.label="Site", 
                  sample.order = "Sample_Type", title="Heatmap, NMDS + Bray-Curtis, 50 most abundant compounds") +
  theme(axis.title.y=element_blank(), axis.title.x =element_blank()) +
  labs(fill = "Area")
print(p)
ggsave(p, file="LVL1_noQC_Heatmap_50mostab.pdf", width = 35, height = 20, units = "cm")


###--------PERMANOVA + post-hoc pairwise---------------------------

pseq <- physeq_noQC
pseq.rel <- microbiome::transform(known_so, "compositional")
otu <- microbiome::abundances(known_so)
meta <- microbiome::meta(known_so)


permanova <- vegan::adonis(t(otu)~ Site,
                           data = meta, permutations=999, method = "bray")

permanova$aov.tab

post_hoc_permanova <- pairwise.adonis(t(otu), meta$Site, sim.function = "vegdist",
                                      sim.method = "bray", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova


###--------DESeq2------------

#Be aware that Deseq2 has a min and a max value. min=0, max= 2.9e9. Which is occasionally exceeded 
# when using chemical data (area). DESeq also gets angry if you log transform. 
#Scale down all counts across entire dataset by 100 or 1000, to remove this error.

head(otu_table(physeq_noQC))
test_phy <- (otu_table(physeq_noQC))/100
head(test_phy)

#DESeq2

temp <- physeq_noQC
otu_table(temp) = test_phy
head(otu_table(temp))

diagdds = phyloseq_to_deseq2(temp, ~Site)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
print(diagdds)

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(temp)[rownames(sigtab), ], "matrix"))
head(sigtab)

# theme_set(theme_bw())

# Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
# Subclass order
x = tapply(sigtab$log2FoldChange, sigtab$Subclass_I, function(x) max(x))
x = sort(x, TRUE)
sigtab$Subclass_I = factor(as.character(sigtab$Subclass_I), levels=names(x))
ggplot(sigtab, aes(x=Subclass_I, y=log2FoldChange, color=Class)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# OTU <- otu_table(physeq_noQC, taxa_are_rows = TRUE)
# OTU <- t(OTU)
# #Merge new phyloseq object
# TAX = as.matrix(tax_table(physeq_noQC))
# met <- as.data.frame(sample_data(physeq_noQC))
# pseq1 <- phyloseq(OTU, TAX)
# physeq2<-merge_phyloseq(pseq1, met)
# rm(OTU, TAX, met, pseq1)
# 
# minTotRelAbun = 0.1
# x = taxa_sums(physeq2)
# keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
# phy3 = prune_taxa(keepTaxa, physeq2)
# phy3 = physeq2
# 
# deseq_sig <- differential_abundance(physeq_noQC, grouping_column = "Site",
#                                     pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = "LVL1_site_DeSeq2")


###---------Heat tree-------------------

# temp <- drugs
temp <- knownmet_so
# temp <- physeq_noQC

tax_table(temp) <- tax_table(temp)[,2:4]
head(tax_table(temp))


# temp <- tax_glom(temp, taxrank = "Subclass_II")


#parse phyloseq object temp

obj_all <- parse_phyloseq(temp)

obj <- obj_all
tissuegroup <- obj$data$sample_data$Site

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", 
                                      cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")


heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "fruchterman-reingold",
                 initial_layout = "reingold-tilford",
                 node_size_axis_label = "Number of compounds",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'LVL1_heattree_knownmet_so.pdf')


# x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
#                    class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
#                    class_regex = "^(.+)__(.+)$")
# 
# x
# # Get per-taxon counts
# x$data$tax_table <- calc_taxon_abund(x, data = "tax_data", cols = hmp_samples$sample_id)
# 
# # Calculate difference between groups
# x$data$diff_table <- calc_diff_abund_deseq2(x, data = "tax_table",
#                                             cols = hmp_samples$sample_id,
#                                             groups = hmp_samples$body_site)

# Remove taxa with only small differences 
per_taxon_fold_changes <- obs(obj, data = 'diff_table', value = 'log2_median_ratio')
per_taxon_fold_changes
per_taxon_max_change <- unlist(lapply(per_taxon_fold_changes, function(tax_changes) max(abs(tax_changes))))
per_taxon_max_change
x <- filter_taxa(obj, per_taxon_max_change > 1, supertaxa = TRUE, reassign_obs = c(diff_table = FALSE))
x

heat_tree_matrix(x,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_size_axis_label = "Number of compounds",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'LVL1_heattree_LfC1.pdf')


# #Test p-values
# # obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
# test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.1] <- 0 )
# obj$data$diff_table <- test
# 
# heat_tree_matrix(obj,
#                  data = "diff_table",
#                  node_size = n_obs,
#                  node_label = taxon_names,
#                  node_color = log2_median_ratio,
#                  node_color_range = diverging_palette(),
#                  node_color_trans = "linear",
#                  node_color_interval = c(-3, 3),
#                  edge_color_interval = c(-3, 3),
#                  layout = "davidson-harel",
#                  initial_layout = "reingold-tilford",
#                  node_size_axis_label = "Number of compounds",
#                  node_color_axis_label = "Log2 ratio median proportions",
#                  output_file = 'LVL1_heattree_pval.pdf')


# ###---------Venn Diagram----------------------------
# 
#Calculate presence in a separate 
# 
# Venn <- read.csv2(file = "venn.csv") 
# colnames(Venn)
# 
# overlap <- calculate.overlap(x = list('A' = set1, 'D' = set2, 'S' = set3))
# area_overlap <- sapply(overlap, length)
# 
# ## areas of each animal
# area1 <- as.numeric(nrow(subset(Venn, A == 1)))
# ## [1] 135
# area2 <- as.numeric(nrow(subset(Venn, D == 1)))
# ## [1] 99
# area3 <- as.numeric(nrow(subset(Venn, S == 1)))
# ## [1] 92
# 
# 
# ## areas of 2-group overlap
# 
# ## A & D
# n12 <- as.numeric(nrow(subset(Venn, A == 1 & D == 1)))
# ## [1] 75
# 
# ## A & S
# n13 <- as.numeric(nrow(subset(Venn, A == 1 & S == 1)))
# ## [1] 52
# 
# ## D & S
# n23 <- as.numeric(nrow(subset(Venn, D == 1 & S == 1)))
# ## [1] 53
# 
# ## 3 group overlap: A, D & S
# n123 <- as.numeric(nrow(subset(Venn, A == 1 & D == 1 & S == 1)))
# ## [1] 47
# 
# draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category =
#                    c("A", "D", "S"), lty = "blank", 
#                  col = rep("black", 3), fill = c("darkblue", "darkorange", "turquoise"))

