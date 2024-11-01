if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
devtools::install_github("RIVM-IIV-Microbiome/biomeUtils")
library(biomeUtils)

BiocManager::install("biomformat")
BiocManager::install("phyloseq",force = TRUE)
library(microbiomeMarker)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(biomformat)
library(phyloseq)
library(dplyr)
library(patchwork)

library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(readxl)
library(arrow)
library(stringr)

# set.seed(0)
# tax_table(kostic_crc) 
# otu_table(kostic_crc)
# # ----import-dada2-------------------------------------------------------------
# seq_tab <- readRDS("./rawdata/HN00199038_ASV/ASV_Results/biom/ASVs.rds")
# tax_tab <- read_biom("./rawdata/HN00199038_ASV/ASV_Results/biom/ASV_table.blast_NCBI_16S.biom")
# sam_tab <- read.table("./rawdata/HN00199038_ASV/ASV_Results/biom/metadata.txt",header = F,
#                       row.names = 1)

library("phyloseq")
PT <- read_excel('./0_PT_PROFILE/PT_CODE.xlsx')
# OTU = otu_table(openxlsx::read.xlsx("./microbiome_marker/otu_table.xlsx"), taxa_are_rows = TRUE) 
OTU = otu_table(openxlsx::read.xlsx("./microbiome_marker/otu_table.xlsx",rowNames = T)%>% t() %>% as.data.frame()%>% rownames_to_column("SampleName") %>% left_join(PT) %>% mutate(SampleName=NULL) %>% column_to_rownames("PTCODE") %>% t(),taxa_are_rows = TRUE) 
TAX = tax_table(as.matrix(openxlsx::read.xlsx("./microbiome_marker/tax_table.xlsx",rowNames = T)))

color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
                   "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                   "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072")

color_palette_spp <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
  "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#8C564B", "#C49C94", "#C5B0D5", "#F7B6D2", "#C7C7C7", "#BCBD22", "#17BECF",
  "#1F77B4", "#AEC7E8", "#FFBB78", "#98DF8A", "#D62728", "#9467BD", "#8C564B",
  "#E377C2", "#7F7F7F", "#BDBDBD", "#17BECF", "#9EDAE5", "#FF69B4", "#8B4513", 
  "#2E8B57", "#FFD700", "#ADFF2F", "#00CED1", "#DC143C", "#8A2BE2", "#FF4500"
)
A <- plot_bar(copd_micro_phyloseq_inclusion_abund_fraction, fill = "Phylum") + 
  geom_bar(aes(fill=Phylum), stat="identity", position="stack")+scale_fill_brewer(palette="Set3") +
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +theme_bw()

B <- plot_bar(copd_micro_phyloseq_inclusion_abund_fraction_top_20, fill = "Genus") + 
  geom_bar(aes(fill=Genus), stat="identity", position="stack")+scale_fill_manual(values = color_palette) +
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme_bw()

C <- plot_bar(copd_micro_phyloseq_inclusion_abund_fraction_top_30_species, fill = "Species") + 
  geom_bar(aes(fill=Species), stat="identity", position="stack")+scale_fill_manual(values = color_palette_spp) +
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme_bw()



jpeg("./Figure_1_copd.jpeg",width = 7000,height = 5000,res=600)
A+B+C+plot_layout(ncol = 2)
dev.off()

jpeg("./Figure_1_2.jpeg",width = 4000,height = 3000,res=600)
C+plot_layout()
dev.off()






#################### By PATIENTS ############
copd_micro_phyloseq_abund_fraction <- copd_micro_phyloseq
sample_names(copd_micro_phyloseq_abund_fraction) <-  as.data.frame(sample_data(copd_micro_phyloseq_abund_fraction))$PTCODE

# copd_micro_phyloseq_abund_fraction <- merge_samples(copd_micro_phyloseq, "COPD_STAGE")
copd_micro_phyloseq_abund_fraction <- transform_sample_counts(copd_micro_phyloseq_abund_fraction, function(x) 100 * x / sum(x))

##GENUS
genus_abundance <- data.frame(FREQ=taxa_sums(copd_micro_phyloseq),ASV=names(taxa_sums(copd_micro_phyloseq))) %>% left_join(data.frame(tax_table(copd_micro_phyloseq),ASV=row.names(data.frame(tax_table(copd_micro_phyloseq)))))
top20_genera <- (genus_abundance %>% group_by(Genus)%>%summarise(Sum_GENUS=sum(FREQ)) %>% arrange(desc(Sum_GENUS)) %>% slice(1:20))$Genus
top30_species <- (genus_abundance %>% group_by(Species)%>%summarise(Sum_Species=sum(FREQ)) %>% arrange(desc(Sum_Species)) %>% slice(1:30))$Species

copd_micro_phyloseq_abund_fraction_top_20 <- copd_micro_phyloseq_abund_fraction
tax_table(copd_micro_phyloseq_abund_fraction_top_20) <- tax_table(copd_micro_phyloseq_abund_fraction_top_20) %>%
  as.data.frame() %>%
  mutate(Genus = ifelse(Genus %in% top20_genera, Genus, "Others"),
         Genus = factor(Genus,levels=c(setdiff(Genus, "Others"),"Others"))) %>%
  as.matrix()


copd_micro_phyloseq_abund_fraction_top_30_species <- copd_micro_phyloseq_abund_fraction
tax_table(copd_micro_phyloseq_abund_fraction_top_30_species) <- tax_table(copd_micro_phyloseq_abund_fraction_top_30_species) %>%
  as.data.frame() %>%
  mutate(Species = ifelse(Species %in% top30_species, Species, "Others"),
         Species = factor(Species,levels=c(setdiff(Species, "Others"),"Others"))) %>%
  as.matrix()

Genus_all <- plot_bar(copd_micro_phyloseq_abund_fraction_top_20, fill = "Genus") + 
  geom_bar(aes(fill = Genus), stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette) +
  ylab("Relative abundance (%)") +
  facet_wrap(~ GROUP_SAMPLE, scales = "free_x", nrow = 1) +
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))

Phylum_all <- plot_bar(copd_micro_phyloseq_abund_fraction, fill = "Phylum") + 
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  scale_fill_brewer(palette="Set3") +
  ylab("Relative abundance (%)") +
  facet_wrap(~ GROUP_SAMPLE, scales = "free_x", nrow = 1) +
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))

Species_all <- plot_bar(copd_micro_phyloseq_abund_fraction_top_30_species, fill = "Species") + 
  geom_bar(aes(fill = Species), stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette_spp) +
  ylab("Relative abundance (%)") +
  facet_wrap(~ GROUP_SAMPLE, scales = "free_x", nrow = 1) +
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))

jpeg("./Figure_1SS.jpeg",width = 4000,height = 4000,res=300)
Phylum_all+Genus_all+Species_all+plot_layout(ncol = 1)
dev.off()


############### Alpha diversity ##########################

                                                              
shannon_data<- 
  estimate_richness(copd_micro_phyloseq, split = TRUE, measures = NULL) %>%
  mutate(SampleName=gsub("X", "",row.names(.)),
         SampleName=gsub("\\.", "-",SampleName)) %>% 
  left_join(merged_data%>%  mutate(SampleName=gsub("_", "-",PTCODE )))

moonBook::mytable(GROUP_SAMPLE~Observed+Shannon+Fisher,data=shannon_data,method=3)
shannon_data <- shannon_data %>% mutate(Antibiotics_Use_2=case_when(PT_NO %in% Ant_USE~"Antibiotics use",TRUE~"Antibiotics non-use"),
Steroids_Use_2=case_when(PT_NO %in% Ste_USE ~1,TRUE~0))
shannon_data <- shannon_data %>% mutate(Season_2=case_when(Season %in% c("spring","summer")~"S/S",TRUE~"F/W"))
shannon_data <- shannon_data %>% mutate(Season_1=factor(Season_2,levels=c("spring","summer","autumn","winter")))
shannon_data <- shannon_data %>% mutate(Season_2=factor(Season_2,levels=c("S/S","F/W")))

shannon_data <- shannon_data %>% mutate(Season_3=case_when(Season %in% c("summer","autumn")~"high humidity",TRUE~"low humidity"))
shannon_data <- shannon_data %>% mutate(Season_3=factor(Season_2,levels=c("low humidity","high humidity")))


moonBook::mytable(AgeGroup~Observed+Shannon+Fisher,data=shannon_data%>% dplyr::filter(GROUP_SAMPLE != "AE-2"),method=3)
moonBook::mytable(COPD_stage_Group~Observed+Shannon+Fisher+Simpson,data=shannon_data%>% dplyr::filter(GROUP_SAMPLE != "AE-2"),method=3)
moonBook::mytable(COPD_STAGE~Observed+Shannon+Fisher,data=shannon_data %>% dplyr::filter(GROUP_SAMPLE != "AE-2"),method=3)
moonBook::mytable(COPD_stage_Group~Observed+Shannon+Fisher+Simpson,data=shannon_data%>% dplyr::filter(GROUP_SAMPLE != "AE-2"),method=3)

moonBook::mytable(Season_2~Observed+Shannon+Fisher,data=shannon_data %>% dplyr::filter(GROUP_SAMPLE != "AE-2"),method=3)

moonBook::mytable(Season_3~Observed+Shannon+Fisher,data=shannon_data %>% dplyr::filter(GROUP_SAMPLE != "AE-2"),method=3)

moonBook::mytable(Season~Observed+Shannon+Fisher,data=shannon_data %>% dplyr::filter(GROUP_SAMPLE != "AE-2",Season%in%c("summer","winter")), method=3)

moonBook::mytable(Antibiotics_Use_2~Observed+Shannon+Fisher,data=shannon_data %>% dplyr::filter(GROUP_SAMPLE == "AE-2"), method=3)

moonBook::mytable(Smoking_Hx_2~Observed+Shannon+Fisher,data=shannon_data %>% dplyr::filter(GROUP_SAMPLE != "AE-2"),method=3)


moonBook::mytable(Antibiotics_Use_2~Observed+Shannon+Fisher,data=shannon_data %>% dplyr::filter(GROUP_SAMPLE == "AE-2"),method=3)
moonBook::mytable(Steroids_Use_2~Observed+Shannon+Fisher,data=shannon_data %>% dplyr::filter(GROUP_SAMPLE != "AE-2"),method=3)

Shannon_plot_Anti <- ggpaired(shannon_data %>% dplyr::filter(GROUP_SAMPLE == "AE-2"), 
  x = "Antibiotics_Use_2", 
  y = "Shannon",
  id = "PT_NO",
  color = "Antibiotics_Use_2",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  # palette = c("#E7B800", "#FC4E07"),
  add="jitter"
  # ,facet.by = "Antibiotics_Use_2"
) + stat_compare_means(comparisons = list(c("Antibiotics non-use", "Antibiotics use")), method = "wilcox.test",paired = F) +ylim(0,6)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Shannon's index")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")













Shannon_plot_C <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="SD",!(연구등록번호 %in% EX_PT)) %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Shannon",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#E7B800", "#FC4E07"),add="jitter"
  ,facet.by = "Antibiotics_Use_2"
) + stat_compare_means(comparisons = list(c("AE-1", "AE-2")), method = "wilcox.test",paired = T) +ylim(0,6)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Shannon's index")+
  theme(
        # axis.title.y = element_blank(),  
        # axis.text.y = element_blank(),  
        # axis.title.x = element_blank(),  
        # axis.text.x = element_blank(),   
        # axis.ticks.x = element_blank(),
        # panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        legend.position = "none")


ASV_plot_C <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="SD",!(연구등록번호 %in% EX_PT)) %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Observed",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#E7B800", "#FC4E07"),add="jitter"
  ,facet.by = "Antibiotics_Use_2"
) + stat_compare_means(comparisons = list(c("AE-1", "AE-2")), method = "wilcox.test",paired = T) +
  ylim(0,200)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Observed features")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
        legend.position = "none")


Fisher_plot_C <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="SD",!(연구등록번호 %in% EX_PT)) %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Fisher",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#E7B800", "#FC4E07"),add="jitter"
  ,facet.by = "Antibiotics_Use_2"
) + stat_compare_means(comparisons = list(c("AE-1", "AE-2")), method = "wilcox.test",paired = T) +
  ylim(0,25)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Fisher's alpha",color="Group")+
  theme(
    # axis.title.y = element_blank(),  
        # axis.text.y = element_blank(),  
        # axis.title.x = element_blank(),  
        # axis.text.x = element_blank(),   
        # axis.ticks.x = element_blank(),
        # panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        legend.position = "right")

ASV_plot_C+Shannon_plot_C+Fisher_plot_C




Shannon_plot_D <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Shannon",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "AgeGroup"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test"
                       # ,paired = T
                       ) +ylim(0,6)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Shannon's index")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")



ASV_plot_D <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Observed",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "AgeGroup"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test") +
  ylim(0,250)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Observed features")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")


Fisher_plot_D <-  ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Fisher",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "AgeGroup"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test") +
  ylim(0,35)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Fisher's alpha",color="Group")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "right")






Shannon_plot_E <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Shannon",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  # line.color = "gray",
  # line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "Season_2"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test"
                       # ,paired = T
) +ylim(0,6)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Shannon's index")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")



ASV_plot_E <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Observed",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "Season_2"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test") +
  ylim(0,250)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Observed features")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")


Fisher_plot_E <-  ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Fisher",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "Season_2"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test") +
  ylim(0,35)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Fisher's alpha",color="Group")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "right")





Shannon_plot_F <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Shannon",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  # line.color = "gray",
  # line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "Season_2"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test"
                       # ,paired = T
) +ylim(0,6)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Shannon's index")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")



ASV_plot_F <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Observed",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "Season_2"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test") +
  ylim(0,250)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Observed features")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")


Fisher_plot_F <-  ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Fisher",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "Season_2"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test") +
  ylim(0,35)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Fisher's alpha",color="Group")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "right")





Shannon_plot_G <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Shannon",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  # line.color = "gray",
  # line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test"
                       # ,paired = T
) +ylim(0,6)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Shannon's index")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")



ASV_plot_G <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Observed",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test") +
  ylim(0,250)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Observed features")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")


Fisher_plot_G <-  ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Fisher",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#00AFBB", "#E7B800"),add="jitter"
  ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("SD", "AE-1")), method = "wilcox.test") +
  ylim(0,35)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Fisher's alpha",color="Group")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "right")



Shannon_plot_H <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "AgeGroup", 
  y = "Shannon",
  id = "PT_NO",
  color = "AgeGroup",
  # line.color = "gray",
  # line.size = 0.4,
  # show.points = T,
  palette = c("#E41A1C", "#377EB8"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("Age < 73", "Age ≥ 73")), method = "wilcox.test"
                       ,paired = F
) +ylim(0,6)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Shannon's index")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")



ASV_plot_H <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "AgeGroup", 
  y = "Observed",
  id = "PT_NO",
  color = "AgeGroup",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#E41A1C", "#377EB8"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("Age < 73", "Age ≥ 73")), method = "wilcox.test") +
  ylim(0,250)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Observed features")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")


Fisher_plot_H <-  ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "AgeGroup", 
  y = "Fisher",
  id = "PT_NO",
  color = "AgeGroup",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#E41A1C", "#377EB8"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("Age < 73", "Age ≥ 73")), method = "wilcox.test") +
  ylim(0,35)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Fisher's alpha",color="Group")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "right")






Shannon_plot_I <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "Smoking_Hx_2", 
  y = "Shannon",
  id = "PT_NO",
  color = "Smoking_Hx_2",
  # line.color = "gray",
  # line.size = 0.4,
  # show.points = T,
  palette = c("#999999", "#984EA3"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("Non-smoker", "Current-smoker")), method = "wilcox.test"
                       ,paired = F
) +ylim(0,6)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Shannon's index")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")



ASV_plot_I <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "Smoking_Hx_2", 
  y = "Observed",
  id = "PT_NO",
  color = "Smoking_Hx_2",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#999999", "#984EA3"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("Non-smoker", "Current-smoker")), method = "wilcox.test") +
  ylim(0,250)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Observed features")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")


Fisher_plot_I <-  ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "Smoking_Hx_2", 
  y = "Fisher",
  id = "PT_NO",
  color = "Smoking_Hx_2",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#999999", "#984EA3"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("Non-smoker", "Current-smoker")), method = "wilcox.test") +
  ylim(0,35)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Fisher's alpha",color="Group")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "right")



Shannon_plot_J <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE=="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "Antibiotics_Use_2", 
  y = "Shannon",
  id = "PT_NO",
  color = "Antibiotics_Use_2",
  # line.color = "gray",
  # line.size = 0.4,
  # show.points = T,
  palette = c(  "#8C564B","#C5B0D5"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("Antibiotics non-use", "Antibiotics use")), method = "wilcox.test"
                       ,paired = F
) +ylim(0,6)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Shannon's index")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")



ASV_plot_J <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE=="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "Antibiotics_Use_2", 
  y = "Observed",
  id = "PT_NO",
  color = "Antibiotics_Use_2",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#8C564B","#C5B0D5"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("Antibiotics non-use", "Antibiotics use")), method = "wilcox.test") +
  ylim(0,250)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Observed features")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")


Fisher_plot_J <-  ggpaired(shannon_data %>% filter(GROUP_SAMPLE=="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "Antibiotics_Use_2", 
  y = "Fisher",
  id = "PT_NO",
  color = "Antibiotics_Use_2",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#8C564B","#C5B0D5"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("Antibiotics non-use", "Antibiotics use")), method = "wilcox.test") +
  ylim(0,35)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Fisher's alpha",color="Group")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "right")








Shannon_plot_K <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "Season_2", 
  y = "Shannon",
  id = "PT_NO",
  color = "Season_2",
  # line.color = "gray",
  # line.size = 0.4,
  show.points = T,
  palette = c( "#66C2A5", "#FC8D62"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("S/S", "F/W")), method = "wilcox.test") +ylim(0,6)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Shannon's index")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")



ASV_plot_K <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "Season_2", 
  y = "Observed",
  id = "PT_NO",
  color = "Season_2",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#66C2A5", "#FC8D62"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("S/S", "F/W")), method = "wilcox.test") +
  ylim(0,250)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Observed features")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")


Fisher_plot_K <-  ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "Season_2", 
  y = "Fisher",
  id = "PT_NO",
  color = "Season_2",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#66C2A5", "#FC8D62"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("S/S", "F/W")), method = "wilcox.test") +
  ylim(0,35)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Fisher's alpha",color="Group")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "right")






Shannon_plot_L <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "COPD_stage_Group", 
  y = "Shannon",
  id = "PT_NO",
  color = "COPD_stage_Group",
  # line.color = "gray",
  # line.size = 0.4,
  show.points = T,
  palette = c(  "#BEBADA", "#FB8072"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("1-2", "3-4")), method = "wilcox.test") +ylim(0,6)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Shannon's index")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")



ASV_plot_L <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "COPD_stage_Group", 
  y = "Observed",
  id = "PT_NO",
  color = "COPD_stage_Group",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c( "#BEBADA", "#FB8072"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("1-2", "3-4")), method = "wilcox.test") +
  ylim(0,250)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Observed features")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")


Fisher_plot_L <-  ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="AE-2") %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "COPD_stage_Group", 
  y = "Fisher",
  id = "PT_NO",
  color = "COPD_stage_Group",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c( "#BEBADA", "#FB8072"),add="jitter"
  # ,facet.by = "COPD_stage_Group"
) + stat_compare_means(comparisons = list(c("1-2", "3-4")), method = "wilcox.test") +
  ylim(0,35)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Fisher's alpha",color="Group")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "right")



jpeg("./Figure_ALPHA_ANALYSIS_.jpeg",width = 8000,height = 6000,res=600)
(ASV_plot_H+Shannon_plot_H+Fisher_plot_H+ASV_plot_J+Shannon_plot_J+Fisher_plot_J+ASV_plot_I+Shannon_plot_I+Fisher_plot_I+ASV_plot_K+Shannon_plot_K+Fisher_plot_K+ASV_plot_L+Shannon_plot_L+Fisher_plot_L)+plot_layout(ncol = 6)
dev.off()







ASV_plot_C <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="SD",!(연구등록번호 %in% EX_PT)) %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Observed",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#E7B800", "#FC4E07"),add="jitter"
  ,facet.by = "Antibiotics_Use_2"
) + stat_compare_means(comparisons = list(c("AE-1", "AE-2")), method = "wilcox.test",paired = T) +
  ylim(0,200)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Observed features")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "none")


Fisher_plot_C <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="SD",!(연구등록번호 %in% EX_PT)) %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Fisher",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#E7B800", "#FC4E07"),add="jitter"
  ,facet.by = "Antibiotics_Use_2"
) + stat_compare_means(comparisons = list(c("AE-1", "AE-2")), method = "wilcox.test",paired = T) +
  ylim(0,25)+
  xlab(NULL) +
  ylab(NULL)  +
  labs(subtitle = "Fisher's alpha",color="Group")+
  theme(
    # axis.title.y = element_blank(),  
    # axis.text.y = element_blank(),  
    # axis.title.x = element_blank(),  
    # axis.text.x = element_blank(),   
    # axis.ticks.x = element_blank(),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "right")






my_comparisons <- list( c("SD", "AE-1"), c("SD", "AE-2"))

ASVs_plot_A <- ggboxplot(shannon_data #%>% filter(GROUP_SAMPLE!="AE-2")
                       , x = "GROUP_SAMPLE", y = "Observed",
                       color = "GROUP_SAMPLE", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                   add = "jitter")+ stat_compare_means(comparisons = my_comparisons, method = "wilcox")+ # Add pairwise comparisons p-value
  ylim(0,300)+
  # stat_compare_means(label.y =300)  +
  labs(subtitle = "Cross-sectional")+
  xlab("Group") +ylab("Observed features") +
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        legend.position = "none")


ASVs_plot_B <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="SD",!(연구등록번호 %in% EX_PT)) %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
                        x = "GROUP_SAMPLE", 
                        y = "Observed",
                        id = "PT_NO",
                        color = "GROUP_SAMPLE",
                        line.color = "gray",
                        line.size = 0.4,
                        show.points = T,
                        palette = c("#E7B800", "#FC4E07"),add="jitter"
                        # ,facet.by = "GROUP_SAMPLE"
                        ) + stat_compare_means(comparisons = list(c("AE-1", "AE-2")), method = "wilcox.test",paired = T) +ylim(0,300)+
  xlab("Group") +ylab("Observed features")  +
  labs(subtitle = "Longitudinal")+
  theme(axis.title.y = element_blank(),  
        axis.text.y = element_blank(),  
        axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        legend.position = "none")

ASVs_plot_A+ASVs_plot_B+plot_layout()


Shannon_plot_A <- ggboxplot(shannon_data #%>% filter(GROUP_SAMPLE!="AE-2")
                         , x = "GROUP_SAMPLE", y = "Shannon",
                         color = "GROUP_SAMPLE", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                         add = "jitter")+ stat_compare_means(comparisons = my_comparisons, method = "wilcox")+ # Add pairwise comparisons p-value
  ylim(0,6)+
  # stat_compare_means(label.y =300)  +
  labs(subtitle = "Cross-sectional")+
  xlab("Group") +ylab("Shannon's index") +
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        legend.position = "none")


Shannon_plot_B <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="SD",!(연구등록번호 %in% EX_PT)) %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Shannon",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#E7B800", "#FC4E07"),add="jitter"
  # ,facet.by = "GROUP_SAMPLE"
) + stat_compare_means(comparisons = list(c("AE-1", "AE-2")), method = "wilcox.test",paired = T) +ylim(0,6)+
  xlab("Group") +ylab("Shannon's index")  +
  labs(subtitle = "Longitudinal")+
  theme(axis.title.y = element_blank(),  
        axis.text.y = element_blank(),  
        axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")


Shannon_plot_A+Shannon_plot_B 




Fisher_plot_A <- ggboxplot(shannon_data #%>% filter(GROUP_SAMPLE!="AE-2")
                            , x = "GROUP_SAMPLE", y = "Fisher",
                            color = "GROUP_SAMPLE", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                            add = "jitter")+ stat_compare_means(comparisons = my_comparisons, method = "wilcox")+ # Add pairwise comparisons p-value
  ylim(0,45)+
  # stat_compare_means(label.y =300)  +
  labs(subtitle = "Cross-sectional")+
  xlab("Group") +ylab("Fisher's alpha") +
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        legend.position = "none")


Fisher_plot_B <- ggpaired(shannon_data %>% filter(GROUP_SAMPLE!="SD",!(연구등록번호 %in% EX_PT)) %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "GROUP_SAMPLE", 
  y = "Fisher",
  id = "PT_NO",
  color = "GROUP_SAMPLE",
  line.color = "gray",
  line.size = 0.4,
  show.points = T,
  palette = c("#E7B800", "#FC4E07"),add="jitter"
  # ,facet.by = "GROUP_SAMPLE"
) + stat_compare_means(comparisons = list(c("AE-1", "AE-2")), method = "wilcox.test",paired = T) +ylim(0,45)+
  xlab("Group") +ylab("Fisher's alpha")  +
  labs(subtitle = "Longitudinal")+
  theme(axis.title.y = element_blank(),  
        axis.text.y = element_blank(),  
        axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")


Fisher_plot_A+Fisher_plot_B


jpeg("./Figure_1_1.jpeg",width = 6000,height = 2250,res=600)
(ASVs_plot_A+ASVs_plot_B+Shannon_plot_A+Shannon_plot_B+Fisher_plot_A+Fisher_plot_B)+plot_layout(nrow = 1)
dev.off()

jpeg("./Figure_ST_ANALYSIS_ANTI.jpeg",width = 6000,height = 3000,res=600)
(ASV_plot_C+Shannon_plot_C+Fisher_plot_C)+plot_layout(nrow = 1)
dev.off()

jpeg("./Figure_ST_ANALYSIS_AGE.jpeg",width = 6000,height = 3000,res=600)
(ASV_plot_D+Shannon_plot_D+Fisher_plot_D)+plot_layout(nrow = 1)
dev.off()

jpeg("./Figure_ST_ANALYSIS_SMK.jpeg",width = 6000,height = 3000,res=600)
(ASV_plot_E+Shannon_plot_E+Fisher_plot_E)+plot_layout(nrow = 1)
dev.off()

jpeg("./Figure_ST_ANALYSIS_SEASON.jpeg",width = 6000,height = 3000,res=600)
(ASV_plot_F+Shannon_plot_F+Fisher_plot_F)+plot_layout(nrow = 1)
dev.off()

jpeg("./Figure_ST_ANALYSIS_GRADE.jpeg",width = 6000,height = 3000,res=600)
(ASV_plot_G+Shannon_plot_G+Fisher_plot_G)+plot_layout(nrow = 1)
dev.off()

pair_data<- shannon_data %>% filter(GROUP_SAMPLE!="SD",!(연구등록번호 %in% EX_PT)) %>% mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO)))

moonBook::mytable(GROUP_SAMPLE~Observed+Shannon+Fisher,data=pair_data,method=2)



###################################
# Paired_BC
library(phangorn)
library(MicrobiomeStat)
library(microbiome)
library(ape)
library(mia)
library(scater)

library("ape")
random_tree = rtree(ntaxa(copd_micro_phyloseq), rooted=TRUE, tip.label=taxa_names(copd_micro_phyloseq))
plot(random_tree)
copd_micro_phyloseq_test <- merge_phyloseq(copd_micro_phyloseq, sam_data(copd_micro_phyloseq), random_tree)
copd_micro_phyloseq_test_2 <- merge_phyloseq(copd_micro_phyloseq, sam_data(copd_micro_phyloseq), random_tree)

copd_micro_phyloseq_AE <- subset_samples(copd_micro_phyloseq, GROUP_SAMPLE !="SD")
copd_micro_phyloseq_AE <- subset_samples(copd_micro_phyloseq_AE, !(연구등록번호 %in% EX_PT))
tse <- makeTreeSEFromPhyloseq(copd_micro_phyloseq_AE)
# assayNames(tse) <- "counts"
tse <- transformAssay(tse, method = "relabundance")
tse <- transformAssay(tse, assay.type = "relabundance",method = "clr", pseudocount = 1)
tse <- runMDS(tse, FUN = vegan::vegdist, method = "bray", name = "PCoA_BC", exprs_values = "counts")




# plot



dis_jac <- vegan::vegdist(t(assays(tse)$counts), method = "jaccard")

# principal coordinate analysis
jaccard_pcoa <- ecodist::pco(dis_jac)

# a data frame from principal coordinates and groupng variable
jaccard_pcoa_df <- data.frame(pcoa1 = jaccard_pcoa$vectors[,1], 
                              pcoa2 = jaccard_pcoa$vectors[,2])
metadata <- colData(tse)[, c("PT_NO", "GROUP_SAMPLE","PTCODE")] %>% as.data.frame() %>% 
  mutate(PT_NO=case_when(
    substr(PTCODE,1,5)=="PT-16"~"6371768_1",
    substr(PTCODE,1,5)=="PT-17"~"6371768_2",
    TRUE~as.character(PT_NO)))

# Combine the coordinates and metadata into a single data frame
plot_data_jac <- cbind(jaccard_pcoa_df, metadata)

# Extract eigenvalues for axis labels
ae1_data_jac <- plot_data_jac %>%  as.data.frame() %>% filter(GROUP_SAMPLE == "AE-1")
ae2_data_jac <- plot_data_jac %>% filter(GROUP_SAMPLE == "AE-2")

jaccard_plot <- ggplot(data = plot_data_jac, aes(x=pcoa1, y=pcoa2, color = GROUP_SAMPLE)) +
  geom_point(size = 3) +  # Plot points
  geom_segment(data = ae1_data_jac,
               aes(x = pcoa1, y = pcoa2, 
                   xend = ae2_data_jac$pcoa1, yend = ae2_data_jac$pcoa2, 
                   group = PT_NO),
               arrow = arrow(type = "closed", length = unit(0.05, "inches")),  # Add arrow
               color = "grey", linetype = "dashed", size = 0.5) +  # Dashed lines with arrows
  labs(x = paste("PCoA 1 (", round(100 * jaccard_pcoa$values[[1]] / sum(jaccard_pcoa$values), 1), "%", ")", sep = ""),
       y = paste("PCoA 2 (", round(100 * jaccard_pcoa$values[[2]] / sum(jaccard_pcoa$values), 1), "%", ")", sep = "")) +  # Adjust color legend title
  stat_ellipse(aes(color = GROUP_SAMPLE)) + 
  scale_color_manual(values = c( "#E7B800", "#FC4E07"))+ theme_bw()+
  theme(legend.position = "none") 

  

# Use geom_segment with arrows to indicate direction from AE-1 to AE-2
p <- ggplot(plot_data, aes(x = V1, y = V2, color = GROUP_SAMPLE)) +
  geom_point(size = 1) +  # Plot points
  geom_segment(data = ae1_data,
               aes(x = V1, y = V2, 
                   xend = ae2_data$V1, yend = ae2_data$V2, 
                   group = PT_NO),
               arrow = arrow(type = "closed", length = unit(0.05, "inches")),  # Add arrow
               color = "grey", linetype = "dashed", size = 0.5) +  # Dashed lines with arrows
  labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
       y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = ""),
       color = "Group") +  # Adjust color legend title
  stat_ellipse(aes(color = GROUP_SAMPLE)) + 
  scale_color_manual(values = c( "#E7B800", "#FC4E07"))+ theme_bw()+
  theme(legend.position = "none") 


dis_bray <- vegan::vegdist(t(assays(tse)$counts), method = "bray")
dis_jaccard <- vegan::vegdist(t(assays(tse)$counts), method = "jaccard")
b_bray <- vegan::betadisper(dis_bray, colData(tse)$GROUP_SAMPLE)
b_jaccard <- vegan::betadisper(dis_jaccard, colData(tse)$GROUP_SAMPLE)

p_1_bray <- ggpaired(cbind(distance = as.numeric(b_bray$distances),
                      cohort = colData(tse)$GROUP_SAMPLE,
                      PT_NO=colData(tse)$PT_NO,
                      PTCODE=colData(tse)$PTCODE) %>% 
                  as_tibble() %>% 
                  mutate(distance = as.numeric(distance)) %>% 
                  mutate(PT_NO=case_when(
                    substr(PTCODE,1,5)=="PT-16"~"6371768_1",
                    substr(PTCODE,1,5)=="PT-17"~"6371768_2",
                    TRUE~as.character(PT_NO))), 
                x = "cohort", 
                y = "distance",
                id = "PT_NO",
                color = "cohort",
                line.color = "gray",
                line.size = 0.4,
                show.points = T,
                palette = c("#E7B800", "#FC4E07"),add="jitter"
                # ,facet.by = "GROUP_SAMPLE"
) + stat_compare_means(comparisons = list(c("AE-1", "AE-2")), method = "wilcox.test",paired = T)+ labs(x = element_blank(), y = "Bray-Curtis distance", color = "Group")+
  theme(legend.position = "right") + theme_bw()

data_dis<- as.data.frame(cbind(distance = as.numeric(b_bray$distances),
      cohort = colData(tse)$GROUP_SAMPLE,
      PT_NO=colData(tse)$PT_NO,
      PTCODE=colData(tse)$PTCODE))
(data_dis %>% filter(cohort=="AE-1"))
(data_dis %>% filter(cohort=="AE-2"))

median((data_dis %>% filter(cohort=="AE-1"))$distance %>% as.numeric())
quantile((data_dis %>% filter(cohort=="AE-1"))$distance %>% as.numeric(),0.25)
quantile((data_dis %>% filter(cohort=="AE-1"))$distance %>% as.numeric(),0.75)

median((data_dis %>% filter(cohort=="AE-2"))$distance %>% as.numeric())
quantile((data_dis %>% filter(cohort=="AE-2"))$distance %>% as.numeric(),0.25)
quantile((data_dis %>% filter(cohort=="AE-2"))$distance %>% as.numeric(),0.75)


wilcox.test((data_dis %>% filter(cohort=="AE-1"))$distance %>% as.numeric(),(data_dis %>% filter(cohort=="AE-2"))$distance %>% as.numeric(),pair=T)

p_1_jaccard <- ggpaired(cbind(distance = as.numeric(b_jaccard$distances),
                           cohort = colData(tse)$GROUP_SAMPLE,
                           PT_NO=colData(tse)$PT_NO,
                           PTCODE=colData(tse)$PTCODE) %>% 
                       as_tibble() %>% 
                       mutate(distance = as.numeric(distance)) %>% 
                       mutate(PT_NO=case_when(
                         substr(PTCODE,1,5)=="PT-16"~"6371768_1",
                         substr(PTCODE,1,5)=="PT-17"~"6371768_2",
                         TRUE~as.character(PT_NO))), 
                     x = "cohort", 
                     y = "distance",
                     id = "PT_NO",
                     color = "cohort",
                     line.color = "gray",
                     line.size = 0.4,
                     show.points = T,
                     palette = c("#E7B800", "#FC4E07"),add="jitter"
                     # ,facet.by = "GROUP_SAMPLE"
) + stat_compare_means(comparisons = list(c("AE-1", "AE-2")), method = "wilcox.test",paired = T)+ labs(x = element_blank(), y = "Bray-Curtis distance", color = "Group")+
  theme(legend.position = "right") + theme_bw()


pcoa_coords <- as.data.frame(reducedDim(tse, "PCoA_BC"))
metadata <- colData(tse)[, c("PT_NO", "GROUP_SAMPLE","PTCODE")] %>% as.data.frame() %>% 
  mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO)))

# Combine the coordinates and metadata into a single data frame
plot_data <- cbind(pcoa_coords, metadata)

# Extract eigenvalues for axis labels
e <- attr(reducedDim(tse, "PCoA_BC"), "eig")
rel_eig <- e / sum(e[e > 0])
ae1_data <- plot_data %>% filter(GROUP_SAMPLE == "AE-1")
ae2_data <- plot_data %>% filter(GROUP_SAMPLE == "AE-2")

# Use geom_segment with arrows to indicate direction from AE-1 to AE-2
p <- ggplot(plot_data, aes(x = V1, y = V2, color = GROUP_SAMPLE)) +
  geom_point(size = 1) +  # Plot points
  geom_segment(data = ae1_data,
               aes(x = V1, y = V2, 
                   xend = ae2_data$V1, yend = ae2_data$V2, 
                   group = PT_NO),
               arrow = arrow(type = "closed", length = unit(0.05, "inches")),  # Add arrow
               color = "darkgray", linetype = "solid", size = 0.5) +  # Dashed lines with arrows
  labs(x = paste("PCoA 1 [", round(100 * rel_eig[[1]], 1), "%", "]", sep = ""),
       y = paste("PCoA 2 [", round(100 * rel_eig[[2]], 1), "%", "]", sep = ""),
       color = "Group") +  # Adjust color legend title
  stat_ellipse(aes(color = GROUP_SAMPLE)) + 
  scale_color_manual(values = c( "#E7B800", "#FC4E07"))+ theme_bw()+
  theme(legend.position = "none") 


dis <- vegan::vegdist(t(assays(tse)$counts), method = "bray")
b <- vegan::betadisper(dis, colData(tse)$GROUP_SAMPLE)
print(anova(b))
b$vectors

p_1 <- ggpaired(cbind(distance = as.numeric(b$distances),
              cohort = colData(tse)$GROUP_SAMPLE,
              PT_NO=colData(tse)$PT_NO,
              PTCODE=colData(tse)$PTCODE) %>% 
           as_tibble() %>% 
           mutate(distance = as.numeric(distance)) %>% 
           mutate(PT_NO=case_when(
  substr(PTCODE,1,5)=="PT-16"~"6371768_1",
  substr(PTCODE,1,5)=="PT-17"~"6371768_2",
  TRUE~as.character(PT_NO))), 
  x = "cohort", 
  y = "distance",
  id = "PT_NO",
  color = "cohort",
  line.color = "darkgray",
  line.size = 0.4,
  show.points = T,
  palette = c("#E7B800", "#FC4E07"),add="jitter"
  # ,facet.by = "GROUP_SAMPLE"
) + ylim(0.62,0.75)+stat_compare_means(comparisons = list(c("AE-1", "AE-2")), method = "wilcox.test",paired = T)+ labs(x = element_blank(), y = "Distance to centroid", color = "Group")+
  theme_bw() + theme(legend.position = "none") 


p+p_1
















tse <- transformSamples(tse, method = "clr", pseudocount = 1)

tse <- runMDS(tse, FUN = vegan::vegdist, name = "MDS_euclidean",
              method = "euclidean", exprs_values = "clr")

#######################################################################################
#######################################################################################
data_AE.obj <- mStat_convert_phyloseq_to_data_obj(copd_micro_phyloseq_AE)
data_AE.obj$meta.dat$PT_NO <- ifelse(substr(data_AE.obj$meta.dat$PTCODE,1,5)=="PT-16","6371768_1",ifelse(substr(data_AE.obj$meta.dat$PTCODE,1,5)=="PT-17","6371768_2",data_AE.obj$meta.dat$연구등록번호))
data_AE.obj$meta.dat$GROUP_SAMPLE <- factor(data_AE.obj$meta.dat$GROUP_SAMPLE ,levels = c("AE-1","AE-2"))
data_AE.obj$meta.dat$TIME_POINT <- ifelse(data_AE.obj$meta.dat$GROUP_SAMPLE == "AE-1",1,2)
data_AE.obj$meta.dat$Season_2 <- factor(data_AE.obj$meta.dat$Season_2 ,levels = c("Non-smoker","Current-smoker"))


# Generate taxa change test pair
test.list <- generate_taxa_change_test_pair(
  data.obj = data_AE.obj,
  subject.var = "PT_NO",
  time.var = "TIME_POINT",
  group.var = "TIME_POINT",
  # adj.vars = c("SEX"),
  change.base = "1",
  feature.change.func = "log fold change",
  feature.level = "Genus",
  prev.filter = 0.01,
  abund.filter = 0.001,
  feature.dat.type = "count"
)

plot.list <- generate_taxa_volcano_single(
  data.obj = data_AE.obj,
  group.var = "GROUP_SAMPLE",
  test.list = test.list,
  feature.sig.level = 0.1,
  feature.mt.method = "none"
)





# Generate the volcano plot
plot.list <- generate_taxa_volcano_single(
  data.obj = peerj32.obj,
  group.var = "group",
  test.list = test.list,
  feature.sig.level = 0.1,
  feature.mt.method = "none"
)


####################################################
library(Maaslin2)
Species_phyloseq <- tax_glom(copd_micro_phyloseq_AE, taxrank = "Species")
genus_phyloseq <- tax_glom(copd_micro_phyloseq_AE, taxrank = "Genus")
genus_data <- as(otu_table(genus_phyloseq), "matrix") 
genus_names <- as.character(tax_table(genus_phyloseq)[, "Genus"])
genus_names <- make.unique(genus_names)  # Ensure unique names to avoid conflicts
# Update the taxa names in the phyloseq object to use genus names
row.names(genus_data) <- genus_names

Species_phyloseq <- tax_glom(copd_micro_phyloseq_AE, taxrank = "Species")
Species_data <- as(otu_table(Species_phyloseq), "matrix") 
Species_names <- as.character(tax_table(Species_phyloseq)[, "Species"])
Species_names <- make.unique(Species_names)  # Ensure unique names to avoid conflicts
# Update the taxa names in the phyloseq object to use Species names
row.names(Species_data) <- Species_names


sample_data_df <- as(sample_data(copd_micro_phyloseq_AE), "data.frame")
sample_data_df$PT_NO <- ifelse(substr(sample_data_df$PTCODE,1,5)=="PT-16","6371768_1",ifelse(substr(sample_data_df$PTCODE,1,5)=="PT-17","6371768_2",sample_data_df$연구등록번호))
# Ensure that ASV data and metadata are formatted correctly for Maaslin2
# Transpose if ASV data is in rows and samples in columns
# if (taxa_are_rows(copd_micro_phyloseq_AE)) {
#   asv_data <- t(asv_data)
# }

# Add sample IDs as a column in metadata
sample_data_df$SampleID <- rownames(sample_data_df)

# Add SampleID as a column in the ASV data
colnames(Species_data) <- rownames(sample_data_df)
colnames(genus_data) <- rownames(sample_data_df)

# Create results directory

dir.create("maaslin2_results", showWarnings = FALSE)
# Run Maaslin2
fit_data <- Maaslin2(
  input_data = Species_data,
  input_metadata = sample_data_df,
  output = "maaslin2_results",
  fixed_effects = "GROUP_SAMPLE",
  random_effects = "PT_NO",
  reference = c("GROUP_SAMPLE", "AE-1"),
  normalization = "CLR", # Centered Log-Ratio normalization if ASV counts are compositional
  analysis_method = "LM" # Linear Model suited for paired comparisons
)

fit_data_gen <- Maaslin2(
  input_data = genus_data,
  input_metadata = sample_data_df,
  output = "maaslin2_results",
  fixed_effects = "GROUP_SAMPLE",
  random_effects = "PT_NO",
  reference = c("GROUP_SAMPLE", "AE-1"),
  normalization = "CLR", # Centered Log-Ratio normalization if ASV counts are compositional
  analysis_method = "LM" # Linear Model suited for paired comparisons
)
# Load the Maaslin2 significant results file
results <- read.delim("./maaslin2_results_GEN/significant_results.tsv")
results_2 <- read.delim("./maaslin2_results_SP/significant_results.tsv")
# Example of a customized bar plot for effect size
PLOT_1 <- ggplot(results, aes(x = reorder(feature, -coef), y = coef, fill = coef > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Genus", y = "Coefficient (reference: AE-1)") + ylim(-2.5,2.5)+
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("#E7B800", "#FC4E07"), name = "Enriched group", labels = c("AE-1", "AE-2"))

PLOT_2 <- ggplot(results_2, aes(x = reorder(feature, -coef), y = coef, fill = coef > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Species", y = "Coefficient (reference: AE-1)") + ylim(-2.5,2.5)+
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#E7B800", "#FC4E07"), name = "Enriched group", labels = c("AE-1", "AE-2"))

(A+PLOT_1)/PLOT_2

jpeg("./MELA_genus_sp.jpeg", height = 2500, width = 4000, res = 600)
PLOT_1/PLOT_2
dev.off()

Species_phyloseq <- tax_glom(copd_micro_phyloseq_AE, taxrank = "Species")
genus_phyloseq <- tax_glom(copd_micro_phyloseq_AE, taxrank = "Genus")
genus_data <- as(otu_table(genus_phyloseq), "matrix") 
genus_names <- as.character(tax_table(genus_phyloseq)[, "Genus"])

# Subset to a specific genus (replace 'YourGenus' with the actual genus name)
genus_phyloseq_abd <- tax_glom(copd_micro_phyloseq_AE, taxrank = "Genus")
genus_phyloseq_abd <- subset_taxa(genus_phyloseq_abd, Genus %in% c("Peptostreptococcus", "Prevotellamassilia","Oribacterium","Filifactor"))

# Extract abundance and metadata
# Make sure to replace 'Genus' with the specific taxonomic rank of interest if necessary
genus_abd <- as.data.frame(otu_table(genus_phyloseq_abd))
sample_data_df <- as(sample_data(genus_phyloseq_abd), "data.frame")
genus_names <- as.character(tax_table(genus_phyloseq_abd)[, "Genus"])
row.names(genus_abd) <- genus_names

sample_data_df$PT_NO <- ifelse(substr(sample_data_df$PTCODE,1,5)=="PT-16","6371768_1",ifelse(substr(sample_data_df$PTCODE,1,5)=="PT-17","6371768_2",sample_data_df$연구등록번호))
sample_data_df$PT_NO <- ifelse(substr(sample_data_df$PTCODE,1,5)=="PT-16","6371768_1",ifelse(substr(sample_data_df$PTCODE,1,5)=="PT-17","6371768_2",sample_data_df$연구등록번호))


genus_abd_long <- genus_abd %>%
  rownames_to_column("Genus") %>%
  mutate(Genus=factor(Genus,levels=c("Peptostreptococcus", "Prevotellamassilia","Oribacterium","Filifactor"))) %>% 
  pivot_longer(cols = -Genus, names_to = "PTCODE", values_to = "Abundance") %>%
  dplyr::left_join(sample_data_df, by = "PTCODE")

A <- ggpaired(
  genus_abd_long,
  x = "GROUP_SAMPLE",
  y = "Abundance",
  color = "GROUP_SAMPLE",
  palette = c("#E7B800", "#FC4E07"),
  line.color = "gray", # Optional: adds lines connecting paired points
  id = "PT_NO", # Column identifying paired samples, such as subject ID
  add = "jitter",
  add.params = list(size = 1, alpha = 0.7)
) +
  facet_wrap(~ Genus, ncol = 2) + # Display each genus in a single column
  labs(x = "Group", y = "Abundance", color = "Group") +
  theme_pubr() +
  theme(legend.position = "none")+
  coord_flip() 
# +
  # stat_compare_means(
  #   method = "t.test",
  #   paired = TRUE,
  #   label = "p.format",
  #   # , # or "p.format" for actual p-value
  #   # label.y = "max" # Adjust y-position of p-value as needed,
  #   comparisons = list(c("AE-1", "AE-2"))
  # )
A+PLOT_1

Species_phyloseq <- tax_glom(copd_micro_phyloseq_AE, taxrank = "Species")
Species_data <- as(otu_table(Species_phyloseq), "matrix") 
Species_names <- as.character(tax_table(Species_phyloseq)[, "Species"])
Species_names <- make.unique(Species_names)  # Ensure unique names to avoid conflicts
# Update the taxa names in the phyloseq object to use Species names
row.names(Species_data) <- Species_names

# Subset to a specific Species (replace 'YourSpecies' with the actual Species name)
Species_phyloseq_abd <- tax_glom(copd_micro_phyloseq_AE, taxrank = "Species")
Species_phyloseq_abd <- subset_taxa(Species_phyloseq_abd, Species %in% c("Leptotrichia massiliensis", "Actinomyces oris","Leptotrichia trevisanii","Capnocytophaga sputigena"))

# Extract abundance and metadata
# Make sure to replace 'Species' with the specific taxonomic rank of interest if necessary
Species_abd <- as.data.frame(otu_table(Species_phyloseq_abd))
sample_data_df <- as(sample_data(Species_phyloseq_abd), "data.frame")
Species_names <- as.character(tax_table(Species_phyloseq_abd)[, "Species"])
row.names(Species_abd) <- Species_names

sample_data_df$PT_NO <- ifelse(substr(sample_data_df$PTCODE,1,5)=="PT-16","6371768_1",ifelse(substr(sample_data_df$PTCODE,1,5)=="PT-17","6371768_2",sample_data_df$연구등록번호))
sample_data_df$PT_NO <- ifelse(substr(sample_data_df$PTCODE,1,5)=="PT-16","6371768_1",ifelse(substr(sample_data_df$PTCODE,1,5)=="PT-17","6371768_2",sample_data_df$연구등록번호))


Species_abd_long <- Species_abd %>%
  rownames_to_column("Species") %>%
  mutate(Species=case_when(Species=="Leptotrichia massiliensis"~"L.massiliensis",
                           Species=="Actinomyces oris"~"A.oris",
                           Species=="Leptotrichia trevisanii"~"L.trevisanii",
                           Species=="Capnocytophaga sputigena"~"C.sputigena",
                           TRUE~"none"),
         Species=factor(Species,levels=c("L.massiliensis", "A.oris","L.trevisanii","C.sputigena"))) %>% 
  pivot_longer(cols = -Species, names_to = "PTCODE", values_to = "Abundance") %>%
  dplyr::left_join(sample_data_df, by = "PTCODE")

B <- ggpaired(
  Species_abd_long,
  x = "GROUP_SAMPLE",
  y = "Abundance",
  color = "GROUP_SAMPLE",
  palette = c("#E7B800", "#FC4E07"),
  line.color = "gray", # Optional: adds lines connecting paired points
  id = "PT_NO", # Column identifying paired samples, such as subject ID
  add = "jitter",
  add.params = list(size = 1, alpha = 0.7)
) +
  facet_wrap(~ Species, ncol = 2
             # ,strip.position = "left"
             ) + # Display each Species in a single column
  labs(x = "Group", y = "Abundance", color = "Group") +
  theme_pubr() +
  theme(legend.position = "none")+
  coord_flip() 


(A+PLOT_1)/(B+PLOT_2)

jpeg("./LMOD_genus_sp_REC.jpeg", height = 3500, width = 7000, res = 600)
A+PLOT_1+B+PLOT_2+patchwork::plot_layout(nrow=2)
dev.off()



print(A)


B <-  plot_abundance(mm_test_Species, group = "GROUP_SAMPLE")+ 
  scale_fill_manual(values = c( "#E7B800", "#00AFBB"))+
  labs(fill = "Group") +
  theme_bw()+
  theme(legend.position = "bottom")




library(aplot)
generate_beta_ordination_single(
  data.obj = data_AE.obj,
  dist.obj = NULL,
  pc.obj = NULL,
  subject.var = "PT_NO",
  time.var = "GROUP_SAMPLE", # Variable representing time points
  t.level = NULL, # Specific time level for subset (if subset desired)
  group.var = "GROUP_SAMPLE",
  strata.var = NULL,
  dist.name = c("BC","Jaccard"),
  base.size = 20,
  theme.choice = "bw",
  custom.theme = NULL,
  palette = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)



generate_beta_ordination_pair(
  data.obj = data_AE.obj,
  dist.obj = NULL,
  time.var = "GROUP_SAMPLE",
  subject.var = "PT_NO",
  group.var = "GROUP_SAMPLE",
  adj.vars = c("GROUP_SAMPLE"),
  change.base = "1",
  dist.name = c('BC', 'Jaccard'),
  base.size = 16,
  theme.choice = "bw",
  custom.theme = NULL,
  palette = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)

generate_beta_pc_boxplot_long(
  data.obj = data_AE.obj,
  dist.obj = NULL,
  pc.obj = NULL,
  subject.var = "PT_NO",
  time.var = "GROUP_SAMPLE",
  t0.level = "AE-1",
  ts.levels = "AE-2",
  # group.var = "GROUP_SAMPLE",
  # strata.var = "sex",
  dist.name = c('BC', 'Jaccard'),
  base.size = 17,
  theme.choice = "bw",
  custom.theme = NULL,
  palette = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)
+ stat_compare_means(method = "wilcox.test")


generate_beta_change_test_pair(
  data.obj = data_AE.obj,
  dist.obj = NULL,
  time.var = "GROUP_SAMPLE",
  subject.var = "PT_NO",
  group.var = "GROUP_SAMPLE",
  adj.vars = c("GROUP_SAMPLE"),
  change.base = "1",
  dist.name = c('BC', 'Jaccard') 
)

bray <- ordinate(
  physeq = copd_micro_phyloseq, #change this to your phyloseq
  method = "PCoA", 
  distance = "bray"
)
sample_metadata <- data.frame(sample_data(copd_micro_phyloseq_test))

TEXT_UNI <- adonis2(phyloseq::distance(copd_micro_phyloseq_test, method="unifrac") ~ GROUP_SAMPLE,data = sample_metadata)

TEXT_UNI_W <- adonis2(phyloseq::distance(copd_micro_phyloseq_test, method="wunifrac") ~ GROUP_SAMPLE,
        data = sample_metadata)
TEXT_BRAY <- adonis2(phyloseq::distance(copd_micro_phyloseq_test, method="bray") ~ GROUP_SAMPLE,
        data = sample_metadata)

TEXT_BRAY_pair <- pairwise.adonis2(phyloseq::distance(copd_micro_phyloseq_test, method="bray") ~ GROUP_SAMPLE,data = sample_metadata)

TEXT_UNI_pair <- pairwise.adonis2(phyloseq::distance(copd_micro_phyloseq_test, method="unifrac") ~ GROUP_SAMPLE,
                    data = sample_metadata)
TEXT_UNI_W_pair <- pairwise.adonis2(phyloseq::distance(copd_micro_phyloseq_test, method="wunifrac") ~ GROUP_SAMPLE,
                      data = sample_metadata)



# Create the label text

bray_plot <- plot_ordination(
  physeq = copd_micro_phyloseq,  # Phyloseq object
  ordination = bray             # Ordination
) + geom_point(aes(color = GROUP_SAMPLE), size = 1.0) + 
  stat_ellipse(aes(color = GROUP_SAMPLE)) + 
  scale_color_manual(values = c( "#E7B800", "#FC4E07","#00AFBB"))+
  theme_bw()+
  theme(legend.position = "right") +
  labs(x = "PCoA 1 [5.2 %]",
       y = "PCoA 2 [5.0 %]",
       color = "Group")
# +
  # annotate("text", x = -Inf, y = Inf, label = label_text_bray, hjust = -0.1, vjust = 1.1, size = 3, color = "black", parse = TRUE)





library("ape")
random_tree = rtree(ntaxa(copd_micro_phyloseq), rooted=TRUE, tip.label=taxa_names(copd_micro_phyloseq))
plot(random_tree)
copd_micro_phyloseq_test <- merge_phyloseq(copd_micro_phyloseq, sam_data(copd_micro_phyloseq), random_tree)
copd_micro_phyloseq_test_2 <- merge_phyloseq(copd_micro_phyloseq, sam_data(copd_micro_phyloseq), random_tree)



UniFrac(copd_micro_phyloseq_test, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
UniFrac(copd_micro_phyloseq_test_2, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)

#ordination
uni <- ordinate(
  physeq = copd_micro_phyloseq_test_2, 
  method = "PCoA", 
  distance = "wunifrac"
)
#summary(distance)
#distance
uni_unweighted <- ordinate(
  physeq = copd_micro_phyloseq_test, 
  method = "PCoA", 
  distance = "unifrac"
)
#summary(distance)
#distance

#plot
label_text_uni_w <- sprintf("R^2==%s~','~%s", round(TEXT_UNI_W$R2[1],3), ifelse(TEXT_UNI_W$P[1] < 0.001, "P<0.001", sprintf("P==%.3f", round(TEXT_UNI_W$P[1], 3))))
label_text_uni <- sprintf("R^2==%s~','~%s", round(TEXT_UNI$R2[1],3), ifelse(TEXT_UNI$P[1] < 0.001, "P<0.001", sprintf("P==%.3f", round(TEXT_UNI$P[1], 3))))

Uni_weight_plot <- plot_ordination(
  physeq = copd_micro_phyloseq_test_2,                                                          #phyloseq object
  ordination = uni)+geom_point(aes(color = GROUP_SAMPLE), size = 1.0) + 
  stat_ellipse(aes(color = GROUP_SAMPLE)) + 
  scale_color_manual(values = c( "#E7B800", "#FC4E07","#00AFBB"))+
  labs(x = "PCoA 1 [13.0 %]",
       y = "PCoA 2 [8.3 %]",
       color = "Group")+ theme_bw()+
  theme(legend.position = "none") 
# +
  # annotate("text", x = -Inf, y = Inf, label = label_text_uni_w, hjust = -0.1, vjust = 1.1, size = 3, color = "black", parse = TRUE)



Uni_weight_plot_un <- plot_ordination(
  physeq = copd_micro_phyloseq_test,                                                          #phyloseq object
  ordination = uni_unweighted)+geom_point(aes(color = GROUP_SAMPLE), size = 1.0) + 
  stat_ellipse(aes(color = GROUP_SAMPLE)) + 
  scale_color_manual(values = c( "#E7B800", "#FC4E07","#00AFBB")) +
  labs(x = "PCoA 1 [6.3 %]",
       y = "PCoA 2 [3.9 %]",
       color = "Group")+ theme_bw()+
  theme(legend.position = "none") 
# +
# +
  # annotate("text", x = -Inf, y = Inf, label = label_text_uni, hjust = -0.1, vjust = 1.1, size = 3, color = "black", parse = TRUE)





jpeg("./Figure_2_1.jpeg",width = 6000,height = 2250,res=600)
(ASVs_plot+shannon_plot+fisher_plot)+plot_layout()
dev.off()





jpeg("./Figure_2_2.jpeg",width = 6000,height = 3600,res=600)
(Uni_weight_plot_un+Uni_weight_plot+bray_plot+p+p_1)+plot_layout(ncol=3)
dev.off()















sample_data(copd_micro_phyloseq_test)
###########################################
##COPD Stage and beta, alpha diversity
TEXT_UNI_COPD <- adonis2(distance(copd_micro_phyloseq_test, method="unifrac") ~ COPD_STAGE,
                    data = sample_data(copd_micro_phyloseq_test))
TEXT_UNI_W_COPD <- adonis2(distance(copd_micro_phyloseq_test, method="wunifrac") ~ COPD_STAGE,
                      data = metadata)
TEXT_BRAY_COPD <- adonis2(distance(copd_micro_phyloseq_test, method="bray") ~ COPD_STAGE,
                     data = metadata)

TEXT_UNI_COPD_pair <- pairwise.adonis2(distance(copd_micro_phyloseq_test, method="unifrac") ~ COPD_STAGE,
                         data = metadata)
TEXT_UNI_W_COPD_pair <- pairwise.adonis2(distance(copd_micro_phyloseq_test, method="wunifrac") ~ COPD_STAGE,
                           data = metadata)
TEXT_BRAY_COPD_pair <- pairwise.adonis2(distance(copd_micro_phyloseq_test, method="bray") ~ COPD_STAGE,
                          data = metadata)


# Create the label text
label_text_bray_copd <- sprintf("R^2==%s~','~%s", round(TEXT_BRAY_COPD$R2[1],3), ifelse(TEXT_BRAY_COPD$P[1] < 0.001, "P<0.001", sprintf("P==%.3f", round(TEXT_BRAY_COPD$P[1], 3))))

bray_plot_copd <- plot_ordination(
  physeq = copd_micro_phyloseq,  # Phyloseq object
  ordination = bray             # Ordination
) + geom_point(aes(color = COPD_STAGE), size = 1.0) + 
  stat_ellipse(aes(color = COPD_STAGE)) + 
  scale_color_manual(values = c( "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"))+
  theme_bw()+
  theme(legend.position = "none")   +
  annotate("text", x = -Inf, y = Inf, label = label_text_bray_copd, hjust = -0.1, vjust = 1.1, size = 3, color = "black", parse = TRUE)


label_text_uni_copd <- sprintf("R^2==%s~','~%s", round(TEXT_UNI_COPD$R2[1],3), ifelse(TEXT_UNI_COPD$P[1] < 0.001, "P<0.001", sprintf("P==%.3f", round(TEXT_UNI_COPD$P[1], 3))))

Uni_weight_plot_copd_un <- plot_ordination(
  physeq = copd_micro_phyloseq_test,                                                          #phyloseq object
  ordination = uni_unweighted)+geom_point(aes(color = COPD_STAGE), size = 1.0) + 
  stat_ellipse(aes(color = COPD_STAGE)) + 
  scale_color_manual(values = c( "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")) + 
  labs(color = "COPD stage")  + theme_bw() +
  annotate("text", x = -Inf, y = Inf, label = label_text_uni_copd, hjust = -0.1, vjust = 1.1, size = 3, color = "black", parse = TRUE)+ 
  theme(legend.position = "none")


label_text_uni_w_copd <- sprintf("R^2==%s~','~%s", round(TEXT_UNI_W_COPD$R2[1],3), ifelse(TEXT_UNI_W_COPD$P[1] < 0.001, "P<0.001", sprintf("P==%.3f", round(TEXT_UNI_W_COPD$P[1], 3))))

Uni_weight_plot_copd <- plot_ordination(
  physeq = copd_micro_phyloseq_test,                                                          #phyloseq object
  ordination = uni)+geom_point(aes(color = COPD_STAGE), size = 1.0) + 
  stat_ellipse(aes(color = COPD_STAGE)) + 
  scale_color_manual(values = c( "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")) + 
  labs(color = "Group") + 
  theme(legend.position = "right") + theme_bw() +
  annotate("text", x = -Inf, y = Inf, label = label_text_uni_w_copd, hjust = -0.1, vjust = 1.1, size = 3, color = "black", parse = TRUE)

jpeg("./Figure_COPD_STAGE.jpeg",width = 6000,height = 1800,res=600)
(bray_plot_copd+Uni_weight_plot_copd_un+Uni_weight_plot)+plot_layout()
dev.off()



##################################################

library(vegan)
library(dplyr)
library(corrplot)
library(Hmisc)

cor_results <- rcorr(as.matrix(shannon_data %>% select(Observed, Shannon, Simpson, Fisher, AGE, BMI,COPD_STAGE,Hemoglobin,WBC,  PLT, EOSINOPHIL,CRP, Antibiotics_Hx, Steroids_Hx,z.score.FEV1,z.score.FVC,z.score.FEV1FVC), type = "spearman"))
cor_matrix <- cor_results$r
p_matrix <- cor_results$P

cor_matrix[p_matrix > 0.05] <- NA
cor_matrix[is.na(cor_matrix)] <- 0
jpeg("./CORPLOT_diversity.jpeg", height = 3000, width = 3000, res = 600)
corrplot(cor_matrix, method = "circle", na.label = " ",
         tl.pos = "n" ,
         # ,order = "hclust", addrect = 4
) # NA를 빈 공간으로 표시
dev.off()

ASVs_plot_anti <- ggboxplot(shannon_data #%>% filter(SampleName!="10045778-2")
                       , x = "COPD_STAGE", y = "Observed",
                       color = "COPD_STAGE",
                       add = "jitter")+ 
  # stat_compare_means(comparisons = my_comparisons,method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y =300)  + xlab("Group") +ylab("Observed features") +theme(legend.position = "none")

shannon_plot_anti <- ggboxplot(shannon_data #%>% filter(SampleName!="10045778-2")
                          , x = "COPD_STAGE", y = "Shannon",
                          color = "COPD_STAGE",
                          add = "jitter")+ 
  # stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y =6)  + xlab("Group") +ylab("Shannon's index")+ labs(color = "Group")+theme(legend.position = "none")

merged_data <- merged_data %>% mutate(GROUP_BRAY=case_when(PT_NO%in%c("6924726","1040632","3587399")~"1",TRUE~"0"))

moonBook::mytable(GROUP_BRAY~AGE+Antibiotics_Use+Steroids_Hx+Hemoglobin+WBC+PLT+CRP+AST+ALT+Fibrosis+Pneumonia+Bronchitis+Emphysema+PUD+EOSINOPHIL_300_GROUP+EOSINOPHIL_100_GROUP+EOSINOPHIL+z.score.FEV1FVC+z.score.FEV1+BMI+z.score.FVC,merged_data,method=3)

merged_data$CRP <- as.numeric(merged_data$CRP )

copd_micro_phyloseq_inclusion_abund <- merge_samples(copd_micro_phyloseq_inclusion_abund, "GROUP_SAMPLE")
plot_bar(copd_micro_phyloseq_inclusion_abund, fill = "Genus") + 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

copd_micro_phyloseq_inclusion_abund

LR.1 <- microbiomeMarker::run_sl(
  ps = copd_micro_phyloseq,
  group = "GROUP_ENG", 
  method = "LR")


microbiomeMarker::plot_sl_roc(LR.1, group = "Genus")  
microbiomeMarker::plot_ef_bar(LR.1)

plot_bar(copd_micro_phyloseq_initial, fill = "Phylum")

copd_micro_phyloseq_initial <- subset_samples(copd_micro_phyloseq, GROUP_SAMPLE !="AE-2")
# copd_micro_phyloseq_initial_abund <- filter_taxa(copd_micro_phyloseq_initial, function(x) sum(x > total*0.20) > 0, TRUE)
mm_test_Genus <- run_lefse( copd_micro_phyloseq_initial,
                      wilcoxon_cutoff = 0.05,
                      norm = "CPM", 
                      group = "GROUP_SAMPLE",
                      taxa_rank = "Genus",
                      kw_cutoff = 0.05,
                      multigrp_strat = TRUE,
                      lda_cutoff = 1,
                      bootstrap_n = 50
)

mm_test_Species <- run_lefse( copd_micro_phyloseq_initial,
                            wilcoxon_cutoff = 0.05,
                            norm = "CPM", 
                            group = "GROUP_SAMPLE",
                            taxa_rank = "Species",
                            kw_cutoff = 0.05,
                            multigrp_strat = TRUE,
                            lda_cutoff = 1,
                            bootstrap_n = 50
)
C <- plot_ef_bar(mm_test_Genus)+ 
  scale_fill_manual(values = c( "#00AFBB","#E7B800"))+
  theme(legend.position = "none")

D <- plot_ef_bar(mm_test_Species)+ 
  scale_fill_manual(values = c( "#00AFBB","#E7B800"))+
  theme(legend.position = "bottom")



A <- plot_abundance(mm_test_Genus, group = "GROUP_SAMPLE") + 
  scale_fill_manual(values = c( "#E7B800", "#00AFBB"))+
  labs(fill = "Group") +
  theme_bw()+
  theme(legend.position = "none")

B <-  plot_abundance(mm_test_Species, group = "GROUP_SAMPLE")+ 
  scale_fill_manual(values = c( "#E7B800", "#00AFBB"))+
  labs(fill = "Group") +
  theme_bw()+
  theme(legend.position = "bottom")


mm_test_Genus@marker_table$feature <- factor(mm_test_Genus@marker_table$feature, levels = unique(mm_test_Genus@marker_table$feature))


jpeg("./LDA_genus_sp.jpeg", height = 7000, width = 3500, res = 600)
A+B+C+D+plot_layout(ncol = 1)
dev.off()






























######################################################################

########## Network analysis


# Install NetCoMi
devtools::install_github("stefpeschel/NetCoMi", 
                         dependencies = c("Depends", "Imports", "LinkingTo"),
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))
devtools::install_github("zdk123/SpiecEasi")
devtools::install_github("GraceYoon/SPRING")

library(SpiecEasi)
library(igraph)
library(NetCoMi)
library(phyloseq)

ps.f = prune_samples(sample_sums(copd_micro_phyloseq)>0, copd_micro_phyloseq) 
ps.f = prune_taxa(rowSums(otu_table(ps.f)) > 0, ps.f)  
ps_genus <- tax_glom(ps.f, taxrank = "Genus")
ps_species <- tax_glom(ps.f, taxrank = "Species")
taxtab <- as(tax_table(ps_genus), "matrix")
taxtab_species <- as(tax_table(ps_species), "matrix")


# necomi 형식으로 변환 및 Genus level 설정 
ps_genus_renamed <- NetCoMi::renameTaxa(ps_genus,
                                        pat = "<name>",
                                        substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus")

ps_species_renamed <- NetCoMi::renameTaxa(ps_species,
                                        pat = "<name>",
                                        substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Species")

ps_genus_renamed_sd <- phyloseq::subset_samples(ps_genus_renamed, 
                                                  GROUP_SAMPLE == "SD")
ps_genus_renamed_ae <- phyloseq::subset_samples(ps_genus_renamed, 
                                                GROUP_SAMPLE == "AE-1")

ps_genus_renamed_ae_1 <- phyloseq::subset_samples(ps_genus_renamed, 
                                                (!(연구등록번호 %in% EX_PT))&(GROUP_SAMPLE == "AE-1"))

ps_genus_renamed_rec <- phyloseq::subset_samples(ps_genus_renamed, 
                                                GROUP_SAMPLE == "AE-2")

ps_species_renamed_sd <- phyloseq::subset_samples(ps_species_renamed, 
                                                GROUP_SAMPLE == "SD")
ps_species_renamed_ae <- phyloseq::subset_samples(ps_species_renamed, 
                                                GROUP_SAMPLE == "AE-1")
ps_species_renamed_ae_1 <- phyloseq::subset_samples(ps_species_renamed, 
                                                  (!(연구등록번호 %in% EX_PT))&(GROUP_SAMPLE == "AE-1"))

ps_species_renamed_rec <- phyloseq::subset_samples(ps_species_renamed, 
                                                 GROUP_SAMPLE == "AE-2")


ps_genus.sparcc_sd_ae <- netConstruct(
  data = ps_genus_renamed_sd,
  data2 = ps_genus_renamed_ae,
  
  filtTax = "highestFreq",               # Top을 뽑는 기준
  filtTaxPar = list(highestFreq  = 100), # Top 100
  taxRank = "Genus",                     # Genus level 
  measure = "sparcc",                    # sparcc
  measurePar = list(nlambda=20, rep.num=10),
  
  normMethod = "clr",                    # transformation
  zeroMethod = "none",                   # zero 값 보정 
  sparsMethod = "threshold", 
  adjust = "adaptBH",                    # p-value 값 보정 
  thresh = 0.3,
  
  dissFunc = "signed",
  verbose = 2,
  seed = 42)

ps_genus.sparcc_sd_ae <- netAnalyze(ps_genus.sparcc_sd_ae, 
                              clustMethod = "cluster_fast_greedy",
                              hubPar = "eigenvector",           # Hub node 판별방법
                              normDeg = FALSE)

summary(ps_genus.sparcc_sd_ae, groupNames = c("SD", "AE-1"))


comparison_results_ae1_ae2 <- list()
sample_data(ps_genus_renamed_ae_1)$SampleName <- substr(sample_data(ps_genus_renamed_ae_1)$PTCODE,4,5)
sample_data(ps_genus_renamed_rec)$SampleName <- substr(sample_data(ps_genus_renamed_rec)$PTCODE,4,5)
PT_ID_AE <- unique(sample_data(ps_genus_renamed_ae_1)$SampleName)

# 각 ID에 대해 네트워크 생성 및 비교
for (id in PT_ID_AE) {
  # 각 ID의 전후 데이터를 이용한 네트워크 생성
  ae_1 <- phyloseq::subset_samples(ps_genus_renamed_ae_1
                                   # , SampleName ==id
                                   )
  ae_2 <- phyloseq::subset_samples(ps_genus_renamed_rec
                                   # , SampleName ==id
                                   )
  
  
  
  ps_genus.sparcc_ae_1 <- netConstruct(
    data = ae_1, 
    filtTax = "highestFreq",               # Top을 뽑는 기준
    filtTaxPar = list(highestFreq  = 100), # Top 100
    taxRank = "Genus",                     # Genus level
    measure = "sparcc",                    # sparcc
    measurePar = list(nlambda=20, rep.num=10),
    normMethod = "clr",                    # transformation
    zeroMethod = "none",                   # zero 값 보정
    sparsMethod = "threshold",
    adjust = "adaptBH",                    # p-value 값 보정
    thresh = 0.3,
    dissFunc = "signed",
    verbose = 2,
    seed = 42)

  
  ps_genus.sparcc_ae_2 <- netConstruct(
    data = ae_2, 
    filtTax = "highestFreq",               # Top을 뽑는 기준
    filtTaxPar = list(highestFreq  = 100), # Top 100
    taxRank = "Genus",                     # Genus level
    measure = "sparcc",                    # sparcc
    measurePar = list(nlambda=20, rep.num=10),
    normMethod = "clr",                    # transformation
    zeroMethod = "none",                   # zero 값 보정
    sparsMethod = "threshold",
    adjust = "adaptBH",                    # p-value 값 보정
    thresh = 0.3,
    dissFunc = "signed",
    verbose = 2,
    seed = 42)
  
  
  
  ps_genus.sparcc_ae_all <- netConstruct(
    data = ae_1, 
    data2 = ae_2, 
    filtTax = "highestFreq",               # Top을 뽑는 기준
    filtTaxPar = list(highestFreq  = 100), # Top 100
    taxRank = "Genus",                     # Genus level
    measure = "sparcc",                    # sparcc
    measurePar = list(nlambda=20, rep.num=10),
    normMethod = "clr",                    # transformation
    zeroMethod = "none",                   # zero 값 보정
    sparsMethod = "threshold",
    adjust = "adaptBH",                    # p-value 값 보정
    thresh = 0.3,
    dissFunc = "signed",
    verbose = 2,
    seed = 42)
  
  
  netAnalyze_ae_1 <- netAnalyze(ps_genus.sparcc_ae_1, 
                                clustMethod = "cluster_fast_greedy",
                                hubPar = "eigenvector",           # Hub node 판별방법
                                normDeg = FALSE)
  netAnalyze_ae_2 <- netAnalyze(ps_genus.sparcc_ae_2, 
                                          clustMethod = "cluster_fast_greedy",
                                          hubPar = "eigenvector",           # Hub node 판별방법
                                          normDeg = FALSE)
  
  netAnalyze_ae_1$centralities
  
  
  netAnalyze_ae_all <- netAnalyze(ps_genus.sparcc_ae_all, 
                                clustMethod = "cluster_fast_greedy",
                                hubPar = "eigenvector",           # Hub node 판별방법
                                normDeg = FALSE)
  
  # comparison <- netCompare(
  #   netAnalyze_ae_all,
  #   permTest = T,
  #   nPerm = 2,
  #   jaccQuant = 0.75,
  #   testRand = TRUE,
  #   nPermRand = 1000L,
  #   fileStoreAssoPerm = "/mnt/ICU_STUDY/COPD_MICRO/network_data/ae_assoPerm/",
  #   gcd = TRUE,
  #   # gcdOrb = c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1),
  #   verbose = TRUE,
  #   # nPerm = 1000L,
  #   adjust = "adaptBH",
  #   trueNullMethod = "convest",
  #   cores = 1
  # )
  # 

  
  summary(ps_genus.sparcc_ae_1_ae_2, groupNames = c("AE-1", "AE-2"))
  summary(comparison, groupNames = c("AE-1", "AE-2"))

  # plot(comparison, 
  #      cexLabels = 2,
  #      cexNodes = 0.7, 
  #      cexLegend = 2.5,
  #      cexTitle = 3,
  #      mar = c(3,2,5,15),
  #      # legendGroupnames = c("Mixed diet", "Vegetarian"),
  #      legendPos = c(1.2,1.5),
  #      legendArgs = list(lwd = 4),
  #      fade = FALSE)
  # 
  
  # 결과 저장
  comparison_results_sd_ae[[id]] <- comparison
}

netAnalyze_ae_all %>% View()















ps_genus.sparcc_ae_rec <- netConstruct(
  data = ps_genus_renamed_ae_1,
  data2 = ps_genus_renamed_rec,
  
  filtTax = "highestFreq",               # Top을 뽑는 기준
  filtTaxPar = list(highestFreq  = 100), # Top 100
  taxRank = "Genus",                     # Genus level 
  measure = "sparcc",                    # sparcc
  measurePar = list(nlambda=20, rep.num=10),
  
  normMethod = "clr",                    # transformation
  zeroMethod = "none",                   # zero 값 보정 
  sparsMethod = "threshold", 
  adjust = "adaptBH",                    # p-value 값 보정 
  thresh = 0.3,
  
  dissFunc = "signed",
  verbose = 2,
  seed = 42)

ps_genus.sparcc_ae_rec <- netAnalyze(ps_genus.sparcc_ae_rec, 
                                    clustMethod = "cluster_fast_greedy",
                                    hubPar = "eigenvector",           # Hub node 판별방법
                                    normDeg = FALSE)

summary(ps_genus.sparcc_ae_rec, groupNames = c("AE-1", "AE-2"))


set.seed(42)

ps_species.sparcc_sd_ae <- netConstruct(
  data = ps_species_renamed_sd,
  data2 = ps_species_renamed_ae,
  
  filtTax = "highestFreq",               # Top을 뽑는 기준
  filtTaxPar = list(highestFreq  = 100), # Top 100
  taxRank = "Species",                     # species level 
  measure = "sparcc",                    # sparcc
  measurePar = list(nlambda=20, rep.num=10),
  
  normMethod = "clr",                    # transformation
  zeroMethod = "none",                   # zero 값 보정 
  sparsMethod = "threshold", 
  adjust = "adaptBH",                    # p-value 값 보정 
  thresh = 0.5,
  
  dissFunc = "signed",
  verbose = 2,
  seed = 42)

ps_species.sparcc_sd_ae <- netAnalyze(ps_species.sparcc_sd_ae, 
                                    clustMethod = "cluster_fast_greedy",
                                    hubPar = "eigenvector",           # Hub node 판별방법
                                    normDeg = FALSE)

summary(ps_species.sparcc_sd_ae, groupNames = c("SD", "AE-1"))












ps_species.sparcc_ae_rec <- netConstruct(
  data = ps_species_renamed_ae,
  data2 = ps_species_renamed_rec,
  
  filtTax = "highestFreq",               # Top을 뽑는 기준
  filtTaxPar = list(highestFreq  = 100), # Top 100
  taxRank = "Species",                     # species level 
  measure = "sparcc",                    # sparcc
  measurePar = list(nlambda=20, rep.num=10),
  
  normMethod = "clr",                    # transformation
  zeroMethod = "none",                   # zero 값 보정 
  sparsMethod = "threshold", 
  adjust = "adaptBH",                    # p-value 값 보정 
  thresh = 0.3,
  
  dissFunc = "signed",
  verbose = 2,
  seed = 42)

ps_species.sparcc_ae_rec <- netAnalyze(ps_species.sparcc_ae_rec, 
                                     clustMethod = "cluster_fast_greedy",
                                     hubPar = "eigenvector",           # Hub node 판별방법
                                     normDeg = FALSE)

summary(ps_species.sparcc_ae_rec, groupNames = c("AE-1", "AE-2"))

#####################
# Paired





# ps_props_sparcc <- netCompare(ps_props_sparcc,
#                                   permTest = TRUE,
#                                   nPerm = 1000,
#                                   cores = 6,
#                                   seed = 13075,
#                                   # storeAssoPerm = TRUE,
#                                   # fileStoreAssoPerm = "general/network_data/spring_assoPerm",
#                                   verbose = TRUE)


# ps_props_sparcc <- diffnet(ps_props_sparcc,
#                           diffMethod = "perm",
#                           # fileLoadAssoPerm = "general/network_data/spring_assoPerm",
#                           adjust = "lfdr")




# Compute layout
ps_graph_sd_ae <- igraph::graph_from_adjacency_matrix(ps_genus.sparcc_sd_ae$adjaMat1,weighted = TRUE)
ps_graph_ae_rec <- igraph::graph_from_adjacency_matrix(ps_genus.sparcc_ae_rec$adjaMat1,weighted = TRUE)
ps_lay_fr_sd_ae <- igraph::layout_with_fr(ps_graph_sd_ae)
ps_lay_fr_ae_rec <- igraph::layout_with_fr(ps_graph_ae_rec)

ps_species_graph_sd_ae <- igraph::graph_from_adjacency_matrix(ps_species.sparcc_sd_ae$adjaMat1,weighted = TRUE)
ps_species_graph_ae_rec <- igraph::graph_from_adjacency_matrix(ps_species.sparcc_ae_rec$adjaMat1,weighted = TRUE)
ps_species_lay_fr_sd_ae <- igraph::layout_with_fr(ps_species_graph_sd_ae)
ps_species_lay_fr_ae_rec <- igraph::layout_with_fr(ps_species_graph_ae_rec)

# Row names of the layout matrix must match the node names
rownames(ps_species_lay_fr_sd_ae) <- rownames(ps_species.sparcc_sd_ae$adjaMat1)
rownames(ps_species_lay_fr_ae_rec) <- rownames(ps_species.sparcc_ae_rec$adjaMat1)

# Get phyla names
taxtab <- as(tax_table(ps_genus_renamed), "matrix")
phyla <- as.factor(taxtab[, "Phylum"])
names(phyla) <- taxtab[, "Genus"]
table(phyla)

taxtab_species <- as(tax_table(ps_species_renamed), "matrix")
phyla_species <- as.factor(taxtab_species[, "Phylum"])
names(phyla_species) <- taxtab_species[, "Species"]
table(phyla_species)


# Define phylum colors
library(RColorBrewer)
phylcol <- c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"))

# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- colToTransp(phylcol, 12)


jpeg("./network_genus_sd_ae.jpeg",res=600,height = 4000,width = 8000)
plot(ps_genus.sparcc_sd_ae,
     layout = "spring",
     nodeSize = "eigenvector",  # "clr",
     nodeColor = "cluster",     #"feature",
     rmSingles = T, 
     cexNodes = 1.3,
     cexHubs = 1.3, 
     cexLabels = 3,
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "#00AFBB",
     negCol = "#E7B800",
     title1 = "Network on genus level with sparcc correlations", 
     showTitle = TRUE,
     groupNames = c("SD","AE-1"))

# legend(-1.2, 10, cex = 1.5, pt.cex = 2.5, title = "Phylum:", 
       # legend=levels(phyla), col = phylcol_transp, bty = "n", pch = 16) 

# legend(0.7, 1.0, cex = 2.2, title = "estimated correlation:",
       # legend = c("+","-"), lty = 1, lwd = 1, col = c("#00AFBB","#E7B800"), 
       # bty = "n", horiz = TRUE)
dev.off()



jpeg("./network_species_sd_ae.jpeg",res=600,height = 4000,width = 8000)
plot(ps_species.sparcc_sd_ae,
     layout = "spring",
     nodeSize = "eigenvector",  # "clr",
     nodeColor = "cluster",     #"feature",
     rmSingles = T, 
     cexNodes = 1.3,
     cexHubs = 1.3, 
     cexLabels = 3,
     featVecCol = phyla_species, 
     colorVec =  phylcol,
     posCol = "#00AFBB",
     negCol = "#E7B800",
     title1 = "Network on species level with sparcc correlations", 
     showTitle = TRUE,
     groupNames = c("SD","AE-1"))

# legend(-1.2, 10, cex = 1.5, pt.cex = 2.5, title = "Phylum:", 
# legend=levels(phyla), col = phylcol_transp, bty = "n", pch = 16) 

# legend(0.7, 1.0, cex = 2.2, title = "estimated correlation:",
# legend = c("+","-"), lty = 1, lwd = 1, col = c("#00AFBB","#E7B800"), 
# bty = "n", horiz = TRUE)
dev.off()

jpeg("./network_genus_ae_rec.jpeg",res=600,height = 4000,width = 8000)
plot(ps_genus.sparcc_ae_1_ae_2,
     layout = "spring",
     nodeSize = "eigenvector",  # "clr",
     nodeColor = "cluster",     #"feature",
     rmSingles = T, 
     cexNodes = 1.3,
     cexHubs = 1.3, 
     cexLabels = 3,
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "#00AFBB",
     negCol = "#E7B800",
     title1 = "Network on genus level with sparcc correlations", 
     showTitle = TRUE,
     groupNames = c("AE-1","AE-2"))

# legend(-1.2, 10, cex = 1.5, pt.cex = 2.5, title = "Phylum:", 
# legend=levels(phyla), col = phylcol_transp, bty = "n", pch = 16) 

# legend(0.7, 1.0, cex = 2.2, title = "estimated correlation:",
# legend = c("+","-"), lty = 1, lwd = 1, col = c("#00AFBB","#E7B800"), 
# bty = "n", horiz = TRUE)
dev.off()
###################### RECOVERY ANALYSIS############
copd_micro_phyloseq_recovery <- subset_samples(copd_micro_phyloseq, GROUP_SAMPLE !="SD")
# copd_micro_phyloseq_initial_abund <- filter_taxa(copd_micro_phyloseq_initial, function(x) sum(x > total*0.20) > 0, TRUE)
mm_test_Genus_recovery <- run_lefse( copd_micro_phyloseq_recovery,
                            wilcoxon_cutoff = 0.05,
                            norm = "CPM", 
                            group = "GROUP_SAMPLE",
                            taxa_rank = "Genus",
                            kw_cutoff = 0.05,
                            multigrp_strat = TRUE,
                            lda_cutoff = 1,
                            bootstrap_n = 50
)

mm_test_Species_recovery <- run_lefse( copd_micro_phyloseq_recovery,
                              wilcoxon_cutoff = 0.05,
                              norm = "CPM", 
                              group = "GROUP_SAMPLE",
                              taxa_rank = "Species",
                              kw_cutoff = 0.05,
                              multigrp_strat = TRUE,
                              lda_cutoff = 1,
                              bootstrap_n = 50
)

C_RE <- plot_ef_bar(mm_test_Genus_recovery)+ 
  scale_fill_manual(values = c( "#E7B800", "#FC4E07"))+
  theme(legend.position = "none")

D_RE <- plot_ef_bar(mm_test_Species_recovery)+ 
  scale_fill_manual(values = c( "#E7B800", "#FC4E07"))+
  theme(legend.position = "bottom")


A_RE <- plot_abundance(mm_test_Genus_recovery, group = "GROUP_SAMPLE") + 
  scale_fill_manual(values = c( "#E7B800", "#FC4E07"))+
  labs(fill = "Group") +
  theme_bw()+
  theme(legend.position = "none")

B_RE <-  plot_abundance(mm_test_Species_recovery, group = "GROUP_SAMPLE")+ 
  scale_fill_manual(values = c( "#E7B800", "#FC4E07"))+
  labs(fill = "Group") +
  theme_bw()+
  theme(legend.position = "bottom")

mm_test_Genus@marker_table$feature <- factor(mm_test_Genus@marker_table$feature, levels = unique(mm_test_Genus@marker_table$feature))


jpeg("./LDA_genus_sp_recovery.jpeg", height = 7000, width = 3500, res = 600)
A_RE+B_RE+C_RE+D_RE+plot_layout(ncol = 1)
dev.off()

mm_test_Genus_recovery@marker_table
mm_test_Species_recovery@marker_table

###############################
mm_test_Genus_ini <- run_lefse( copd_micro_phyloseq_initial,
                                     wilcoxon_cutoff = 0.05,
                                     norm = "CPM", 
                                     group = "GROUP_SAMPLE",
                                     taxa_rank = "all",
                                     kw_cutoff = 0.05,
                                     multigrp_strat = TRUE,
                                     lda_cutoff = 1,
                                     bootstrap_n = 50
)

plot_ef_bar(mm_test_Genus_ini)
C_PLOT_A <- plot_cladogram(mm_test_Genus_ini,color = c(SD="#00AFBB",`AE-1`="#E7B800"))


mm_test_Genus_rec <- run_lefse( copd_micro_phyloseq_recovery,
                                wilcoxon_cutoff = 0.05,
                                norm = "CPM", 
                                group = "GROUP_SAMPLE",
                                taxa_rank = "all",
                                kw_cutoff = 0.05,
                                multigrp_strat = TRUE,
                                lda_cutoff = 1,
                                bootstrap_n = 50
)

plot_ef_bar(mm_test_Genus_rec)
C_PLOT_B <- plot_cladogram(mm_test_Genus_rec,color = c(`AE-1`="#E7B800",`AE-2`="#FC4E07"))

jpeg("./C_plot_genus.jpeg", height = 6000, width = 3500, res = 300)
(C_PLOT_A/C_PLOT_B)+plot_layout(ncol = 1)
dev.off()





# recovery network

ps.f = prune_samples(sample_sums(copd_micro_phyloseq_recovery)>0, copd_micro_phyloseq_recovery) 
ps.f = prune_taxa(rowSums(otu_table(ps.f)) > 0, ps.f)  
ps_genus <- tax_glom(ps.f, taxrank = "Genus")
taxtab <- as(tax_table(ps_genus), "matrix")
# necomi 형식으로 변환 및 Genus level 설정 
ps_genus_renamed <- NetCoMi::renameTaxa(ps_genus,
                                        pat = "<name>",
                                        substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus")

ps_genus_renamed_ae <- phyloseq::subset_samples(ps_genus_renamed, 
                                                GROUP_SAMPLE == "AE-1")
ps_genus_renamed_rec <- phyloseq::subset_samples(ps_genus_renamed, 
                                                 GROUP_SAMPLE == "AE-2")

ps_genus.sparcc <- netConstruct(
  data = ps_genus_renamed_ae,
  data2 = ps_genus_renamed_rec,
  filtTax = "highestFreq",               # Top을 뽑는 기준
  filtTaxPar = list(highestFreq  = 100), # Top 100
  taxRank = "Genus",                     # Genus level 
  measure = "sparcc",                    # sparcc
  measurePar = list(nlambda=20, rep.num=10),
  
  normMethod = "clr",                    # transformation
  zeroMethod = "none",                   # zero 값 보정 
  sparsMethod = "threshold", 
  adjust = "adaptBH",                    # p-value 값 보정 
  thresh = 0.3,
  
  dissFunc = "signed",
  verbose = 2,
  seed = 42)

ps_props_sparcc <- netAnalyze(ps_genus.sparcc, 
                              clustMethod = "cluster_fast_greedy",
                              hubPar = "eigenvector",           # Hub node 판별방법
                              normDeg = FALSE)

summary(ps_props_sparcc, groupNames = c("AE-1", "AE-2"))
# 
# ps_props_sparcc <- netCompare(ps_props_sparcc,
#                               permTest = TRUE,
#                               nPerm = 1000,
#                               cores = 6,
#                               seed = 13075,
#                               # storeAssoPerm = TRUE,
#                               # fileStoreAssoPerm = "general/network_data/spring_assoPerm",
#                               verbose = TRUE)

# 
# ps_props_sparcc <- diffnet(ps_props_sparcc,
#                            diffMethod = "perm",
#                            # fileLoadAssoPerm = "general/network_data/spring_assoPerm",
#                            adjust = "lfdr")
# 



# Compute layout
ps_graph <- igraph::graph_from_adjacency_matrix(ps_genus.sparcc$adjaMat1,
                                                weighted = TRUE)

ps_lay_fr <- igraph::layout_with_fr(ps_graph)

# Row names of the layout matrix must match the node names
rownames(ps_lay_fr) <- rownames(ps_genus.sparcc$adjaMat1)


# Get phyla names
taxtab <- as(tax_table(ps_genus_renamed), "matrix")
phyla <- as.factor(taxtab[, "Phylum"])
names(phyla) <- taxtab[, "Genus"]
table(phyla)

# Define phylum colors
library(RColorBrewer)
phylcol <- c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"))

# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- colToTransp(phylcol, 12)


jpeg("./network_genus_rec.jpeg",res=600,height = 4000,width = 8000)
plot(ps_props_sparcc,
     layout = "spring",
     nodeSize = "eigenvector",  # "clr",
     nodeColor = "cluster",     #"feature",
     rmSingles = T, 
     cexNodes = 1.3,
     cexHubs = 1.3, 
     cexLabels = 3,
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "#00AFBB",
     negCol = "#E7B800",
     title1 = "Network on genus level with sparcc correlations", 
     showTitle = TRUE,
     groupNames = c("AE-1","AE-2"))

# legend(-1.2, 10, cex = 1.5, pt.cex = 2.5, title = "Phylum:", 
# legend=levels(phyla), col = phylcol_transp, bty = "n", pch = 16) 

# legend(0.7, 1.0, cex = 2.2, title = "estimated correlation:",
# legend = c("+","-"), lty = 1, lwd = 1, col = c("#00AFBB","#E7B800"), 
# bty = "n", horiz = TRUE)
dev.off()

















SPP_B <- ggboxplot(SPP_data, x = "GROUP_SAMPLE", y = "Spodiobacter_cordis",
                   color = "GROUP_SAMPLE", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                   add = "jitter", shape = "GROUP_SAMPLE")+ stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y =60)  + xlab("Group") +ylab("H.parainfluenzae (%)")+theme(legend.position = "none")






plot_net(copd_micro_phyloseq_initial, maxdist = 0.3, point_label = "SAMPLE_NO",color = "GROUP_SAMPLE")


ig <- make_network(mm_test_Genus, max.dist=0.3)

plot_network(ig, enterotype)





plot_cladogram(mm_test_Species, color = c(SD = "darkgreen")) +
  theme(plot.margin = margin(0, 0, 0, 0))












pht <- run_posthoc_test(copd_micro_phyloseq, group = "GROUP_SAMPLE")
plot_postHocTest(pht, feature = "k__Bacteria|p__Actinomycetota")
merged_data






copd_micro_phyloseq_stage <- subset_samples(copd_micro_phyloseq, COPD_STAGE !="1")
mm_test <- run_lefse( copd_micro_phyloseq_stage,
                      wilcoxon_cutoff = 0.01,
                      norm = "CPM", 
                      group = "COPD_STAGE",
                      taxa_rank = "Species",
                      kw_cutoff = 0.01,
                      multigrp_strat = TRUE,
                      lda_cutoff = 2,
                      bootstrap_n = 50
)
plot_abundance(mm_test, group = "COPD_STAGE")

marker_table(mm_test)
plot_ef_bar(mm_test)
plot_cladogram(mm_test, color = c(`SD` = "darkgreen",`AE-1`="red"))




















####################### Community analysis
threshold <- 0.05 / 100

# Function to filter taxa based on the relative abundance condition
filter_taxa_by_abundance <- function(physeq, threshold, min_samples) {
  # Convert counts to relative abundances
  relative_abundance <- transform_sample_counts(physeq, function(x) x / sum(x))
  
  # Filter taxa based on the condition
  keep_taxa <- apply(otu_table(relative_abundance), 1, function(x) sum(x >= threshold) >= min_samples)
  
  # Prune taxa to keep only those that meet the condition
  pruned_physeq <- prune_taxa(keep_taxa, physeq)
  
  return(pruned_physeq)
}

# Apply the function to the phyloseq object
copd_micro_phyloseq_filtered <- filter_taxa_by_abundance(copd_micro_phyloseq, threshold, 3)

# Check the filtered phyloseq object
copd_micro_phyloseq_filtered
