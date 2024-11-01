# Divertisty analysis

# copd_micro_phyloseq - phyloseq package object

shannon_data<- 
  estimate_richness(copd_micro_phyloseq, split = TRUE, measures = NULL) %>%
  left_join(metadata)

# graph for each indices (Fisher's plot in Figure )

Fisher_plot_D <-  ggpaired(shannon_data, 
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

