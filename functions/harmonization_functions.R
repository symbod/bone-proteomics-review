
# Function to create Venn diagrams showing harmonization purpose
create_harmonization_venn_diagrams <- function(meta_info){
  data <- meta_info[c("Tissue","Species"),]
  data <- t(data)
  studies <- row.names(data)
  data <- as.data.table(data)
  data$Study <- studies
  data$GeneName <- "Yes"
  data[data$Study == "9_Human"]$GeneName <- "No"  # TODO: change this to make it more flexible
  data[data$Study == "30_Rat"]$GeneName <- "ENSEMBL"
  data$UniprotID <- "Yes"
  data[data$Tissue == "ECM"]$UniprotID <- "No"
  data[data$Study == "30_Rat"]$UniprotID <- "ENSEMBL"
  

  species_colors <- colors_for_levels[unique(data$Species)]
  
  # Availability of Gene Names
  tmp <- data %>% group_by(Species, GeneName) %>% count()
  
  
  # Basic piechart
  tmp$Species <- factor(tmp$Species, levels = unique(data$Species))
  gn_plot <- ggplot(as.data.table(tmp), aes(x=1, y=n, fill = Species, pattern = GeneName)) +
    geom_col_pattern(
      colour = "black",
      pattern_fill = "grey50",
      pattern_density = 0.35,
      pattern_spacing = 0.03)  + 
    coord_polar("y") +
    theme_void() +
    theme(legend.key.size = unit(1, 'cm')) +
    scale_pattern_manual(
      name = "Availability",
      values = c(
        "Yes" = "none",
        "No" = "stripe",
        "ENSEMBL" = "circle"
      ),
      guide = guide_legend(override.aes = list(fill = "white", pattern = c("none", "stripe", "circle"))),
      labels = c("Available", "Not available", "Only ENSEMBL")
    )  +
    scale_fill_manual(
      name = "Species",
      values = species_colors,
      guide = guide_legend(override.aes = list(pattern = "none")))+
    ggtitle("Using Gene Symbols as Input") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  # Availability of Uniprot IDs
  tmp <- data %>% group_by(Species, UniprotID) %>% count()
  
  # Basic piechart
  tmp$Species <- factor(tmp$Species, levels = unique(data$Species))
  pi_plot <- ggplot(as.data.table(tmp), aes(x=1, y=n, fill = Species, pattern = UniprotID)) +
    geom_col_pattern(
      colour = "black",
      pattern_fill = "grey50",
      pattern_density = 0.35,
      pattern_spacing = 0.03)  + 
    coord_polar("y") +
    theme_void() +
    theme(legend.key.size = unit(1, 'cm')) +
    scale_pattern_manual(
      name = "Availability",
      values = c(
        "Yes" = "none",
        "No" = "stripe",
        "ENSEMBL" = "circle"
      ),
      guide = guide_legend(override.aes = list(fill = "white", pattern = c("none", "stripe", "circle"))),
      labels = c("Available", "Not available", "Only ENSEMBL")
    )  +
    scale_fill_manual(
      name = "Species",
      values = species_colors,
      guide = guide_legend(override.aes = list(pattern = "none"))) +
    ggtitle("Using Protein IDs as Input") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # ProHarMeD Plot
  tmp <- data.table(Species = c("Human"), "GeneName" = c("Yes"), "n" = c(31))
  proharmed_plot <- ggplot(as.data.table(tmp), aes(x=1, y=n, pattern = GeneName, fill = Species)) +
    geom_col_pattern(
      colour = "black",
      pattern_fill = "grey50",
      pattern_density = 0.5,
      pattern_key_scale_factor = 0.3)  + 
    coord_polar("y") +
    theme_void() +
    theme(legend.position = "none") + scale_fill_manual(values = species_colors) +
    scale_pattern_manual(
      name = "Availability of Protein IDs",
      values = c(
        "Yes" = "none"),
      guide = guide_legend(override.aes = list(fill = "white", pattern = c("none"))),
      labels = c("Available")
    ) +
    ggtitle("After Harmonization with ProHarMeD \n on Gene Symbol Basis") +
    theme(plot.title = element_text(hjust = 0.5))
  
  legend_theme <-theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14))
  gn_plot <- gn_plot + legend_theme
  pie_plots <- ggarrange(pi_plot, gn_plot, proharmed_plot, ncol = 3, common.legend = TRUE, legend = "right", labels =  c("A", "B", "C"), font.label = list(size = 17))
  return(pie_plots)

}
