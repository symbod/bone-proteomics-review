
generate_sankey_plot <- function(data_info, colors_for_levels){
  ## First get necessary information in the right format
  input <- t(data_info[c(1,2,4,5),])
  input <- as.data.frame(input)
  colnames(input)[1] <- "Approach"
  
  
  # study sizes
  sizes <- input[, c("Approach", "Species", "Tissue")]
  sizes <- sizes %>% mutate(row_id = row_number()) %>%
    pivot_longer(-row_id, names_to = "variable", values_to = "value") %>%
    group_by(variable, value) %>%
    summarize(count = n(), .groups = "drop") #%>%
  sizes <- sizes[, c("value", "count")] %>%
    mutate(across(value, ~recode(.,
                                 "Cells" = "Cells",
                                 "Bone" = "Bone",
                                 "Serum-Plasma" = "Liquid-Biopsy",
                                 "Plasma-EVs" = "Liquid-Biopsy",
                                 "In vivo" = "In vivo")))
  sizes <- sizes %>% group_by(value) %>% summarize(count = sum(count), .groups = "drop")
  sizes$labels <- paste0(sizes$value, " (", sizes$count, ")")
  size_vector <- setNames(sizes$labels, sizes$value)
  
  
  df <- input %>% make_long(Species, Approach, Tissue)
  # TODO: change this
  df <- df %>%
    mutate(across(everything(), ~recode(.,
                                        "Cells" = "Cells",
                                        "Bone" = "Bone",
                                        "Serum-Plasma" = "Liquid-Biopsy",
                                        "Plasma-EVs" = "Liquid-Biopsy",
                                        "In vivo" = "In vivo")))
  
  levels <- c("-", "PXD", 
              "Rat", "Rabbit", "Mouse", "Human", 
              "In vivo", "In vitro", 
              "Liquid-Biopsy", "Bone", "ECM","Cells", "EVs")
  df$node <- factor(df$node, levels = levels)
  df$next_node <- factor(df$next_node, levels = levels)
  
  df <- df %>% mutate(across(everything(), ~recode(., !!!size_vector)))
  
  colors_for_levels_sizes <- setNames(colors_for_levels, size_vector[names(colors_for_levels)])
  
  ## generate plot
  sankey <- ggplot(df, aes(x = x, 
                           next_x = next_x, 
                           node = node, 
                           next_node = next_node,
                           fill = node,
                           label = node
  )) +
    geom_sankey(flow.alpha = .6, node.color = "black") +
    geom_sankey_label( size = 3, color = "black", fill="white") +
    scale_fill_manual(values = colors_for_levels_sizes) +
    theme_sankey(base_size = 10, base_family = "Arial")  +
    scale_x_discrete(position = "top") +
    theme(legend.position = "none",
          plot.title = element_text(hjust=.5),
          axis.text.x = element_text(colour = "black", size = 11, face = "bold")) +
    labs(x=NULL)
  return(sankey)
}

plot_boxplot_DAPS <- function(intersections_all, data_info, tissue_list, colors_for_levels){
  # Assuming intersections_all is already defined and contains 'Gene' and 'Datasets' columns
  intersections_all_exploded <- intersections_all %>%
    select(Gene, Datasets) %>%
    separate_rows(Datasets, sep = ";")
  
  # Count occurrences of each unique element in "Datasets"
  dataset_counts <- table(intersections_all_exploded$Datasets)
  
  # Convert the counts to a dataframe
  dataset_counts_df <- data.frame(Datasets = names(dataset_counts), Count = as.numeric(dataset_counts))
  dataset_counts_df$Datasets <- factor(dataset_counts_df$Datasets, levels = names(data_info))
  dataset_counts_df$Tissues <- unlist(c(data_info["Tissue", ])[dataset_counts_df$Datasets])
  dataset_counts_df$Tissues <- dataset_counts_df$Tissues
  
  # Create a mapping from original tissue names to new names using tissue_list
  tissue_mapping <- unlist(lapply(names(tissue_list), function(tissue) setNames(rep(tissue, length(tissue_list[[tissue]]$label)), tissue_list[[tissue]]$label)))
  
  # Replace tissue names using the mapping
  dataset_counts_df$Tissues  <- sapply(dataset_counts_df$Tissues, function(x) tissue_mapping[[x]])
  dataset_counts_df$Assay <- factor(sapply(dataset_counts_df$Tissues, function(x) tissue_list[[x]]$assay))
  dataset_counts_df$Tissue  <- factor(dataset_counts_df$Tissue, levels = c("Cells", "EVs", "ECM","Bone", "Liquid-Biopsy"))
  
  
  # Create the boxplot with log-transformed y-axis
  boxplot <- ggplot(dataset_counts_df, aes(x = Tissues, y = Count, fill = Tissues)) +
    geom_boxplot(linewidth = 0.3, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.3) +  # Add jittered points
    scale_fill_manual(values = colors_for_levels) +  # Custom colors
    scale_y_log10() +  # Log-transform the y-axis
    facet_grid(. ~ Assay, scales = "free_x", space = "free_x") +  # Separate 'In vitro' and 'In vivo'
    labs(x = "",
         y = "Nr. of Proteins (Log Scale)",
         fill = "Tissue") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          #panel.spacing = unit(2, "lines"),
          strip.text.x = element_text(size = 11, face = "bold"),
          panel.border = element_rect(color = "black", linewidth = 0.3, fill = NA)) 
  return(boxplot)
}

plot_barplot_CAPS <- function(intersections_all, filtered_intersections, tissue_list, colors_for_levels){
  # Your existing pre-filtered data preparation
  genes_per_tissue <- intersections_all %>%
    rename(Tissue = Type) %>%
    group_by(Tissue) %>%
    summarise(SumUniqueGenes = n_distinct(Gene)) %>%
    mutate(DataType = "Pre-Filtered")
  
  # Prepare post-filtered data
  filtered_genes_per_tissue <- filtered_intersections %>%
    rename(Tissue = Type) %>%
    group_by(Tissue) %>%
    summarise(SumUniqueGenes = n_distinct(Gene)) %>%
    mutate(DataType = "Post-Filtered")
  
  # Combine both datasets
  combined_genes_data <- rbind(genes_per_tissue, filtered_genes_per_tissue)
  
  combined_genes_data$DataType <- factor(combined_genes_data$DataType, levels = c("Pre-Filtered", "Post-Filtered"))
  
  combined_genes_data$Assay <- factor(sapply(combined_genes_data$Tissue, function(x) tissue_list[[x]]$assay))
  
  combined_genes_data$Tissue  <- factor(combined_genes_data$Tissue, levels = c("Cells", "EVs", "ECM","Bone", "Liquid-Biopsy"))
  
  
  # Create the bar plot
  bar_plot <- ggplot(combined_genes_data, aes(x = Tissue, y = SumUniqueGenes, fill = Tissue, alpha = DataType)) +
    facet_grid(.~Assay, scales = "free_x", space = "free") +
    geom_bar(stat = "identity", position = position_dodge(), color = "black", linewidth = 0.3) +
    #geom_text(aes(label = SumUniqueGenes), 
    #          position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
    scale_y_log10() + # Log-transform the y-axis
    scale_fill_manual(values = colors_for_levels) +
    scale_alpha_manual(values = c(0.4, 1)) + # Pre-filtered bars are 50% transparent
    labs(
      x = "",
      y = "Nr. of Proteins (Log Scale)") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 11, face = "bold"), 
          legend.position = "bottom",
          panel.border = element_rect(color = "black", linewidth = 0.3, fill = NA),
          legend.spacing.x = unit(0.1, "cm"),
          legend.text = element_text(size = 11),
          legend.title = element_text(face = "bold")) + 
    guides(fill = guide_legend(title = "Tissue", order = 1, byrow = TRUE),
           alpha = guide_legend(title = "WSI", order = 2,byrow = TRUE))
  return(bar_plot)
}

