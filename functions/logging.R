
plot_proharmed_logging <- function(logging_path){
  ordered_studies <- c("1_Human", "2_Human", "3_Human", "4_Human", "5_Human", "6_Human", "7_Human", "8_Human", "9_Human", "10_Human", "11_Human", "12_Human",
                       "13_Human", "14_Human", "15_Human", "16_Human", "17_Human", "18_Human", "19_Human", "20_Human", "21_Human", "22_Human", "23_Human", "24_Human", "25_Mouse", "26_Mouse",
                       "27_Mouse", "28_Rat", "29_Rat","30_Rat" ,"31_Rabbit")
  
  # ==== Filter IDs =====
  # Load the data from the first file
  data <- read.csv(paste0(logging_path, "filter_ids_overview.csv"))
  data$X <- NULL
  
  # Group the data by 'Study' and sum the 'Nr Filtered IDs' and 'Nr Removed IDs', reordering according to the original data
  grouped_data <- data %>% group_by(Study) %>% summarise(Nr.Filtered.IDs = sum(Nr.Filtered.IDs), Nr.Removed.IDs = sum(Nr.Removed.IDs))
  colnames(grouped_data)[colnames(grouped_data) == "Nr.Filtered.IDs"] <- "Nr.Kept.IDs"
  
  # add studies without loss
  studies_without_loss <- ordered_studies[!ordered_studies %in% unique(grouped_data$Study)]
  studies_without_loss <- data.frame(Study = studies_without_loss,
                                     Nr.Kept.IDs = rep(0, length(studies_without_loss)),
                                     Nr.Removed.IDs = rep(0, length(studies_without_loss)))
  grouped_data <- rbind(grouped_data, studies_without_loss)
  grouped_data$Study <- factor(grouped_data$Study, level = ordered_studies)
  grouped_data <- grouped_data[order(grouped_data$Study),]
  
  # Load the data from the second file
  data_detailed <- read.csv(paste0(logging_path, "filter_ids_detailed.csv"))
  data_detailed$X <- NULL
  
  # Counting 'Obsolete ID' and 'Wrong Organism' for each study
  data_detailed$Obsolete.ID <- data_detailed$Organism == "Not found"
  data_detailed$Wrong.Organism <- ! data_detailed$Obsolete.ID
  grouped_reasons <- data_detailed %>% group_by(Study) %>% summarise(Obsolete.ID = sum(Obsolete.ID), Wrong.Organism = sum(Wrong.Organism))
  
  # add studies without loss
  studies_without_loss <- ordered_studies[!ordered_studies %in% unique(grouped_reasons$Study)]
  studies_without_loss <- data.frame(Study = studies_without_loss,
                                     Obsolete.ID = rep(0, length(studies_without_loss)),
                                     Wrong.Organism = rep(0, length(studies_without_loss)))
  
  grouped_reasons <- rbind(grouped_reasons, studies_without_loss)
  grouped_reasons$Study <- factor(grouped_reasons$Study, level = ordered_studies)
  grouped_reasons <- grouped_reasons[order(grouped_reasons$Study),]
  
  
  # Plotting
  
  total <- grouped_data$Nr.Kept.IDs + grouped_data$Nr.Removed.IDs
  colnames(grouped_data) <- c("Study", "Kept ID", "Removed ID")
  long_grouped_data <- melt(grouped_data, id.vars = "Study", value.name = "Nr.IDs", variable.name = "Type")
  long_grouped_data$Perc <- round(((long_grouped_data$Nr.IDs / total) * 100), 1)
  long_grouped_data$Label <- paste0(long_grouped_data$Nr.IDs, " (", long_grouped_data$Perc, "%)")
  long_grouped_data$Label[long_grouped_data$Nr.IDs == 0] <- ""
  long_grouped_data$Type <- factor(long_grouped_data$Type, levels = c("Removed ID", "Kept ID"))
  
  p1 <- ggplot(long_grouped_data, aes( x = Study, y = Perc, fill = Type)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("plum", "lightsteelblue")) + 
    labs(x="Study", y = "Number of IDs (%)") + 
    theme_bw() + 
    scale_y_continuous(limits = c(0, 105)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom") +
    geom_text(aes(label = Label),
              position = position_stack(vjust = 0.3),
              size = 3,
              color = "black",
              angle = 90)
  
  
  
  
  total <- grouped_reasons$Obsolete.ID + grouped_reasons$Wrong.Organism
  colnames(grouped_reasons) <- c("Study", "Obsolete ID", "Wrong Organism")
  long_grouped_reasons_data <- melt(grouped_reasons, id.vars = "Study", value.name = "Nr.IDs", variable.name = "Type")
  long_grouped_reasons_data$Perc <- round(((long_grouped_reasons_data$Nr.IDs / total) * 100), 1)
  long_grouped_reasons_data$Label <- paste0(long_grouped_reasons_data$Nr.IDs, " (", long_grouped_reasons_data$Perc, "%)")
  long_grouped_reasons_data$Label[long_grouped_reasons_data$Nr.IDs == 0] <- ""
  long_grouped_reasons_data$Type <- factor(long_grouped_reasons_data$Type, levels = c("Wrong Organism", "Obsolete ID"))
  
  
  p2 <- ggplot(long_grouped_reasons_data, aes( x = Study, y = Perc, fill = Type)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("darkseagreen", "lightgray")) + 
    labs(x="Study", y = "Number of IDs (%)") + 
    scale_y_continuous(limits = c(0, 105)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom") +
    geom_text(aes(label = Label),
              position = position_stack(vjust = 0.5),
              size = 3,
              color = "black",
              angle = 90)
  
  
  # ==== Filter IDs =====
  
  filter_grouped <- data %>% group_by(Study) %>% summarize(Rows.Initial = sum(Nr.Filtered.IDs > 0))
  
  # add studies without loss
  studies_without_loss <- ordered_studies[!ordered_studies %in% unique(filter_grouped$Study)]
  studies_without_loss <- data.frame(Study = studies_without_loss,
                                     Rows.Initial = rep(0, length(studies_without_loss)))
  
  filter_grouped <- rbind(filter_grouped, studies_without_loss)
  
  
  # ==== Remap Symbols =====
  
  remap_data <- read.csv(paste0(logging_path, "remap_genenames_overview.csv"))
  remap_data$X <- NULL
  remap_grouped <- remap_data %>% group_by(Study) %>% summarize(Rows.After.Remapping = sum(Nr.Remapped.Gene.Names > 0))
  
  # add studies without loss
  studies_without_loss <- ordered_studies[!ordered_studies %in% unique(remap_grouped$Study)]
  studies_without_loss <- data.frame(Study = studies_without_loss,
                                     Rows.After.Remapping = rep(0, length(studies_without_loss)))
  
  remap_grouped <- rbind(remap_grouped, studies_without_loss)
  
  
  # ==== Reduction Symbols =====
  reduced_data <- read.csv(paste0(logging_path, "reduce_genenames_overview.csv"))
  reduced_data$X <- NULL
  reduced_grouped <- reduced_data %>% group_by(Study) %>% summarize(Rows.After.Reduction = sum(Nr.Reduced.Gene.Names > 0), Rows.Before.Reduction = sum(Nr.Gene.Names > 0))
  
  # add studies without loss
  studies_without_loss <- ordered_studies[!ordered_studies %in% unique(reduced_data$Study)]
  studies_without_loss <- data.frame(Study = studies_without_loss,
                                     Rows.After.Reduction = rep(0, length(studies_without_loss)),
                                     Rows.Before.Reduction = rep(0, length(studies_without_loss)))
  
  reduced_grouped <- rbind(reduced_grouped, studies_without_loss)
  
  
  
  # ==== Ortholog Symbols =====
  ortholog_data <- read.csv(paste0(logging_path,"map_orthologs_overview.csv"))
  ortholog_data$X <- NULL
  
  ortholog_grouped <- ortholog_data %>% group_by(Study) %>% summarize(Rows.After.Ortholog = sum(Nr.Ortholog.Gene.Names > 0 ))
  # add studies without loss
  studies_without_loss <- ordered_studies[!ordered_studies %in% unique(ortholog_grouped$Study)]
  studies_without_loss <- data.frame(Study = studies_without_loss,
                                     Rows.After.Ortholog = rep(0, length(studies_without_loss)))
  
  ortholog_grouped <- rbind(ortholog_grouped, studies_without_loss)
  
  
  ## Merge all 
  
  all_grouped <- merge(reduced_grouped, ortholog_grouped, by = "Study")
  all_grouped <- merge(remap_grouped, all_grouped, by = "Study")
  all_grouped <- merge(filter_grouped, all_grouped, by = "Study")
  
  # Change initial and final gene list size based on study
  
  all_grouped$Rows.Final <- all_grouped$Rows.After.Reduction
  
  not_human_studies <- c("25_Mouse", "26_Mouse","28_Rat", "29_Rat","30_Rat","31_Rabbit") # 27_Human not included because human proteins reported
  # final list of not humans studies from ortholog data
  all_grouped$Rows.Final[all_grouped$Study %in% not_human_studies] <- all_grouped$Rows.After.Ortholog[all_grouped$Study %in% not_human_studies]
  
  # initial list of 25_mouse, 26_mouse, and 30_rat
  initial_gene_symbol_studies <- c("25_Mouse", "26_Mouse", "30_Rat")
  all_grouped$Rows.Initial[all_grouped$Study %in% initial_gene_symbol_studies] <- all_grouped$Rows.Before.Reduction[all_grouped$Study %in% initial_gene_symbol_studies]
  
  # calculate loss of all mappings of all human studies
  all_grouped_humans <- all_grouped[!all_grouped$Study %in% not_human_studies,]
  all_grouped_humans$Loss.Remapping <- all_grouped_humans$Rows.Initial - all_grouped_humans$Rows.After.Remapping
  all_grouped_humans$Loss.Reduction <- all_grouped_humans$Rows.After.Remapping - all_grouped_humans$Rows.After.Reduction
  all_grouped_humans$Loss.Ortholog <- 0
  
  # calculate loss of all mapping of not human studies
  all_grouped_not_humans <- all_grouped[all_grouped$Study %in% not_human_studies,]
  all_grouped_not_humans$Loss.Remapping <- 0
  all_grouped_not_humans$Loss.Reduction <- all_grouped_not_humans$Rows.Initial - all_grouped_not_humans$Rows.After.Reduction
  all_grouped_not_humans$Loss.Ortholog <- all_grouped_not_humans$Rows.After.Reduction - all_grouped_not_humans$Rows.After.Ortholog
  
  all_grouped_final <- rbind(all_grouped_humans, all_grouped_not_humans)
  
  all_grouped_final$Study <- factor(all_grouped_final$Study, level = ordered_studies)
  all_grouped_final <- all_grouped_final[order(all_grouped_final$Study),]
  
  total <- all_grouped_final$Rows.Initial
  all_grouped_final <- all_grouped_final[, c("Study", "Loss.Remapping", "Loss.Reduction", "Loss.Ortholog", "Rows.Final")]
  colnames(all_grouped_final) <- c("Study", "Loss of Remapping", "Loss of Reduction", "Loss of Ortholog Mapping", "Final List")
  
  
  long_all_grouped_final <- melt(all_grouped_final, id.vars = "Study", value.name = "Nr.IDs", variable.name = "Type")
  long_all_grouped_final$Perc <- round(((long_all_grouped_final$Nr.IDs / total) * 100), 1)
  long_all_grouped_final$Label <- paste0(long_all_grouped_final$Nr.IDs, " (", long_all_grouped_final$Perc, "%)")
  long_all_grouped_final$Type <- factor(long_all_grouped_final$Type, levels = c("Final List", "Loss of Remapping", "Loss of Reduction", "Loss of Ortholog Mapping"))
  long_all_grouped_final$Label[long_all_grouped_final$Nr.IDs == 0] <- ""
  
  
  long_all_grouped_final$Study <- factor(long_all_grouped_final$Study, level = ordered_studies)
  long_all_grouped_final <- long_all_grouped_final[order(long_all_grouped_final$Study),]
  
  
  p3 <- ggplot(long_all_grouped_final, aes( x = Study, y = Perc, fill = Type)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("lightgray","darkseagreen", "plum", "lightsteelblue")) + 
    labs(x="Study", y = "Number of Gene Symbols (%)") + 
    theme_bw() +
    scale_y_continuous(limits = c(0,110)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom") 
  geom_richtext(data = long_all_grouped_final[long_all_grouped_final$Perc > 10, ], aes(label = Label),
                position = position_stack(vjust = 0.5),
                size = 3,
                color = "black",
                angle = 90)
  
  
  plot <- ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 1)
  return(plot)
}

