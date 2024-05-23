#' Explode Data Frame Columns
#'
#' This function takes a data frame and "explodes" each column by a specified delimiter (';'). 
#' It creates a binary matrix indicating the presence of unique elements in each exploded column.
#'
#' @param data A data frame with columns that need to be exploded.
#'
#' @return A data frame where each column represents an exploded version of the original columns.
#' Each row indicates the presence (1) or absence (0) of each unique element in the original column.
#'
#' @examples
#' df <- data.frame(col1 = c("A;B", "B;C"), col2 = c("1;2", "2;3"))
#' exploded_df <- explode_data(df)
#' exploded_df
#'
#' @importFrom tidyr separate_rows
#' @export
explode_data <- function(data) {
  # Initialize a list to store exploded data
  exploded_list <- list()
  
  # Explode each column by ';' and store in exploded_list
  for (col_name in colnames(data)) {
    exploded <- data.frame(c = data %>% select(all_of(col_name)))
    names(exploded) <- col_name
    exploded <- separate_rows(exploded, col_name, sep = ";")
    exploded <- na.omit(exploded)
    exploded_list[[col_name]] <- exploded
  }
  
  # Create a binary matrix indicating presence of unique elements
  exploded_matrix <- lapply(exploded_list, function(exploded) {
    sapply(exploded, function(x) {
      ifelse(unique(unlist(exploded_list)) %in% x, 1, 0)
    })
  })
  
  # Convert the list of matrices into a data frame
  exploded_data <- do.call(cbind, exploded_matrix)
  rownames(exploded_data) <- unique(unlist(exploded_list))
  
  return(exploded_data)
}


#' Process Data Intersection
#'
#' This function filters data based on assay, label, and value, calculates dataset sizes,
#' and determines intersection sizes and datasets that intersect.
#'
#' @param data The data frame to be processed.
#' @param data_info Data frame containing assay information.
#' @param label Column name used for labeling in data_info.
#' @param value Values to filter on in the label column of data_info.
#' @param assay Assay type to filter the data.
#' @return A list containing the processed data frame and set numbers.
process_data_intersection <- function(data, data_info, label, value, assay) {
  # Filter data based on assay, label, and value
  data_info_assay <- data_info[,data_info["Assay",] == assay]
  set_nr <- colnames(data_info_assay[data_info_assay[label, ] %in% value])
  dt <- data[, set_nr]

  # Calculate the size of each dataset
  sizes <- colSums(dt)
  
  # Convert to numeric and divide by dataset sizes
  dt <- sweep(as.data.frame(dt), 2, sizes, "/")
  
  # Add intersection size and datasets that intersect
  dt$Intersected <- rowSums(dt != 0)
  dt$Datasets <- apply(dt[, !colnames(dt) %in% "Intersected"], 1, function(row) paste(names(row)[row != 0], collapse = ";"))
  
  return(list(dt=dt, set_nr=set_nr))
}


#' Get Intersection Data Table
#'
#' This function processes the data to get intersection information, adds weights to each protein,
#' and organizes the data into a final data table.
#'
#' @param data The data frame to be processed.
#' @param data_info Data frame containing assay information.
#' @param label Column name used for labeling in data_info.
#' @param value Values to filter on in the label column of data_info.
#' @param assay Assay type to filter the data.
#' @param type Type of analysis or data category.
#' @return A data frame with intersection information and weights for each protein.
get_intersection_dt <- function(data, data_info, label, value, assay, type) {
  # Call the subfunction for common operations
  processed_data <- process_data_intersection(data, data_info, label, value, assay)
  dt <- processed_data$dt
  
  # Add weight of each protein
  dt$Weight <- rowSums(dt[,1:(ncol(dt)-2)])
  
  # Subset and reorder columns
  dt <- dt[dt$Intersected != 0, c("Intersected", "Weight", "Datasets")]
  
  # Add extra info
  dt[, c("Assay", "Type", "Size", "Gene")] <- list(assay, type, length(processed_data$set_nr), rownames(dt))
  dt <- dt[, c("Gene", "Datasets", "Intersected", "Size", "Assay", "Type", "Weight")]
  return(dt)
}

#' Get Intersection Data Table with Double Weighting
#'
#' This function extends the get_intersection_dt functionality by applying double weighting
#' based on type sizes and organizing the data into a final data table.
#'
#' @param data The data frame to be processed.
#' @param data_info Data frame containing assay information.
#' @param label Column name used for labeling in data_info.
#' @param value Values to filter on in the label column of data_info.
#' @param assay Assay type to filter the data.
#' @param type Type of analysis or data category.
#' @return A data frame with intersection information, double weighting, and additional data.
get_intersection_dt_double_weighting <- function(data, data_info, label, value, assay, type) {
  # Call the subfunction for common operations
  processed_data <- process_data_intersection(data, data_info, label, value, assay)
  dt <- processed_data$dt
  
  # Calculate type sizes for double weighting
  studies <- sapply(value, function(t) {
    length(colnames(processed_data$dt[data_info[data_info[label, ] == t, ] %in% value]))
  })
  names(studies) <- value
  
  # Apply double weighting
  dt <- sweep(dt, 2, studies[value], '/')
  
  # Add protein weight and select columns
  dt$Weight <- rowSums(dt[, 1:(ncol(dt) - 3)])
  dt <- dt[dt$Intersected != 0, c("Intersected", "Weight", "Datasets")]
  
  # Add extra info and reorder columns
  dt[, c("Assay", "Type", "Size", "Gene")] <- list(assay, type, length(processed_data$set_nr), rownames(dt))
  dt <- dt[, c("Gene", "Datasets", "Intersected", "Size", "Assay", "Type", "Weight")]
  return(dt)
}


#' Get Intersection Table for All Tissues
#'
#' @param data The data frame to be processed.
#' @param data_info Data frame containing assay information.
#' @param tissue_list 
#' @param weighting 
#'
#' @return data table of intersections of all tissues
#' @export
#'
calculate_intersections <- function(data, data_info, tissue_list, weighting){
  # Initialize a list to store intermediate data frames
  intersections_list <- list()
  
  for(tissue in names(tissue_list)) {
    if(weighting == "single") {
      intersections_list[[tissue]] <- get_intersection_dt(data = data, data_info = data_info, label = "Tissue", 
                                                          value = c(tissue_list[[tissue]]$label), 
                                                          assay = tissue_list[[tissue]]$assay, 
                                                          type = tissue)
    } else if(weighting == "double") {
      intersections_list[[tissue]] <- get_intersection_dt_double_weighting(data = data, data_info = data_info, label = "Tissue", 
                                                                           value = c(tissue_list[[tissue]]$label), 
                                                                           assay = tissue_list[[tissue]]$assay, 
                                                                           type = tissue)
    } else {
      stop("Exiting RMarkdown due to wrong/missing value for 'weighting' inside config.R.")
    }
  }
  
  # Combine all data frames into one
  intersections <- do.call(rbind, intersections_list)
  return(intersections)
}


#' Select Intersections
#'
#' This function selects intersections based on a given threshold and weighting criterion from a subset study dataframe.
#'
#' @param subset_study_df A dataframe containing the study data.
#' @param intersection_thr The threshold for the normal intersection, default is 3.
#' @param weighting_thr The threshold type for weighting, can be "normal" or "max", default is "normal".
#' @return A dataframe with selected intersections.
#' @examples
#' select_intersections(my_study_df, 3, "normal")
#' select_intersections(my_study_df, 3, "max")
#' @export
select_intersections <- function(subset_study_df, intersection_thr = 3, weighting_thr = "normal") {
  # Calculate the normal intersection threshold
  normal_thresh <- floor(max(subset_study_df$Size) / intersection_thr)
  
  # Subset for normal intersections
  normal_intersection <- subset_study_df %>%
    filter(Intersected >= normal_thresh) %>%
    mutate(Intersection = "normal")
  
  # Determine the minimum weight based on weighting threshold
  min_weight <- if (weighting_thr == "normal") {
    min(normal_intersection$Weight)
  } else if (weighting_thr == "max") {
    min(subset_study_df %>% filter(Intersected == max(Intersected))$Weight)
  } else {
    stop("Invalid weighting threshold. Choose 'normal' or 'max'.")
  }
  
  # Subset for weighted intersections
  weighted_intersection <- subset_study_df %>%
    filter(Weight >= min_weight) %>%
    mutate(Intersection = "weighted")
  
  # Combine and mark genes present in both intersections
  intersection <- rbind(normal_intersection, weighted_intersection) %>%
    arrange(Gene) %>%
    mutate(Intersection = ifelse(duplicated(Gene) | duplicated(Gene, fromLast = TRUE), "both", Intersection)) %>%
    distinct()
  
  return(intersection)
}

filter_intersections <- function(intersections, intersection_thr, weighting_thr){
  # Group by 'Type' and apply select_intersections() function to each group
  resulting_data_frames <- intersections %>%
    group_by(Type) %>%
    do(select_intersections(., intersection_thr=intersection_thr, weighting_thr=weighting_thr)) %>%
    ungroup()
  
  # Combine resulting data frames into one
  filtered_intersections <- bind_rows(resulting_data_frames)
  return(filtered_intersections)
}

### Plots


#' Plot Weighted Intersection
#'
#' This function creates a bar plot for weighted intersections from a given dataset.
#'
#' @param weighted_dt A dataframe with weighted intersection data.
#' @param min_intersected Minimum number of intersections.
#' @param min_weight Minimum weight threshold.
#' @return A list containing the plot and the modified dataframe.
#' @examples
#' plot_weighted_intersection2(weighted_data, 3, 0.5)
#' @export

plot_weighted_intersection <- function(weighted_dt, min_intersected, min_weight) {
  # Remove intersections 0 and 1
  dt <- weighted_dt#[!weighted_dt$Intersected %in% c(0,1),]
  
  # Generate plot per intersection size
  dt <- dt %>% 
    group_by(Intersected, Weight) %>% 
    summarise(Nr_Genes = length(unique(Gene)), Gene = paste(unique(Gene), collapse = "; "), .groups = 'drop')
  
  # Sort to give pseudonames
  dt <- dt[order(dt$Intersected, dt$Weight),]
  
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  dt$Labels <- LETTERS702[1:nrow(dt)]
  dt$Labels <- paste(dt$Labels, " (", dt$Nr_Genes, ")", sep = "")
  
  dt$Labels <- factor(dt$Labels, levels = dt$Labels)
  
  # Define colors based on conditions
  dt$FillColor <- ifelse(dt$Intersected == 1, "Below Min", ifelse(dt$Intersected >= min_intersected, "Above Min Intersected",
                                                                  ifelse(dt$Weight > min_weight, "Above Min Weight", "Below Min")))
  
  # Start the plot
  p <- ggplot(dt, aes(x = Labels, y = Weight, fill = FillColor)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = c("Below Min" = "grey", "Above Min Intersected" = "darkblue", "Above Min Weight" = "orange"),
                      name = "Status") + 
    geom_hline(aes(yintercept = min_weight, linetype = "Minimum Weight"), 
               color = "orange", linetype = "dotted",linewidth = 1, show.legend = FALSE) +
    scale_linetype_manual(values = c("Minimum Weight" = "dotted"),
                          name = "Min Weight Line") +
    facet_wrap(~Intersected, ncol = max(dt$Intersected), scales = "free_x") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.box = "vertical") +
    guides(fill = guide_legend(order = 1, title.position = "top"), 
           linetype = guide_legend(order = 2, title.position = "bottom")) + # Control the order of the legends
    ggtitle(paste0(subset_df$Type[[1]], " (Size:", subset_df$Size[[1]], ")")) 
  
  return(list(Plot = p, Data = dt))
}

plot_original_vs_harmonized_intersections <- function(harmonized_intersections, harmonized_intersections_filtered, original_intersections){
  harmonized_dt <- harmonized_intersections[, c("Gene", "Intersected", "Type")]
  colnames(harmonized_dt) <- c("Gene", "ProHarMeD", "Type")
  original_dt <- original_intersections[, c("Gene", "Intersected", "Type")]
  colnames(original_dt) <- c("Gene", "Original", "Type")
  
  dt <- merge(harmonized_dt, original_dt, by = c("Gene", "Type"), all.x = TRUE, all.y = TRUE)
  dt$ProHarMeD[is.na(dt$ProHarMeD)] <- 0
  dt$Original[is.na(dt$Original)] <- 0
  
  # Define base colors for each Type
  base_colors <- c(
    "EVs" = "#A3A500",
    "Cells" = "#F8766D",
    "Liquid-Biopsy" = "#E76BF3",
    "Bone" = "#00B0F6",
    "ECM" = "#00BF7D"
  )
  
  # Pie Chart of All Intersections
  tmp <- dt %>% group_by(ProHarMeD, Original, Type) %>% summarize(n())
  tmp <- tmp[tmp$ProHarMeD != 0 & tmp$Original != 0,]
  colnames(tmp) <- c("ProHarMeD","Original", "Type", "Intersection")
  wide_tmp <- dcast(tmp, ProHarMeD+Original~Type, value.var = "Intersection")
  wide_tmp[is.na(wide_tmp)] <- 0
  
  summary_tmp <- tmp %>% group_by(ProHarMeD, Original) %>% summarize("Intersection" = sum(Intersection))
  comparison_plot <- ggplot() + geom_abline() + geom_scatterpie(aes(x=Original, y=ProHarMeD), data=wide_tmp,
                                                                cols=names(base_colors), pie_scale = 2) + coord_equal() +
    scale_x_continuous(breaks = c(1,2,3,4,5,6,7), limits = c(0.5,7.5)) + scale_y_continuous(breaks = c(1,2,3,4,5,6,7), limits = c(0.5,7.5)) +  geom_text(data = summary_tmp, aes(x = Original, y = ProHarMeD -0.5, label = Intersection)) +
    scale_fill_manual(name = "Tissue", values = base_colors)  + labs(y = "ProHarMeD") + theme_bw() + theme(axis.text = element_text(size = 12)) +
    labs(x = "Intersection Size in Original Data", y = "Intersection Size in Harmonized Data")

  # Pie Chart of Filtered Intersections
  # check with tissues
  dt <- merge(dt, harmonized_intersections_filtered[, c("Gene", "Type")], by = c("Gene", "Type"), all.y = TRUE, all.x = FALSE)
  tmp <- dt %>% group_by(ProHarMeD, Original, Type) %>% summarize(n())
  tmp <- tmp[tmp$ProHarMeD != 0 & tmp$Original != 0,]
  colnames(tmp) <- c("ProHarMeD","Original", "Type", "Intersection")
  wide_tmp <- dcast(tmp, ProHarMeD+Original~Type, value.var = "Intersection")
  wide_tmp[is.na(wide_tmp)] <- 0
  
  summary_tmp <- tmp %>% group_by(ProHarMeD, Original) %>% summarize("Intersection" = sum(Intersection))
  comparison_filtered_plot <- ggplot() + geom_abline() + geom_scatterpie(aes(x=Original, y=ProHarMeD), data=wide_tmp,
                                                                         cols=names(base_colors), pie_scale = 2) + coord_equal() +
    scale_x_continuous(breaks = c(1,2,3,4,5,6,7), limits = c(0.5,7.5)) + scale_y_continuous(breaks = c(1,2,3,4,5,6,7), limits = c(0.5,7.5)) + 
    geom_text(data = summary_tmp, aes(x = Original, y = ProHarMeD -0.5, label = Intersection)) +
    scale_fill_manual(name = "Tissue", values = base_colors)  + labs(y = "ProHarMeD") + theme_bw() + theme(axis.text = element_text(size = 12)) +
    labs(x = "Intersection Size in Original Data", y = "Intersection Size in Harmonized Data")
  plot <- ggarrange(comparison_plot, comparison_filtered_plot, ncol = 2, common.legend = TRUE, legend = "right", labels = c("A", "B"), vjust = 2)
  return(plot)
}
