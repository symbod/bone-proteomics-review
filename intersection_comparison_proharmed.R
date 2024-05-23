########## Install packages ##########

required_packages <- c("rmarkdown", "gprofiler2", "data.table", "ggplot2","tidyr", "dplyr", "readxl", "scatterpie", "ggpubr")
for(package in required_packages){
  if(!require(package,character.only = TRUE, quietly = TRUE)) install.packages(package, dependencies = TRUE, quietly = TRUE)
  library(package, character.only = TRUE, quietly = TRUE)
}

source("config.R")
source("functions/scoringFunctions.R")
source("functions/IntersectionAnalysis.R")

#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

########## Read ProHarMeD Data Table ##########

# read data
harmonized_data <- fread("data/cleaned_protein_list_uhh.csv", sep=",", header = TRUE, na.strings="")
harmonized_data <- as.data.frame(harmonized_data)
rownames(harmonized_data) <- harmonized_data$V1
harmonized_data$V1 <- NULL

# save data info
data_info <- harmonized_data[1:12,]
harmonized_data <- harmonized_data[13:nrow(harmonized_data),]

# Explode data (only 1 gene name per line, and make data numeric) (1 if gene name in data set, 0 if not, rownames = all possible gene names):
harmonized_data_exploded <- explode_data(harmonized_data)

########## Read Original Data Table ##########

# read data
excel_file <- "data/Protein_lists_UHH_v6.xlsx"

# Load the "Gene names" sheet into a data frame with specified column names
original_data <- read_excel(excel_file, sheet = "Gene names", col_names = c("Info",names(harmonized_data)))
original_data <- original_data[,-1]
original_data <- original_data[15:nrow(original_data),]

original_data_exploded <- explode_data(original_data)                            

########## Perform Intersections ##########

perform_intersection <- function(data){
  intersections <- data.frame()
  for(tissue in names(tissue_list)){
    if(weighting == "single") {
      intersections <- rbind(intersections, get_intersection_dt(data = data, data_info = data_info, label = "Tissue", 
                                                                value = c(tissue_list[[tissue]]$label), 
                                                                assay = tissue_list[[tissue]]$assay, 
                                                                type = tissue))
    } else if(weighting == "double") {
      intersections <- rbind(intersections, get_intersection_dt_double_weighting(data = data, data_info = data_info, label = "Tissue", 
                                                                                 value = c(tissue_list[[tissue]]$label), 
                                                                                 assay = tissue_list[[tissue]]$assay, 
                                                                                 type = tissue))
    } 
  }
  return(intersections)
}

harmonized_intersections <- perform_intersection(harmonized_data_exploded)
original_intersections <- perform_intersection(original_data_exploded)
# remove 9_human -> no gene names provided intersection
original_intersections <- original_intersections[original_intersections$Gene != "no Gene names provided",]

########## Filter Intersections ##########

# Group by 'Type' and apply select_intersections() function to each group
intersections <- harmonized_intersections[harmonized_intersections$Intersected != 1,]

resulting_data_frames <- intersections %>%
  group_by(Type) %>%
  do(select_intersections(., intersection_thr=intersection_thr, weighting_thr=weighting_thr)) %>%
  ungroup()

# Combine resulting data frames into one
harmonized_intersections_filtered <- bind_rows(resulting_data_frames)


########## Compare Number of Intersected Studies ##########

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

#ggplot(dt, aes(x = Original, y = ProHarMeD, color = Type)) + geom_jitter() + 
#  scale_x_continuous(breaks = seq(0,10)) + scale_y_continuous(breaks = seq(0,10)) + 
#  scale_color_manual(name = "Tissue", values = base_colors) + geom_abline() + 
#  geom_rect(aes(xmin = 1 - 0.5, xmax = 1 + 0.5, ymin = 1 - 0.5, ymax = 1 + 0.5),color = "black", alpha = 0.2, fill = NA) +
#  geom_rect(aes(xmin = 2 - 0.5, xmax = 2 + 0.5, ymin = 2 - 0.5, ymax = 2 + 0.5),color = "black", alpha = 0.2, fill = NA) +
#  geom_rect(aes(xmin = 3 - 0.5, xmax = 3 + 0.5, ymin = 3 - 0.5, ymax = 3 + 0.5),color = "black", alpha = 0.2, fill = NA) +
#  geom_rect(aes(xmin = 4 - 0.5, xmax = 4 + 0.5, ymin = 4 - 0.5, ymax = 4 + 0.5),color = "black", alpha = 0.2, fill = NA) +
#  geom_rect(aes(xmin = 5 - 0.5, xmax = 5 + 0.5, ymin = 5 - 0.5, ymax = 5 + 0.5),color = "black", alpha = 0.2, fill = NA) +
#  geom_rect(aes(xmin = 6 - 0.5, xmax = 6 + 0.5, ymin = 6 - 0.5, ymax = 6 + 0.5),color = "black", alpha = 0.2, fill = NA) +
#  geom_rect(aes(xmin = 7 - 0.5, xmax = 7 + 0.5, ymin = 7 - 0.5, ymax = 7 + 0.5),color = "black", alpha = 0.2, fill = NA) #

# x <- 1,2,3,4, ..
# y <- 1,2,3,4, ..
tmp <- dt %>% group_by(ProHarMeD, Original) %>% summarize(n())
colnames(tmp) <- c("ProHarMeD", "Original", "GeneCount")
ggplot(tmp, aes(x = Original, y = ProHarMeD, fill = GeneCount, label = GeneCount)) + geom_tile(color = "white") +
  scale_x_continuous(breaks = seq(0,10)) + scale_y_continuous(breaks = seq(0,10)) + geom_text(color = "white")

ggplot(tmp[!(tmp$Original == 1 & tmp$ProHarMeD == 1),], aes(x = Original, y = ProHarMeD, fill = GeneCount, label = GeneCount)) + geom_tile(color = "white") +
  scale_x_continuous(breaks = seq(0,10)) + scale_y_continuous(breaks = seq(0,10)) + geom_text(color = "white") 

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
ggsave("/Users/lisiarend/Desktop/comparison_plot.png", width = 8, height = 6)

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
ggsave("/Users/lisiarend/Desktop/comparison_plot_filtered.png", width = 8, height = 6)

ggarrange(comparison_plot, comparison_filtered_plot, ncol = 2, common.legend = TRUE, legend = "right", labels = c("A", "B"), vjust = 2)
ggsave("/Users/lisiarend/Desktop/comparison_plot_combined.png", width = 8, height = 4)


  
