# load libraries
library(readxl)
library(ggplot2)
library(ggpubr)
library(dendsort)
library(ggdendro)
library(ggsci)
library(data.table)
library(viridis)
library(ComplexHeatmap)
library(circlize)

# set graphic theme
theme.cust = theme_bw() +
  theme(axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 1))

theme.cust2 = theme_bw() +
  theme(plot.title = element_text(color = "black", size = 15, face = "bold", hjust = 0.5), 
        axis.title = element_text(color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black", size = 16))

# intersecting terms
## import data
### for dendrogram
int.terms.dend <- read_xlsx("data/term_intersection.xlsx",
                            sheet = "R3") # load from Excel sheet

int.terms.dend.rn = int.terms.dend$term_name  # extract term names for row names
int.terms.dend = int.terms.dend[, 2:6]          # tuncate to relevant data
int.terms.dend <- as.data.frame(int.terms.dend)
rownames(int.terms.dend) = int.terms.dend.rn   # set row names

### main information
int.terms.main <- read_xlsx("data/term_intersection.xlsx",
                            sheet = "R2") # load from Excel sheet

int.terms.main$term_name = factor(int.terms.main$term_name, levels = rev(unique(int.terms.main$term_name))) # make term names factors and reverse for later matching of dendogram in plot
int.terms.main$source = factor(int.terms.main$source, levels = unique(int.terms.main$source)) # make factors to control order of plotting
int.terms.main$type = factor(int.terms.main$type, levels = unique(int.terms.main$type)) # make factors to control order of plotting

### annotations
int.terms.anno <- read_xlsx("data/term_intersection.xlsx",
                            sheet = "R4") # load from Excel sheet

int.terms.anno$term_name = factor(int.terms.anno$term_name, levels = rev(unique(int.terms.anno$term_name))) # make term names factors and reverse for later matching of dendogram in plot
int.terms.anno$source = factor(int.terms.anno$source, levels = unique(int.terms.anno$source)) # make factors to control order of plotting

## plot
### preliminary plot for sanity checks
#### dendogram
#int.terms.dend.dist = dendsort(hclust(dist(int.terms.dend, method = "binary"), method = "ward.D2")) # sorted hierachical Jaccard clustering based on re-identification of terms in tissue types 
#int.terms.dend.dist.od = int.terms.dend.dist$order # extract order of clustered terms for later plotting

#ggdendrogram(int.terms.dend.dist, rotate = TRUE, size = 4) # plot clustering for manual check

### ComplexHeatmap

# Extract data for heatmap matrix
data <- int.terms.main
data$prop <- data$prop * 100
heatmap_data <- dcast(as.data.table(data), term_name~type, value.var = "prop") %>% as.data.frame()
row.names(heatmap_data) <- heatmap_data$term_name
heatmap_data$term_name <- NULL
colnames(heatmap_data) <- c("Cells", "EVs", "ECM", "Bone", "Liquid-Biopsy")

# Clustering
int.terms.dend <- int.terms.dend[rownames(heatmap_data),]
dist_cluster <- dist(int.terms.dend, method = "binary")
clustered_terms <- hclust(dist_cluster, method = "ward.D2")
sorted_clustered_dendro <- dendsort(clustered_terms)

ggdendrogram(sorted_clustered_dendro)
plot(sorted_clustered_dendro)

# Extract data for circles
circle_data <- dcast(as.data.table(data), term_name~type, value.var = "count") %>% as.data.frame()
row.names(circle_data) <- circle_data$term_name
circle_data$term_name <- NULL
colnames(circle_data) <- c("Cells", "EVs", "ECM", "Bone", "Liquid-Biopsy")

circle_data <- circle_data[rownames(heatmap_data),]

# Database source annotation
sources <- int.terms.anno[, c("term_name", "source")]
sources_vec <- sources$source
names(sources_vec) <- sources$term_name

sources_vec <- sources_vec[rownames(heatmap_data)]

sources_colors <- pal_rickandmorty(palette = c("schwifty"), alpha = 1)(length(unique(sources$source)))
names(sources_colors) <- unique(sources$source)

sources_anno <- rowAnnotation(Sources = sources_vec, 
                              col = list(Sources = sources_colors), 
                              gp = gpar(col = "black"), 
                              show_annotation_name = FALSE,
                              annotation_legend_param = list(Sources = list(title = "Database Source")),
                              show_legend = FALSE)

# Total no. proteins annotation
total_proteins <- unique(int.terms.anno[, c("term_name", "count")])
total_proteins_vec <- total_proteins$count
names(total_proteins_vec) <- total_proteins$term_name

total_proteins_vec <- total_proteins_vec[rownames(heatmap_data)]

col_fun_total_proteins <- colorRamp2(c(0, max(total_proteins_vec)), hcl_palette = "YlOrBr", reverse = TRUE)

total_proteins_anno <- rowAnnotation(Total = total_proteins_vec, 
                                     gp = gpar(col = "black"), 
                                     col = list(Total = col_fun_total_proteins), 
                                     show_annotation_name = FALSE,
                                     annotation_legend_param = list(Total = list(title = "Total Nr. Proteins")),
                                     show_legend = FALSE)

# In vitro, in vivo annotation 

assay_type <- sapply(tissue_list, function(x) x$assay)
assay_type[assay_type == "in vivo"] <- "In vivo"
assay_type <- assay_type[colnames(heatmap_data)]

assay_colors <- colors_for_levels[unique(assay_type)]

assay_anno <- HeatmapAnnotation(Assay = assay_type, 
                                col = list(Assay = assay_colors), 
                                gp = gpar(col = "black"),
                                show_annotation_name = FALSE,
                                annotation_legend_param = list(Assay = list(title = "Testing Approach")),
                                show_legend = FALSE)



# Color circles inside cells

# Define custom cell function
custom_cell_fun <- function(j, i, x, y, width, height, fill) {
  radius <- circle_data[i, j] / 2000 # Adjust circle size
  grid.circle(x = x, y = y, r = unit(radius, "npc"), gp = gpar(col = "black", fill = fill, lwd = 2))
}

# Legend for circles
#sizes <- c(0.1, 0.3, 0.6)
# Custom annotation function for the legend
# draw_custom_legend <- function() {
#   grid.newpage()
#   layout_heights <- unit(rep(2.5, length(sizes)), "lines")
#   pushViewport(viewport(layout = grid.layout(length(sizes), 2, heights = layout_heights)))
#   for (i in seq_along(sizes)) {
#     pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
#     grid.circle(x = unit(0.8, "npc"), y = unit(0.5, "npc"), 
#                 r = unit(sizes[i], "cm"), gp = gpar(fill = "white", col = "black", lwd = 2))
#     upViewport()
#     pushViewport(viewport(layout.pos.row = i, layout.pos.col = 2))
#     grid.text(sizes[i], x = unit(0, "npc"), just = "left")
#     upViewport()
#   }
# }

col_fun_heatmap <- colorRamp2(c(0, max(heatmap_data, na.rm = TRUE)), hcl_palette = "viridis", reverse = TRUE)

heatmap <- ComplexHeatmap::Heatmap(as.matrix(heatmap_data),
                        name = "Shared Proteins [%]",
                        col = col_fun_heatmap,
                        right_annotation = total_proteins_anno,
                        left_annotation = sources_anno,
                        bottom_annotation = assay_anno,
                        cluster_rows = clustered_terms,
                        cluster_columns = FALSE,
                        show_row_dend = TRUE,
                        row_dend_reorder = TRUE,
                        na_col = "white",
                        row_names_side = "left",
                        border = TRUE,
                        row_names_max_width = unit(10,"cm"),
                        border_gp = gpar(col = "black"),
                        rect_gp = gpar(col = "black"),
                        cell_fun = custom_cell_fun,
                        show_heatmap_legend = FALSE
                        )

assay_legend <- Legend(labels = unique(assay_type), title = "Approach", legend_gp = gpar(fill = assay_colors, type = "grid"))
sources_legend <- Legend(labels = levels(sources_vec), title = "Database Sources", legend_gp = gpar(fill = sources_colors, type = "grid"))
total_proteins_legend <- Legend(title = "Total Nr. of Proteins", col_fun = col_fun_total_proteins)
heatmap_legend <- Legend(title = "Shared Proteins [%]", col_fun = col_fun_heatmap)


png(paste0(out_dir, "enrichment_heatmap.png"), width=9,height=8.5,units="in",res=1200)

draw(heatmap, annotation_legend_list = list(heatmap_legend, sources_legend, total_proteins_legend, assay_legend), annotation_legend_side = "right")

dev.off()


lgd <- Legend(title = "Nr. Assigned Proteins", 
            labels = c("  20", "  40", "  60"),
            type = "points",
            legend_gp = gpar(col = "black", fill = "white", lwd = 2),
            size = unit(c((20/2000) * 2, (40/2000) * 2, (60/2000)*2), "npc"),
            background = "white", row_gap = unit(0.3, "cm"))

