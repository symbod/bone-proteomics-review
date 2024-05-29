#' Score Gene Network
#'
#' Computes various centrality measures and their ranks for nodes in a network graph. 
#' It calculates closeness, betweenness, degree, harmonic centrality, and PageRank scores. Optionally, renames IDs based on a specified vertex attribute.
#'
#' @param g An igraph object representing the network graph.
#' @param attribute_name Optional name of the vertex attribute to rename the IDs in the output. If NULL, default vertex names are used.
#' @return A data.table containing network scores and ranks for each node.
#' @examples
#' # Assuming 'my_graph' is an igraph object
#' score_gene_network(my_graph, attribute_name = "name")
#' @importFrom data.table setDT
#' @export
score_gene_network <- function(g, attribute_name = NULL) {
  # Calculating centrality measures
  centrality_measures <- list(
    closeness = igraph::closeness(g),
    betweenness = igraph::betweenness(g),
    degree = igraph::degree(g),
    harmonic_centrality = igraph::harmonic_centrality(g),
    page_rank = igraph::page_rank(g)$vector
  )
  
  # Convert to data frame and add IDs
  network_scores <- as.data.frame(centrality_measures)
  network_scores$IDs <- if (is.null(attribute_name)) rownames(network_scores) else vertex_attr(g, attribute_name)
  
  # Add ranks for each centrality measure
  sources <- names(centrality_measures)
  for (source in sources) {
    network_scores[[paste0(source, "_Rank")]] <- rank(-network_scores[[source]], ties.method = "average")
  }
  
  return(setDT(network_scores))
}

#' Scoring of gene names based on enrichment terms
#'
#' @param genes character vector of gene names (gene symbols)
#' @param organism organism name (for gProfiler)
#' @param sources vector of data sources to sue (per default: GO:BP, GO:CC, GO:MF, KEGG, REAC)
#'
#' @return data table with scores and ranks for each gene name
#' @export
#'
#' @examples
score_genes_enrichment <- function(genes, organism = "hsapiens", sources = c("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC")){
  enrich <- gost(query = genes, organism = organism, ordered_query = FALSE, significant = TRUE, evcodes=TRUE)$result
  scoring <- NULL
  # scoring sources 
  for(source in sources){
    res <- enrich[enrich$source == source,]
    res <- res[, c("term_id", "term_name", "intersection")]
    res <- separate_rows(res, "intersection", sep=",")
    scores <- res %>% group_by(intersection) %>% count()
    # add missing gene names
    missing <- genes[!genes %in% scores$intersection]
    missing <- data.table(intersection = missing, n = rep(0, length(missing)))
    scores <- rbind(scores, missing)
    scores$rank <- rank(-scores$n, ties.method = "average")
    colnames(scores) <- c("IDs", source, paste0(source, "_Rank"))
    if(is.null(scoring)){
      scoring <- scores
    } else {
      scoring <- merge(scoring, scores, by = "IDs", all = TRUE)
    }
  }
  return(setDT(scoring))
}


get_top_drugs <- function(genes, number){
  parameters = list("target"= "drug",
                    "algorithm"= "trustrank",
                    "pdiDataset"= 'nedrex',
                    "includeNonApprovedDrugs" = TRUE,
                    "resultSize"= number)
  result <- ds$new_task(seeds = unique(genes), parameters = parameters)
  result <- result$get_result()
  drugs <- as.data.frame(do.call(rbind, result$get_drugs()))
  drugs$hasEdgesTo <- lapply(drugs$hasEdgesTo, FUN=function(x) paste(x, collapse = ","))
  return(drugs)
}

#' Extract Union of Top Genes and Their Rank Data
#'
#' Identifies the top 10 genes based on each ranking metric in a scoring dataframe and creates a union of these genes.
#' Additionally, extracts the rank data for these top genes.
#'
#' @param scoring_df A dataframe with gene IDs and their ranks across various metrics.
#' @return A list containing a union of top genes and their corresponding rank data.
#' @examples
#' network_scoring <- score_gene_network(g = core_ami_ds$ds_graph, attribute_name = "name")
#' top_genes_info <- extract_top_genes_union(network_scoring)
#' @export
extract_top_genes_union <- function(scoring_df) {
  # Identify top 10 genes for each rank
  top_genes_list <- list()
  ranking_metrics <- grep("_Rank$", names(scoring_df), value = TRUE)
  for (metric in ranking_metrics) {
    top_genes_this_metric <- scoring_df[order(scoring_df[[metric]]), "IDs"][1:10]
    top_genes_list[[metric]] <- top_genes_this_metric
  }
  
  # Create a union of these genes
  top_genes_union <- unique(unlist(top_genes_list))
  
  # Extract rank data for these genes
  use_columns <- c("IDs", ranking_metrics)
  top_genes_data <- scoring_df[scoring_df$IDs %in% top_genes_union, ..use_columns]
  
  return(list(top_genes_union = top_genes_union, top_genes_data = top_genes_data))
}

calculate_drug_connections <- function(drugs, gene_ids) {
  # Split the 'hasEdgesTo' column and count occurrences
  gene_counts <- strsplit(as.character(drugs$hasEdgesTo), ",") %>%
    unlist() %>%
    table() %>%
    as.data.frame()
  
  # Rename columns for clarity
  names(gene_counts) <- c("IDs", "drugs")
  
  # Create a complete list of gene_ids with 0 counts for those not in gene_counts
  complete_gene_list <- data.frame(IDs = gene_ids, drugs = 0)
  gene_counts <- merge(complete_gene_list, gene_counts, by = "IDs", all.x = TRUE)
  gene_counts$drugs <- rowSums(gene_counts[, c("drugs.x", "drugs.y")], na.rm = TRUE)
  
  # Drop the extra columns after merging
  gene_counts <- gene_counts[, c("IDs", "drugs")]
  
  # Add ranking based on the number of drugs connected
  gene_counts$drugs_Rank <- rank(-gene_counts$drugs, ties.method = "average")
  
  return(gene_counts)
}

check_vertex <- function(vertex, core_graph, combined_graph) {
  if (V(core_graph)[vertex]$node_attributes == "new") {
    return("new")
  } else if (vertex %in% V(combined_graph)$name) {
    return(V(combined_graph)[vertex]$source)
  } else {
    return("problem")
  }
}

get_gene_drug_clustering <- function(drugs, all_ranks){
  # 1. Parse the 'drugs' DataFrame
  drugs_list <- strsplit(as.character(drugs$hasEdgesTo), ",")
  names(drugs_list) <- drugs$label
  
  # 2. Create Binary Matrix for Gene-Drug Associations
  gene_drug_matrix <- do.call(cbind, lapply(drugs_list, function(drug_genes) {
    +(all_ranks$IDs %in% drug_genes)
  }))
  rownames(gene_drug_matrix) <- all_ranks$IDs
  colnames(gene_drug_matrix) <- names(drugs_list)
  
  # 3. Calculate Distance Matrix
  # Here, Jaccard distance is used as an example. You can choose other metrics as appropriate.
  gene_distance_matrix <- as.dist(1 - proxy::simil(gene_drug_matrix, method = "jaccard"))
  
  # 4. Cluster the Genes
  gene_clustering <- hclust(gene_distance_matrix)
  return(gene_clustering)
}

get_median_rank_data <- function(scoring_dt, rank_string = "Rank"){
  rank_colnames <- colnames(scoring_dt)[str_detect(colnames(scoring_dt), "_Rank")]
  ranks <- scoring_dt[, c("IDs", rank_colnames), with = FALSE]
  ranks_long <- melt(ranks, id.vars = "IDs", value.name = "Rank", variable.name = "Measure")
  ranks_median <- ranks_long %>% group_by(IDs) %>% summarize(Rank= median(Rank)) %>% as.data.table()
  colnames(ranks_median) <- c("IDs", rank_string)
  return(ranks_median)
}

generate_combined_ranking_dt <- function(network_scoring, enrich_scoring, access_scoring, boneabund_scoring){
  # get median over all network ranks
  network_medians <- get_median_rank_data(network_scoring, rank_string = "Network_Rank")
  enrich_medians <- get_median_rank_data(enrich_scoring, "Enrich_Rank")
  access_medians <- get_median_rank_data(as.data.table(access_scoring), "Accessibility_Rank")
  boneabund_medians <- get_median_rank_data(as.data.table(boneabund_scoring), "BoneAbundance_Rank")
  
  # merge scores
  all_ranks <- merge(enrich_medians, network_medians, by = "IDs") %>% as.data.frame()
  all_ranks <- merge(all_ranks, access_medians, by = "IDs") %>% as.data.frame()
  all_ranks <- merge(all_ranks, boneabund_medians, by = "IDs") %>% as.data.frame()
  
  return(all_ranks)
}

generate_tissue_assay_file_anno_vec <- function(tissue_annotation, tissue_list, ids_order, colors_for_levels, out_dir){
  # create temp directory for saving pie chart
  out_tmp_tissue <- paste0(out_dir, "/scoring/temp_pie_charts_tissue")
  out_tmp_assay <- paste0(out_dir, "/scoring/temp_pie_charts_assay")
  dir.create(out_tmp_tissue)
  dir.create(out_tmp_assay)
  
  sapply(names(tissue_annotation[!tissue_annotation %in% c("new", "problem")]), function(gene){
    tissues <- strsplit(tissue_annotation[gene], ";")[[1]]
    data <- data.frame(tissue = tissues, value = rep(1, length(tissues)))
    data$assay <- sapply(data$tissue, function(x) tissue_list[[x]]$assay)
    # tissue piechart
    ggplot(data, aes(x = "", y = value, fill = tissue)) + geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) + scale_fill_manual(values = colors_for_levels[tissues]) + theme_void() +   theme(legend.position = "none")
    ggsave(paste0(out_tmp_tissue, "/", gene, ".eps"), width = 1, height = 1, device = "eps")
    
    # assay piechart
    ggplot(unique(data[, c("value", "assay")]), aes(x = "", y = value, fill = assay)) + geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) + scale_fill_manual(values = colors_for_levels[unique(data$assay)]) + theme_void() +   theme(legend.position = "none")
    ggsave(paste0(out_tmp_assay, "/", gene, ".eps"), width = 1, height = 1, device = "eps")
  })
  
  pie_charts_tissue_files <- list.files(out_tmp_tissue, full.names = TRUE)
  pie_charts_assay_files <- list.files(out_tmp_assay, full.names = TRUE)
  
  # rearrange vector for problem and new 
  names(pie_charts_tissue_files) <- sapply(list.files(out_tmp_tissue), function(x) strsplit(x, ".eps")[[1]])
  names_new_problem_genes <- ids_order[! ids_order %in% names(tissue_annotation)]
  new_problem_genes <- rep("", length(names_new_problem_genes))
  names(new_problem_genes) <- names_new_problem_genes
  
  tissue_anno_vec <- c(pie_charts_tissue_files, new_problem_genes)
  tissue_anno_vec <- tissue_anno_vec[ids_order]
  
  # assay annotation
  names(pie_charts_assay_files) <- sapply(list.files(out_tmp_assay), function(x) strsplit(x, ".eps")[[1]])
  assay_anno_vec <- c(pie_charts_assay_files, new_problem_genes)
  assay_anno_vec <- assay_anno_vec[ids_order]
  return(list("tissue" = tissue_anno_vec, "assay" = assay_anno_vec))
}

