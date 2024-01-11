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


#' Create Top Genes Heatmap
#'
#' Identifies the top 10 genes based on each ranking metric in a scoring dataframe, creates a union of these genes, 
#' extracts their ranking data, and generates a heatmap displaying the ranks.
#'
#' @param scoring_df A dataframe with gene IDs and their ranks across various metrics.
#' @param output_file Optional file path to save the heatmap.
#' @return A ggplot object representing the heatmap.
#' @importFrom ggplot2 ggplot geom_tile geom_text scale_fill_gradient theme_minimal labs theme
#' @importFrom reshape2 melt
#' @examples
#' network_scoring <- score_gene_network(g = core_ami_ds$ds_graph, attribute_name = "name")
#' enrich_scoring <- score_genes_enrichment(genes = core_ami_ds$ds_genes)
#' create_top_genes_heatmap(network_scoring)
#' @export
create_top_genes_heatmap <- function(scoring_df, output_file = NULL) {
  
  top_genes_data <- extract_top_genes_union(scoring_df = scoring_df)$top_genes_data
  
  # Melt the data for heatmap
  melted_data <- melt(top_genes_data, id.vars = "IDs")
  
  # Calculate the sum of ranks for each gene
  sum_ranks <- aggregate(. ~ IDs, melted_data, function(x) sum(x, na.rm = TRUE))
  
  # Order the genes by sum of ranks
  ordered_genes <- sum_ranks[order(sum_ranks$value), "IDs"]
  
  # Reorder the melted data based on the sum of ranks
  melted_data$IDs <- factor(melted_data$IDs, levels = ordered_genes)
  
  # Create the heatmap with values
  p <- ggplot(melted_data, aes(x = IDs, y = variable, fill = value)) + 
    geom_tile() +
    geom_text(aes(label = round(value, 2)), size = 3, vjust = 1) + # Add text labels
    scale_fill_gradient(low = "darkgreen", high = "grey") +
    theme_minimal() + 
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white")
    ) +
    labs(title = "Top Gene Ranks Heatmap", x = "Gene ID", y = "Rank Metric", fill = "Rank") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 0, hjust = 1))
  
  # Optionally save the heatmap
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 10, height = 6, units = "in", dpi = 300)
  }
  
  return(p)
}

#' Create Combined Heatmap from Network and Enrichment Scorings
#'
#' Combines top genes from network and enrichment scorings into a single heatmap. 
#' This function extracts top genes from both network and enrichment scoring dataframes, 
#' merges their data, and then creates a heatmap showing the ranks of these genes.
#'
#' @param network_scoring A scoring dataframe obtained from network analysis.
#' @param enrich_scoring A scoring dataframe obtained from enrichment analysis.
#' @param output_file Optional file path to save the heatmap.
#' @return A ggplot object representing the combined heatmap.
#' @importFrom ggplot2 ggplot geom_tile geom_text scale_fill_gradient theme_minimal labs theme
#' @importFrom reshape2 melt
#' @examples
#' network_scoring <- score_gene_network(g = core_ami_ds$ds_graph, attribute_name = "name")
#' enrich_scoring <- score_genes_enrichment(genes = core_ami_ds$ds_genes)
#' create_combined_heatmap(network_scoring, enrich_scoring)
#' @export
create_combined_heatmap <- function(network_scoring, enrich_scoring, output_file = NULL) {
  # Extract top genes from both scorings
  get_top_genes <- function(scoring_df) {
    top_genes <- lapply(grep("_Rank$", names(scoring_df), value = TRUE), function(metric) {
      head(sort(scoring_df[[metric]]), 10)
    })
    unique(unlist(top_genes))
  }
  
  top_genes_network <- extract_top_genes_union(scoring_df = network_scoring)
  top_genes_enrich <- extract_top_genes_union(scoring_df = enrich_scoring)
  
  # Combine and deduplicate the list of genes
  all_top_genes <- unique(c(top_genes_network$top_genes_union, top_genes_enrich$top_genes_union))
  
  # Merge the scoring data for these genes
  merged_data <- cbind(
    network_scoring[network_scoring$IDs %in% all_top_genes, ],
    enrich_scoring[enrich_scoring$IDs %in% all_top_genes, ] )
  
  cols <- union(names(top_genes_network$top_genes_data), names(top_genes_enrich$top_genes_data))
  merged_data <- merged_data[, ..cols]
  
  # Melt the data for the heatmap
  melted_data <- melt(merged_data, id.vars = "IDs")
  
  # Calculate the sum of ranks for each gene
  sum_ranks <- aggregate(. ~ IDs, melted_data, function(x) sum(x, na.rm = TRUE))
  
  # Order the genes by sum of ranks
  ordered_genes <- sum_ranks[order(sum_ranks$value), "IDs"]
  
  # Reorder the melted data based on the sum of ranks
  melted_data$IDs <- factor(melted_data$IDs, levels = ordered_genes)
  
  # Create the heatmap
  p <- ggplot(melted_data, aes(x = IDs, y = variable, fill = value)) + 
    geom_tile() +
    geom_text(aes(label = round(value, 2)), size = 3, vjust = 1) + 
    scale_fill_gradient(low = "darkgreen", high = "grey") +
    theme_minimal() + 
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white")
    ) +
    labs(title = "Combined Top Gene Ranks Heatmap", x = "Gene ID", y = "Rank Metric", fill = "Rank") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 0, hjust = 1))
  
  # Optionally save the heatmap
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 10, height = 6, units = "in", dpi = 300)
  }
  
  return(p)
}
