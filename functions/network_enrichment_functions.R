#' Run DrugStone AMI Analysis
#'
#' This function runs the DrugStone analysis with the AMI algorithm and returns genes, edges, and a graph object.
#'
#' @param genes Vector of genes for analysis.
#' @return A list containing genes, graph object, and edges from the DrugStone AMI analysis.
#' @examples
#' run_drugstone_ami(c("gene1", "gene2"))
#' @export
#' 
run_drugstone_ami <- function(genes){
  parameters = list("algorithm" = "multisteiner", "target" = "drug-target", "tolerance" = 5, "hubPenalty" = 0.5, "num_trees" = 5, "ppiDataset"= 'nedrex')
  ds_out <- ds$new_task(genes, parameters)
  ds_results <- ds_out$get_result()
  ds_genes <- ds_results$get_genes()
  ds_edges <- do.call(rbind, lapply(ds_results$get_edges(), function(x) { data.frame(from = x$from, to = x$to) }))
  ds_graph <- graph_from_edgelist(as.matrix(ds_edges), directed = FALSE)
  return(list("ds_genes" = names(ds_genes), "ds_graph" = ds_graph, "ds_edges" = ds_edges))
}

#' Run DrugStone Drug Analysis
#'
#' This function runs the DrugStone analysis for drugs and returns drugs, genes, and a graph object.
#'
#' @param genes Vector of genes for analysis.
#' @param edges Dataframe of edges for graph construction.
#' @return A list containing drugs, genes, and a graph object from the DrugStone analysis.
#' @examples
#' run_drugstone_drugs(c("gene1", "gene2"), my_edges)
#' @export
run_drugstone_drugs <- function(genes, edges){
  ds_out_2 <- ds$new_task(genes, list("algorithm" = "trustrank", "target" = "drug","pdiDataset"= 'nedrex', "includeNonApprovedDrugs" = TRUE))
  ds_results_2 <- ds_out_2$get_result()
  ds_drugs_2 <- ds_results_2$get_drugs()
  ds_genes <- ds_results_2$get_genes()
  ds_edges_2 <- do.call(rbind, lapply(ds_results_2$get_edges(),  function(x) { data.frame(from = x$from, to = x$to) }))
  ds_edges_2 <- rbind(ds_edges_2, edges)
  ds_graph_2 <- graph_from_edgelist(as.matrix(ds_edges_2), directed = FALSE)
  return(list("ds_drugs" = ds_drugs_2, "ds_gene" = names(ds_genes), "ds_graph" = ds_graph_2))
}

#' Save DrugStone Analysis Results
#'
#' Saves the results of DrugStone analysis including drugs, gene lists, and graph data to files.
#'
#' @param ds_result The result object from DrugStone analysis.
#' @param filename The base filename for the output files.
#' @param out_dir Directory where output files will be saved.
#' @examples
#' save_results(ds_result, "my_analysis")
#' @export
save_results <- function(ds_result, filename){
  if ( "ds_drugs" %in% names(ds_result) ) {
    ds_drugs <- do.call(rbind, lapply(ds_result$ds_drugs, as.data.frame))
    ds_drugs <- ds_drugs %>% group_by(across(-hasEdgesTo)) %>%
      summarize(hasEdgesTo = paste(hasEdgesTo, collapse = ";"))
    write.table(ds_drugs, paste0(out_dir, filename, "_drugs.csv"), quote=FALSE, sep=",", row.names = FALSE)
  }
  if ( "ds_gene" %in% names(ds_result) ) {
    ds_gene <- ds_result$ds_gene
    write.table(ds_gene, paste0(out_dir, filename, "_genelist.csv"), quote=FALSE, sep=",", row.names = FALSE)
  }
  if ( "ds_graph" %in% names(ds_result) ) {
    ds_graph <- ds_result$ds_graph
    # save as graphml
    write.graph(ds_graph, paste0(out_dir, filename, "_graph.graphml"), format = "graphml")
    write.table(get.data.frame(ds_graph), paste0(out_dir, filename, "_graph.csv"), quote=FALSE, sep=",", row.names = FALSE)
  }
}

