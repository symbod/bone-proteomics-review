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
save_results <- function(ds_result, filename, out_dir){
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

#' Extract Gene Categories
#'
#' This function extracts various categories of genes based on operations like set differences and intersections.
#'
#' @param seeds Vector containing seeds.
#' @param ds_results Dataframe containing core drugs and associated genes.
#' @param filtered_intersections Dataframe containing filtered intersections.
#' @param all_genes Vector with all listed genes..
#' @param nedrex_database (Optional) Dataframe containing Nedrex database information.
#' @return A list containing different categories of genes.
#' @examples
#' extract_gene_categories(seeds, ds_results, filtered_intersections, all_genes, nedrex_database)
#' @export

extract_gene_categories <- function(seeds, ds_results, filtered_intersections, all_genes, nedrex_database = NULL) {
  drugs <- names(ds_results$ds_drugs)
  exception_genes <- setdiff(setdiff(ds_results$ds_gene, seeds), drugs)
  
  # Uncomment and use these lines if nedrex_database and type are provided and needed in your analysis
  # not_in_network <- if (!is.null(nedrex_database)) setdiff(seeds, nedrex_database) else NULL
  # lost_in_final_network <- if (!is.null(type)) setdiff(setdiff(seeds, ds_results$ds_gene), not_in_network) else NULL
  
  new_genes <- setdiff(exception_genes, all_genes)
  exception_genes_in_filtered_intersection <- intersect(exception_genes, filtered_intersections$Gene)
  exception_genes_filtered_out <- setdiff(setdiff(exception_genes, exception_genes_in_filtered_intersection), new_genes)
  
  return(list(
    seeds = seeds,
    drugs = drugs,
    exception_genes = exception_genes,
    # not_in_network = not_in_network,  # Uncomment if needed
    # lost_in_final_network = lost_in_final_network,  # Uncomment if needed
    new_genes = new_genes,
    exception_genes_in_filtered_intersection = exception_genes_in_filtered_intersection,
    exception_genes_filtered_out = exception_genes_filtered_out
  ))
}

#' Prepare Node Attributes and Types for Graph
#'
#' This function prepares node attributes and types for a given graph. It categorizes nodes based on
#' various gene categories obtained from 'extract_gene_categories' function.
#'
#' @param graph An igraph object representing the network graph.
#' @param gene_categories A list of gene categories obtained from 'extract_gene_categories' function.
#' @return A list containing two elements: 'attributes' which are the node attributes and 'types' which are the node types for each vertex in the graph.
#' @examples
#' # Assuming 'my_graph' is an igraph object and 'gene_cats' is the output from 'extract_gene_categories'
#' prepare_node_attributes(my_graph, gene_cats)
#' @export
prepare_node_attributes <- function(graph, gene_categories) {
  vertex_names <- igraph::V(graph)$name
  seeds <- gene_categories$seeds
  drugs <- gene_categories$drugs 
  new_genes <- gene_categories$new_genes
  exception_genes_in_filtered_intersection <- gene_categories$exception_genes_in_filtered_intersection
  
  seed_att <- ifelse(vertex_names %in% seeds, "seed",
                     ifelse(vertex_names %in% drugs, "drug",
                            ifelse(vertex_names %in% new_genes, "new",
                                   ifelse(vertex_names %in% exception_genes_in_filtered_intersection, 
                                          "significant only in single tissue", "filtered out for single tissue"))))
  
  node_type <- ifelse(vertex_names %in% seeds, "circle",
                      ifelse(vertex_names %in% drugs, "diamond", "triangle"))
  
  list(attributes = seed_att, types = node_type)
}


draw_pie_chart_network <- function(graph, attributes, all_pie_values, pie_column, pie_colors, label, seed = 123, legend_position = "topleft",  vertex.size=3, vertex.label.size = 1.5, vertex.label.dist=0.35){
  # make sure that in attributes only genes that are also in graph
  attributes <- attributes[attributes$name %in% V(graph)$name,]
  # make sure graph nodes and attribute data frame are in the same order
  attributes  <- attributes[match(V(graph)$name, attributes$name),]
  
  # construct value list for pie-charts
  values <- lapply(attributes$name, function(gene){
    pie_values <- attributes[attributes$name == gene,][[pie_column]]
    pie_values <- strsplit(pie_values, ";")[[1]]
    pie_values <- as.numeric(all_pie_values %in% pie_values)
    return(pie_values)
  })
  
  # set seed
  # plot with ggnet2
  set.seed(seed)
  p <- ggnet2(graph, size=3,label=TRUE, label.size = 2) 
  
  # extract coordinates of nodes
  coordinates <- p$data
  coordinates <- coordinates[, c("label", "x", "y")]
  # order
  coordinates <- coordinates[match(attributes$name, coordinates$label),]
  
  # set pie colors
  V(graph)$pie.color=list(pie_colors)
  
  # built layout
  layout <- as.matrix(coordinates[, c("x", "y")])
  
  # see how it works
  par(mar = c(0, 0, 0, 0)) 
  
  # labeling
  if(label){
    label <- V(graph)$name
  } else {
    label <- NA
  }
  
  # I think this tells igraph to normalize the coordinates of the 
  # layout relative to the space you're working with
  # see how it works
  plot(graph,
       vertex.shape="pie", 
       vertex.pie=values,
       vertex.label=label,
       vertex.size=vertex.size,#3
       vertex.label.family = "Helvetica",
       vertex.label.font = 1,
       vertex.label.dist = vertex.label.dist, #0.35,
       vertex.label.degree = pi/2,
       vertex.label.color = "black",
       vertex.label.cex=vertex.label.size,#1.5,
       edge.width = 2,
       rescale = TRUE,
       asp = 0,
       layout = layout
  )
  legend(legend_position, legend = all_pie_values, pch = 16, col = pie_colors, bty = "n", cex = 3)
}