
read_meta_study_data <- function(file_path, sheet){
  
  # Read data as data.table
  original_data <- read_excel(file_path, sheet = sheet, col_names = FALSE) 
  original_data <- as.data.table(original_data)
  
  # Extract meta information of the studies
  data_info <- as.data.frame(original_data[1:12,]) # TODO: automatically detection
  
  # Construct naming vector for studies
  study_species <- as.character(as.vector(data_info[data_info$...1 == "Species",-1]))
  study_names <- paste0(seq(1, length(study_species)), "_",study_species)
  
  colnames(data_info) <- c("Info", study_names)
  rownames(data_info) <- data_info$Info
  data_info$Info <- NULL
  
  # Separate data from meta-info
  original_data <- original_data[13:nrow(original_data),-1]
  colnames(original_data) <- study_names
  
  return(list("data" = original_data, "meta-info" = data_info))
}
