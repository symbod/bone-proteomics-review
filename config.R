# Configuration for Bone Biomarker Review

# Setup file names
original_study_file <- "data/Protein_lists_UHH_v6.xlsx" # TODO: remove two first rows


# Setup overall

tissue_list <- list("EVs"= list("label" = "EVs",
                                "assay" = "In vitro"),
                    "Cells"= list("label" = "Cells",
                                  "assay" = "In vitro"), 
                    "ECM"= list("label" = "ECM",
                                "assay" = "In vitro"), 
                    "Bone"= list("label" = "Bone",
                                 "assay" = "In vivo"), 
                    "Liquid-Biopsy" = list("label" = c("Serum-Plasma", "Plasma-EVs"),
                                              "assay" = "In vivo"))

out_dir <- "/home/kikky/Projects/symbod/bone-proteomics-review//data/results/"

colors_for_levels <- c(
  "Human" = "#66A61E",
  "Rat" = "#B2DF8A",
  "Mouse" = "#A6CEE3",
  "Rabbit" = "#1F78B4",
  "In vitro" = "#FFD92F",
  "In vivo" = "#A6761D",
  "EVs" = "#A3A500",
  "Cells" = "#F8766D",
  "Liquid-Biopsy" = "#E76BF3",
  "Bone" = "#00B0F6",
  "ECM" = "#00BF7D"
)


# Setup for Intersection Analysis 

IA_rerun = TRUE # Set on false, if previously calculated results should be taken

weighting = "single" # alternative options: "double"

intersection_thr=2
weighting_thr="normal"

# Setup for Network Analysis 

input_network = "nedrex" # alternative options: string, iid, biogrid

