library(readr)
sc_markers_6oct25 <- read_csv("~/scRNA/Diane/sc_markers_6oct25.csv")
library(Biobase)
sc_markers = sc_markers_6oct25
e_this <- readRDS("~/scRNA/Diane/e_this_7oct25.rds")
sc_genes <- sc_markers$names
bulk_genes <- featureNames(e_this)
matching_genes <- intersect(sc_genes, bulk_genes)
print(matching_genes)
filtered_sc_markers <- sc_markers[sc_markers$names %in% matching_genes, ]
filtered_e.this <- e_this[matching_genes, ]
R2cutoff <- 0.1
specificity_cutoff <- 0.1
markers <- list()
celltypes <- levels(as.factor(filtered_sc_markers$group))
for (i in seq_along(celltypes)) {
markers[[i]] <- subset(
filtered_sc_markers,
group == celltypes[i] &
pvals_adj > specificity_cutoff
)$names
}
names(markers) <- celltypes
save_path <- "~/scRNA/Diane/markers_7oct25.rds"
saveRDS(markers, file = save_path)
View(markers)
library(dplyr)
max_length <- max(sapply(markers, length))
# Create a data frame by padding the shorter vectors with NA
markers_df <- data.frame(lapply(markers, function(x) {
length(x) <- max_length
return(x)
}), check.names = FALSE)
View(markers_df)
write.csv(markers_df, file = "/home/brassl1/scRNA/Diane/markers_7oct25.csv", row.names = FALSE)
all_genes <- unlist(markers)
unique_genes <- unique(all_genes)
print(unique_genes)
write.csv(unique_genes, file = "/home/brassl1/scRNA/Diane/unique_genes_7oct25.csv", row.names = FALSE)