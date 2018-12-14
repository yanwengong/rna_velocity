## taste loomR

## loomR is R package developed to manipulate loom file in R

## This script use "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/mm10_10x_1010.loom"
## to experiance loomR

# Install devtools from CRAN
install.packages("devtools")
# Use devtools to install hdf5r and loomR from GitHub
## installing loomR will run line 13 to install hdf5r automatically if not already 
devtools::install_github(repo = "hhoeflin/hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop", force=TRUE)

library(loomR)

input_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/mm10_10x_1010.loom"

## connect to the loom file in read/write mode
mm10_loomR <- connect(filename = input_path, mode = "r+")


## view
mm10_loomR[["matrix"]]
mm10_loomR$row.attrs$Chromosome
mm10_loomR$col.attrs$CellID
## the one below is not working for some reason
mm10_loomR[["col_attres/CellID"]]

## access
mm10_loomR[["matrix"]][1:5, 1:5]
full.matrix <- mm10_loomR$matrix[, ]
dim(x = full.matrix)
## access cell names
cell.names <- mm10_loomR$col.attrs$CellID


# Print the number of genes
mm10_loomR[["col_attrs/CellID"]]$dims
# Is the number of genes the same as the second dimension (typically
# columns) for the matrix?
mm10_loomR[["col_attrs/CellID"]]$dims == mm10_loomR[["matrix"]]$dims[1]
# For the sake of consistency within the single-cell community, we've
# reversed the dimensions for the `shape` field.  As such, the number of
# genes is stored in `lfile$shape[1]`; the number of cells is stored in the
# second field
mm10_loomR[["col_attrs/CellID"]]$dims == mm10_loomR$shape[2]


# Pull gene expression data for the gene MS4A1. Note that we're using the
# column position for genes
data.cell <- mm10_loomR[["matrix"]][, mm10_loomR$col.attrs$CellID[] == "mm10_outs:AAAGATGTCTGTCAAGx"]
head(x = data.cell)
str(data.cell)           
unique(data.cell)
## why it's 5372x6 ?


data.gene <- mm10_loomR[["matrix"]][, mm10_loomR$row.attrs$Accession[] == "Xkr4"]
head(x=data.gene)


