library(TED)
library(data.table)

set.seed(42)
###----------------------------------------------
# Paths to inputs

cfrna_counts = "./misc_cfrna_ftcounts.tsv"

ref_counts = "./tabsap_subCells-allGenes_thresh-100000.COUNTS.csv"
ref_meta = "./tabsap_subCells-allGenes_thresh-100000.META.csv"


###----------------------------------------------
# Read in data
cat("-----> READING IN DATA \n")
X = read.delim(cfrna_counts,row.names=c(1)) #%>% t() %>% data.frame()
X = X[!grepl("PAR_Y",rownames(X)),1:5]
rownames(X) <- gsub("\\_.*","",rownames(X))
X = t(X)

ref.dat = as.matrix(fread(ref_counts),rownames=1)

meta = fread(ref_meta)
    
###----------------------------------------------
# Clean up reference matrix
cat("-----> CLEANING UP REFERENCE \n")
ref.dat.filtered <- cleanup.genes(ref.dat = ref.dat,
                                   species="hs",
                                   gene.type=c("RB","chrM","chrX","chrY"),
                                   input.type="scRNA",
                                   exp.cells=5)

rm("ref.dat")                  
###----------------------------------------------
# Run deconvolution                          
cat("-----> PERFORMING DECONVOLUTION \n")                 
prevail_bayes_prism <- run.Ted(ref.dat=ref.dat.filtered,
                    X=X,
                    cell.type.labels=meta$cell_type_collapsed,
                    input.type="scRNA",
                    n.cores=15,
                    pdf.name="misc.full")
                  
###----------------------------------------------
# Save Outputs  
cat("-----> SAVING OUTPUT \n")                    
# saveRDS(prevail_bayes_prism,"prevail_bayes-prism.rds")
                  
results = prevail_bayes_prism$res$final.gibbs.theta
write.csv(results,"./prevail_bayes-prism.csv")

