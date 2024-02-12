library(tidyverse)
library(data.table)
source("../0_support-files//theme_CRP-MISC.R")

annotation <- fread(file="../0_support-files/gencode.biotype.name.key.tsv")

gene_conv <- function(gene,annotation){
    
    if (substr(gene,1,4) == "ENSG"){
        return(annotation[which(annotation$gene_id == gene),]$gene_name)
    }else{
        return(annotation[which(annotation$gene_name == gene),]$gene_id)
    }
}


#------------------------
## META DATA

wb_meta <- read.csv("../1_sample-data/STable7_wbrna-samples.csv") %>% select(PTID,wbrna_sample_id,Diagnosis,timepoint) %>%
    mutate(PTID_TIME = paste0(PTID,"_",timepoint))
cf_meta <- read.csv("../1_sample-data/STable6_cfrna-samples.csv") %>% select(PTID,cfrna_sample_id,Diagnosis,timepoint) %>%
    mutate(PTID_TIME = paste0(PTID,"_",timepoint))

all_meta <- merge(wb_meta,cf_meta, by = c("PTID_TIME"))

#------------------------
## COUNT MATRICES

wb_ftcounts <- read.delim("../1_sample-data/wbrna_ftcounts.txt") %>% 
    rename(gene_id = X) %>% 
    mutate(gene_id = gsub("\\_.*","",gene_id)) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    column_to_rownames(var = "gene_id")

cf_ftcounts <- read.delim("../1_sample-data/cfrna_ftcounts.txt") %>%
    rename(gene_id = Geneid) %>% 
    mutate(gene_id = gsub("\\_.*","",gene_id)) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    column_to_rownames(var = "gene_id")


wb_ftcounts <- wb_ftcounts[,all_meta$wbrna_sample_id]
cf_ftcounts <- cf_ftcounts[,all_meta$cfrna_sample_id]

#------------------------
## MATCH ROW ORDER

wb_ftcounts <- wb_ftcounts[match(rownames(cf_ftcounts),rownames(wb_ftcounts)),]

#------------------------
## FILTER

gene.list <- read.delim("../0_support-files/genelist.hs.tsv",col.names = c("type,","ENSMBL","gene_symbol"))

gene.ids <- gsub("\\..*","",rownames(wb_ftcounts))

exclude.idx <- gene.ids %in% gene.list[,2]

wb_ftcounts = wb_ftcounts[!exclude.idx,] 
cf_ftcounts = cf_ftcounts[!exclude.idx,] 

#------------------------
## NORMALIZE

wb_ftcounts <- edgeR::cpm(wb_ftcounts)
cf_ftcounts <- edgeR::cpm(cf_ftcounts)

#------------------------
## TRANSPOSE

wb_ftcounts <- t(wb_ftcounts)
cf_ftcounts <- t(cf_ftcounts)

#------------------------
## CHECKS

length(colnames(wb_ftcounts)) == length(colnames(cf_ftcounts))
length(rownames(wb_ftcounts)) == length(rownames(cf_ftcounts))

all(colnames(wb_ftcounts) == colnames(cf_ftcounts))


#------------------------
## FUNCTION

calculate_cors <- function(wb_ftcounts,cf_ftcounts){
    
    cors <- list()

    for (i in 1:ncol(wb_ftcounts)){
#     for (i in 1:800){

        df <- data.frame(cf = log(cf_ftcounts[,i]+1),
                        wb = log(wb_ftcounts[,i]+1))
        
        if( (mean(df$cf) > log(10)) & (mean(df$wb) > log(10)) ){ #& (nrow(df) > 50)){
            
            df$keep <- 1
            df$keep <- ifelse( abs(df$cf - mean(df$cf)) > (3*sd(df$cf)),0,1)
            df$keep <- ifelse( abs(df$wb - mean(df$wb)) > (3*sd(df$wb)),0,1) 

            df <- df %>% filter(keep == 1) #%>% filter(!(cf == log(1) & wb == log(1)))
            
            cor_out <- cor.test(df$cf,df$wb, method = "pearson")
            gene = colnames(cf_ftcounts)[i]
            cors[[gene]] <- c(cor_out$estimate, cor_out$p.value)
            
        }
        }
    

    cor_df <- data.frame(do.call("rbind",cors))
    colnames(cor_df) <- c("pearson","p.value")

    cor_df$adj.p <- p.adjust(cor_df$p.value,method = "BH")


    cor_df$gene_name <- lapply(rownames(cor_df),gene_conv,annotation)
    
    return(cor_df)

}

#--------------------------
## CORRELATE ALL

cor_df <- calculate_cors(wb_ftcounts,cf_ftcounts)

keepers <- cor_df %>% rownames()

wb_ftcounts <- wb_ftcounts[,keepers]
cf_ftcounts <- cf_ftcounts[,keepers]


#--------------------------
## SIMULATION by genes
cat("starting simulation 1\n")
set.seed(42)

N = 1000
# N = 10

num_sig <- c()
num_up <- c()
num_down <- c()
total_genes <- c()


for (i in 1:N){

    wb_ftcounts_scram <- wb_ftcounts[sample(1:nrow(wb_ftcounts)),]
    cf_ftcounts_scram <- cf_ftcounts[sample(1:nrow(cf_ftcounts)),]
    
    output <- calculate_cors(wb_ftcounts_scram,cf_ftcounts_scram)
    
    sig <- output %>% filter(adj.p < 0.05)
    up <- output %>% filter(adj.p < 0.05) %>% filter(pearson>0)
    down <- output %>% filter(adj.p < 0.05) %>% filter(pearson<0)
    
    num_sig <- c(num_sig,nrow(sig))
    num_up <- c(num_up,nrow(up))
    num_down <- c(num_down,nrow(down))
    total_genes <- c(total_genes,nrow(output))


}

df = data.frame(num_sig,num_up,num_down,total_genes)
write.table(df,"./simulation_output.rows.tsv",sep="\t",row.names=FALSE)

summary(num_sig)
summary(num_up)
summary(num_down)

#--------------------------
## SIMULATION by samples
cat("starting simulation 2\n")
set.seed(42)

N = 1000
# N = 10

num_sig <- c()
num_up <- c()
num_down <- c()
total_genes <- c()

for (i in 1:N){

    wb_ftcounts_scram <- wb_ftcounts[,sample(1:ncol(wb_ftcounts))]
    cf_ftcounts_scram <- cf_ftcounts[,sample(1:ncol(cf_ftcounts))]
    
    output <- calculate_cors(wb_ftcounts_scram,cf_ftcounts_scram)
    
    sig <- output %>% filter(adj.p < 0.05)
    up <- output %>% filter(adj.p < 0.05) %>% filter(pearson>0)
    down <- output %>% filter(adj.p < 0.05) %>% filter(pearson<0)
    
    num_sig <- c(num_sig,nrow(sig))
    num_up <- c(num_up,nrow(up))
    num_down <- c(num_down,nrow(down))
    total_genes <- c(total_genes,nrow(output))
    
}

df = data.frame(num_sig,num_up,num_down,total_genes)
write.table(df,"./simulation_output.cols.tsv",sep="\t",row.names=FALSE)

summary(num_sig)
summary(num_up)
summary(num_down)