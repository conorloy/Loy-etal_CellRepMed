suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(dendextend))
suppressMessages(library(reshape2))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gridExtra))
source("../0_support-files//theme_CRP-MISC.R")

###----------------------------------

PVAL_THRESH = 0.05
LOGFC_THRESH = 0.01

###----------------------------------

annotation <- fread(file="../0_support-files/gencode.biotype.name.key.tsv")

#----------------------------------------------------------------------------------------------
# FUNCTIONS 
#----------------------------------------------------------------------------------------------

colLab<<-function(n,metadata){
    if(is.leaf(n)){
        
#         print(attributes(n)$label)
        
        #I take the current attributes
        a=attributes(n)
        
        #I deduce the line in the original data, and so the treatment and the specie.
        ligne=match(attributes(n)$label,rownames(metadata))
        Diagnosis=metadata[ligne,]$Diagnosis;
            if(Diagnosis=="COVID-19"){col_diag="#F0484E"};if(Diagnosis=="MIS-C"){col_diag="#5CB2EB"};if(Diagnosis=="Control_Non-inflammatory"){col_diag="#FBE77C"};if(Diagnosis=="Control_Inflammatory"){col_diag="purple"}
        #Modification of leaf attribute
        attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=15,col=col_diag,lab.col=col_diag,lab.font=1,lab.cex=1))
        }
    return(n)
}


get_gene_id <- function(X,annotation){
    return(gsub(".*_","",X))
    # return(annotation[which(annotation$gene_name == X),]$gene_id)
    }

get_dend <- function(IDS,sig_genes,meta,ftcounts){
    
    metadata <- meta 
    rownames(metadata) <- metadata$PTID
    metadata <- metadata[IDS,]
    metadata <- metadata %>% select(Diagnosis, severity)

    ###----------------------------------
    ## Subset count matrix
    mat <- data.frame(ftcounts)[,IDS] %>% 
            filter(row.names(ftcounts) %in% all_of(sig_genes)) %>% 
            as.matrix()

    mat <- t(scale(t(mat)))

    ###----------------------------------
    ## Perform clustering

    h <- hclust(as.dist(1 - cor(mat, method = "pearson", use = 'pairwise.complete.obs')))
#     h <- hclust(as.dist(mat))

    dend <- as.dendrogram(h)

    dendL <- dendrapply(dend, colLab,metadata)

    return(dendL)
}


make_dend_list <- function(cf_tmp_rds, wb_tmp_res, 
                           cf_meta, wb_meta,
                           cf_ftcounts,wb_ftcounts,
                           PVAL_THRESH = PVAL_THRESH, LOGFC_THRESH = LOGFC_THRESH,N=0){
    
    #----------------------------------------
    ## Cell-free
    IDS_cf <- colnames(cf_tmp_rds[['dds']])
    if (N>0){IDS_cf <- sample(IDS_cf,N)}

    PTIDs <- cf_meta %>% filter(cfrna_sample_id %in% all_of(IDS_cf)) %>% pull(PTID)

    sig_genes_cf <- data.frame(cf_tmp_rds[['res']]) %>% filter(padj < PVAL_THRESH & abs(log2FoldChange) > LOGFC_THRESH) %>% rownames()
    cf_dendL = get_dend(PTIDs,sig_genes_cf,cf_meta,cf_ftcounts)

    #----------------------------------------
    ## Whole blood
    wb_tmp_res$ensmbl <- gsub("_.*","",wb_tmp_res$GeneID)

    sig_genes_wb <- wb_tmp_res %>% filter(padj < PVAL_THRESH & abs(log2FoldChange) > LOGFC_THRESH) %>% pull(ensmbl)
    wb_dendL = get_dend(PTIDs,sig_genes_wb,wb_meta,wb_ftcounts)

    return(dendlist(cf_dendL,wb_dendL))
}


get_pval <- function(dendlist,the_cor,R){
    
    set.seed(42)
    
    cor_bakers_gamma_results <- numeric(R)
    
    d1 <- dendlist[[1]]
    d2 <- dendlist[[2]]
    
    for(i in 1:R) {
       dend_mixed <- sample.dendrogram(d2, replace = FALSE)
       cor_bakers_gamma_results[i] <- cor_bakers_gamma(d1, dend_mixed)
    }
    
    if(sum(the_cor < cor_bakers_gamma_results)==0){return("ALL BELOW")}
    
    pval <- sum(the_cor < cor_bakers_gamma_results)/ R
                  
    return(pval)

}

get_conf <- function(dendlist,the_cor,R){

    set.seed(42)
    
    dend1 <- dendlist[[1]]
    dend2 <- dendlist[[2]]

    dend1_labels <- labels(dend1)
    dend2_labels <- labels(dend2)
    cor_bakers_gamma_results <- numeric(R)
    
    for(i in 1:R) {
       sampled_labels <- sample(dend1_labels, replace = TRUE)
       # members needs to be fixed since it will be later used in nleaves
       dend_mixed1 <- sample.dendrogram(dend1, 
                                        dend_labels=dend1_labels,
                                        fix_members=TRUE,fix_order=TRUE,fix_midpoint=FALSE,
                                        replace = TRUE, sampled_labels=sampled_labels
                                          )
       dend_mixed2 <- sample.dendrogram(dend2, dend_labels=dend2_labels,
                                        fix_members=TRUE,fix_order=TRUE,fix_midpoint=FALSE,
                                        replace = TRUE, sampled_labels=sampled_labels
                                          )                                    
       cor_bakers_gamma_results[i] <- cor_bakers_gamma(dend_mixed1, dend_mixed2, warn = FALSE)
    }


    CI95 <- quantile(cor_bakers_gamma_results, probs=c(.025,.975))
    
    return(CI95)
}

#----------------------------------------------------------------------------------------------
## LOAD OUTPUTS
#----------------------------------------------------------------------------------------------

cf_covid_cntrl <- readRDS("./subset_CF/daa_output/covid-control_paired.rds")
cf_misc_cntrl <- readRDS("./subset_CF/daa_output/misc-control_paired.rds")
cf_misc_covid <- readRDS("./subset_CF/daa_output/misc-covic_paired.rds")

wb_covid_cntrl <- read.delim("./subset_WB/tables/wb_covid-cntrl_paired_DESeq.tsv")
wb_misc_cntrl <- read.delim("./subset_WB/tables/wb_misc-cntrl_paired_DESeq.tsv")
wb_misc_covid <- read.delim("./subset_WB/tables/wb_misc-covid_paired_DESeq.tsv")


## PAIRED SAMPLES

paired <- read.delim("../1_sample-data/paired_sample_key.tsv") %>%
    filter(!is.na(cfrna_sample_id) & !is.na(wbrna_sample_id))


## FILTER

globin = c('HBA1','HBA2','HBB','HBBP1','HBD',
     'HBE1','HBG1','HBG2','HBM','HBQ1',
     'HBZ','HBZP1')
globin = annotation %>% filter(gene_name %in% all_of(globin)) %>% pull(gene_id) 

gene.list <- read.delim("../0_support-files/genelist.hs.tsv",col.names = c("type,","ENSMBL","gene_symbol"))

#----------------------------------------------------------------------------------------------
## WHOLE BLOOD
#----------------------------------------------------------------------------------------------

#------------------------
## READ OUTPUTS

wb_meta <- read.csv("../1_sample-data/STable7_wbrna-samples.csv") %>% 
    filter(wbrna_sample_id %in% all_of(paired$wbrna_sample_id))

wb_meta = wb_meta[!duplicated(wb_meta$PTID),]

wb_ftcounts <- read.delim("../1_sample-data/wbrna_ftcounts.txt") %>% 
    rename(gene_id = X) %>% 
    mutate(gene_id = gsub("\\_.*","",gene_id)) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    column_to_rownames(var = "gene_id")

wb_ftcounts <- wb_ftcounts[,wb_meta$wbrna_sample_id]
colnames(wb_ftcounts) <- wb_meta$PTID

#------------------------
## FILTER

gene.ids_wb <- gsub("\\..*","",rownames(wb_ftcounts))
exclude.idx_wb <- gene.ids_wb %in% c(gene.list[,2], globin)
wb_ftcounts = wb_ftcounts[!exclude.idx_wb,] 
wb_ftcounts <- wb_ftcounts[,colSums(wb_ftcounts) > 0]

#------------------------
## NORMALIZE

wb_ftcounts <- edgeR::cpm(wb_ftcounts)


#----------------------------------------------------------------------------------------------
## CELL FREE
#----------------------------------------------------------------------------------------------

#------------------------
## READ OUTPUTS

cf_meta <-read.csv("../1_sample-data/STable6_cfrna-samples.csv") %>%
    filter(cfrna_sample_id %in% all_of(paired$cfrna_sample_id))

cf_meta = cf_meta[!duplicated(cf_meta$PTID),]


cf_ftcounts <- read.delim("../1_sample-data/cfrna_ftcounts.txt") %>%
    rename(gene_id = Geneid) %>% 
    mutate(gene_id = gsub("\\_.*","",gene_id)) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    column_to_rownames(var = "gene_id")

cf_ftcounts <- cf_ftcounts[,cf_meta$cfrna_sample_id]
colnames(cf_ftcounts) <- cf_meta$PTID

#------------------------
## FILTER

gene.ids_cf <- gsub("\\..*","",rownames(cf_ftcounts))
exclude.idx_cf <- gene.ids_cf %in% c(gene.list[,2], globin)
cf_ftcounts = cf_ftcounts[!exclude.idx_cf,] 


#------------------------
## NORMALIZE
cf_ftcounts <- edgeR::cpm(cf_ftcounts)



#----------------------------------------------------------------------------------------------
## SIMULATION
#----------------------------------------------------------------------------------------------

output_misc_covid <- c()
output_misc_cntrl <- c()


for (n in c(1:1000)){

    NUM_SAMP = 26

    ## MISC v COVID
    output <- make_dend_list(cf_misc_covid, wb_misc_covid, 
                           cf_meta, wb_meta,
                           cf_ftcounts,wb_ftcounts,
                           PVAL_THRESH = PVAL_THRESH, LOGFC_THRESH = LOGFC_THRESH, N=NUM_SAMP)
    output_misc_covid <- c(output_misc_covid,cor_bakers_gamma(output[[1]],output[[2]]))
    
    
    ## MISC v CONTROL
    output <- make_dend_list(cf_misc_cntrl, wb_misc_cntrl, 
                           cf_meta, wb_meta,
                           cf_ftcounts,wb_ftcounts,
                           PVAL_THRESH = PVAL_THRESH, LOGFC_THRESH = LOGFC_THRESH, N=NUM_SAMP)
    
    output_misc_cntrl <- c(output_misc_cntrl,cor_bakers_gamma(output[[1]],output[[2]]))

}

df <- data.frame(output_misc_covid, output_misc_cntrl)
write.table(df, "./simulation_bakers_gamma.tsv",sep="\t",row.names=FALSE)