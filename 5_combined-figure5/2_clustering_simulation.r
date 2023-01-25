library(tidyverse)
library(data.table)
library(dendextend)
library(reshape2)
suppressMessages(library(pheatmap))
library(RColorBrewer)
library(gridExtra)

set.seed(42)

###----------------------------------

PVAL_THRESH = 0.05
LOGFC_THRESH = 1


###----------------------------------

annotation <- fread(file="../../0_metadata/gencode.biotype.name.key.tsv")

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
            if(Diagnosis=="COVID-19"){col_diag="#740C33"};if(Diagnosis=="MIS-C"){col_diag="#0234C1"};if(Diagnosis=="Control_Non-inflammatory"){col_diag="#4b7a47"};if(Diagnosis=="Control_Inflammatory"){col_diag="purple"}
        #Modification of leaf attribute
        attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=15,col=col_diag,lab.col=col_diag,lab.font=1,lab.cex=1))
        }
    return(n)
}


get_gene_id <- function(X,annotation){
    return(annotation[which(annotation$gene_name == X),]$gene_id)
    }


get_dend <- function(IDS,sig_genes,meta,ftcounts){
    
    ANNOTATIONS <- c("Diagnosis","severity",'pcr_positive_MC_1yes','antibody_positive_MC_1yes','kd_like_MC_1yes','req_vasopressors_inotropes_MC_1yes')

    if (is.null(meta$cfrna_file_id)) {meta$cfrna_file_id <- meta$SEQ_ID..UCSFonly.}
    
    metadata <- meta %>% 
    filter(cfrna_file_id %in% all_of(IDS)) %>%                                 ###
    mutate(PTID_DATE = gsub("-","\\.",PTID_DATE)) %>%
    column_to_rownames(var = "PTID_DATE") %>%  
    mutate(severity = ifelse(grepl("ontrol",severity),NA,severity)) %>%
    select(all_of(ANNOTATIONS))

    colnames(metadata) <- c("Diagnosis","Severity","PCR","Antibody")

    ###----------------------------------
    ## Subset count matrix
    mat <- data.frame(ftcounts) %>% 
            filter(row.names(ftcounts) %in% all_of(sig_genes)) %>% 
            select(all_of(rownames(metadata))) %>% 
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

    sig_genes_cf <- data.frame(cf_tmp_rds[['res']]) %>% filter(padj < PVAL_THRESH & abs(log2FoldChange) > LOGFC_THRESH) %>% rownames()
    cf_dendL = get_dend(IDS_cf,sig_genes_cf,cf_meta,cf_ftcounts)

    #----------------------------------------
    ## Whole blood
    wb_tmp_res$ensmbl <- lapply(wb_tmp_res$gene,get_gene_id,annotation)

    cf_pt_dates <- cf_meta %>% filter(cfrna_sample_id %in% all_of(IDS_cf)) %>% pull(PTID_DATE)

    IDS_wb <- wb_meta %>% 
        filter(PTID_DATE %in% all_of(cf_pt_dates)) %>% 
        pull(SEQ_ID..UCSFonly.)

    sig_genes_wb <- wb_tmp_res %>% filter(padj < PVAL_THRESH & abs(log2FoldChange) > LOGFC_THRESH) %>% pull(gene)
    wb_dendL = get_dend(IDS_wb,sig_genes_wb,wb_meta,wb_ftcounts)
    
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

# cf_covid_cntrl <- readRDS("../../1_cfRNA/DAA/comparisons_to_controls_ftcount/2_ALL-results.rds")
# cf_misc_cntrl <- readRDS("../../1_cfRNA/DAA/comparisons_to_controls_ftcount/3_ALL-results.rds")
# cf_misc_covid <- readRDS("../../1_cfRNA/DAA/comparisons_to_controls_ftcount/5_ALL-results.rds")


cf_covid_cntrl <- readRDS("./subset_CF/2_ALL-results.rds")
cf_misc_cntrl <- readRDS("./subset_CF/3_ALL-results.rds")
cf_misc_covid <- readRDS("./subset_CF/5_ALL-results.rds")

wb_covid_cntrl <- read.delim("./subset_WB/paired_wb_covid_ctrls.txt")
wb_misc_cntrl <- read.delim("./subset_WB/paired_wb_misc_ctrls.txt")
wb_misc_covid <- read.delim("./subset_WB/paired_wb_misc_covid.txt")

PatientDate_df <- read.delim("/workdir/cjl332/prevail/meta_data/patient-date_df.tsv") %>% mutate(PTID_DATE = gsub("\\-",".",paste0(PTID,"_",Date)))


#----------------------------------------------------------------------------------------------
## WHOLE BLOOD
#----------------------------------------------------------------------------------------------

#------------------------
## READ OUTPUTS

wb_meta <- read.delim("../../0_metadata/wbrna.tsv") %>% 
    filter(sample_group_matched != -1)
#     filter(Coverage..10...0.no..1.yes..NA.not.applicable.or.sample.not.available. == '1') %>%
#     mutate(PTID_DATE = paste0(PTID,"_",Date))

# wb_meta <- wb_meta[match(unique(wb_meta$PTID_DATE), wb_meta$PTID_DATE),]

wb_ftcounts <- read.delim("../../1_sample-output/wbrna_ftcounts.txt") %>% 
# wb_ftcounts <- read.delim("../../3_wbRNA/misc72_cnh19_covid37_ctrls26_counts.txt") %>% 
    rename(gene_id = X) %>% 
    mutate(gene_id = gsub("\\_.*","",gene_id)) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    column_to_rownames(var = "gene_id")


# wb_meta <- wb_meta[which(wb_meta$SEQ_ID..UCSFonly. %in% colnames(wb_ftcounts)),]
#------------------------
## FILTER

globin = c('HBA1','HBA2','HBB','HBBP1','HBD',
     'HBE1','HBG1','HBG2','HBM','HBQ1',
     'HBZ','HBZP1')
globin = annotation %>% filter(gene_name %in% all_of(globin)) %>% pull(gene_id) 

gene.list <- read.delim("../../0_metadata/genelist.hs.tsv",col.names = c("type,","ENSMBL","gene_symbol"))
gene.ids_wb <- gsub("\\..*","",rownames(wb_ftcounts))
exclude.idx_wb <- gene.ids_wb %in% c(gene.list[,2], globin)
wb_ftcounts = wb_ftcounts[!exclude.idx_wb,] 

wb_ftcounts <- wb_ftcounts[,colSums(wb_ftcounts) > 0]
#------------------------
## NORMALIZE

wb_ftcounts <- edgeR::cpm(wb_ftcounts)



#------------------------
## RENAME

wb_ftcounts <- wb_ftcounts[,colnames(wb_ftcounts) %in% wb_meta$SEQ_ID..UCSFonly.]
wb_ids <- colnames(wb_ftcounts)
wb_ptid_dates <- wb_meta %>% filter(SEQ_ID..UCSFonly. %in% wb_ids) %>% arrange(factor(SEQ_ID..UCSFonly., levels = wb_ids)) %>% pull(PTID_DATE)
colnames(wb_ftcounts) <- wb_ptid_dates



#----------------------------------------------------------------------------------------------
## CELL FREE
#----------------------------------------------------------------------------------------------

#------------------------
## READ OUTPUTS

cf_meta <-read.delim("../../0_metadata/cfrna.tsv") %>% filter(passQC_1yes_.1no == 1) %>% 
    filter(sample_group_matched != -1)

cf_ftcounts <- read.delim("../../1_sample-output/cfrna_ftcounts.txt") %>%
    rename(gene_id = Geneid) %>% 
    mutate(gene_id = gsub("\\_.*","",gene_id)) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    column_to_rownames(var = "gene_id")

# cf_meta <- cf_meta[which(cf_meta$cfrna_file_id %in% colnames(cf_ftcounts)),]

#------------------------
## FILTER

globin = c('HBA1','HBA2','HBB','HBBP1','HBD',
     'HBE1','HBG1','HBG2','HBM','HBQ1',
     'HBZ','HBZP1')
globin = annotation %>% filter(gene_name %in% all_of(globin)) %>% pull(gene_id) 

gene.list <- read.delim("../../0_metadata/genelist.hs.tsv",col.names = c("type,","ENSMBL","gene_symbol"))

gene.ids_cf <- gsub("\\..*","",rownames(cf_ftcounts))
exclude.idx_cf <- gene.ids_cf %in% c(gene.list[,2], globin)
cf_ftcounts = cf_ftcounts[!exclude.idx_cf,] 


#------------------------
## NORMALIZE
cf_ftcounts <- edgeR::cpm(cf_ftcounts)


#------------------------
## RENAME

cf_ftcounts <- cf_ftcounts[,colnames(cf_ftcounts) %in% cf_meta$cfrna_file_id]
cf_ids <- colnames(cf_ftcounts)
cf_ptid_dates <- cf_meta %>% filter(cfrna_file_id %in% cf_ids) %>% arrange(factor(cfrna_file_id, levels = cf_ids)) %>% pull(PTID_DATE)
colnames(cf_ftcounts) <- cf_ptid_dates



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