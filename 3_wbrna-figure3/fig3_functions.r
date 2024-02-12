topDiffGenes <- function(allScore) {  # function will be called in goData object to filter genes
  return(allScore < 0.05
      )
}

topGO <- function(res,SIG_THRESH = 0.01,sign){

    res_df = data.frame(res)
    res_df$geneid <- gsub("\\.[0-9]*","",rownames(res_df))                            # save just ensembl ID
    res_df = unique(res_df)  

    ##------------------------------------
    # Filter results and create list
    
    if (sign == "+"){
        res_df <- res_df[(res_df$log2FoldChange > 0) & (!is.na(res_df$padj)),]
        }else{
        res_df <- res_df[(res_df$log2FoldChange < 0) & (!is.na(res_df$padj)),]
    }
    
    
    
    geneList <- res_df$padj
    names(geneList) <- res_df$geneid


    ##------------------------------------
    # Create GOdata object

    # topGO accepts ("entrez", "genbank", "alias", "ensembl", "symbol", "genename", "unigene")
    #ontology = "BP" "CC", "MF"
    GOdata <- new("topGOdata",
                  ontology = "BP",
                  allGenes = geneList,
                  geneSel = topDiffGenes,
                  annot = annFUN.org,
                  mapping = "org.Hs.eg.db",
                  ID = "ensembl")


    ##------------------------------------
    # Test for Signficance

    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    allRes <- GenTable(GOdata,  classicFisher = resultFisher, classicKS = resultKS, KS.elim = resultKS.elim,
                       orderBy = "classicFisher",
                       ranksOf = "classicFisher", topNodes = 50, numChar = 100)


    allRes$genes <- sapply(allRes$GO.ID, function(x)
    {
      genes<-genesInTerm(GOdata, x)
      genes[[1]][genes[[1]] %in% names(geneList)] # myGenes is the queried gene list
    })
    allRes$genes[which(allRes$resultKS<0.05)] # print those only with p-value < 0.05

    # add the genes to the results table
    allRes$genes <-vapply(allRes$genes, paste, collapse = ",", character(1L))
    
    return(allRes)
    }

create_heatmap <- function(RDS,cpm,meta_data) {
    
    ANNOTATIONS = c("Diagnosis","severity",'group','ivig_rel_samp')
    temp_rds <- RDS

    ###----------------------------------
    ## Read in DESeq2 outputs
    dds <- temp_rds[['dds']]
    res <- temp_rds[['res']]
    
    SAMP_IDS = colnames(dds)
    
    if ("MIS-C" %in% dds$Diagnosis){
        
        NEW_IDS = meta_data %>% filter(group == "validation") %>% pull(cfrna_sample_id)
        SAMP_IDS = c(colnames(dds), NEW_IDS)
    
    }

    ###----------------------------------
    ## Extract significant genes
    sig_genes <- data.frame(res) %>% filter(padj < SIG_THRESH & abs(log2FoldChange) > FC_THRESH) %>% rownames()

    ###----------------------------------
    ## Subset metadata 
    metadata <- meta_data %>% 
        filter(cfrna_sample_id %in% all_of(SAMP_IDS)) %>% 
        column_to_rownames(var = "cfrna_sample_id") %>%  
        mutate(severity = ifelse(grepl("ontrol",severity),NA,severity)) %>%
        select(all_of(ANNOTATIONS))

    colnames(metadata) <- c("Diagnosis","Severity","Group","IVIG")

    ###----------------------------------
    ## Subset count matrix
    mat <- data.frame(cpm) %>% 
            filter(row.names(cpm) %in% all_of(sig_genes)) %>% 
            select(all_of(rownames(metadata))) %>% 
            as.matrix()

    mat <- t(scale(t(mat)))

    ###----------------------------------
    ## Colors
    color = colorRampPalette(c("blue","yellow"))(50)
    breaksList = seq(-2, +2, length = 51)

    #color_groups = c('COVID-19\ndiscovery' = '#c1272d', 'MIS-C\ndiscovery' = '#0000a7', 'Control_Non-inflammatory\ndiscovery' = '#eecc16',"MIS-C\nvalidation" = '#008176')
    
    my_colour <- list(
    Diagnosis = c('COVID-19' = '#F0484E', 'MIS-C' = '#5CB2EB', 'Control_Non-inflammatory' = '#FBE77C'), #"MISC_acute_validation" = '#00FFFF'),
    Severity = c("-1" = "white", "0" = "white", "2" = "#efe5d7", "3" = "#bc8e52"),
    Group = c("validation" = "orange", "discovery" = "light blue"),
    IVIG = c("after" = "maroon", "before" = "darkseagreen2", "concurrent with" = "dark blue","noivig"="white") #,
    # PCR = c("1" = "green", "2" = "blue"),
    # Antibody =c("1" = "green", "2" = "blue")
    )

    ###----------------------------------
    ## Plot

    # pdf("./tmp.pdf")
    heatmap_plt <- pheatmap(mat,

             # Colors
             col=color,
             breaks=breaksList, 
             annotation_col=metadata,
             annotation_colors=my_colour,
             na_col = "#FFFFFF",

             # Fonts
             show_colnames=F,
             show_rownames=F,
             fontsize=12,
             fontsize_col=3,
             annotation_names_col=F,
             annotation_names_row=F,

             # Clustering

             clustering_distance_cols="correlation",
             clustering_distance_rows="correlation", 
    #          clustering_distance_cols="euclidean",
    #          clustering_distance_rows="euclidean",
    #          cluster_cols=hc,
    #          cluster_rows=hr,
    #          clustering_distance_rows="euclidean",
    #          clustering_method="complete",
             treeheight_row=0,
            treeheight_col= 15,

             # Misc.
             border_color=NA,
            legend=FALSE,
            annotation_legend=FALSE
            ) 
    # dev.off()




    return(heatmap_plt)

        }