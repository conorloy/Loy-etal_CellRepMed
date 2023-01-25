scaleFUN <- function(x) sprintf("%.2f", x)

get_all_plots <- function(df_melt,val_vector,maxY,severity){
    
    labels <- names(val_vector)
    
    df_melt <- df_melt[which(df_melt$celltype %in% names(val_vector)),]

##---------------------------
# MIS-C
##---------------------------
##---------------------------

    diagnosis = "MIS-C"

    # Create master dataframe
    diag_df <- df_melt[which((df_melt$Diagnosis == diagnosis) & (df_melt$group == "discovery") ),]
    
    if (severity != FALSE){
    diag_df <- diag_df[which(diag_df$severity == severity),]}

    ##---------------------------
    # Summarize
    diag_sum <- diag_df %>% group_by(plottingName,celltype) %>% summarize(mean_fraction = mean(fraction))
    
    ##---------------------------
    # Plot area plot

    misc_area_plot <- diag_sum %>% 
    mutate(celltype= factor(celltype, levels = labels)) %>%
    ggplot(aes(x=plottingName,y=mean_fraction)) + 
    geom_area(aes(fill=celltype,group=celltype), colour="black",size = 0.15)+
#     theme_minimal()+
    theme_prevail()+
    theme(legend.position = "none",
          plot.margin=grid::unit(c(1,0,0,-0.9), "mm"))+
    scale_fill_manual(values=val_vector)+
    scale_x_discrete(expand=c(0,0.05))+
    scale_y_continuous(position = "left")+ 
    coord_cartesian(ylim = c(0,maxY))+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())
##---------------------------
# COVID-19
##---------------------------

    diagnosis = "COVID-19"

    # Create master dataframe
    diag_df <- df_melt[which((df_melt$Diagnosis == diagnosis) & (df_melt$group == "discovery") ),]
    
    if (severity != FALSE){
    diag_df <- diag_df[which(diag_df$severity == severity),]}

    ##---------------------------
    # Summarize
    diag_sum <- diag_df %>% group_by(plottingName,celltype) %>% summarize(mean_fraction = mean(fraction))
    
    ##---------------------------
    # Plot area plot

    covid_area_plot <- diag_sum %>% 
    mutate(celltype= factor(celltype, levels = labels)) %>%
    ggplot(aes(x=plottingName,y=mean_fraction)) + 
    geom_area(aes(fill=celltype,group=celltype), colour="black",size = 0.15)+
#     theme_minimal()+
    theme_prevail()+
    theme(legend.position = "none",
          plot.margin=grid::unit(c(1,0,0,-0.9), "mm"))+
    scale_fill_manual(values=val_vector)+
    scale_x_discrete(expand=c(0,0.05))+
    scale_y_continuous(position = "left")+ 
    coord_cartesian(ylim = c(0,maxY))+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())
    

##---------------------------
# MIS-C validation
##---------------------------

        ##---------------------------
        # Create master dataframe
        cntrl_df <- df_melt[which(df_melt$group == "validation"),]

        ##---------------------------
        # Summarize
        cntrl_sum <- cntrl_df %>% group_by(plottingName,celltype) %>% summarize(mean_fraction = mean(fraction))

        ##---------------------------
        # Plot area plot

        validation_plot <- cntrl_sum %>%
        mutate(celltype=factor(celltype,levels=labels)) %>%
        ggplot(aes(x=plottingName,y=mean_fraction,fill=celltype))+
        geom_bar(position="stack",stat="identity", colour="black",size = 0.15)+
#         theme_minimal()+
    theme_prevail()+
        theme(plot.margin=grid::unit(c(1,0,0,0), "mm"),
             legend.position = "none")+
        scale_fill_manual(values=val_vector)+
        scale_y_continuous(position = "left",labels=scaleFUN)+ 
        coord_cartesian(ylim = c(0,maxY))+
        theme(axis.title.y = element_blank(),
             axis.title.x = element_blank())

##---------------------------
# Vaccine Controls
##---------------------------

    ##---------------------------
    # Create master dataframe
    cntrl_df <- df_melt[which(df_melt$Diagnosis == "Control_Non-inflammatory"),]
    
    ##---------------------------
    # Summarize
    cntrl_sum <- cntrl_df %>% group_by(plottingName,celltype) %>% summarize(mean_fraction = mean(fraction))

    ##---------------------------
    # Plot area plot

    control_plot <- cntrl_sum %>%
    mutate(celltype=factor(celltype,levels=labels)) %>%
    ggplot(aes(x=plottingName,y=mean_fraction,fill=celltype))+
    geom_bar(position="stack",stat="identity", colour="black",size = 0.15)+
#     theme_minimal()+
    theme_prevail()+
    theme(plot.margin=grid::unit(c(1,0,0,0), "mm"),
         legend.position = "top")+
    scale_fill_manual(values=val_vector)+
    scale_y_continuous(position = "right",labels=scaleFUN)+ 
    coord_cartesian(ylim = c(0,maxY))+
    theme(axis.title.y=element_blank(),
         axis.title.x = element_blank())+ 
    theme(legend.title = element_blank(),
          legend.text = element_text(family="Helvetica", size=6, color='black'),
         legend.key.size = unit(0.2, "cm"),,
        legend.spacing.x = unit(0.01, 'cm'))

##---------------------------
# Legend
##---------------------------
    
    ##---------------------------
    # Extract legend from control
    
    legend <- cowplot::get_legend(control_plot)

    control_plot <- control_plot + theme(legend.position="none")

    return(list(validation_plot,misc_area_plot,covid_area_plot,control_plot,legend))
}



make_boxplot <- function(diag_df,CELLTYPE,YLIM,expGroupPalette){

    plot_df <- diag_df %>% filter(celltype == CELLTYPE) %>%
    mutate(expGroup = factor(expGroup,levels = c('Control_Non-inflammatory Not-hospitalized discovery','COVID-19 acute discovery','MIS-C acute discovery','MIS-C acute validation')))

    
    bxplt <- plot_df %>%
        ggplot() + 
        geom_boxplot(aes(x=expGroup,y=fraction,fill=expGroup),outlier.shape = NA,size = 0.2, color = "black")+ #, color = "black"
        geom_jitter(aes(x=expGroup,y=fraction),width=.2,height=0,size = 0.01)+
#         geom_point(aes(colour = severity,shape = severity,group=plottingName),position = position_jitterdodge(jitter.width=0.1,seed=100,jitter.height=0),size = 0.5,alpha= 0.5)+       # ,pch=21 ,color = "black"
#         annotate("text", x = 1, y = YLIM - 0.01, label = CELLTYPE,size = 2.5,family = "Helvetica")+
        theme_prevail()+
#         scale_color_manual(values=c('1' = "black",
#                                        '2' = "red")) + 
#         scale_shape_manual(values=c('1'= 16,
#                                        '2'= 4))+
        scale_fill_manual(values=expGroupPalette) + 
        theme(panel.grid.minor = element_blank())


    ##---------------------------------------
    # ADD ARROWS    

    outliers <- plot_df%>% 
            filter(fraction > YLIM)

    subset <- unique(plot_df$expGroup)

    y_stop = YLIM 
    y_start= YLIM - (YLIM*0.05)

    XORDER <- list('Control_Non-inflammatory discovery'= 1,
                'COVID-19 acute discovery'= 2,
                'MIS-C acute discovery'= 3,
                'MIS-C acute validation'= 4)

    ffset_x = 1  # offset for different condition / plotting names
    offset_y = .20     # offset from top for different outliers in same condition




    text_offset = 0.09
    text_size = 1.5
    text_lineheight = .75

    arrow_size = 0.5

    for (i  in 1:length(subset)){
        group_name <- subset[i]
        group_outliers <- outliers %>% filter(expGroup == group_name)

        ## ADD ARROW
        if (nrow(group_outliers) > 0){

            outlier_values <- round(group_outliers[, "fraction",drop=TRUE],2)
            outlier_values <- as.character(outlier_values[order(outlier_values,decreasing=TRUE)])

            text_center <- y_start * (1 - (.01*(length(outlier_values)-1)))

            XVAL <- XORDER[[group_name]]        

            bxplt <- bxplt + annotate("segment", x = XVAL, xend = XVAL,
                                              y = YLIM-0.007, yend = YLIM,
                                              size = arrow_size, lineend="butt", linejoin="mitre", arrow=arrow(length=unit(.06,"npc")))



            for (ii in 1:length(outlier_values)){

                text_color = "black"
                if (group_outliers[ii, "severity"] == 2){
                    text_color = "red"} 

                if (ii == 1){YVAL = YLIM
                            }else { 
                    YVAL = YLIM- ((ii-1)*0.006)
                }

                VAL <- outlier_values[ii]

                bxplt <- bxplt + annotate("text", x = XVAL, y = YVAL, hjust=-0.2,vjust=0.5,
                                              label = VAL,
                                              size = text_size, family = "Helvetica", lineheight = text_lineheight,color = text_color)

            }}}

    ##---------------------------------------
    # ADD PVALS
    
    #YLIM
    stat.test <- data.frame(plot_df) %>%
        wilcox_test(fraction ~ expGroup, paired = FALSE) 
    
    
    stat.test <- stat.test[c(1,2,4,6),] %>% 
        adjust_pvalue(method = "BH") %>% 
        add_significance("p.adj") %>% 
        add_xy_position(x = "expGroup") 
        print(stat.test)


    
#     stat.test <- stat.test %>% filter( (group1 == "Control_Non-inflammatory Not-hospitalized discovery" & group2 == "COVID-19 acute discovery" )|
#                                       (group1 == "Control_Non-inflammatory Not-hospitalized discovery" & group2 == "MIS-C acute discovery" )|
#                                       (group1 == "COVID-19 acute discovery" & group2 == "MIS-C acute discovery" )|
#                                       (group1 == "MIS-C acute discovery" & group2 == "MIS-C acute validation")
#                                       )
    
#     stat.test <- stat.test %>% 
#         adjust_pvalue(method = "BH") %>% 
#         add_significance("p.adj") %>% 
#         add_xy_position(x = "expGroup") 
#     print(data.frame(stat.test))
    
    
    
#     # rescale
    max.y <- max(stat.test$y.position) 
    stat.test$y.position <- YLIM - (max.y - stat.test$y.position) - (.1*YLIM)
    
        
    bxplt <- bxplt +  stat_pvalue_manual(stat.test, label = "p.adj.signif")
    
    return(bxplt)
    }



make_severity_plt <- function(diag_df, CELLTYPE){
    bxplt <- diag_df %>%
        filter(celltype == CELLTYPE) %>%
        ggplot() + 
        geom_boxplot(aes(x=group,y=fraction,fill=group),outlier.shape = NA,size = 0.2, color = "black")+ #, color = "black"
        geom_jitter(aes(x=group,y=fraction),width=.2,height=0,size = 0.01)+
#         geom_point(aes(colour = severity,shape = severity,group=plottingName),position = position_jitterdodge(jitter.width=0.1,seed=100,jitter.height=0),size = 0.5,alpha= 0.5)+       # ,pch=21 ,color = "black"
#         annotate("text", x = 1, y = YLIM - 0.01, label = CELLTYPE,size = 2.5,family = "Helvetica")+
        theme_prevail()+
#         scale_color_manual(values=c('1' = "black",
#                                        '2' = "red")) + 
#         scale_shape_manual(values=c('1'= 16,
#                                        '2'= 4))+
        labs(title = CELLTYPE)+
        scale_fill_manual(values=expGroupPalette) + 
        theme(panel.grid.minor = element_blank(),
              axis.title = element_blank(),
             title = element_blank())
    
    return(bxplt)

}



make_alpha_plot <- function(alpha,SAMP_LEVELS,METHOD){
alpha <- alpha %>% filter(expGroup %in% all_of(SAMP_LEVELS)) %>% mutate(expGroup = factor(expGroup, levels = all_of(SAMP_LEVELS)))

PLOT <- alpha %>%
    ggplot()+
    geom_boxplot(aes(x=expGroup, y=.data[[method]], fill=expGroup), outlier.shape=NA,size = .2,fatten = 0.7, color = "black" )+ #, color = "black"
    geom_jitter(aes(x=expGroup, y=.data[[method]]),width=.2,height=0,size = 0.1) #, colour = factor(severity),shape=factor(severity)


ADD_ON <- list(theme_prevail(),
                theme(legend.position = "none",
                    axis.title.y = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.margin=grid::unit(c(0.04,.04,0,0), "in")),
#                 scale_color_manual(values=c('0' = "black",'1' = "black",'2' = "black",'3' = "black")),
#                 scale_shape_manual(values=c('0'= 1,'1'= 1, '2'= 4, '3'= 4)),
#                 scale_x_discrete(guide = guide_axis(n.dodge = 2)),
               scale_fill_manual(values=expGroupPalette) ,
#                scale_y_continuous(breaks = seq(0,1,length.out  = 4)),
               coord_cartesian(ylim = c(0,1))
              )

alpha$method <- alpha[[METHOD]]
stat.test <- alpha %>%
    wilcox_test(method ~ expGroup, ref.group = "Control Not-hospitalized discovery") %>% 
    adjust_pvalue(method = "BH") %>% 
    add_significance("p.adj") %>% 
    add_xy_position(x = "expGroup"
                   
)

PLOT <- PLOT + 
        ADD_ON + 
        stat_pvalue_manual(stat.test, label = "p.adj.signif") + 
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
    
    

return(PLOT) 
}


make_bc_plot <- function(bc_melt,SAMP_LEVELS){
    
    bc_melt <- bc_melt %>%
        mutate(Group2 = factor(Group2,levels = all_of(SAMP_LEVELS))) %>%
        filter(!is.na(Group2)) %>%
        mutate(expGroup = sub(".*?\\+","",Group2))
    
#     bc_melt_mean <- bc_melt %>%
#         mutate(Group2 = factor(Group2,levels = all_of(SAMP_LEVELS))) %>%
#         filter(!is.na(Group2)) %>%
#         mutate(expGroup = sub(".*?\\+","",Group2)) %>%
#         group_by(Sam.ID_1,expGroup) %>% summarize(mean_BC = mean(BC)) %>%
#         mutate(Group2 = paste0("Control Not-hospitalized discovery+",expGroup))%>%
#         mutate(Group2 = factor(Group2,levels = all_of(SAMP_LEVELS)))%>%
#         filter(!is.na(Group2))
        
    
    PLOT <- bc_melt %>% ggplot()+
                        geom_boxplot(aes(x=Group2,y=BC,fill = expGroup), outlier.shape=NA,size = .2,fatten = 0.7, color = "black")+
#                         geom_boxplot(data=bc_melt_mean,aes(x=Group2,y=mean_BC,fill = expGroup), outlier.shape=NA,size = .4,fatten = 0.7, color = "black")+
                        geom_jitter(aes(x=Group2,y=BC,fill = expGroup), width=.2,height=0,size = 0.01) #+
#                         geom_jitter(data=bc_melt_mean,aes(x=Group2,y=mean_BC),width=.2,height=0,size = 0.01)

    ADD_ON <- list(
                theme_prevail(),
                theme(legend.position = "none",
                    axis.title.y = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.margin=grid::unit(c(0.75,.25,0,0), "mm")),
                scale_color_manual(values=c('1' = "black",'2' = "red")),
                scale_shape_manual(values=c('1'= 16, '2'= 4)),
                scale_x_discrete(guide = guide_axis(n.dodge = 2)),
               scale_fill_manual(values=expGroupPalette),
#                scale_y_continuous(breaks = seq(0,1,length.out  = 4)),
               coord_cartesian(ylim = c(0,1.0))
              )
    
    

    stat.test <- bc_melt %>%
        wilcox_test(BC ~ Group2, ref.group = "Control Not-hospitalized discovery+Control Not-hospitalized discovery") %>% 
#     stat.test <- bc_melt_mean %>%
#         wilcox_test(mean_BC ~ Group2, ref.group = "Control Not-hospitalized discovery+Control Not-hospitalized discovery") %>% 
        adjust_pvalue(method = "BH") %>% 
        add_significance("p.adj") %>% 
        add_xy_position(x = "Group2")
    
    print(stat.test)

    PLOT <- PLOT + 
            ADD_ON + 
            stat_pvalue_manual(stat.test, label = "p.adj.signif") #+ 
#             scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
    
    return(PLOT)
}



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

get_exp_plot <- function(GENE,meta_data,count_cpm,gene_lookup= annotation){
    
    samp_df <- data.frame("cfrna_sample_id" = meta_data$cfrna_sample_id,
                      "Diagnosis" = meta_data$Diagnosis,
                     "severity" = meta_data$severity,
                     "group" = meta_data$group)
    
    GENE <- gene_lookup[which(gene_lookup$gene_name == GENE),]$gene_id

    samp_df <- merge(samp_df,count_cpm[,GENE],by.x="cfrna_sample_id",by.y=0)

    samp_df$plotting_name <- paste(samp_df$Diagnosis,samp_df$group,sep="\n")

    samp_df$plotting_name <- factor(samp_df$plotting_name, levels = c("Control_Non-inflammatory\ndiscovery","COVID-19\ndiscovery","MIS-C\ndiscovery","MIS-C\nvalidation"))

    samp_df %>%
    ggplot(aes(x=plotting_name,y=y,fill=plotting_name))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(height=0,width=0.2,size = 1)+
    annotate(geom="text",label=gsub(".*_","",GENE),x=-Inf,y=Inf,hjust=-.1,vjust=1.25)+
    scale_fill_manual(values=color_groups)+
    theme_prevail()+
    theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5))+
    labs(y="CPM")+
#     labs(y="TPM",title=gsub(".*\\_","",GENE))+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))
}


get_fig_plot <- function(GENE,meta_data,count_cpm,YLIM, gene_lookup= annotation, expGroupPalette=expGroupPalette){
    
    ##---------------------------------------
    # PREPARE DATA
    samp_df <- data.frame("cfrna_sample_id" = meta_data$cfrna_sample_id,
                      "Diagnosis" = meta_data$Diagnosis,
                     "severity" = meta_data$severity,
                     "group" = meta_data$group)
    GENE_N <- GENE
    GENE <- gene_lookup[which(gene_lookup$gene_name == GENE),]$gene_id

    samp_df <- merge(samp_df,count_cpm[,GENE],by.x="cfrna_sample_id",by.y=0)

    samp_df$plotting_name <- paste(samp_df$Diagnosis,samp_df$group,sep="\n")

    samp_df$plotting_name <- factor(samp_df$plotting_name, levels = c("Control_Non-inflammatory\ndiscovery","COVID-19\ndiscovery","MIS-C\ndiscovery","MIS-C\nvalidation"))
    
    ##---------------------------------------
    # MAKE PLOT
    PLOT <- samp_df %>%
    ggplot(aes(x=plotting_name,y=y,color=plotting_name))+
#     geom_jitter(height=0,width=0.2,size = .65, color="black")+
    geom_point(  size = .75, position = position_jitter(seed = 42,height=0,width=0.2))+ #colour="black",pch=21,aes(fill=plotting_name),
#     annotate(geom="text",label=gsub(".*_","",GENE_N),x=-Inf,y=Inf,hjust=-.1,vjust=1.25,size = 2)+
    scale_color_manual(values=color_groups)+
    theme_prevail()+
    theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
         plot.margin=grid::unit(c(0.02,0.02,0,0), "in"))+
    labs(y="CPM")+
#     labs(y="TPM",title=gsub(".*\\_","",GENE))+
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) #+
#     scale_color_manual(values=expGroupPalette)
    
    ##---------------------------------------
    # ADD ARROWS
    outliers <- samp_df%>% 
            filter(y > YLIM)
    
    print("num outliers:")
    print(nrow(outliers))

    y_stop = YLIM 
    y_start= YLIM - (YLIM*0.05)

    XORDER <- list('Control_Non-inflammatory\ndiscovery'= 1,
                'COVID-19\ndiscovery'= 2,
                'MIS-C\ndiscovery'= 3,
                'MIS-C\nvalidation'= 4)

    offset_x = 1  # offset for different condition / plotting names
    offset_y = .20     # offset from top for different outliers in same condition


    text_offset = 0.09
    text_size = 1.5
    text_lineheight = .75

    arrow_size = 0.5

    for (i  in 1:length(XORDER)){
        group_name <- names(XORDER[i])
        group_outliers <- outliers %>% filter(plotting_name == group_name)

        ## ADD ARROW
        if (nrow(group_outliers) > 0){

            outlier_values <- round(group_outliers[, "y",drop=TRUE],2)
            outlier_values <- as.character(outlier_values[order(outlier_values,decreasing=TRUE)])

            text_center <- y_start * (1 - (.01*(length(outlier_values)-1)))

            XVAL <- XORDER[[group_name]]        

            PLOT <- PLOT + annotate("segment", x = XVAL, xend = XVAL,
                                              y = YLIM-0.007, yend = YLIM,
                                              size = arrow_size, lineend="butt", linejoin="mitre", arrow=arrow(length=unit(.06,"npc")))



            for (ii in 1:length(outlier_values)){

                text_color = "black"
                if (group_outliers[ii, "severity"] == 2){
                    text_color = "red"} 

                if (ii == 1){YVAL = YLIM
                            }else { 
                    YVAL = YLIM- ((ii-1)*0.006)
                }

                VAL <- outlier_values[ii]

                PLOT <- PLOT + annotate("text", x = XVAL, y = YVAL, hjust=-0.2,vjust=0.5,
                                              label = VAL,
                                              size = text_size, family = "Helvetica", lineheight = text_lineheight,color = text_color)

            }}}
    
    ##---------------------------------------
    # ADD SIG BARS
    stat.test <- samp_df %>%
    wilcox_test(y ~ plotting_name) %>%  #, ref.group = "MIS-C\ndiscovery"
    adjust_pvalue(method = "BH") %>% 
    add_significance("p.adj") %>% 
    add_xy_position(x = "plotting_name")
    
#     print(stat.test)
    
    stat.test <- stat.test %>% arrange(desc(y.position))
    stat.test$rank <- c(1:nrow(stat.test))
    stat.test$y.position <- YLIM - ((0.05*stat.test$rank)*YLIM)
    
    
    PLOT <- PLOT + stat_pvalue_manual(stat.test, label = "p.adj.signif",tip.length=0) 
    
    return(PLOT)
}