
`%nin%` <- Negate(`%in%`)

findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

#### Access sample data ####
read_cf_samples <- function(data_tsv, excl_list) {
  exclude_samples <- excl_list
  read_tsv(data_tsv, col_names = T) %>%
    dplyr::filter(Stage == "Larva" & Sample %nin% exclude_samples) %>%
    dplyr::select(Run, Sample, Name = `Submitted Name`, Host_form, CF, Treatment, Family, Stress, Perf_plasticity, CH_plasticity, CO_plasticity) %>%
    mutate(Sample = str_sub(Sample, 1,10),
           Treatment = factor(Treatment, levels = c("H", "HH", "HO", "O", "OO", "OH")),
           Switch = NA,
           Switch = if_else(Treatment == "OH", "OH", 
                            if_else(Treatment == "HO", "HO", NA)))
}

#### Check mapping rates ####

load_mapping_rates <- function(data, samples) {
  read_delim(delim = " ", file = data, col_names = c("Sample", "total_reads", "perc_map")) %>% 
    right_join(samples)
}

plot_mapping_rates <- function(data, pal) {
  
  outliers <- findoutlier(data$perc_map)
  
  ggplot(data, aes(x = Host_form, y = perc_map)) +
    geom_jitter(aes(fill = Treatment, shape = Treatment), size = 2, height = 0.01, width = 0.05) +
    scale_fill_manual(values = pal) +
    scale_shape_manual(values = c(22,22,22,21,21,21)) +
    geom_boxplot(fill = NA, width = 0.5) +
    stat_compare_means(method = "t.test", hjust =-0.5) +
    scale_x_discrete(labels = c("CH", "CO")) +
    labs(x = "Host race", y = "Quantification rate (% of reads)") + 
    # geom_text(data = data[outliers,], aes(label=Sample), na.rm=TRUE, hjust=0, size = 3) +
    theme_classic()  +
    theme(axis.text.x = ggtext::element_markdown())
}

#### Compile trimming, quantification statistics ####

comb_multiqc <- function(raw, trim, samples){
  cols_raw <- c("Read","raw_percent_duplicates","raw_percent_gc",
                "raw_avg_sequence_length", "raw_percent_fails", "raw_total_sequences")
  cols_trim <- c("Read","percent_trimmed", "trim_percent_duplicates","trim_percent_gc",
                 "trim_avg_sequence_length", "trim_percent_fails", "trim_total_sequences")
  raw_mqc <- read_tsv(raw, col_names = cols_raw, skip = 1) %>% 
    mutate(Run = str_split_i(Read, "_R", i = 1))
  trim_mqc <- read_tsv(trim, col_names = cols_trim, skip = 1) %>% 
    mutate(Read = gsub("_val_1||_val_2", "", Read),
           Run = str_split_i(Read, "_R", i = 1) )
  
  comb <- dplyr::select(samples, Run, Sample, Family, Host_form, Treatment) %>% 
    left_join(raw_mqc) %>% 
    left_join(dplyr::select(trim_mqc, Read, Run, percent_trimmed) %>% na.omit()) %>% 
    left_join(dplyr::select(trim_mqc, -percent_trimmed) %>% na.omit()) 
}


#### Access count data ####

load_dds <- function(dds) {
  load(dds)
  ddsTxi
}

subset_by_samples <- function(dds, samples){
  sample_ids <- samples$Sample
  colnames(dds) <- str_sub(colnames(dds), 1, 10)
  dds[,colnames(dds) %in% sample_ids]
  
}

filter_by_counts <- function(dds) {
  keep.genes <- rowSums(assay(dds) >= 5 ) >= 6
  dds[keep.genes]
}


#### Unsupervised clustering: PCA ####

rlog2pca <- function(mat, n) {
  # identify and filter to keep exons with highest variation
  top <-names(sort(apply(mat, 1, var), decreasing = T))[1:n]
  mat_top <- t(mat[top,])
  mat_top %>% prcomp() 
}

prcompVarSumm <- function(prcomp_out) {
  as.data.frame(t(summary(prcomp_out)$importance)) %>% 
    rownames_to_column(var = "Principal Component")
} 

plot_var_comp <- function (prcomp_summ) {
  ggplot(prcomp_summ[1:6,], aes( x= `Principal Component`, y = `Proportion of Variance`)) +
    geom_bar(stat = "identity") +
    labs(x = "Principal component", "Prop. of variance") +
    geom_label(stat = "identity", aes(label = round(`Proportion of Variance`, 3))) +
    theme_classic(base_size = 15)
}

prcomp2df <- function(prcomp_out, samples) {
  as.data.frame(prcomp_out$x) %>%
    rownames_to_column(var = "Sample") %>% 
    left_join(dplyr::select(samples, Sample, Treatment, Host_form)) 
}

plot_pca_list <- function(var_sum, prcompdf, pal){
  plotlist <- list()
  for (i in 1:4){
    PC_axes <- c(i,i+1)
    PC_columns <- c(paste0("PC",PC_axes[1]), paste0("PC",PC_axes[2]))
    PC_var <- c(paste0(PC_columns[1], " (", round(var_sum[PC_axes[1],3]*100, 1), "%)"),
                paste0(PC_columns[2], " (", round(var_sum[PC_axes[2],3]*100, 1), "%)"))
    
    plotlist[[i]] <- ggplot(prcompdf, aes(.data[[PC_columns[1]]],.data[[PC_columns[2]]])) +
      geom_point( aes(fill = Treatment, shape = Host_form), size = 2)+ #size = 2, stroke = 1.25,
      scale_fill_manual(values =pal, name = "Condition", guide = "none")+
      scale_shape_manual(values =c(22,21), guide = "none")+
      labs(x = PC_var[1], y =PC_var[2]) +
      theme_classic() 
  }
  plotlist
}


#### Unsupervised clustering: MDS ####
mds2df <- function(mds, samples) {
  data.frame(Sample = rownames(mds$distance.matrix.squared),
             mds.x = mds$x,
             mds.y = mds$y) %>% 
    left_join(dplyr::select(samples, Sample, Treatment, Host_form)) 
}

plot_mds <- function(mds, mds_df, pal){
  
  mds_varExplained <- as.numeric(mds$var.explained) 
  mds_axes <- c(1,2)
  mds_columns <- c(paste0("Axis ",mds_axes[1]), paste0("Axis",mds_axes[2]))
  mds_var <- c(paste0(mds_columns[1], " (", signif(mds_varExplained[1], 3) *100, "%)"),
               paste0(mds_columns[2], " (", signif(mds_varExplained[2], 3) *100, "%)"))
  
  ggplot(mds_df, aes(x = mds.x, y = mds.y)) +
    geom_point(aes(fill = Treatment, shape = Host_form), size = 2)+ #size = 2, stroke = 1.25,
    scale_fill_manual(values =pal, name = "Condition", guide = "none")+
    scale_shape_manual(values =c(22,21), guide = "none")+
    labs(x = mds_var[1], y =mds_var[2]) +
    theme_classic()
}


#### Unsupervised clustering: heatmap ####

rlog2heatmap <- function(mat, n, samples, pal){
  top <-names(sort(apply(mat, 1, var), decreasing = T))[1:n]
  mat_top <- mat[top,]
  
  colfun <- colorRamp2(breaks = c(-2,9,20), colors = c("#1f449c","white", "#f05039"))
  
  ht_trt <- samples$Treatment
  names(ht_trt) <- samples$Sample
  col_trt <- pal
  names(col_trt) <- levels(samples$Treatment)
  
  Heatmap(mat_top, 
          top_annotation = HeatmapAnnotation(
            Trt = ht_trt,
            col = list(Trt = col_trt), 
            show_annotation_name = TRUE
          ),
          clustering_method_rows = "ward.D2",
          clustering_method_columns = "ward.D2",
          show_row_names = F, 
          column_names_gp = gpar(fontsize = 6),
          # column_labels = modulelabel,
          column_names_side = "top",
          col = colfun,
          na_col = "grey",
          # cluster_rows = F, 
          show_heatmap_legend = T, 
          row_names_max_width = unit(10, units = "cm"))
}


#### DGE ####

run_deseq_cond_mat <- function(dds, cond_mat) {
  
  cond_names <- paste0(cond_mat[,1],"_vs_", cond_mat[,2])
  
  deseq_results_list <- list()
  deseq_results_list$results <- list()
  deseq_results_list$shrink <- list()
  deseq_results_list$up <- list()
  deseq_results_list$down <- list()
  
  for (i in 1:dim(cond_mat)[1]){
    # conditions are too heterogeneous to test using the full design matrix
    # recommendation from M. Love (https://support.bioconductor.org/p/9150170/)
    dds_sub <- dds[,dds$Treatment %in% cond_mat[i,]]
    dds_sub$Treatment <- droplevels(dds_sub$Treatment)
    dds_sub$Treatment <- relevel(dds_sub$Treatment, ref = cond_mat[i,2])
    dds_sub <- DESeq(dds_sub, minReplicatesForReplace=8)
    
    deseq_results_list$results[[i]] <- results(dds_sub, 
                                               name=paste0("Treatment_", cond_names[i]))
    
    deseq_results_list$shrink[[i]] <- lfcShrink(dds_sub, 
                                                coef=paste0("Treatment_",cond_names[i]), 
                                                type="apeglm",
                                                svalue = T)
    
    deseq_results_list$up[[i]] <- deseq_results_list$shrink[[i]] %>% 
      as.data.frame() %>% 
      na.omit() %>% 
      rownames_to_column("geneID") %>% 
      # dplyr::filter(padj < 0.05 & log2FoldChange > 0.1 ) %>% 
      dplyr::filter(svalue < 0.001 & log2FoldChange > 0 ) %>%
      dplyr::select(geneID)
    
    deseq_results_list$down[[i]] <- deseq_results_list$shrink[[i]] %>% 
      as.data.frame() %>% 
      na.omit() %>% 
      rownames_to_column("geneID") %>% 
      # dplyr::filter(padj < 0.05 & log2FoldChange < -0.1 ) %>% 
      dplyr::filter(svalue < 0.001 & log2FoldChange < 0 ) %>%
      dplyr::select(geneID)
    
    names(deseq_results_list$results)[i] <- paste0(cond_mat[i,2],"_vs_", cond_mat[i,1])
    names(deseq_results_list$shrink)[i] <- paste0(cond_mat[i,2],"_vs_", cond_mat[i,1])
    names(deseq_results_list$up)[i] <- paste0(cond_mat[i,2],"_vs_", cond_mat[i,1])
    names(deseq_results_list$down)[i] <- paste0(cond_mat[i,2],"_vs_", cond_mat[i,1])
  }
  deseq_results_list
}


run_deseq_cond_mat_fullMat <- function(dds, cond_mat) {
  
  cond_names <- paste0(cond_mat[,1],"_vs_", cond_mat[,2])
  
  deseq_results_list <- list()
  deseq_results_list$results <- list()
  deseq_results_list$shrink <- list()
  deseq_results_list$up <- list()
  deseq_results_list$down <- list()
  
  dds$Treatment <- relevel(dds$Treatment, ref = cond_mat[1,2])
  dds <- DESeq(dds, minReplicatesForReplace=8)
  
  for (i in 1:dim(cond_mat)[1]){
    if (levels(dds$Treatment)[1] != cond_mat[i,2]) {
      dds$Treatment <- relevel(dds$Treatment, ref = cond_mat[i,2])
      dds <- DESeq(dds, minReplicatesForReplace=8)
    }
    
    deseq_results_list$results[[i]] <- results(dds, 
                                               name=paste0("Treatment_", cond_names[i]))
    
    deseq_results_list$shrink[[i]] <- lfcShrink(dds, 
                                                coef=paste0("Treatment_",cond_names[i]), 
                                                type="apeglm",
                                                svalue = T)
    
    deseq_results_list$up[[i]] <- deseq_results_list$shrink[[i]] %>% 
      as.data.frame() %>% 
      na.omit() %>% 
      rownames_to_column("geneID") %>% 
      # dplyr::filter(padj < 0.05 & log2FoldChange > 0.1 ) %>% 
      dplyr::filter(svalue < 0.001 & log2FoldChange > 0 ) %>%
      dplyr::select(geneID)
    
    deseq_results_list$down[[i]] <- deseq_results_list$shrink[[i]] %>% 
      as.data.frame() %>% 
      na.omit() %>% 
      rownames_to_column("geneID") %>% 
      # dplyr::filter(padj < 0.05 & log2FoldChange < -0.1 ) %>% 
      dplyr::filter(svalue < 0.001 & log2FoldChange < 0 ) %>%
      dplyr::select(geneID)
    
    names(deseq_results_list$results)[i] <- paste0(cond_mat[i,2],"_vs_", cond_mat[i,1])
    names(deseq_results_list$shrink)[i] <- paste0(cond_mat[i,2],"_vs_", cond_mat[i,1])
    names(deseq_results_list$up)[i] <- paste0(cond_mat[i,2],"_vs_", cond_mat[i,1])
    names(deseq_results_list$down)[i] <- paste0(cond_mat[i,2],"_vs_", cond_mat[i,1])
  }
  deseq_results_list
}


plot_pvals <- function(deseq_results_list) {
  plot_list <- list()
  # create for loop to produce ggplot2 graphs 
  for (i in 1:length(deseq_results_list)) { 
    
    df <- deseq_results_list$results[[i]] %>% 
      as.data.frame() %>% 
      rownames_to_column("geneID") %>% 
      na.omit()
    
    plot_list[[i]] <- ggplot(df, aes(x = pvalue)) +
      geom_histogram(binwidth = 0.05) +
      theme_classic()
  }
  names(plot_list) <- names(deseq_results_list)
  plot_list
}

plot_volcano <- function(deseq_results_list,pal){
  names(pal) <- c("H", "HH", "HO", "O", "OO", "OH")
  plot_list <- list()
  # create for loop to produce ggplot2 graphs 
  for (i in 1:length(deseq_results_list$shrink)) { 
    trt_names <- str_split(names(deseq_results_list$shrink)[i], pattern = "_",simplify = T)[c(1,3)]
    shap <- c(22, 21)
    names(shap) <- c("H", "O")
    
    down <- unname(pal[names(pal) %in% trt_names[1]])
    up <- unname(pal[names(pal) %in% trt_names[2]])
    
    down_shape <- unname(shap[names(shap) %in% str_sub(trt_names[1], 1,1)])
    up_shape <- unname(shap[names(shap) %in% str_sub(trt_names[2], 1,1)])
    
    down_trt <- paste0("Up in ", trt_names[1])
    up_trt <- paste0("Up in ", trt_names[2])
    
    df <- deseq_results_list$shrink[[i]] %>% 
      as.data.frame() %>% 
      rownames_to_column("geneID") %>% 
      na.omit() %>% 
      filter(abs(log2FoldChange) <= 10) %>%
      mutate(DE = factor(
        if_else(geneID %in% deseq_results_list$down[[i]]$geneID, down_trt, 
                if_else(geneID %in% deseq_results_list$up[[i]]$geneID, 
                        up_trt, 
                        "Not DE")), 
        levels = c("Not DE", 
                   down_trt,
                   up_trt))) %>% 
      # arrange(desc(padj))
      arrange(desc(svalue))
    
    labels <- df %>% 
      filter(DE != "Not DE") %>% 
      # arrange(desc(pvalue)) %>% 
      arrange(desc(svalue)) %>% 
      tail(10)
    
    plot_list[[i]] <- ggplot(data = df, 
                             # aes(x = log2FoldChange, y = (-log10(pvalue)))) + 
                             aes(x = log2FoldChange, y = (-log10(svalue)))) +
      geom_point(aes(fill = DE, color = DE, size = DE, shape = DE), alpha = 0.75) +
      # geom_text_repel(data = labels, aes(x = log2FoldChange, y = (-log10(svalue)), label = geneID)) +
      scale_fill_manual(values = c("grey", down, up),
                        breaks = c("Not DE",down_trt, up_trt),
                        name = "Diff. expr.") +
      scale_color_manual(values = c("grey", "black", "black"), 
                         breaks = c("Not DE",down_trt, up_trt),
                         name = "Diff. expr.") +
      scale_size_manual(values = c(0.75, 1.5, 1.5),
                        breaks = c("Not DE",down_trt, up_trt),
                        name = "Diff. expr.") +
      scale_shape_manual(values = c(1, down_shape, up_shape),
                        breaks = c("Not DE",down_trt, up_trt),
                        name = "Diff. expr.") +
      annotate(geom = "text", x = -Inf, y = Inf, label = down_trt, hjust = -0.1, vjust = 1.1) + 
      annotate(geom = "text", x = Inf, y = Inf, label = up_trt, hjust = 1.1, vjust = 1.1) + 
      labs(x = expression(log[2]~"(Fold change)"), 
           y = expression(-log[10]~"(s-value)"),
           title = paste(trt_names, collapse =" vs. ")) +
      xlim(c(-10,10))+
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "none")
    
    names(plot_list)[i] <- names(deseq_results_list$shrink)[i]
  }
  plot_list
}


make_results_df <- function(deseq_results_list) {
  for (i in 1:length(deseq_results_list$shrink)){
    deseq_results_list$shrink[[i]] <- deseq_results_list$shrink[[i]] %>% 
      as.data.frame() %>% 
      rownames_to_column("geneID")
  }
  
  deseq_results_df <- dplyr::bind_rows( deseq_results_list$shrink, .id = "comp" )
  deseq_results_df
}

plot_MA <- function(deseq_results_list, pal) {
  names(pal) <- c("H", "HH", "HO", "O", "OO", "OH")
  plot_list <- list()
  # create for loop to produce ggplot2 graphs 
  for (i in 1:length(deseq_results_list$shrink)) { 
    trt_names <- str_split(names(deseq_results_list$shrink)[i], pattern = "_",simplify = T)[c(1,3)]
    
    shap <- c(22, 21)
    names(shap) <- c("H", "O")
    
    down <- unname(pal[names(pal) %in% trt_names[1]])
    up <- unname(pal[names(pal) %in% trt_names[2]])
    
    down_shape <- unname(shap[names(shap) %in% str_sub(trt_names[1], 1,1)])
    up_shape <- unname(shap[names(shap) %in% str_sub(trt_names[2], 1,1)])
    
    down_trt <- paste0("Up in ", trt_names[1])
    up_trt <- paste0("Up in ", trt_names[2])
    
    df <- deseq_results_list$shrink[[i]] %>% 
      as.data.frame() %>% 
      rownames_to_column("geneID") %>% 
      na.omit() %>% 
      # filter(abs(log2FoldChange) <= 10) %>%
      mutate(DE = factor(
        if_else(geneID %in% deseq_results_list$down[[i]]$geneID, down_trt, 
                if_else(geneID %in% deseq_results_list$up[[i]]$geneID, 
                        up_trt, 
                        "Not DE")), 
        levels = c("Not DE", 
                   down_trt,
                   up_trt)))
    
    plot_list[[i]] <- ggplot(data = df %>% 
                               # arrange(desc(padj)), 
                               arrange(desc(svalue)), 
                             aes(x = log10(baseMean), y = log2FoldChange)) + 
      geom_point(aes(fill = DE, color = DE, size = DE, shape = DE), alpha = 0.75) +
      scale_fill_manual(values = c("grey", down, up),
                        breaks = c("Not DE",down_trt, up_trt),
                        name = "Diff. expr.") +
      scale_color_manual(values = c("grey", "black", "black"), 
                         breaks = c("Not DE",down_trt, up_trt),
                         name = "Diff. expr.") +
      scale_size_manual(values = c(0.75, 1.5, 1.5),
                        breaks = c("Not DE",down_trt, up_trt),
                        name = "Diff. expr.") +
      scale_shape_manual(values = c(1, down_shape, up_shape),
                        breaks = c("Not DE",down_trt, up_trt),
                        name = "Diff. expr.") +
      labs(y = expression(log[2](Fold~change)), 
           x = expression(log[10](Base~mean)),
           title = paste(trt_names, collapse =" vs. ")) +
      # xlim(c(-10,10))+
      theme_bw() +
      theme(panel.grid = element_blank())
    
    names(plot_list)[i] <- names(deseq_results_list$shrink)[i]
  }
  panel <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 3, labels = "AUTO")
  panel
}


    
plot_scatter <- function(deseq_results_df, pal, comparisons, rev_x = FALSE)   {
  names(pal) <- c("H", "HH", "HO", "O", "OO", "OH")
  
  comps <- str_split(comparisons, pattern = "_", simplify = T)[,c(1,3)]
  
  comp_levels <- c( "Both",
                    paste0("Up in ", comps[1,1], " (vs ", comps[1,2], ") only"), # y-axis, 
                    paste0("Up in ", comps[1,2], " (vs ", comps[1,1], ") only"), # y-axis, 
                    paste0("Up in ", comps[2,1], " (vs ", comps[2,2], ") only"), # x-axis, 
                    paste0("Up in ", comps[2,2], " (vs ", comps[2,1], ") only"), # x-axis, 
                    "Not DE")
  
  pal1 <- c( "gold", 
             pal[names(pal) == comps[1,1]],
             pal[names(pal) == comps[1,2]],
             pal[names(pal) == comps[2,1]],
             pal[names(pal) == comps[2,2]],
             "grey")
  
  names(pal1) <- NULL
  
  shap <- c(22, 21)
  names(shap) <- c("H", "O")
  shap1 <- c( 23,
              shap[names(shap) == str_sub(comps[1,1], 1, 1)],
              shap[names(shap) == str_sub(comps[1,2], 1, 1)],
              shap[names(shap) == str_sub(comps[2,1], 1, 1)],
              shap[names(shap) == str_sub(comps[2,2], 1, 1)],
             1)
  
  names(shap1) <- NULL
  
  
  # pvals <- paste0("padj_", comparisons)
  sig <- paste0("svalue_", comparisons)
  lfcs <- paste0("log2FoldChange_", comparisons)
  
  df <- deseq_results_df %>% 
    dplyr::filter(comp %in% comparisons) %>% 
    # dplyr::select(geneID, comp, log2FoldChange, padj) %>% 
    dplyr::select(geneID, comp, log2FoldChange, svalue) %>%
    pivot_wider(id_cols = geneID, 
                names_from = comp, 
                # values_from = c(log2FoldChange, padj)) %>% 
                values_from = c(log2FoldChange, svalue)) %>% 
    mutate(DE = factor(
      case_when(
        .data[[sig[1]]] < 0.001  & .data[[sig[2]]] < 0.001  ~ "Both", 
        (.data[[sig[1]]] < 0.001 & .data[[lfcs[1]]] < 0) & .data[[sig[2]]] >= 0.01  ~ comp_levels[2], # y-axis
        (.data[[sig[1]]] < 0.001 & .data[[lfcs[1]]] > 0) & .data[[sig[2]]] >= 0.01  ~ comp_levels[3], # y-axis
        (.data[[sig[2]]] < 0.001 & .data[[lfcs[2]]] < 0) & .data[[sig[1]]] >= 0.01  ~ comp_levels[4], # x-axis
        (.data[[sig[2]]] < 0.001 & .data[[lfcs[2]]] > 0) & .data[[sig[1]]] >= 0.01  ~ comp_levels[5], # x-axis
        TRUE ~ "Not DE"), 
      levels = comp_levels))
  
  labs <- c(paste0(gsub("_", "~", comparisons[1]),"~log[2](fold~change)"), # x-axis
            paste0(gsub("_", "~", comparisons[2]),"~log[2](fold~change)") # y-axis
  )
  
  if(rev_x == TRUE) {
    df[,colnames(df) %in% lfcs[1]] <-  df[,colnames(df) %in% lfcs[1]] * -1
    labs[1] <- paste0(gsub("_", "~", comparisons[1]),"~log[2](fold~change)")
  }
  df %>% arrange(desc(DE)) %>% 
    filter(DE == "Both") %>% 
    write_tsv(paste0("00_data/02_DEgenes/Tcon_", comparisons[1], "_", comparisons[2], "_both_DE.tsv"))
  
  ggplot(df %>% arrange(desc(DE)), aes(.data[[lfcs[1]]],.data[[lfcs[2]]])) +
    geom_point(aes(fill = DE, color = DE, size = DE, shape = DE)) +
    scale_fill_manual(values = scales::alpha(pal1, 0.75), 
                      breaks = comp_levels)+
    scale_color_manual(values = c("black","black","black","black","black", "grey"), 
                       breaks = comp_levels) +
    scale_size_manual(values = c(2,2,2,2,2,0.5), 
                      breaks = comp_levels) +
    scale_shape_manual(values = shap1, 
                      breaks = comp_levels) +
    # xlim(c(-10, 10)) + ylim(c(-10, 10)) +
    labs(x =parse(text = labs[1]), 
         y = parse(text = labs[2])) +
    theme_bw() +
    theme(panel.grid = element_blank()
    )
}

#### Functional enrichment ####

topgo_genelist <- function(deseq_results_list, GO_path){
  
  # analyzing GO terms with at least 5 members,
  # as this yield more stable results.
  node_size= 5
  
  # use topGO to read the functional annotation
  geneID2GO <- readMappings(GO_path)
  
  # define the gene universe
  geneUniverse <- names(geneID2GO)
  
  temp0 <- list()
  dirs <- c("up", "down")
  for (k in 1:2){
    
    direction <-dirs[[k]]
    temp1 <- list()
    for (j in 1:length(deseq_results_list[[k+2]])){
      
      contrast=names(deseq_results_list[[k+2]])[j]
      
      genesOfInterest.bv <- deseq_results_list[[k+2]][[j]]$geneID
      
      
      geneList.bv <- factor(as.integer(geneUniverse %in% genesOfInterest.bv))
      
      if(sum(as.numeric(as.character(geneList.bv))) < 5) {
        temp1[[j]] <- data.frame(
          comp = factor(contrast, levels = names(deseq_results_list[[k+2]])),
          direction = factor(direction, levels = dirs))
      } else {
        names(geneList.bv) <- geneUniverse
    
        categories <- c("BP", "MF")
        temp2 <- list()
        for (i in 1:2) {
          GO_category=categories[i]
          # Create a new topGO GSEA objected
          myGOdata.bv <- new("topGOdata",
                             description="Candidate genes", 
                             ontology=GO_category, 
                             allGenes=geneList.bv, 
                             annot = annFUN.gene2GO, 
                             gene2GO = geneID2GO, 
                             nodeSize = node_size, 
          )  
          resultClassic <- runTest(myGOdata.bv, algorithm="classic", statistic="fisher")
          
          # Fisher with parent-child algorithm
          resultParentchild <- runTest(myGOdata.bv, algorithm="parentchild", statistic="fisher")
          
          # Create a GSEA results table
          # GenTable() is  annoying and requires a number of terms we want reported (default is 10)
          ## so we leverage our topGOresult object
          top_nodes <-  sum(resultParentchild@score < 0.05)
          
          temp2[[i]] <-  GenTable(
            myGOdata.bv,
            classicFisher = resultClassic,
            parentchildFisher = resultParentchild,
            orderBy = "parentchildFisher",
            topNodes = top_nodes + 1 ) %>% 
            mutate(cat = factor(GO_category, levels = categories),
                   comp = factor(contrast, levels = names(deseq_results_list[[k+2]])),
                   direction = factor(direction, levels = dirs))
        }
        temp1[[j]] <- compact(temp2) %>% purrr::reduce(full_join)
      }
    }
    temp0[[k]] <- compact(temp1) %>% purrr::reduce(full_join)
  }
  GO_df <- temp0 %>% 
    purrr::reduce(full_join)
  GO_df
}

plot_topgo <- function(GO_df, pal){
  dirs = c("up", "down")
  comps = levels(GO_df$comp)
  categories = c("BP", "MF")
  grid <- expand.grid(categories = categories, dirs =  dirs, comps = comps)
  
  names(pal) <- c("H", "HH", "HO", "O", "OO", "OH")
  shap <- c(22, 21)
  names(shap) <- c("H", "O")
  
  plot_list <- list()
  for (i in 1:dim(grid)[1]){
    
    trt_names <- str_split(grid$comps[i], pattern = "_",simplify = T)[c(1,3)]
    if(grid$dirs[i] == "up"){
      go_fill <- pal[names(pal) %in% trt_names[2]]
      go_shape <- shap[names(shap) %in% str_sub(trt_names[2],1,1)]
    } else if (grid$dirs[i] == "down") {
      go_fill <- pal[names(pal) %in% trt_names[1]]
      go_shape <- shap[names(shap) %in% str_sub(trt_names[1],1,1)]
    } else{
      go_fill <- "black"
      go_shape <- 1
    }
    go_fill <- unname(go_fill)
    go_shape <- unname(go_shape)
    
    go_analysis <-  GO_df %>% 
      dplyr::filter(comp == grid$comp[i], 
                    direction == grid$dirs[i],
                    cat == grid$categories[i]) %>% 
      mutate(GO_term = paste(GO.ID, Term, sep = "-"),
             parentchildFisher = as.numeric(gsub("< ", "", parentchildFisher)),
             sig= if_else(parentchildFisher < 0.01, "p < 0.01", 
                          if_else(parentchildFisher <= 0.05, "p < 0.05", "Not sig.")))
    
    plot_title <- paste(go_analysis$comp[1], go_analysis$direction[1], go_analysis$cat[1], sep = "_")
    
    if(dim(go_analysis)[1] < 20){
      max_rank <- dim(go_analysis)[1] 
    } else {
      max_rank <- 20
    }
    
    plot_list[[i]] <- ggplot(
      go_analysis[1:max_rank,], 
      aes( x= -log10(parentchildFisher), 
           y = forcats::fct_reorder(GO_term, -log10(parentchildFisher)), 
           fill = sig)) +
      geom_point(aes(size = Annotated), shape = go_shape, alpha = 0.75) +
      # geom_vline(xintercept = c(-log10(0.05), -log10(0.01)),
      #            linetype = c("dotted", "dashed")) +
      labs(x = expression(-log[10](adj.~p-value)), 
           y = "GO term and function", title = plot_title ) +
      scale_fill_manual(values = c(go_fill, "grey55", "grey"), 
                        breaks = c("p < 0.01", "p < 0.05", "Not sig."),
                        guide = "none") +
      theme_classic() +
      theme(legend.position = c(1,0), 
            legend.justification = c("right", "bottom"))
    
    fig_height  <- 2 + 0.15 * max_rank
    
    ggsave(paste("04_targetPlots/topgo_plots/topgo", plot_title,"plot.pdf", sep = "_"), height = fig_height, width = 7 )
  }
  plot_list
}

topgo_genelist_both <- function(deseq_results_list, GO_path){
  
  # analyzing GO terms with at least 5 members,
  # as this yield more stable results.
  node_size= 5
  
  # use topGO to read the functional annotation
  geneID2GO <- readMappings(GO_path)
  
  # define the gene universe
  geneUniverse <- names(geneID2GO)
  
  temp1 <- list()

  for (j in 1:length(deseq_results_list[[1]])){
    # for (j in 1:2){
      
      contrast=names(deseq_results_list[[1]])[j]
      
      genesOfInterest.bv <- unique(c(deseq_results_list[[3]][[j]]$geneID, 
                                     deseq_results_list[[4]][[j]]$geneID))
      
      
      geneList.bv <- factor(as.integer(geneUniverse %in% genesOfInterest.bv))
      
      if(sum(as.numeric(as.character(geneList.bv))) < 5) {
        temp1[[j]] <- data.frame(
          comp = factor(contrast, levels = names(deseq_results_list[[1]])))
      } else {
        names(geneList.bv) <- geneUniverse
        
        categories <- c("BP", "MF")
        temp2 <- list()
        for (i in 1:2) {
          GO_category=categories[i]
          # Create a new topGO GSEA objected
          myGOdata.bv <- new("topGOdata",
                             description="Candidate genes", 
                             ontology=GO_category, 
                             allGenes=geneList.bv, 
                             annot = annFUN.gene2GO, 
                             gene2GO = geneID2GO, 
                             nodeSize = node_size, 
          )  
          resultClassic <- runTest(myGOdata.bv, algorithm="classic", statistic="fisher")
          
          # Fisher with parent-child algorithm
          resultParentchild <- runTest(myGOdata.bv, algorithm="parentchild", statistic="fisher")
          
          # Create a GSEA results table
          # GenTable() is  annoying and requires a number of terms we want reported (default is 10)
          ## so we leverage our topGOresult object
          top_nodes <-  sum(resultParentchild@score < 0.05)
          
          temp2[[i]] <-  GenTable(
            myGOdata.bv,
            classicFisher = resultClassic,
            parentchildFisher = resultParentchild,
            orderBy = "parentchildFisher",
            topNodes = top_nodes + 1 ) %>% 
            mutate(cat = factor(GO_category, levels = categories),
                   comp = factor(contrast, levels = names(deseq_results_list[[1]])))
        }
        temp1[[j]] <- compact(temp2) %>% purrr::reduce(full_join)
      }
    }
  GO_df <- temp1%>% 
    purrr::reduce(full_join)
  GO_df
}

plot_topgo_both <- function(GO_df, pal){
  comps = levels(GO_df$comp)
  categories = c("BP", "MF")
  grid <- expand.grid(categories = categories, comps = comps)
  
  names(pal) <- c("H", "HH", "HO", "O", "OO", "OH")
  
  plot_list <- list()
  for (i in 1:dim(grid)[1]){
    
    trt_names <- str_split(grid$comps[i], pattern = "_",simplify = T)[c(1,3)]
    go_fill <- "black"
    
    go_analysis <-  GO_df %>% 
      dplyr::filter(comp == grid$comp[i], 
                    cat == grid$categories[i]) %>% 
      mutate(GO_term = paste(GO.ID, Term, sep = "-"),
             parentchildFisher = as.numeric(gsub("< ", "", parentchildFisher)),
             sig= if_else(parentchildFisher < 0.01, "p < 0.01", 
                          if_else(parentchildFisher <= 0.05, "p < 0.05", "Not sig.")))
    
    plot_title <- paste(go_analysis$comp[1], go_analysis$cat[1], sep = "_")
    
    if(dim(go_analysis)[1] < 20){
      max_rank <- dim(go_analysis)[1] 
    } else {
      max_rank <- 20
    }
    
    plot_list[[i]] <- ggplot(
      go_analysis[1:max_rank,], 
      aes( x= -log10(parentchildFisher), 
           y = forcats::fct_reorder(GO_term, -log10(parentchildFisher)), 
           fill = sig)) +
      geom_point(aes(size = Annotated), shape = 23, alpha = 0.75) +
      # geom_vline(xintercept = c(-log10(0.05), -log10(0.01)),
      #            linetype = c("dotted", "dashed")) +
      labs(x = expression(-log[10](adj.~p-value)), 
           y = "GO term and function", title = plot_title ) +
      scale_fill_manual(values = c(go_fill, "grey55", "grey"), 
                        breaks = c("p < 0.01", "p < 0.05", "Not sig."),
                        guide = "none") +
      theme_classic() +
      theme(legend.position = c(1,0), 
            legend.justification = c("right", "bottom"))
    
    fig_height  <- 2 + 0.15 * max_rank
    
    ggsave(paste("04_targetPlots/topgo_plots/topgo", plot_title,"plot.pdf", sep = "_"), height = fig_height, width = 7 )
  }
  plot_list
}

topgo_targeted <- function(genelist, GO_path, name){
  # analyzing GO terms with at least 5 members,
  # as this yield more stable results.
  node_size= 5
  # use topGO to read the functional annotation
  geneID2GO <- readMappings(GO_path)
  
  # define the gene universe
  geneUniverse <- names(geneID2GO)
  
  genesOfInterest.bv <- genelist$geneID
  geneList.bv <- factor(as.integer(geneUniverse %in% genesOfInterest.bv))
  
  names(geneList.bv) <- geneUniverse
  
  categories <- c("BP", "MF")
  temp2 <- list()
  for (i in 1:2) {
    GO_category=categories[i]
    # Create a new topGO GSEA objected
    myGOdata.bv <- new("topGOdata",
                       description="Candidate genes", 
                       ontology=GO_category, 
                       allGenes=geneList.bv, 
                       annot = annFUN.gene2GO, 
                       gene2GO = geneID2GO, 
                       nodeSize = node_size, 
    )  
    resultClassic <- runTest(myGOdata.bv, algorithm="classic", statistic="fisher")
    
    # Fisher with parent-child algorithm
    resultParentchild <- runTest(myGOdata.bv, algorithm="parentchild", statistic="fisher")
    
    # Create a GSEA results table
    # GenTable() is  annoying and requires a number of terms we want reported (default is 10)
    ## so we leverage our topGOresult object
    top_nodes <-  sum(resultParentchild@score < 0.05)
    
    temp2[[i]] <-  GenTable(
      myGOdata.bv,
      classicFisher = resultClassic,
      parentchildFisher = resultParentchild,
      orderBy = "parentchildFisher",
      topNodes = top_nodes + 1 ) %>% 
      mutate(name = name,
             cat = factor(GO_category, levels = categories))
  }
  GO_df <- temp2 %>% 
    purrr::reduce(full_join)
  GO_df
}

plot_topgo_targeted <- function(GO_df, color, name){
  categories <- c("BP", "MF")
  plot_list <- list()
  for (i in 1:2){
    go_fill <- color
    
    go_analysis <-  GO_df %>% 
      dplyr::filter(cat == categories[i]) %>% 
      mutate(GO_term = paste(GO.ID, Term, sep = "-"),
             parentchildFisher = as.numeric(gsub("< ", "", parentchildFisher)),
             sig= if_else(parentchildFisher < 0.01, "p < 0.01", 
                          if_else(parentchildFisher <= 0.05, "p < 0.05", "Not sig.")))
    
    plot_title <- paste(name, go_analysis$cat[1], sep = "_")
    
    if(dim(go_analysis)[1] < 20){
      max_rank <- dim(go_analysis)[1] 
    } else {
      max_rank <- 20
    }
    
    plot_list[[i]] <- ggplot(
      go_analysis[1:max_rank,], 
      aes( x= -log10(parentchildFisher), 
           y = forcats::fct_reorder(GO_term, -log10(parentchildFisher)), 
           fill = sig)) +
      geom_point(aes(size = Annotated), shape = 23, alpha = 0.75) +
      # geom_vline(xintercept = c(-log10(0.05), -log10(0.01)),
      #            linetype = c("dotted", "dashed")) +
      labs(x = expression(-log[10](adj.~p-value)), 
           y = "GO term and function", title = plot_title ) +
      scale_fill_manual(values = c(go_fill, "grey55", "grey"), 
                        breaks = c("p < 0.01", "p < 0.05", "Not sig."),
                        guide = "none") +
      theme_classic() +
      theme(legend.position = c(1,0.01), 
            legend.justification = c("right", "bottom"),
            axis.text.y = element_text(size = 6))
    
    fig_height  <- 2 + 0.15 * max_rank
    
    ggsave(paste("04_targetPlots/topgo_plots/topgo", plot_title,"plot.pdf", sep = "_"), height = fig_height, width = 7 )
  }
  plot_list
}


plot_semsim <- function(GO_df, pal){
  TcGO_bp <- godata('org.Tconura.eg.db', keytype = "REFSEQ", ont = "BP")
  
  dirs = c("up", "down")
  comps = levels(GO_df$comp)
  categories = c("BP", "MF")
  grid <- expand.grid(categories = categories, dirs =  dirs, comps = comps)
  
  names(pal) <- c("H", "HH", "HO", "O", "OO", "OH")
  
  
  
  plot_list <- list()
  for (i in 1:dim(grid)[1]){
    
    trt_names <- str_split(grid$comps[i], pattern = "_",simplify = T)[c(1,3)]
    if(grid$dirs[i] == "up"){
      go_fill <- pal[names(pal) %in% trt_names[2]]
      go_shape <- shap[names(shap) %in% str_sub(trt_names[2],1,1)]
    } else if (grid$dirs[i] == "down") {
      go_fill <- pal[names(pal) %in% trt_names[1]]
      go_shape <- shap[names(shap) %in% str_sub(trt_names[2],1,1)]
    } else{
      go_fill <- "black"
      go_shape <- 1
    }
    go_fill <- unname(go_fill)
    go_shape <- unname(go_shape)
    
    
    go_analysis <-  GO_df %>% 
      dplyr::filter(comp == grid$comp[i], 
                    direction == grid$dirs[i],
                    cat == grid$categories[i]) %>% 
      mutate(GO_term = paste(GO.ID, Term, sep = "-"),
             parentchildFisher = as.numeric(gsub("< ", "", parentchildFisher)),
             sig= if_else(parentchildFisher < 0.01, "p < 0.01", 
                          if_else(parentchildFisher <= 0.05, "p < 0.05", "Not sig.")))
    
    plot_title <- paste(go_analysis$comp[1], go_analysis$direction[1], go_analysis$cat[1])
    
    if(dim(go_analysis)[1] < 20){
      max_rank <- dim(go_analysis)[1] 
    } else {
      max_rank <- 20
    }
    
    plot_list[[i]] <- ggplot(
      go_analysis[1:max_rank,], 
      aes( x= -log10(parentchildFisher), 
           y = forcats::fct_reorder(GO_term, -log10(parentchildFisher)), 
           fill = sig)) +
      geom_point(aes(size = Annotated), shape = go_shape, alpha = 0.75) +
      # geom_vline(xintercept = c(-log10(0.05), -log10(0.01)),
      #            linetype = c("dotted", "dashed")) +
      labs(x = expression(-log[10](adj.~p-value)), 
           y = "GO term and function", title = plot_title ) +
      scale_fill_manual(values = c(go_fill, "grey55", "grey"), 
                        breaks = c("p < 0.01", "p < 0.05", "Not sig."),
                        guide = "none") +
      theme_classic() +
      theme(legend.position = c(1,0), 
            legend.justification = c("right", "bottom"))
    
    fig_height  <- 2 + 0.15 * max_rank
    
    ggsave(paste("04_targetPlots/topgo_plots/topgo_", plot_title,"plot.pdf", collapse = "_"), height = fig_height, width = 7 )
  }
  plot_list
}

save_genesets <- function(deseq_results_list){
  
  for(i in 1:length(deseq_results_list$up)){
    trt_names <- str_split(names(deseq_results_list$up)[i], pattern = "_",simplify = T)[c(1,3)]
    
    sub <- deseq_results_list$up[[i]]
    write_delim(sub, 
                file = paste0("00_data/02_DEgenes/Tcon_" , names(deseq_results_list$up)[i], 
                              "_up_in_", trt_names[2],".txt")) 
    
    sub <- deseq_results_list$down[[i]]
    write_delim(sub, 
                file = paste0("00_data/02_DEgenes/Tcon_STRINGinput_", names(deseq_results_list$up)[i],
                              "_up_in_", trt_names[1],".txt")) 
  }
}

#### Save plots ####

pca_panel <- function(plotlist, range) {
  
  plots <- plotlist[range]
  ncol <- length(range)
  panel <- ggarrange(plotlist = plots, ncol = ncol, 
                     nrow = 1, 
                     legend = "bottom", 
                     common.legend = T, 
                     labels = "AUTO")
}

save_volc_list <- function(list, path) {
  for(i in 1:length(list)){
    plot <- list[[i]] + theme(legend.position = "none")
    name <- names(list)[i]
    ggsave(plot = plot, file = paste0(path, "tcon_volc_", name, ".pdf"), width = 2.75, height = 3)
  }
}

#### WGCNA #### 
standardiseGenes <- function(dat){
  tmp <- tempfile()
  write.table(dat,file=tmp, sep='\t', quote = F, col.names=NA)
  
  #read it back in as an expression set
  dat <- table2eset(file=tmp)
  
  dat <- standardise(dat)
  
  scaledata <- dat@assayData$exprs
  scaledata
}

rename_object <- function(dat){
 dat
}


get_datTraits <- function(mat, samples){
  traitData = samples
  # Form a data frame analogous to expression data that will hold the clinical traits.
  Samples = rownames(t(mat));
  traitRows = match(Samples, traitData$Sample);
  datTraits = traitData[traitRows, ];
  rownames(datTraits) = rownames(traitData)
  datTraits
}

# Re-cluster samples
plot_dendro <- function(dat, datTraits, pal) {
  # cluster samples
  sampleTree2 = hclust(dist(t(dat)), method = "average")
  
  # Form a data frame analogous to expression data that will hold the traits.
  pal_75 <- scales::alpha(pal, 0.75)
  
  traitColors = data.frame(
    Condition = pal_75[as.numeric(datTraits$Treatment)],
    `Host race` = pal_75[c(1,4)][as.numeric(as.factor(datTraits$Host_form))])
  
  plot <- WGCNA::plotDendroAndColors(sampleTree2, traitColors,
                                     groupLabels = c("Condition", "Host race"),
                                     main = "Sample dendrogram and trait heatmap")
  
  print_plots <- function(x) { 
    plts <- grepl("\\.plt", names(x))
  }
  print_plots(plot)

}

pick_threshold <- function(dat){
  set.seed(2023)
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(t(dat), 
                          networkType = "signed", 
                          powerVector = powers,
                          verbose = 5,
                          RsquaredCut = 0.8)
  sft
  # sft should maximize corr, minimize connectivity
}

plot_soft_threshold <- function(sftThreshold) {
  set.seed(2023)
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  
  # Plot the results:
  p1 <- ggplot(sftThreshold$fitIndices, aes(x = Power, y = -sign(slope)* SFT.R.sq) ) +
    geom_hline(yintercept = 0.8, color = "red") +
    geom_text(aes(label = Power)) + 
    theme_classic() +
    labs(x="Soft Threshold (power)",
         y=expression("Scale Free Topology Model Fit,signed"~R^2),
         title ="Scale independence") 
  
  p2 <- ggplot(sftThreshold$fitIndices, aes(x = Power, y = mean.k. )) +
    geom_text(aes(label = Power)) + 
    theme_classic() +
    labs(x="Soft Threshold (power)",
         y="Mean connectivity",
         title ="Mean connectivity" )
  
  panel <- ggarrange(p1, p2, ncol = 2, nrow = 1, labels = "AUTO")
  panel
}

run_modules <- function(dat, sftThreshold) {
  threshold=sftThreshold$powerEstimate
  net =  blockwiseModules(t(dat),
                          randomSeed =123,
                          power = threshold,
                          corType = "bicor",
                          networkType = "signed", 
                          TOMType = "signed",
                          minModuleSize = 30,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.3,
                          numericLabels = TRUE,
                          pamRespectsDendro = FALSE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "00_data/Tcon_transcriptomeTOM", 
                          verbose = 3)
  net
}
load_modules <- function(path){
  load(path)
  tcon_network
}

plot_network <- function(network){
  # open a graphics window
  par(mfrow = c(3,1));
  sizeGrWindow(9, 12)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(network$colors)
  
  # Plot the dendrogram and the module colors underneat
  pdf("04_targetPlots/tcon_networkModules_plot.pdf", height = 6, width = 8)

  plotDendroAndColors(network$dendrograms[[1]], mergedColors[network$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  plotDendroAndColors(network$dendrograms[[2]], mergedColors[network$blockGenes[[2]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  plotDendroAndColors(network$dendrograms[[3]], mergedColors[network$blockGenes[[3]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}

pull_modules <- function(network, data, datTraits){
  # module <- list()
  # module$nSamples <- nrow(t(data))
  # module$labels <- network$colors
  # module$colors<- labels2colors(network$colors)
  # module$MEs = network$MEs;
  # # Recalculate MEs with color labels
  # module$MEs0 <- moduleEigengenes(t(data),module$colors)$eigengenes
  # module$MEs <- orderMEs(module$MEs0)
  # module$labels_numeric <- paste0("mod", seq(1,dim(module$MEs)[2]))
  # names(module$labels_numeric) <- names(module$MEs)
  module <- list()
  module$nSamples <- nrow(t(data))
  module$genes <- names(network$colors)
  module$labels <- network$colors
  module$colors<- labels2colors(network$colors)
  module$labels[str_count(module$labels) == 1] <- 
    paste0("0",  module$labels[str_count(module$labels) == 1] )
  module$labels <- paste("mod", module$labels, sep = "")
  
  # generate color labels
  module$colors <- labels2colors(network$colors)
  
  # pull eigengenes
  module$MEs <- orderMEs(
    moduleEigengenes(t(data),network$colors)$eigengenes)
  colnames(module$MEs)[str_count(colnames(module$MEs)) == 4] <- 
    paste("mod",
          gsub("ME", "", colnames(module$MEs)[str_count(colnames(module$MEs)) == 4]), 
          sep = "")
  colnames(module$MEs)[str_count(colnames(module$MEs)) == 3] <- 
    paste("mod",
          gsub("ME", "",  colnames(module$MEs)[str_count(colnames(module$MEs)) == 3]), 
          sep = "0")
  module$freq_df <- data.frame(mod = module$labels, 
                               color = module$colors) %>% 
    group_by(mod,color) %>% 
    summarise(freq = n())
  
  module$moduleHostCor <- WGCNA::cor(module$MEs, as.numeric(as.factor(datTraits$Host_form=="het")), use = "p")
  module$moduleHostPvalue = corPvalueStudent(module$moduleHostCor, module$nSamples)
  module$moduleHostPvalueAdj = matrix(p.adjust(module$moduleHostPvalue, method = "BH"), 
                                      nrow = dim(module$moduleHostPvalue)[1], ncol = 1)
  rownames(module$moduleHostPvalueAdj) = names(module$MEs)
  colnames(module$moduleHostPvalueAdj) = names(datTraits$Treatment)
  module$moduleStressCor <- WGCNA::cor(module$MEs, datTraits$Stress, use = "p")
  module$moduleStressPvalue = corPvalueStudent(module$moduleStressCor , module$nSamples)
  module$moduleStressPvalueAdj = matrix(p.adjust(module$moduleStressPvalue, method = "BH"), 
                                        nrow = dim(module$moduleStressPvalue)[1], ncol = 1)
  rownames(module$moduleStressPvalueAdj) = names(module$MEs)
  colnames(module$moduleStressPvalueAdj) = names(datTraits$Treatment)
  module$modulePlasticityCor <- WGCNA::cor(module$MEs, datTraits$Perf_plasticity, use = "p")
  module$modulePlasticityPvalue = corPvalueStudent(module$modulePlasticityCor , module$nSamples)
  module$modulePlasticityPvalueAdj = matrix(p.adjust(module$modulePlasticityPvalue, method = "BH"), 
                                            nrow = dim(module$modulePlasticityPvalue)[1], ncol = 1)
  rownames(module$modulePlasticityPvalueAdj) = names(module$MEs)
  colnames(module$modulePlasticityPvalueAdj) = names(datTraits$Treatment)
  module$moduleCHplasticityCor <- WGCNA::cor(module$MEs, datTraits$CH_plasticity, use = "p")
  module$moduleCHplasticityPvalue = corPvalueStudent(module$moduleCHplasticityCor , module$nSamples)
  module$moduleCHplasticityPvalueAdj = matrix(p.adjust(module$moduleCHplasticityPvalue, method = "BH"), 
                                              nrow = dim(module$moduleCHplasticityPvalue)[1], ncol = 1)
  rownames(module$moduleCHplasticityPvalueAdj) = names(module$MEs)
  colnames(module$moduleCHplasticityPvalueAdj) = names(datTraits$Treatment)
  module$moduleCOplasticityCor <- WGCNA::cor(module$MEs, datTraits$CO_plasticity, use = "p")
  module$moduleCOplasticityPvalue = corPvalueStudent(module$moduleCOplasticityCor , module$nSamples)
  module$moduleCOplasticityPvalueAdj = matrix(p.adjust(module$moduleCOplasticityPvalue, method = "BH"), 
                                              nrow = dim(module$moduleCOplasticityPvalue)[1], ncol = 1)
  rownames(module$moduleCOplasticityPvalueAdj) = names(module$MEs)
  colnames(module$moduleCOplasticityPvalueAdj) = names(datTraits$Treatment)
  module
}

draw_cor_heatmap <- function(module, name, pal) {
  Cor=paste0("module", name, "Cor")
  Padj=paste0("module", name, "PvalueAdj")
  
  testMatrix =  paste(signif(module[[Cor]], 2), " (",
                      signif(module[[Padj]], 1), ")", sep = "");
  dim(testMatrix) = dim(module[[Cor]])
  rownames(testMatrix) = names(module$MEs)
  colnames(testMatrix) = name
  
  if(name == "Host") {
    col_fun = colorRamp2(c(-0.75, 0, 0.75), c(pal[4],"white", pal[1]))
  } else if (name == "Stress") {
    col_fun = colorRamp2(c(-0.75, 0, 0.75), c("#f05039","white", "#f05039"))
  } else if (name == "Plasticity") {
    col_fun = colorRamp2(c(-0.75, 0, 0.75), c("goldenrod","white", "goldenrod"))
  } else if (name == "CHplasticity") {
    col_fun = colorRamp2(c(-0.75, 0, 0.75), c(pal[3],"white", pal[3]))
  } else if (name == "COplasticity") {
    col_fun = colorRamp2(c(-0.75, 0, 0.75), c(pal[6],"white", pal[6]))
  }
  

  ht1 <- ComplexHeatmap::Heatmap(module[[Cor]], 
                 col = col_fun, 
                 show_row_dend = FALSE,
                 show_column_names = TRUE, 
                 column_names_side = "top",
                 show_row_names = T,
                 row_names_gp = gpar(fontsize = 9),
                 row_names_side = "left",
                 row_labels = rownames(module[[Cor]]),
                 name = paste0("Correlation\nwith ", name),
                 cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                   grid::grid.text(testMatrix[i, j], x, y, gp = grid::gpar(fontsize = 7))
                 },
                 width = 2) 
  ht1
}


draw_module_heatmap <- function(module, datTraits, pal) {
  ht_ecotype <- datTraits$Host_form
  names(ht_ecotype) <- datTraits$Sample
  col_ecotype<- pal[c(4,1)]
  names(col_ecotype) <- unique(datTraits$Host_form)
  
  ht_condition <- as.factor(datTraits$Treatment)
  names(ht_condition) <- datTraits$Sample
  col_condition <- pal
  names(col_condition) <- levels(datTraits$Treatment)
  
  ht2 <- Heatmap(t(module$MEs),
                 top_annotation = HeatmapAnnotation(
                   Treatment = ht_condition, 
                   `Host race` = ht_ecotype,
                   col = list(`Host race` = col_ecotype, 
                              Treatment = col_condition), 
                   show_annotation_name = T
                 ),
                 name = "Module\neigengenes", column_dend_height = unit(2, "cm"),
                 show_row_dend = FALSE,
                 show_column_names = T,
                 row_labels = rownames(t(module$MEs)),
                 row_names_side = "left",
                 row_names_gp = gpar(fontsize = 9),
                 # show_row_names = FALSE,
                 column_title = NULL,
                 clustering_method_columns = "ward.D2",
                 column_names_gp = gpar(fontsize = 7),
                 width = 10)
  
  ht2
}

draw_cor_module_combo_plot <- function(ht1, ht2, ht3, ht4, ht5, ht6) {
  ht_list = ht1 + ht2 + ht3 + ht4 +ht5+ ht6
  draw(ht_list,  row_km = 1, cluster_rows = TRUE, show_row_dend = FALSE)
}

get_moduleMembership <- function(scaledata, module, datTraits){
  clusters <- data.frame(cluster = module$freq_df$mod,
                         clusterCol = module$freq_df$color)  
  
  geneModuleMembership = as.data.frame(cor(t(scaledata), module$MEs, use = "p"));
  
  WGCNA_membership <- geneModuleMembership %>% rownames_to_column("gene") %>% 
    gather(-c(gene), key = cluster, value = membership) %>% 
    left_join(clusters) %>% 
    mutate(clusterCol = gsub("ME", "", clusterCol))

  data <-  data.frame(gene = module$genes,
                      cluster = module$labels,
                      clusterCol =  module$colors) %>%
    left_join(rownames_to_column(as.data.frame(scaledata), "gene")) %>% 
    gather(-c(gene, clusterCol, cluster), key = Sample, value = value) %>% 
    right_join(WGCNA_membership) %>%
    na.omit() %>%
    # filter(membership > 0.3) %>%
    group_by(clusterCol,Sample) %>% 
    mutate(centers = median(value, na.rm = T)) %>% 
    left_join(datTraits) 
  data
}

plot_moduleMembership <- function(data, pal){
  ggplot(data, aes(x=Sample, y=value)) + 
    facet_grid(cluster~Treatment, scales = "free", space = "free") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_line(aes(group=gene, alpha = membership/3, color = Treatment), show.legend = FALSE) +
    # scale_colour_gradientn(colours=c('blue1','red2')) +
    #this adds the core 
    geom_point(aes(Sample,centers, group=cluster, shape = Host_form), 
               color="black",fill = "white", inherit.aes=FALSE, size = 1) +
    xlab("Sample") +
    ylab("Normalized expression") +
    labs(color = "Score") +
    scale_color_manual(values = pal)+
    scale_shape_manual(values = c(22, 21), guide = "none") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5), 
          panel.grid = element_blank() , 
          strip.text.x = element_text(face = "bold"),
          strip.background = element_blank(), 
          strip.placement = "outside", 
          strip.text = element_text(face = "bold", angle = 0, vjust = 1)) 
}

summarize_moduleMembership <- function(data){
  dat <- data %>% 
    dplyr::select(Sample, clusterCol, cluster, Host_form, Treatment, Stress, Perf_plasticity, centers) %>% 
    group_by(Sample, clusterCol, cluster, Host_form, Treatment, Stress, Perf_plasticity) %>% 
    mutate(n_genes = n()) %>% 
    unique()
  dat
}

plot_moduleMembershipSumm <- function(dat, pal){
 
  mods <- levels(as.factor(dat$cluster))
  plotlist <- list()
  for (i in 1:length(mods)){
    plotlist[[i]] <- ggplot(dat[dat$cluster == mods[i],], aes(x=Treatment, y=centers)) + 
      # geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      geom_point(aes(color = Treatment, fill = Treatment, shape = Treatment), size = 1, alpha = 0.5, guide = "none") +
      stat_summary(aes(fill = Treatment, shape = Treatment), fun.data = "mean_cl_boot", colour = "black", linewidth = 1, size = 0.5) +
      labs(x = "Treatment", y = "Norm. Expr.", title = mods[i]) + 
      labs(color = "Score") +
      scale_fill_manual(values = pal)+
      scale_color_manual(values = pal, guide = "none")+
      scale_shape_manual(values = c(22,22,22, 21,21,21))+
      scale_y_continuous(breaks = seq(from = 1.5, to = 10, by = 0.5)) +
      theme_bw() +
      theme(
        panel.grid = element_blank() , 
        axis.title = element_blank(),
        plot.title = element_text(face = "bold", angle = 0, vjust = 1, size = 10))
  }
  plotlist
}

plot_moduleEigengenes <- function(datTraits, module, pal){
  dat <- t(module$MEs) %>% as.data.frame() %>% 
    rownames_to_column("mod") %>% 
    left_join(module$freq_df) %>% 
    pivot_longer (-c(mod, color, freq), names_to = "Sample", values_to = "ME") %>% 
    left_join(datTraits)
  mods <- sort(unique(dat$mod))
  plotlist <- list()
  for (i in 1:length(mods)){
    plotlist[[i]] <- ggplot(dat[dat$mod == mods[i],], aes(x=Treatment, y=ME)) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      geom_point(aes(color = Treatment, fill = Treatment, shape = Treatment), size = 1, alpha = 0.5) +
      stat_summary(aes(fill = Treatment, shape = Treatment), fun.data = "mean_cl_boot", colour = "black", linewidth = 1, size = 0.5) +
      labs(x = "Treatment", y = "Eigengene", title = mods[i]) + 
      labs(color = "Score") +
      scale_fill_manual(values = pal)+
      scale_color_manual(values = pal, guide = "none")+
      scale_shape_manual(values = c(22,22,22, 21,21,21))+
      scale_y_continuous(breaks = seq(from = -0.5, to = 0.7, by = 0.1)) +
      theme_bw() +
      theme(legend.position = "none", 
        panel.grid = element_blank() , 
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", angle = 0, vjust = 1, size = 10))
  }
  names(plotlist) <- mods
  plotlist
}

extract_mod_genes <- function(dat, modList){
  mod_genes <- list()
  for (i in 1: length(modList)){
    mod_genes[[i]] <- unique(dplyr::filter(dat, cluster == modList[i]) %>% 
                               ungroup() %>% 
                             dplyr::select(geneID = gene))
  }
  names(mod_genes) <- modList
  mod_genes
}

topgo_modList <- function(geneList, path, modList){
  topGO_list <- list()
  for(i in 1:length(geneList)){
    topGO_list[[i]] <- topgo_targeted(geneList[[i]], path, modList[i])
  }
 names(topGO_list) <- modList
 topGO_list
}

plot_topgo_modList <- function(topgoList, modList, palList){
  plot_list <- list()
  for(i in 1:length(topgoList)){
    plot_list[[i]] <- plot_topgo_targeted(topgoList[[i]], palList[i], modList[i])
  }
  names(plot_list) <- modList
  plot_list
}

#### Population genomics ####

get_scaffold_master <- function(file){
  file %>% 
    group_by(scaffold, ncontigs_scaffold, nbp_scaffold, color_scaffold) %>% 
    summarise(start = min(contig_start_continuous), 
              end = max(contig_end_continuous)) %>% 
    arrange(start) %>% 
    rownames_to_column("scaffold_order") %>% 
    mutate(color_scaffold = if_else(nbp_scaffold < 1e7, "white", color_scaffold))
}

get_fst_list <- function(file_list, master) {
  window_size = 5e4
  fst_cols <- c("region",	
                "chr",	
                "midPos",	
                "Nsites",	
                "Fst")
  
  fst_read <- function(x) {
    read_tsv(x, col_names = fst_cols, skip = 1) %>% 
      mutate(chr = factor(chr)) %>% 
      arrange(chr, midPos) %>% 
      left_join(master) %>% 
      arrange(contig_order, midPos) %>% 
      mutate(midPos_scaffold = if_else(contig_dir == "+", midPos, contig_length - midPos),
             midPos_contigContinuous = midPos_scaffold + contig_start_continuous,
             uniq = paste(contig_order, midPos, sep = "_")) %>% 
      filter(Nsites >= 0.2* window_size)
  }
  list <- lapply(file_list, fst_read)
  names(list) <- str_sub( gsub("05_popGen/", "", file_list), 1,9)
  list
}


plot_fst <- function(list,  s) {
  plotlist <- list()
  for (x in 1:length(list)) {
    outlier_cutoff <- mean(list[[x]]$Fst[list[[x]]$keepMicrobial != F], na.rm = T) + 
      sd(list[[x]]$Fst[list[[x]]$keepMicrobial != F], na.rm = T)*3
    
    plotlist[[x]] <-  ggplot() +
      geom_rect(data = s, alpha = 0.5, 
                aes(xmin = start, xmax = (end - 1), fill = scaffold), ymin = -12, ymax = 12) + 
      scale_fill_manual(breaks = s$scaffold, values = s$color_scaffold, guide = "none") +
      geom_hline(yintercept =  outlier_cutoff, linetype = "dashed") +
      geom_point(data = list[[x]][list[[x]]$keepMicrobial != F & 
                                    list[[x]]$Fst < outlier_cutoff,], 
                 aes(x = midPos_contigContinuous, y = Fst, size = Nsites), 
                 color = "grey55", fill = "grey55", shape = 21, alpha = 0.5) +
      geom_point(data = list[[x]][list[[x]]$keepMicrobial != F & 
                                    list[[x]]$Fst >= outlier_cutoff,], 
                 aes(x = midPos_contigContinuous, y = Fst, size = Nsites), 
                 color = "black", fill = "black", shape = 21, alpha = 0.5)+ 
      annotate("text", x = Inf, y = Inf,  vjust = 1.1, hjust = 1.1, label = "Unplaced\nscaffolds") +
      scale_size_continuous(range = c(0, 2), guide = "none" ) +
      labs(x = "SNP position", y = "Fst") +
      theme_bw() +theme( panel.grid = element_blank())+
      # ylim(c(-0.05, 1))+
      scale_x_continuous( expand = c(0, 0))
  }
  names(plotlist) <- names(list)
  plotlist
}

get_dxy_list <- function(file_list, master) {
  window_size = 5e4
  dxy_cols <- c("chr",	
                "start",	
                "end",
                "Nsites",
                "nVariants",	
                "dxy")
  
  dxy_read <- function(x) {
    read_tsv(x, col_names = dxy_cols, skip = 1) %>% 
      mutate(chr = factor(chr),
             midPos = start + window_size/2) %>% 
      arrange(chr, midPos) %>% 
      left_join(master) %>% 
      arrange(contig_order, midPos) %>% 
      mutate(midPos_scaffold = if_else(contig_dir == "+", midPos, contig_length - midPos),
             midPos_contigContinuous = midPos_scaffold + contig_start_continuous,
             uniq = paste(contig_order, midPos, sep = "_")) %>% 
      filter(Nsites >= 0.2* window_size)
  }
  list <- lapply(file_list, dxy_read)
  names(list) <- str_sub( gsub("05_popGen/", "", file_list), 1,9)
  list
}


plot_dxy <- function(list,  s) {
  plotlist <- list()
  for (x in 1:length(list)) {
    outlier_cutoff <- mean(list[[x]]$dxy[list[[x]]$keepMicrobial != F], na.rm = T) + 
      sd(list[[x]]$dxy[list[[x]]$keepMicrobial != F], na.rm = T)*3
    
    plotlist[[x]] <-  ggplot() +
      geom_rect(data = s, alpha = 0.5, 
                aes(xmin = start, xmax = (end - 1), fill = scaffold), ymin = -12, ymax = 12) + 
      scale_fill_manual(breaks = s$scaffold, values = s$color_scaffold, guide = "none") +
      geom_hline(yintercept =  outlier_cutoff, linetype = "dashed") +
      geom_point(data = list[[x]][list[[x]]$keepMicrobial != F & 
                                    list[[x]]$dxy < outlier_cutoff,], 
                 aes(x = midPos_contigContinuous, y = dxy, size = Nsites), 
                 color = "grey55", fill = "grey55", shape = 21, alpha = 0.5) +
      geom_point(data = list[[x]][list[[x]]$keepMicrobial != F & 
                                    list[[x]]$dxy >= outlier_cutoff,], 
                 aes(x = midPos_contigContinuous, y = dxy, size = Nsites), 
                 color = "black", fill = "black", shape = 21, alpha = 0.5)+ 
      # annotate("text", x = Inf, y = Inf,  vjust = 1.1, hjust = 1.1, label = "Unplaced\nscaffolds") +
      scale_size_continuous(range = c(0, 2), guide = "none" ) +
      labs(x = "SNP position", y = "dxy") +
      theme_bw() +theme( panel.grid = element_blank())+
      # ylim(c(-0.05, 1))+
      scale_x_continuous( expand = c(0, 0))
  }
  names(plotlist) <- names(list)
  plotlist
}

get_theta_list <- function(file_list, master) {
  window_size = 5e4
  theta_cols <- c("region",	
                  "chr",                                                                         
                  "midPos",
                  "tW" ,
                  "tP",
                  "tF",
                  "tH",
                  "tL",
                  "Tajima",
                  "fuf",
                  "fud", 
                  "fayh",
                  "zeng",
                  "Nsites")
  
  theta_read <- function(x) {
    read_tsv(x, col_names = theta_cols, skip = 1) %>% 
      filter(midPos %% 2.5e4 == 0) %>% 
      mutate(chr = factor(chr)) %>% 
      arrange(chr, midPos) %>% 
      left_join(master) %>% 
      arrange(contig_order, midPos) %>% 
      mutate(midPos_scaffold = if_else(contig_dir == "+", 
                                       midPos, 
                                       contig_length - midPos),
             midPos_contigContinuous = midPos_scaffold + contig_start_continuous,
             uniq = paste(contig_order, midPos+1, sep = "_"),
             nucl_div = tP/Nsites) %>% 
      filter(Nsites >= 0.2* window_size)
  }
  
  list <- lapply(file_list, theta_read)
  names(list) <- str_sub( gsub("05_popGen/", "", file_list), 1,4)
  list
}



plot_pi <- function(list,  s) {
  plotlist <- list()
  for (x in 1:length(list)) {
    outlier_cutoff <- mean(list[[x]]$nucl_div[list[[x]]$keepMicrobial != F], na.rm = T) + 
      sd(list[[x]]$nucl_div[list[[x]]$keepMicrobial != F], na.rm = T)*3
    
    plotlist[[x]] <-  ggplot() +
      geom_rect(data = s, alpha = 0.5, 
                aes(xmin = start, xmax = (end - 1), fill = scaffold), ymin = -12, ymax = 12) + 
      scale_fill_manual(breaks = s$scaffold, values = s$color_scaffold, guide = "none") +
      geom_hline(yintercept =  outlier_cutoff, linetype = "dashed") +
      geom_point(data = list[[x]][list[[x]]$keepMicrobial != F & 
                                    list[[x]]$nucl_div < outlier_cutoff,], 
                 aes(x = midPos_contigContinuous, y = nucl_div, size = Nsites), 
                 color = "grey55",fill = "grey55" ,shape = 21, alpha = 0.5) +
      geom_point(data = list[[x]][list[[x]]$keepMicrobial != F & 
                                    list[[x]]$nucl_div >= outlier_cutoff,], 
                 aes(x = midPos_contigContinuous, y = nucl_div, size = Nsites), 
                 color = "black",fill = "black", shape = 21, alpha = 0.5)+ 
      annotate("text", x = -Inf, y = c(Inf, -Inf),  vjust = 1.1, hjust = 1.1, label = c("higher in CHST", "higher in COSK")) +
      scale_size_continuous(range = c(0, 2), guide = "none" ) +
      labs(x = "SNP position", y = expression(pi)) +
      theme_bw() +theme( panel.grid = element_blank())+
      # ylim(c(-0.05, 1))+
      scale_x_continuous( expand = c(0, 0))
  }
  names(plotlist) <- names(list)
  plotlist
}

plot_tajD <- function(list,  s) {
  plotlist <- list()
  for (x in 1:length(list)) {
    outlier_cutoff <- c(mean(list[[x]]$Tajima[list[[x]]$keepMicrobial != F], na.rm = T) -
                          sd(list[[x]]$Tajima[list[[x]]$keepMicrobial != F], na.rm = T)*3,
                        mean(list[[x]]$Tajima[list[[x]]$keepMicrobial != F], na.rm = T) + 
                          sd(list[[x]]$Tajima[list[[x]]$keepMicrobial != F], na.rm = T)*3)
    
    plotlist[[x]] <-  ggplot() +
      geom_rect(data = s, alpha = 0.5, 
                aes(xmin = start, xmax = (end - 1), fill = scaffold), ymin = -12, ymax = 12) + 
      scale_fill_manual(breaks = s$scaffold, values = s$color_scaffold, guide = "none") +
      geom_hline(yintercept =  outlier_cutoff, linetype = "dashed") +
      geom_point(data = list[[x]][list[[x]]$keepMicrobial != F & 
                                    (list[[x]]$Tajima > outlier_cutoff[1] | 
                                       list[[x]]$Tajima < outlier_cutoff[2]) ,], 
                 aes(x = midPos_contigContinuous, y = Tajima, size = Nsites), 
                 color = "grey55", shape = 21, alpha = 0.5) +
      geom_point(data = list[[x]][list[[x]]$keepMicrobial != F & 
                                    (list[[x]]$Tajima <= outlier_cutoff[1] | 
                                       list[[x]]$Tajima >= outlier_cutoff[2]),], 
                 aes(x = midPos_contigContinuous, y = Tajima, size = Nsites), 
                 color = "black", shape = 21, alpha = 0.5)+ 
      # annotate("text", x = -Inf, y = c(Inf, -Inf),  vjust = c(1.1, -0.1), hjust =-0.1, label = c("higher in CHST", "higher in COSK")) +
      scale_size_continuous(range = c(0, 2), guide = "none" ) +
      labs(x = "SNP position", y = "Tajima's D") +
      theme_bw() +theme( panel.grid = element_blank())+
      scale_x_continuous( expand = c(0, 0))
  }
  names(plotlist) <- names(list)
  plotlist
}

comb_theta <- function(list) {
  poplist <- c("CHST", "COSK")
  
  for (i in 1:length(list)) {
    list[[i]]$pop <- names(list)[i]
  }
  purrr::reduce(list, full_join) %>% 
    mutate(transect = factor( if_else(pop %in% poplist[1:4], "West", "East"), levels = c("West", "East")),
           host = factor( str_sub(pop, 1,2), levels = c("CH", "CO")),
           pop = factor(pop, levels = poplist))
}

stat_sum_df <- function(fun, geom="errorbar", ...) {
  stat_summary(fun.data = fun, geom = geom, width = 0.2, ...)
}

plot_theta_dist <- function(df, pal, var) {
  if(var == "nucl_div") {
    ylab = expression("Nucleotide diversity"~(pi))
  } else if (var == "Tajima") {
    ylab = "Tajima's D"
  }
  
  ggplot(df, aes(y = .data[[var]], x = pop, fill = pop, group = pop)) +
    geom_violin(aes(color = pop), draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5, width = 0.75 ) + 
    stat_sum_df("mean_cl_boot", mapping = aes(group = pop, shape = transect, fill = pop), alpha = 0.75) +
    # stat_summary(geom = "point", fun = "mean", aes(shape = transect, fill = pop), alpha = 0.75, size = 2) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_color_manual(values = pal, guide = "none") +
    scale_shape_manual(values = c(21, 23), guide = "none") +
    facet_grid(~transect, scales = "free", space = "free") +
    labs( x = "Population", y = ylab) +
    theme_bw() +
    # coord_cartesian(ylim = c(-0.03, 0.1)) +
    theme( panel.grid = element_blank(), 
           strip.background = element_blank(),
           axis.text.x = element_text(angle = 35, hjust = 0.9))
}

get_diff <- function(combined_data, poplist){
  nucl_div_vars <- paste("nucl_div", c(poplist[1],poplist[2]), sep = "_")
  Tajima_vars <- paste("Tajima",  c(poplist[1],poplist[2]), sep = "_")
  Nsites_vars <- paste("Nsites",  c(poplist[1],poplist[2]), sep = "_")
  
  combined_data %>% na.omit() %>% 
    dplyr::filter(pop == poplist[1] | pop == poplist[2]) %>% 
    dplyr::select(pop, chr, midPos_contigContinuous, keepMicrobial, Nsites, nucl_div, Tajima, uniq) %>% 
    pivot_wider(id_cols = c(chr, midPos_contigContinuous, keepMicrobial,uniq), names_from = pop, values_from = c(Nsites, nucl_div, Tajima)) %>% 
    mutate(nucl_div_diff = .data[[nucl_div_vars[1]]] - .data[[nucl_div_vars[2]]],
           Tajima_diff = .data[[Tajima_vars[1]]] - .data[[Tajima_vars[2]]],
           Nsites_avg = (.data[[Nsites_vars[1]]]+ .data[[Nsites_vars[2]]])/2,
           comp = paste(poplist, collapse = ".")) %>% 
    dplyr::select(comp, chr, midPos_contigContinuous, keepMicrobial, Nsites_avg, nucl_div_diff, Tajima_diff,uniq) 
  
}


plot_diff <- function(data, s, pal, var) {
  
  if(var == "nucl_div_diff") {
    ylab = expression(Delta~pi)
  } else if (var == "Tajima_diff") {
    ylab = expression(Delta~"D")
  }
  
  outlier_cutoff <- c(mean(data[[var]][data$keepMicrobial != F], na.rm = T) -
                        sd(data[[var]][data$keepMicrobial != F], na.rm = T)*3,
                      mean(data[[var]][data$keepMicrobial != F], na.rm = T) + 
                        sd(data[[var]][data$keepMicrobial != F], na.rm = T)*3)
  
  plotname <- unique(data$comp)
  plot <- ggplot() +
    geom_rect(data = s, alpha = 0.5, 
              aes(xmin = start, xmax = (end - 1), fill = scaffold), ymin = -12, ymax = 12) + 
    scale_fill_manual(breaks = s$scaffold, values = s$color_scaffold, guide = "none") +
    geom_hline(yintercept =  outlier_cutoff, linetype = "dashed") +
    geom_point(data = data[data$keepMicrobial != F & 
                             (data[[var]]> outlier_cutoff[1] | 
                                data[[var]] < outlier_cutoff[2]) ,], 
               aes(x = midPos_contigContinuous, y = .data[[var]], size = Nsites_avg), 
               color = "grey55", shape = 1, alpha = 0.5) +
    geom_point(data = data[data$keepMicrobial != F & 
                             (data[[var]] >= outlier_cutoff[2]),], 
               aes(x = midPos_contigContinuous, y = .data[[var]], size = Nsites_avg), 
               color = pal[1],fill = pal[1], shape = 22, alpha = 0.5) + 
    geom_point(data = data[data$keepMicrobial != F & 
                             (data[[var]] <= outlier_cutoff[1]),], 
               aes(x = midPos_contigContinuous, y = .data[[var]], size = Nsites_avg), 
               color = pal[2],fill = pal[2], shape = 21, alpha = 0.5) + 
    # annotate("text", x = -Inf, y = Inf,  vjust = 1.1, hjust = -0.1, label = plotname) +
    # annotate("text", x = Inf, y = Inf,  vjust = 1.1, hjust = 1.1, label = "Unplaced\nscaffolds") +
    scale_size_continuous(range = c(0, 2), guide = "none" ) +
    labs(x = "SNP position", y = ylab) +
    theme_bw() +
    theme( panel.grid = element_blank(), 
           # axis.title = element_blank(), 
           plot.margin = margin(0,0,0,0, unit = "pt"))+
    scale_x_continuous( expand = c(0, 0))
  
  plot
}


#### Outlier windows ####
get_outliers_1way <- function(data, var) {
  outlier_cutoff <- c(mean(data[[var]][data$keepMicrobial != F], na.rm = T) + 
                        sd(data[[var]][data$keepMicrobial != F], na.rm = T)*3)
  outlier_list <- data %>% 
    dplyr::filter(.data[[var]] > outlier_cutoff) %>% 
    dplyr::select(uniq)
  
  outlier_list$uniq
}

get_outliers_2way <- function(data, var) {
  outlier_cutoff <- c(mean(data[[var]][data$keepMicrobial != F], na.rm = T) -
                        sd(data[[var]][data$keepMicrobial != F], na.rm = T)*3,
                      mean(data[[var]][data$keepMicrobial != F], na.rm = T) + 
                        sd(data[[var]][data$keepMicrobial != F], na.rm = T)*3)
  outlier_list <- data %>% 
    dplyr::filter(.data[[var]] < outlier_cutoff[1] |.data[[var]] > outlier_cutoff[2]) %>% 
    dplyr::select(uniq)
  outlier_list$uniq
}

upset_outliers <- function(comb_mat){
  outlier_colors <- c("grey55", "grey55", "black", "black")[comb_degree(comb_mat)]
  
  UpSet(comb_mat, top_annotation =
          HeatmapAnnotation(
            "Intersection size" = anno_barplot(
              comb_size(comb_mat), 
              gp = gpar(fill = outlier_colors, 
                        col = outlier_colors), 
              axis_param = list(side = "left"),
              add_numbers = TRUE,
              border = FALSE, 
              numbers_rot = 0, 
              height = unit(4, "cm")),
            annotation_name_rot = 90, 
            annotation_name_side = "left"),
        comb_col = outlier_colors,
        right_annotation = 
          upset_right_annotation(comb_mat, add_numbers = TRUE,  
                                 axis_param = list(
                                   labels_rot = 0)))
}

get_bed <- function(bedPath){
  bed_cols <- c("chr", "start", "end", "slopStart", "slopEnd", "geneID")
  bed <- read_tsv(bedPath, col_names = bed_cols) %>% 
    mutate(midPos = start + 2.5e4 + 1)
  bed
}


get_outlier_genes <- function(outliers, bed, master) {
  bed %>% 
    dplyr::select(chr, midPos, geneID) %>% 
    arrange(chr, midPos) %>% 
    left_join(master) %>% unique() %>% 
    arrange(contig_order, midPos) %>% 
    mutate(midPos_scaffold = if_else(contig_dir == "+", midPos, contig_length - midPos),
           midPos_contigContinuous = midPos_scaffold + contig_start_continuous,
           uniq = paste(contig_order, midPos, sep = "_")) %>% 
    dplyr::filter(uniq %in% outliers)
}

get_outlier_genes_comb <- function(comb_mat, bed, master) {
  comb_names <- comb_name(comb_mat)
  
  outlier_genes_list <- list()
  for (i in 1:length(comb_names)){
    outliers <- extract_comb(comb_mat, comb_names[i])
    outlier_genes_list[[i]] <- bed %>% 
      dplyr::select(chr, midPos, geneID) %>% 
      arrange(chr, midPos) %>% 
      left_join(master) %>% unique() %>% 
      arrange(contig_order, midPos) %>% 
      mutate(midPos_scaffold = if_else(contig_dir == "+", midPos, contig_length - midPos),
             midPos_contigContinuous = midPos_scaffold + contig_start_continuous,
             uniq = paste(contig_order, midPos, sep = "_")) %>% 
      dplyr::filter(uniq %in% outliers)
    names(outlier_genes_list)[i] <- comb_names[i]
  } 
  outlier_genes_list
}

topgo_outliers <- function(outliers, GO_path){
  outlier_df <- outliers %>% 
    purrr::reduce(full_join) %>% 
    dplyr::select(geneID) %>% 
    unique()
  
  # analyzing GO terms with at least 5 members,
  # as this yield more stable results.
  node_size= 5
  
  # use topGO to read the functional annotation
  geneID2GO <- readMappings(GO_path)
  
  # define the gene universe
  geneUniverse <- names(geneID2GO)
  
  genesOfInterest.bv <- outlier_df$geneID
  
  geneList.bv <- factor(as.integer(geneUniverse %in% genesOfInterest.bv))
  names(geneList.bv) <- geneUniverse
  
  categories <- c("BP", "MF")
  GO_list <- list()
  
  for (i in 1:2) {
    tryCatch({
      GO_category=categories[i]
      # Create a new topGO GSEA objected
      myGOdata.bv <- new("topGOdata",
                         description="Candidate genes", 
                         ontology=GO_category, 
                         allGenes=geneList.bv, 
                         annot = annFUN.gene2GO, 
                         gene2GO = geneID2GO, 
                         nodeSize = node_size)  
      resultClassic <- runTest(myGOdata.bv, algorithm="classic", statistic="fisher")
      
      # Fisher with parent-child algorithm
      resultParentchild <- runTest(myGOdata.bv, algorithm="parentchild", statistic="fisher")
      
      # Create a GSEA results table
      # GenTable() is  annoying and requires a number of terms we want reported (default is 10)
      ## so we leverage our topGOresult object
      top_nodes <-  sum(resultParentchild@score < 0.05)
      
      GO_list[[i]] <-  GenTable(
        myGOdata.bv,
        classicFisher = resultClassic,
        parentchildFisher = resultParentchild,
        orderBy = "parentchildFisher",
        topNodes = top_nodes) %>% 
        mutate(cat = GO_category)
    }, error=function(e){})
  }
  GO_list
}


plot_topgo_outliers <- function(GO_list){
  plot_list <- list()
  for (i in 1:2) {
    tryCatch({
      go_analysis <- GO_list[[i]] %>% 
        mutate(GO_term = paste(GO.ID, Term, sep = "-"),
               parentchildFisher = as.numeric(parentchildFisher), 
               sig= if_else(parentchildFisher < 0.01, "p < 0.01", 
                            if_else(parentchildFisher <= 0.05, "p < 0.05", "Not sig.")))
      
      plot_title <- paste(go_analysis[1]$comp, go_analysis[1]$dir, go_analysis[1]$cat)
      
      if(dim(go_analysis)[1] < 30){
        max_rank <- dim(go_analysis)[1] } else {
          max_rank <- 30
        }
      
      plot_list[[i]] <- ggplot(
        go_analysis[1:max_rank,], 
        aes( x= -log10(parentchildFisher), 
             y = forcats::fct_reorder(GO_term, -log10(parentchildFisher)),
             color = sig)) +
        geom_point(aes(size = Annotated)) +
        scale_color_manual(values = c("black", "grey55", "grey"), 
                           breaks = c("p < 0.01", "p < 0.05", "Not sig."),
                           guide = "none") +
        # geom_vline(xintercept = c(-log10(0.05), -log10(0.01)),
        #            linetype = c("dotted", "dashed")) +
        labs(x = expression(-log[10](adj.~p-value)), 
             y = "GO term and function", title = plot_title ) +
        theme_classic() +
        theme(legend.position = c(1,0), 
              legend.justification = c("right", "bottom"))
    }, error=function(e){})
  }
  plot_list
}

#### DEG vs. outlier comparisons ####

locate_deg <- function(deseq_results_list, comp, bed, master) {
    trt_names <- str_split(comp, pattern = "_",simplify = T)[c(1,3)]
    
    down_trt <- paste0("Up in ", trt_names[1])
    up_trt <- paste0("Up in ", trt_names[2])
    
    down_genes <- data.frame(geneID =deseq_results_list$down[names(deseq_results_list$down) == comp][[1]]$geneID) %>% 
      mutate(DE = down_trt, 
             comp = factor(comp, levels = names(deseq_results_list$down)))
    up_genes <- data.frame(geneID = deseq_results_list$up[names(deseq_results_list$up) == comp][[1]]$geneID) %>% 
      mutate(DE = up_trt, 
             comp = factor(comp, levels = names(deseq_results_list$down)))
    rbind(down_genes, up_genes) %>% 
      left_join(bed) %>%
      left_join(master) %>% 
      mutate(slop_midPos_continuous = (slopEnd - slopStart) + contig_start_continuous,
             contig_midPos_continuous = contig_start_continuous + midPos) %>% 
      dplyr::select(comp, geneID, DE, chr, slopStart, slopEnd, scaffold, slop_midPos_continuous, contig_midPos_continuous, X_linked_enriched, keepMicrobial) %>% 
      unique()
}


locate_modgenes <- function(module_genes, bed, master) {
  module_genes %>%
    ungroup %>% 
    dplyr::select(geneID) %>% 
    unique() %>% 
    left_join(bed) %>%
    left_join(master) %>% 
    mutate(slop_midPos_continuous = (slopEnd - slopStart) + contig_start_continuous,
           contig_midPos_continuous = contig_start_continuous + midPos) %>% 
    dplyr::select(geneID, chr, slopStart, slopEnd, scaffold, slop_midPos_continuous, contig_midPos_continuous, X_linked_enriched, keepMicrobial) %>% 
    unique()
}


deg_rug <- function(deg_locations, comp, pal) {
  
  names(pal) <- c("H", "HH", "HO", "O", "OO", "OH")
  
  trt_names <- str_split(comp, pattern = "_",simplify = T)[c(1,3)]
  
  down <- unname(pal[names(pal) %in% trt_names[1]])
  up <- unname(pal[names(pal) %in% trt_names[2]])
  
  down_trt <- paste0("Up in ", trt_names[1])
  up_trt <- paste0("Up in ", trt_names[2])
  
  down_genes <- data.frame(geneID =deseq_results_list$down[names(deseq_results_list$down) == comp][[1]]$geneID) %>% 
    mutate(DE = down_trt) 
  up_genes <- data.frame(geneID = deseq_results_list$up[names(deseq_results_list$up) == comp][[1]]$geneID) %>% 
    mutate(DE = up_trt) 
}

fit_euler_DE_all <- function(up, down, outliers){
  comps <- names(up)
  
  if("list" %in% class(outliers)){
  outlier_df <- outliers %>% 
    purrr::reduce(full_join) %>% 
    dplyr::select(geneID) %>% 
    unique()} else {
      outlier_df <- outliers %>% 
        dplyr::select(geneID) %>% 
        unique()
    }
  
  fit_list <- list()
  # create for loop to produce ggplot2 graphs 
  for (i in 1:length(up)) { 
    if(length(down[[i]]$geneID) >= 1 & 
       length(up[[i]]$geneID) >= 1) {
      fit_list[[i]] <- eulerr::euler( list(up = up[[i]]$geneID, down = down[[i]]$geneID, outliers = outlier_df$geneID), shape = "ellipse")
    } else if (length(down[[i]]$geneID) < 1 & 
               length(up[[i]]$geneID) >= 1) {
      fit_list[[i]] <- eulerr::euler( list(up = up[[i]]$geneID, outliers = outlier_df$geneID), shape = "ellipse")
    } else if (length(down[[i]]$geneID) >= 1 & 
               length(up[[i]]$geneID) < 1) {
      fit_list[[i]] <- eulerr::euler( list(down = down[[i]]$geneID, outliers = outlier_df$geneID), shape = "ellipse")
    }
    names(fit_list)[[i]] <- comps[i]
  }
  
  fit_list
}

fit_euler_DE_expr <- function(up, down, outliers, expr_genes){
  comps <- names(up)
  
  if("list" %in% class(outliers)){
    outlier_df <- outliers %>% 
      purrr::reduce(full_join) %>% 
      dplyr::select(geneID) %>% 
      unique()} else {
        outlier_df <- outliers %>% 
          dplyr::select(geneID) %>% 
          unique()
      }
  
  outlier_df <- outlier_df[outlier_df$geneID %in% expr_genes,]
  
  fit_list <- list()
  # create for loop to produce ggplot2 graphs 
  for (i in 1:length(up)) { 
    if(length(down[[i]]$geneID) >= 1 & 
       length(up[[i]]$geneID) >= 1) {
      fit_list[[i]] <- eulerr::euler( list(up = up[[i]]$geneID, down = down[[i]]$geneID, outliers = outlier_df$geneID), shape = "ellipse")
    } else if (length(down[[i]]$geneID) < 1 & 
               length(up[[i]]$geneID) >= 1) {
      fit_list[[i]] <- eulerr::euler( list(up = up[[i]]$geneID, outliers = outlier_df$geneID), shape = "ellipse")
    } else if (length(down[[i]]$geneID) >= 1 & 
               length(up[[i]]$geneID) < 1) {
      fit_list[[i]] <- eulerr::euler( list(down = down[[i]]$geneID, outliers = outlier_df$geneID), shape = "ellipse")
    }
    names(fit_list)[[i]] <- comps[i]
  }
  
  fit_list
}

plot_euler_DE <- function(fit_list, pal){
  names(pal) <- c("H", "HH", "HO", "O", "OO", "OH")
  comps <- names(fit_list)
  
  plot_list <- list()
  # create for loop to produce ggplot2 graphs 
  for (i in 1:length(comps)) { 
    trt_names <- str_split(comps[i], pattern = "_",simplify = T)[c(1,3)]
    
    down_pal <- scales::alpha(unname(pal[names(pal) %in% trt_names[2]]), alpha = 0.75)
    up_pal <- scales::alpha(unname(pal[names(pal) %in% trt_names[1]]), alpha = 0.75)
    outlier_pal <- scales::alpha("black", alpha = 0.5)
    
    down_trt <- paste0("Up in ", trt_names[2])
    up_trt <- paste0("Up in ", trt_names[1])
    outlier_trt <- "Outliers"
    
    if(sum(names(fit_list[[i]]$original.values)[1:2] == c("up", "down")) == 2) {
      plot_list[[i]] <- plot(fit_list[[i]],
                             quantities = TRUE,
                             fill = c(up_pal, down_pal, outlier_pal),
                             legend = list(labels = c(up_trt, down_trt, outlier_trt)))
    } else if (sum(names(fit_list[[i]]$original.values)[1:2] == c("up", "outliers")) == 2) {
      plot_list[[i]] <- plot(fit_list[[i]],
                             quantities = TRUE,
                             fill = c(up_pal, outlier_pal),
                             legend = list(labels = c(up_trt, outlier_trt)))
    } else if (sum(names(fit_list[[i]]$original.values)[1:2] == c("down", "outliers")) == 2) {
      fit <- eulerr::euler( list(down = down[[i]]$geneID, outliers = outlier_df$geneID), shape = "ellipse")
      plot_list[[i]] <- plot(fit_list[[i]],
                             quantities = TRUE,
                             fill = c(down_pal, outlier_pal),
                             legend = list(labels = c(down_trt, outlier_trt)))
    }
    # names(plot_list)[i] <- comps[i]
  }
  panel <- ggpubr::ggarrange(plotlist = plot_list, ncol = 3, nrow = 3, labels = comps)
  ggsave(plot = panel, filename =  "04_targetPlots/euler_DE_outliers.pdf", height = 10, width = 10)
  panel
}



fit_euler_module_all <- function(module_genes_list, outliers){
  
  if("list" %in% class(outliers)){
    outlier_df <- outliers %>% 
      purrr::reduce(full_join) %>% 
      dplyr::select(geneID) %>% 
      unique()} else {
        outlier_df <- outliers %>% 
          dplyr::select(geneID) %>% 
          unique()
      }
  
  fit_list <- list()
  for (i in 1:length(module_genes_list)){
    mod_genes <- module_genes_list[[i]]$geneID %>% unique()
    fit_list[[i]] <- eulerr::euler(list(mod = mod_genes,  outliers = outlier_df$geneID), shape = "ellipse")
  }
  names(fit_list) <- names(module_genes_list)
  fit_list
}


fit_euler_module_expr <- function(module_genes_list, outliers, expr_genes){
  
  if("list" %in% class(outliers)){
    outlier_df <- outliers %>% 
      purrr::reduce(full_join) %>% 
      dplyr::select(geneID) %>% 
      unique()} else {
        outlier_df <- outliers %>% 
          dplyr::select(geneID) %>% 
          unique()
      }
  
  outlier_df <- outlier_df[outlier_df$geneID %in% expr_genes,]
  
  fit_list <- list()
  for (i in 1:length(module_genes_list)){
    mod_genes <- module_genes_list[[i]]$geneID %>% unique()
    fit_list[[i]] <- eulerr::euler(list(mod = mod_genes,  outliers = outlier_df$geneID), shape = "ellipse")
  }
  names(fit_list) <- names(module_genes_list)
  fit_list
}


hypergeom_test_DE <- function(fit_list, n_genes) {
  comps <- names(fit_list)
  fit_list_wide <- list()
  for (i in 1:length(comps)){
    fit_list_wide[[i]] <- data.frame(values = fit_list[[i]]$original.values) %>% 
      rownames_to_column("overlap") %>% 
      mutate(comp = comps[i],
             overlap = gsub("&","_", overlap)) %>% 
      pivot_wider(id_cols = comp, names_from = overlap, values_from = values)
  }
  df <- fit_list_wide %>% 
    purrr::reduce(full_join) %>% 
    group_by(comp) %>% 
    replace(is.na(.), 0) %>% 
    mutate(A_DE_outliers = sum(up_outliers, down_outliers),
           B_DE_notOutliers = sum(up,down), 
           C_notDE_outliers = outliers, 
           D_notDE_notOutliers = n_genes - sum(A_DE_outliers, B_DE_notOutliers, C_notDE_outliers, na.rm = T),
           OR = (A_DE_outliers*D_notDE_notOutliers)/ (B_DE_notOutliers*C_notDE_outliers),
           q_DE_outliers = A_DE_outliers,
           m_tot_outliers = A_DE_outliers + C_notDE_outliers,
           n_not_outliers = n_genes - m_tot_outliers,
           k_tot_DE = A_DE_outliers + B_DE_notOutliers, 
           phyper = phyper(q_DE_outliers, m_tot_outliers, n_not_outliers,k_tot_DE, lower.tail = F)[[1]])
  
  df
}


plot_euler_module <- function(fit_list, names){
  pal <- scales::alpha(c("white","black"), alpha = 0.5)
  
  fit_plot_list <- list()
  for (i in 1:length(fit_list)){
  plot(fit,
       quantities = TRUE,
       fill =outlier_pal,
       legend = list(labels = c(name, "Outliers")))
  }
  fit_plot_list
}
  
hypergeom_test_module <- function(fit_list, n_genes) {
  comps <- names(fit_list)
  fit_list_wide <- list()
  for (i in 1:length(comps)){
    fit_list_wide[[i]] <- data.frame(values = fit_list[[i]]$original.values) %>% 
      rownames_to_column("overlap") %>% 
      mutate(comp = comps[i],
             overlap = gsub("&","_", overlap)) %>% 
      pivot_wider(id_cols = comp, names_from = overlap, values_from = values)
  }
  df <- fit_list_wide %>% 
    purrr::reduce(full_join) %>% 
    group_by(comp) %>% 
    replace(is.na(.), 0) %>% 
        mutate(A_mod_outliers = mod_outliers,
           B_mod_notOutliers = mod, 
           C_notmod_outliers = outliers, 
           D_notmod_notOutliers = n_genes - sum(A_mod_outliers, B_mod_notOutliers, C_notmod_outliers, na.rm = T),
           OR = (A_mod_outliers*D_notmod_notOutliers)/ (B_mod_notOutliers*C_notmod_outliers),
           q_mod_outliers = A_mod_outliers,
           m_tot_outliers = A_mod_outliers + C_notmod_outliers,
           n_not_outliers = n_genes - m_tot_outliers,
           k_tot_mod = A_mod_outliers + B_mod_notOutliers, 
           phyper = phyper(q_mod_outliers, m_tot_outliers, n_not_outliers,k_tot_mod, lower.tail = F)[[1]],
           phyper_adj = p.adjust(phyper, n = length(comps), method = "BH"))
  
  df
}


add_popgen_to_bed <- function(bed, contig_master, fst_df, dxy_df, theta_diff_df, theta_list) {
  theta_df_ch <- theta_list[[1]] %>% 
    mutate(midPos = midPos + 1) %>% 
    dplyr::select(chr, midPos, nucl_div_ch = nucl_div, tajima_ch = Tajima)
  theta_df_co <- theta_list[[2]] %>% 
    mutate(midPos = midPos + 1) %>% 
    dplyr::select(chr, midPos, nucl_div_co = nucl_div, tajima_co = Tajima)
  
  bed %>% full_join(contig_master %>% dplyr::select(chr, contig_order)) %>% 
    full_join(dplyr::select(fst_df, chr, midPos, Fst, uniq)) %>% 
    full_join(dplyr::select(dxy_df, chr, midPos, dxy, uniq)) %>% 
    full_join(dplyr::select(theta_diff_df, chr, nucl_div_diff, Tajima_diff, uniq)) %>% 
    full_join(theta_df_ch) %>% 
    full_join(theta_df_co)
}

add_popgen_to_modules <- function(moduleMembership, modules, bed_popgen_df){
  module_coefs <- data.frame(cluster = rownames(modules$moduleHostCor), 
                             moduleHostCor = modules$moduleHostCor,
                             moduleStressCor = modules$moduleStressCor, 
                             modulePlasticityCor = modules$modulePlasticityCor,
                             moduleCHplasticityCor =  modules$moduleCHplasticityCor,
                             moduleCOplasticityCor = modules$moduleCOplasticityCor)
  
  df <- bed_popgen_df %>% 
    left_join(moduleMembership %>% 
                mutate(n_genes = n()) %>% 
                ungroup() %>% 
                dplyr::select(geneID = gene, clusterCol, cluster, n_genes) %>% 
                unique())  %>% 
    left_join(module_coefs) %>% 
    mutate(cluster = factor(if_else(is.na(cluster), "Not expr.", cluster)),
          clusterCol = if_else(is.na(clusterCol), "white", clusterCol))
  
  df
}

plot_popgen_modules <- function(df, pal){
  metrics <- c("Fst", "dxy", "nucl_div_diff", "Tajima_diff", "nucl_div_ch", "nucl_div_co", "tajima_ch", "tajima_co")
  col_fun_eco = circlize::colorRamp2(c(-0.75, 0, 0.75), c(pal[4],"white", pal[1]))
  col_fun_stress = circlize::colorRamp2(c(-0.75, 0, 0.75), c("#f05039","white", "#f05039"))
  col_fun_CH = circlize::colorRamp2(c(-0.75, 0, 0.75), c(pal[3],"white", pal[3]))

    df$cluster <- fct_reorder(df$cluster, df$moduleHostCor)
  # df$cluster <- factor(df$cluster, levels = c("none", levels(df$cluster)))
    df$moduleHostCor[df$cluster == "Not expr."] <- 0
    df <- df %>% mutate(plot_col = case_when(cluster == "mod16" ~ col_fun_CH(moduleCHplasticityCor), 
                                             abs(moduleStressCor) > 0.5 ~ col_fun_stress(moduleStressCor), 
                                             abs(moduleHostCor) > 0.5 ~ col_fun_eco(moduleHostCor),
                                             .default = "white"))
    
    
    col_sum <- df %>% 
      dplyr::select(cluster, moduleHostCor, plot_col) %>% 
      unique() %>% 
      arrange(cluster) %>% 
      mutate(plot_col, as.factor(plot_col))
      
    
  module_plot_list <- list()
  dunn_list <- list()
  for (i in 1:8){
    
    temp <- df %>% 
      dplyr::select(geneID, cluster,moduleHostCor, any_of(metrics[i])) %>% 
      na.omit() %>% 
      unique() 
    temp_kruskal <-  temp %>% rstatix::kruskal_test(as.formula(paste(metrics[i], '~ cluster')))
    dunn_list[[i]] <- temp %>% rstatix::dunn_test(as.formula(paste(metrics[i], '~ cluster'))) %>% arrange(p.adj)
    
    module_plot_list[[i]] <- df %>% 
      dplyr::select(geneID, cluster, moduleHostCor, plot_col, any_of(metrics[i])) %>% 
      na.omit() %>% 
      unique() %>% 
      # filter(cluster != "mod20") %>% 
      ggplot(aes(x = cluster, y = .data[[metrics[i]]])) +
      # geom_jitter(aes(color = cluster), alpha = 0.5, size = 1, width = 0.25) +
      # stat_summary(aes(fill = cluster), fun.data = "mean_cl_boot", colour = "black", linewidth = 1, size = 0.5, shape = 21) +
       geom_boxplot(aes(fill = plot_col), outlier.shape = 21, alpha = 0.75) +
      labs(x = "Modules", y = metrics[i]) + 
      # scale_fill_gradient2(low = pal[4],mid = "white", high = pal[1], midpoint = 0, guide = "none", na.value = "white")+
      scale_fill_manual(values = col_sum$plot_col, breaks = col_sum$plot_col,  guide = "none") +
      # scale_y_continuous(breaks = seq(from = -0.5, to = 0.7, by = 0.1)) +
      theme_bw() +
      theme(
        panel.grid = element_blank() , 
        axis.text.x = element_text(angle = 45, hjust = 0.9, size = 7),
        plot.title = element_text(face = "bold", angle = 0, vjust = 1, size = 10)) +
      annotate( geom = "text", label = paste0("Kruskal-Wallis test, p=", temp_kruskal$p), x = 0, y = Inf, vjust = 1.1, hjust = -0.1)
 
  }
  dunn_df <- dunn_list %>% purrr::reduce(full_join)
  write_tsv(dunn_df, file = paste0("00_data/03_popgen_expression/popgen_modules_dunns_test_results.tsv"))
  module_plot_list
}

plot_popgen_DE <- function(bed_popgen_df, deseq_results_list, pal, comp, sig_test_df){
  
  trt_names <- str_split(comp, pattern = "_",simplify = T)[c(1,3)]
  
  names(pal) <- c("H", "HH", "HO", "O", "OO", "OH")
  pal <- c("grey40", pal[names(pal) %in% trt_names])
  
  down_genes <- deseq_results_list$down[names(deseq_results_list$down) == comp][[1]]
  up_genes <- deseq_results_list$up[names(deseq_results_list$up) == comp][[1]]
  
  metrics <- c("Fst", "dxy", "nucl_div_diff", "Tajima_diff",
               "nucl_div_ch", "nucl_div_co", "tajima_ch", "tajima_co")
  
  DE_popgen_plotlist <- list()
  for (i in 1:8){
    df_sub <- bed_popgen_df %>% 
      left_join(deseq_results_list$shrink[names(deseq_results_list$shrink) == comp][[1]] %>% 
                  as.data.frame() %>% 
                  rownames_to_column("geneID")) %>% 
      # mutate(DE = if_else(geneID %in% down_genes$geneID, paste0("Up in ", trt_names[1]), 
      #                     if_else(geneID %in% up_genes$geneID, paste0("Up in ", trt_names[2]),
      #                             "not DE")),
      mutate(DE = factor(case_when(geneID %in% c(down_genes$geneID, up_genes$geneID) ~ "DE",
                          is.na(baseMean) ~ "Not expr.", 
                          .default = "Not DE"), 
                         levels = c("DE", "Not DE", "Not expr.")),
             
             Inversion = if_else(contig_order > 611 & contig_order < 699, "Inside", "Outside")) %>% 
      dplyr::select(geneID, Inversion, DE, metric = any_of(metrics[i])) %>% 
      filter(DE != "Not expr.") %>% 
      na.omit() %>% 
      unique() 
    
    sig_test_df_sub <- sig_test_df %>% 
      mutate(
        p.signif = if_else(padj >= 0.05, "KW, p > 0.05", 
                           paste("KW, p =", format(padj, scientific = TRUE, digits = 2)))) %>% 
      dplyr::select(Inversion, metric = .y., p.signif, dataset) %>% 
      dplyr::filter(metric == metrics[i] & dataset == "all") 
    
    DE_popgen_plotlist[[i]] <- df_sub %>% 
      ggplot(aes(x = DE, y = metric) )+
      facet_grid(~Inversion) +
      ggbeeswarm::geom_quasirandom(aes(fill = DE), color = "grey50", width = 0.2, shape = 21) +
      geom_boxplot(fill = NA, width = 0.5, outlier.shape = NA) +
      geom_text(data = sig_test_df_sub,aes( label =  p.signif), x = 0.5, y = Inf, hjust = 0, vjust = 1.4)   +
      # scale_fill_manual(values = unname(pal)) +
      # scale_color_manual(values = unname(pal)) +
      scale_fill_brewer(palette = "Greys", direction = -1) +
      labs(y = metrics[i]) +
      theme_bw() +
      theme(
        panel.grid = element_blank(), 
        strip.background = element_blank(), 
        legend.position = "none")
  }
  DE_popgen_plotlist
}

popgen_DE_kw_matched <- function(bed_popgen_df, deseq_results_list, comp){
  
  down_genes <- deseq_results_list$down[names(deseq_results_list$down) == comp][[1]]
  up_genes <- deseq_results_list$up[names(deseq_results_list$up) == comp][[1]]
  
  df_mod <- bed_popgen_df %>% 
    dplyr::left_join(
      deseq_results_list$shrink[names(deseq_results_list$shrink) == comp][[1]] %>% 
        as.data.frame() %>% 
        rownames_to_column("geneID")) %>% 
    # mutate(DE = if_else(geneID %in% down_genes$geneID, paste0("Up in ", trt_names[1]), 
    #                     if_else(geneID %in% up_genes$geneID, paste0("Up in ", trt_names[2]),
    #                             "not DE")),
    dplyr::mutate(
      DE = factor(case_when(geneID %in% c(down_genes$geneID, up_genes$geneID) ~ "DE",
                            is.na(baseMean) ~ "Not expr.", 
                            .default = "Not DE"), 
                  levels = c("DE", "Not DE", "Not expr.")),
      Inversion = if_else(contig_order > 611 & contig_order < 699, "Inside", "Outside"),
      gene_length_slop = slopEnd - slopStart, 
      contig = chr) %>% 
    dplyr::rename(
      windowStart = start, 
      windowEnd = end) %>% 
    filter(DE != "Not expr.") 
  
  kw_test <- data.frame()
  
  metrics <- c("Fst", "dxy", "nucl_div_diff", "Tajima_diff","nucl_div_ch", "nucl_div_co", "tajima_ch", "tajima_co")
  for (i in 1:8){
    temp1 <- df_mod %>% 
      dplyr::select(
        chr, slopStart, slopEnd, geneID, gene_length_slop, Inversion,
        contig, DE, metric = any_of(metrics[i])) %>% 
      na.omit() %>% 
      unique() 
    
    kw <-  temp1 %>%
      group_by(Inversion) %>% 
      rstatix::kruskal_test(metric ~ DE) %>%
      dplyr::select(-method) %>% 
      mutate(.y. = metrics[i], 
             dataset = "all") 
    
    es <- temp1 %>% 
      group_by(Inversion) %>% 
      rstatix::kruskal_effsize(metric ~ DE) %>%
      dplyr::select(-method) %>% 
      mutate(.y. = metrics[i],
             dataset = "all") 
    
    kw_test <-  rbind(kw_test,
                      kw %>% 
                        left_join(es, by = join_by(Inversion, .y., n, dataset))) 
    
    df_granges <- GenomicRanges::makeGRangesFromDataFrame(
      temp1, start.field = "slopStart", end.field = "slopEnd", 
      na.rm = T, keep.extra.columns=T)
    # df_granges@ranges@NAMES<- df_granges$geneID
    # df_granges@elementMetadata@rownames<- df_granges$geneID
    
    set.seed(123)
  
  focal = df_granges[df_granges$DE == "DE"]
  pool = df_granges[df_granges$DE != "DE"]
  mgr <- nullranges::matchRanges(focal = focal,
                          pool = pool,
                          covar = ~ gene_length_slop + Inversion,
                          method = "stratified",
                          replace = F)
  temp2 <- rbind(as.data.frame(focal),
                as.data.frame(mgr)) %>% 
    ungroup()
  
  kw <-  temp2 %>%
    group_by(Inversion) %>% 
    rstatix::kruskal_test(metric ~ DE) %>%
    dplyr::select(-method) %>% 
        mutate(.y. = metrics[i],
           dataset = "matched") 

  es <- temp2 %>%
    group_by(Inversion) %>% 
    rstatix::kruskal_effsize(metric ~ DE) %>%
    dplyr::select(-method) %>% 
        mutate(.y. = metrics[i],
           dataset = "matched") 
  
  kw_test <-  rbind(kw_test,
                    kw %>% 
                      left_join(es, by = join_by(Inversion, .y., n, dataset))) 
  }
  
  kw_test_results <- kw_test %>% 
    dplyr::mutate(padj = p.adjust(p))
  write_tsv(kw_test_results, file = paste("00_data/03_popgen_expression/popgen_DE_kw_results_", comp, ".tsv"))
  kw_test_results
}
  
  
