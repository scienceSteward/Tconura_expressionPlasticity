# _targets.R file
library(targets)
source("Tcon_CF_DGE_continental_funcs.R")
tar_option_set(
  packages = c("tidyverse", 
               "DESeq2", 
               "ggpubr", 
               "ggrepel",
               "ggtext",
               "grid",
               "topGO",
               "GOSemSim",
               "Mfuzz",
               "WGCNA", 
               "circlize",
               "ComplexHeatmap", 
               "ggbeeswarm", 
               "rstatix", 
               "nullranges"))


list(
  tar_target(cf_pal, c("#54278F", "#9E9AC8", "#DADAEB", "#006D2C", "#74C476", "#C7E9C0")),
  
  #### Access sample data ####
  tar_target(sample_data, "00_data/Tconura_RNAsampleMetadata.txt"), 
  tar_target(cf_samples, read_cf_samples(sample_data, c("P18653_176"))),
  
  #### Check mapping rates ####
  tar_target(mapRates, "03_pseudoCounts/mapping_rate_braker.txt"),
  tar_target(mapRates_df, load_mapping_rates(mapRates, cf_samples)),
  tar_target(mapRates_plot, plot_mapping_rates(mapRates_df, cf_pal)),
  
  #### Compile trimming, quantification statistics ####
  tar_target(pretrim_mqc, "00_data/01_RNA/multiqc_general_stats_raw.txt"),
  tar_target(posttrim_mqc, "00_data/01_RNA/multiqc_general_stats_trimmed.txt"),
  tar_target(mqc, comb_multiqc(pretrim_mqc, pretrim_mqc, cf_samples)),
  tar_target(save_mqc, write_tsv(mqc, file = "00_data/01_RNA/multiqc_stats.tsv")),
  
  #### Load dds object ####
  tar_target(dds_rdata, "03_pseudoCounts/salmon_braker_txi_dds_cfLarvae.Rdata"),
  tar_target(tcon_dds, load_dds(dds_rdata)),
  tar_target(tcon_dds_sub, subset_by_samples(tcon_dds, cf_samples)),
  tar_target(tcon_dds_filt, filter_by_counts(tcon_dds_sub)),
  tar_target(tcon_expr_genes, rownames(tcon_dds_filt)),
  
  #### Unsupervised clustering: PCA ####
  tar_target(tcon_rlog, assay(tcon_dds_filt %>% rlog())),
  tar_target(tcon_pca, rlog2pca(tcon_rlog, 5000)),
  tar_target(tcon_pca_var, prcompVarSumm(tcon_pca)),
  tar_target(tcon_pca_var_plot, plot_var_comp(tcon_pca_var)),
  tar_target(tcon_pca_df, prcomp2df(tcon_pca, cf_samples)),
  tar_target(tcon_pca_plotlist, plot_pca_list(tcon_pca_var,tcon_pca_df, cf_pal)),
  
  #### Unsupervised clustering: MDS ####
  tar_target(tcon_mds, edgeR::plotMDS.DGEList(tcon_dds_filt, top = 5000)),
  tar_target(tcon_mds_df, mds2df(tcon_mds, cf_samples)), 
  tar_target(tcon_mds_plot, plot_mds(tcon_mds, tcon_mds_df, cf_pal)),
  
  #### Unsupervised clustering: heatmap ####
  tar_target(tcon_heatmap, rlog2heatmap(tcon_rlog, 5000, cf_samples, cf_pal)),
  
  #### DGE ####
  tar_target(cond_mat, read.table("00_data/conditionMatrix.txt")),
  # tar_target(tcon_deseq_results, run_deseq_cond_mat(tcon_dds_filt, cond_mat)),
  tar_target(tcon_deseq_results, run_deseq_cond_mat_fullMat(tcon_dds_filt, cond_mat)),
  tar_target(tcon_deseq_pval_plot, plot_pvals(tcon_deseq_results)),
  tar_target(tcon_deseq_volc_plot, plot_volcano(tcon_deseq_results, cf_pal)),
  tar_target(tcon_deseq_MA_plot, plot_MA(tcon_deseq_results, cf_pal)),
  tar_target(tcon_deseq_df, make_results_df(tcon_deseq_results)),
  tar_target(HH_HO_vs_HH_OO_scatter, plot_scatter(tcon_deseq_df, cf_pal, 
                                                  c("HH_vs_OO", "HH_vs_HO"), FALSE)),
  tar_target(OO_OH_vs_OO_HH_scatter, plot_scatter(tcon_deseq_df, cf_pal, 
                                                  c("HH_vs_OO", "OO_vs_OH"), TRUE)),
  tar_target(H_HO_vs_H_O_scatter, plot_scatter(tcon_deseq_df, cf_pal, 
                                               c("H_vs_O", "H_vs_HO"), FALSE)),
  tar_target(O_OH_vs_O_H_scatter, plot_scatter(tcon_deseq_df, cf_pal, 
                                               c("H_vs_O", "O_vs_OH"), TRUE)),
  tar_target(OH_HO_vs_HH_OO_scatter, plot_scatter(tcon_deseq_df, cf_pal,
                                                  c("HH_vs_OO", "HO_vs_OH"), FALSE)),
  tar_target(OH_HO_vs_H_O_scatter, plot_scatter(tcon_deseq_df, cf_pal,
                                                  c("H_vs_O", "HO_vs_OH"), FALSE)),
  tar_target(HH_OO_vs_H_O_scatter, plot_scatter(tcon_deseq_df, c(cf_pal), 
                                               c("H_vs_O", "HH_vs_OO"), FALSE)),
  
  #### Unsupervised clustering of DE genes: heatmap ####
  # tar_target(tcon_heatmap_DE, rlog2heatmap(tcon_rlog[rownames(tcon_rlog) %in% unique(tcon_deseq_df$geneID[tcon_deseq_df$padj < 0.05 & abs(tcon_deseq_df$log2FoldChange) > 0.1]),], 2280, cf_samples, cf_pal)),
  
  #### Functional annotation ####
  tar_target(GO_path, "00_data/00_refs/braker.adj.UTR.mod.longIso_2line.default.emapper.annotations.GO.mod.tsv"),
  tar_target(tcon_topgo_list, topgo_genelist(tcon_deseq_results, GO_path)),
  tar_target(tcon_topgo_list_save, write_tsv(tcon_topgo_list, "00_data/02_DEgenes/Tcon_topGO_list.tsv")),
  tar_target(tcon_topgo_plots, plot_topgo(tcon_topgo_list, cf_pal)),
  tar_target(tcon_topgo_list_both, topgo_genelist_both(tcon_deseq_results, GO_path)),
  tar_target(tcon_topgo_plots_both, plot_topgo_both(tcon_topgo_list_both, cf_pal)),
  tar_target(HH_OO_vs_H_O_both, read_tsv("00_data/02_DEgenes/Tcon_H_vs_O_HH_vs_OO_both_DE.tsv")),
  tar_target(tcon_topgo_HH_OO_vs_H_O_both, topgo_targeted(HH_OO_vs_H_O_both, GO_path, "HH_OO_vs_H_O_both")),
  tar_target(tcon_topgo_plots_HH_OO_vs_H_O_both, plot_topgo_targeted(tcon_topgo_HH_OO_vs_H_O_both, "gold", "HH_OO_vs_H_O_both")),
  
  tar_target(HO_OH_vs_H_O_both, read_tsv("00_data/02_DEgenes/Tcon_H_vs_O_HO_vs_OH_both_DE.tsv")),
  tar_target(tcon_topgo_HO_OH_vs_H_O_both, topgo_targeted(HO_OH_vs_H_O_both, GO_path, "HO_OH_vs_H_O_both")),
  tar_target(tcon_topgo_plots_HO_OH_vs_H_O_both, plot_topgo_targeted(tcon_topgo_HO_OH_vs_H_O_both, "gold", "HO_OH_vs_H_O_both")),
  
  tar_target(HO_OH_vs_HH_OO_both, read_tsv("00_data/02_DEgenes/Tcon_HH_vs_OO_HO_vs_OH_both_DE.tsv")),
  tar_target(tcon_topgo_HO_OH_vs_HH_OO_both, topgo_targeted(HO_OH_vs_HH_OO_both, GO_path, "HO_OH_vs_HH_OO_both")),
  tar_target(tcon_topgo_plots_HO_OH_vs_HH_OO_both, plot_topgo_targeted(tcon_topgo_HO_OH_vs_HH_OO_both, "gold", "HO_OH_vs_HH_OO_both")),
  
  tar_target(H_OH_vs_H_O_both, read_tsv("00_data/02_DEgenes/Tcon_H_vs_O_H_vs_HO_both_DE.tsv")),
  tar_target(tcon_topgo_H_OH_vs_H_O_both, topgo_targeted(H_OH_vs_H_O_both, GO_path, "H_OH_vs_H_O_both")),
  tar_target(tcon_topgo_plots_H_OH_vs_H_O_both, plot_topgo_targeted(tcon_topgo_H_OH_vs_H_O_both, "gold", "H_OH_vs_H_O_both")),
  
  tar_target(stringInput, save_genesets(tcon_deseq_results)),
  
  #### Save plots ####
  tar_target(tcon_pca_var_plot_saved, ggsave("04_targetPlots/tcon_pca_var_plot.png", tcon_pca_var_plot, height = 4, width = 5), format = "file"), 
  tar_target(tcon_pca_panel_1_2, pca_panel(tcon_pca_plotlist, c(1,2))),
  tar_target(tcon_pca_panel_1_2_saved, ggsave("04_targetPlots/tcon_pca_panel_1_2.png", tcon_pca_panel_1_2, height = 2.5, width = 5.5), format = "file"), 
  tar_target(tcon_pca_panel_2_3, pca_panel(tcon_pca_plotlist, c(2,3))),
  tar_target(tcon_pca_panel_2_3_saved, ggsave("04_targetPlots/tcon_pca_panel_2_3.png", tcon_pca_panel_2_3, height = 2.5, width = 5.5), format = "file"),
  tar_target(tcon_pca_panel_3_4, pca_panel(tcon_pca_plotlist, c(3,4))),
  tar_target(tcon_pca_panel_3_4_saved, ggsave("04_targetPlots/tcon_pca_panel_3_4.png", tcon_pca_panel_3_4, height = 2.5, width = 5.5), format = "file"),
  tar_target(tcon_volc_save, save_volc_list(tcon_deseq_volc_plot, "04_targetPlots/")),
  tar_target(tcon_ma_panel_save, ggsave("04_targetPlots/tcon_MA_panel.png", tcon_deseq_MA_plot, height = 8, width = 10), format = "file"),
  tar_target(HH_HO_vs_HH_OO_scatter_save, ggsave("04_targetPlots/HH_HO_vs_HH_OO_scatter.pdf", HH_HO_vs_HH_OO_scatter, height = 3, width = 4.5), format = "file"),
  tar_target(OO_OH_vs_OO_HH_scatter_save, ggsave("04_targetPlots/OO_OH_vs_OO_HH_scatter.pdf", OO_OH_vs_OO_HH_scatter, height = 3.5, width = 4.5), format = "file"),
  tar_target(HH_OO_vs_H_O_scatter_save, ggsave("04_targetPlots/HH_OO_vs_H_O_scatter.pdf", HH_OO_vs_H_O_scatter, height = 3, width = 4.5), format = "file"),
  
  #### WGCNA ####
  # tar_target(scaledata, standardiseGenes(tcon_rlog)),
  tar_target(scaledata, rename_object(tcon_rlog)),
  tar_target(datTraits, get_datTraits(scaledata, cf_samples)),
  tar_target(gsg, goodSamplesGenes(t(scaledata), verbose = 3)),
  tar_target(tcon_dendro, plot_dendro(scaledata, datTraits, cf_pal)),
  tar_target(sft, pick_threshold(scaledata)),
  tar_target(sft_plot, plot_soft_threshold(sft)),
  tar_target(sft_plot_save, ggsave(plot = sft_plot, filename = "04_targetPlots/tcon_networkSoftThreshold.pdf", height = 4, width = 8), format = "file"),
  # tar_target(net, run_modules(scaledata, sft)),
  tar_target(net, load_modules("00_data/Tcon_network.Rdata")),
   tar_target(net_plot, plot_network(net)), #update path in function
  tar_target(modules, pull_modules(net, scaledata, datTraits)),
  tar_target(ht_host, draw_cor_heatmap(modules, "Host", cf_pal)),
  tar_target(ht_stress, draw_cor_heatmap(modules, "Stress", cf_pal)),
  tar_target(ht_plasticity, draw_cor_heatmap(modules, "Plasticity", cf_pal)), 
  tar_target(ht_CH_plasticity, draw_cor_heatmap(modules, "CHplasticity", cf_pal)), 
  tar_target(ht_CO_plasticity, draw_cor_heatmap(modules, "COplasticity", cf_pal)), 
  tar_target(ht_modules, draw_module_heatmap(modules, datTraits, cf_pal)), 
  tar_target(ht_combo, draw_cor_module_combo_plot(ht_host, 
                                                  ht_stress, 
                                                  ht_plasticity, 
                                                  ht_CH_plasticity, 
                                                  ht_CO_plasticity, 
                                                  ht_modules)),
  tar_target(moduleMembership, get_moduleMembership(scaledata, modules, datTraits)),
  tar_target(moduleMembership_plot, plot_moduleMembership(moduleMembership, cf_pal)),
  tar_target(moduleMembership_plot_save, ggsave(plot = moduleMembership_plot, 
                                                filename = "04_targetPlots/moduleMembership_plot.pdf", 
                                                height = 15, width = 7)),
  tar_target(moduleMembershipSumm, summarize_moduleMembership(moduleMembership)),
  tar_target(moduleMembershipSumm_plots, plot_moduleMembershipSumm(moduleMembershipSumm,
                                                                   cf_pal)),
  tar_target(moduleEigengenes_plots, plot_moduleEigengenes(datTraits, modules,
                                                                   cf_pal)),
  tar_target(targetMods, c("mod00", "mod01", "mod06", "mod16", "mod15", "mod18", "mod19", "mod20", "mod04")),
  tar_target(targetMods_pal, c(cf_pal[1], cf_pal[1], cf_pal[1], cf_pal[3], cf_pal[1], "#f05039", "#f05039", cf_pal[4], cf_pal[4])),
  tar_target(targetMod_genes, extract_mod_genes(moduleMembership, targetMods)), 
  tar_target(targetMod_topgo, topgo_modList(targetMod_genes, GO_path, targetMods)),
  tar_target(targetMod_topgo_plots, plot_topgo_modList(targetMod_topgo, targetMods, targetMods_pal)),
  
  
  #### Population genomics ####
  tar_target(contig_master, read_tsv("00_data/00_refs/pt_042_hifiasm20201214.primary.yahsHiCnoBreak.contigMaster.txt")),
  tar_target(scaffold_master, get_scaffold_master(contig_master)),
  tar_target(fst_files, list.files(path= "05_popGen", pattern = "fst.windowed", full.names = T)),
  tar_target(fst_list, get_fst_list(fst_files, contig_master)),
  tar_target(fst_plot_list, plot_fst(fst_list, scaffold_master)),
  tar_target(dxy_files, list.files("05_popGen", pattern = "windowed.dxy", full.names = T)),
  tar_target(dxy_list, get_dxy_list(dxy_files, contig_master)),
  tar_target(dxy_plot_list, plot_dxy(dxy_list, scaffold_master)),
  tar_target(theta_files, list.files(path = "05_popGen", pattern = "pestPG", full.names = T)),
  tar_target(theta_list, get_theta_list(theta_files, contig_master)),
  tar_target(pi_plot_list, plot_pi(theta_list, scaffold_master)),
  tar_target(tajD_plot_list, plot_tajD(theta_list, scaffold_master)),
  tar_target(theta_combined, comb_theta(theta_list)),
  tar_target(pi_dist_plot, plot_theta_dist(theta_combined, cf_pal[c(1,4)], "nucl_div" )),
  tar_target(tajD_dist_plot, plot_theta_dist(theta_combined, cf_pal[c(1,4)], "Tajima" )),
  tar_target(theta_diff_list, get_diff(theta_combined, c("CHST", "COSK"))),
  tar_target(pi_plot_diff_list, plot_diff(theta_diff_list, scaffold_master,  cf_pal[c(1,4)], "nucl_div_diff")),
  tar_target(tajD_plot_diff_list, plot_diff(theta_diff_list, scaffold_master, cf_pal[c(1,4)], "Tajima_diff")),
  tar_target(pi_D_panel,
             annotate_figure(ggarrange(pi_plot_list[[1]] + 
                                         theme(axis.title.x = element_blank()), 
                                       pi_plot_list[[2]]+ 
                                         theme(axis.title.x = element_blank()),
                                       tajD_plot_list[[1]] + 
                                         theme(axis.title.x = element_blank()), 
                                       tajD_plot_list[[2]]+ 
                                         theme(axis.title.x = element_blank()), 
                                       ncol = 1, labels = "AUTO", align = "hv"),
                             bottom = ggpubr::text_grob("SNP position"))),
  tar_target(manhattan_panel,
             annotate_figure(ggarrange(fst_plot_list[[1]] + 
                                         theme(axis.title.x = element_blank()), 
                                       dxy_plot_list[[1]] +
                                         theme(axis.title.x = element_blank()),
                                       pi_plot_diff_list + 
                                         theme(axis.title.x = element_blank()), 
                                       tajD_plot_diff_list + 
                                         theme(axis.title.x = element_blank()), 
                                       ncol = 1, labels = "AUTO", align = "hv"),
                             bottom = ggpubr::text_grob("SNP position"))),
  
  #### Outlier windows and Genes ####
  tar_target(fst_outliers, get_outliers_1way(fst_list[[1]], "Fst")),
  tar_target(dxy_outliers, get_outliers_1way(dxy_list[[1]], "dxy")),
  tar_target(pi_diff_outliers, get_outliers_2way(theta_diff_list, "nucl_div_diff")),
  tar_target(tajD_diff_outliers, get_outliers_2way(theta_diff_list, "Tajima_diff")),
  tar_target(outlier_comb_mat, ComplexHeatmap::make_comb_mat(list(Fst=fst_outliers, 
                                                                  dxy=dxy_outliers,
                                                     "\u0394\u03C0" = pi_diff_outliers, 
                                                     "\u0394D" = tajD_diff_outliers))),
  tar_target(outlier_upset_plot, upset_outliers(outlier_comb_mat)),
  tar_target(bed_path, "00_data/00_refs/braker.adj.UTR.mod.longIso.mod.slop2kb.50kb.50kb.mod.bed"),
  tar_target(bed, get_bed(bed_path)),  
  tar_target(gene_path, "00_data/00_refs/braker.adj.UTR.mod.longIso_2line.default.emapper.annotations.gene_function.gene_name.tsv"),
  tar_target(gene_names, read_tsv(gene_path, col_names = c("geneID", "geneFunction", "geneName"))),
  tar_target(outlier_genes_comb_list, get_outlier_genes_comb(outlier_comb_mat, bed, contig_master)),
  tar_target(tcon_topgo_outliers, topgo_outliers(outlier_genes_comb_list[1:4], GO_path)),
  tar_target(tcon_topgo_outlier_plots, plot_topgo_outliers(tcon_topgo_outliers)),
  
  tar_target(inversion_genes, contig_master %>% filter(contig_order <=699 & contig_order >=611) %>% left_join(bed) %>% dplyr::select(geneID)),
  
  #### DEG vs. outlier comparisons ####
  tar_target(DE_inversion_all_euler_fits, fit_euler_DE_all(tcon_deseq_results$up, tcon_deseq_results$down, inversion_genes)),
  tar_target(DE_inversion_expr_euler_fits, fit_euler_DE_expr(tcon_deseq_results$up, tcon_deseq_results$down, inversion_genes, tcon_expr_genes)),
  tar_target(DE_inversion_all_hypergeom_test, hypergeom_test_DE(DE_inversion_all_euler_fits, 25175)),
  tar_target(DE_inversion_expr_hypergeom_test, hypergeom_test_DE(DE_inversion_expr_euler_fits, 11701)),
  
  
  tar_target(DE_outliers_all_euler_fits, fit_euler_DE_all(tcon_deseq_results$up, tcon_deseq_results$down, outlier_genes_comb_list[1:4])),
  tar_target(DE_outliers_expr_euler_fits, fit_euler_DE_expr(tcon_deseq_results$up, tcon_deseq_results$down, outlier_genes_comb_list[1:4], tcon_expr_genes)),
  tar_target(DE_outliers_all_hypergeom_test, hypergeom_test_DE(DE_outliers_all_euler_fits, 25175)),
  tar_target(DE_outliers_expr_hypergeom_test, hypergeom_test_DE(DE_outliers_expr_euler_fits, 11701)),
  
  # tar_target(DE_euler_plots, plot_euler_DE(DE_outliers_all_euler_fits, cf_pal), format = "file"),
  tar_target(module_outliers_all_euler_fits, fit_euler_module_all(targetMod_genes, 
                                                 outlier_genes_comb_list[1:4])),
  tar_target(module_outliers_expr_euler_fits, fit_euler_module_expr(targetMod_genes, 
                                                              outlier_genes_comb_list[1:4], 
                                                              tcon_expr_genes)),
  tar_target(module_outliers_all_hypergeom_test, hypergeom_test_module(module_outliers_all_euler_fits, 25175)),
  tar_target(module_outliers_expr_hypergeom_test, hypergeom_test_module(module_outliers_expr_euler_fits, 11701)),
  
  #### DEG and Pop Gen metrics ####
  tar_target(H_vs_O_locs, locate_deg(tcon_deseq_results, "H_vs_O", bed, contig_master)),
  tar_target(HH_vs_OO_locs, locate_deg(tcon_deseq_results, "HH_vs_OO", bed, contig_master)),
  tar_target(HH_vs_HO_locs, locate_deg(tcon_deseq_results, "HH_vs_HO", bed, contig_master)),
  tar_target(OO_vs_OH_locs, locate_deg(tcon_deseq_results, "OO_vs_OH", bed, contig_master)),
  tar_target(mod16_loc, locate_modgenes(targetMod_genes[[4]], bed, contig_master)), 
  tar_target(bed_popgen_df, add_popgen_to_bed(bed, contig_master, fst_list[[1]], dxy_list[[1]], theta_diff_list, theta_list)), 
  tar_target(modules_popgen_df, add_popgen_to_modules(moduleMembership, modules, bed_popgen_df)),
  tar_target(modules_popgen_plots, plot_popgen_modules(modules_popgen_df, cf_pal)),
  tar_target(HH_vs_OO_popgen_DE_kw_test, popgen_DE_kw_matched(bed_popgen_df, tcon_deseq_results, "HH_vs_OO")),
  tar_target(H_vs_O_popgen_DE_kw_test, popgen_DE_kw_matched(bed_popgen_df, tcon_deseq_results, "H_vs_O")),
  tar_target(HH_vs_OO_popgen_plots, plot_popgen_DE(bed_popgen_df, tcon_deseq_results, cf_pal, "HH_vs_OO", HH_vs_OO_popgen_DE_kw_test)),
  tar_target(H_vs_O_popgen_plots, plot_popgen_DE(bed_popgen_df, tcon_deseq_results, cf_pal, "H_vs_O", H_vs_O_popgen_DE_kw_test))
)
  

