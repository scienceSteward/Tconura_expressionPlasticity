# _targets.R file
library(targets)
source("R/01_samples_functions.R")
source("R/02_popStructure_functions.R")
source("R/03_windowBased_functions.R")
source("UpSet_mod.R")
tar_option_set(packages = c("tidyverse", 
                            "Hmisc",
                            "ggpubr", 
                            "ggbeeswarm",
                            "ComplexHeatmap", 
                            "circlize", 
                            "ggrastr",
                            "ggnewscale", 
                            "ggforce", 
                            "topGO", 
                            "eulerr", 
                            "ggtree", 
                            "treeio"))


list(
  tar_target(sample_data, "Data/metadata/Tcon_sample_IDs.txt", format = "file"),
  tar_target(sample_df, get_samples(sample_data)),
  tar_target(ecotype_pal_100,  c(  "#9E9AC8", "black", "#74C476")),
  tar_target(pal, "Data/metadata/tcon_pop_pal_100.txt", format = "file"),
  tar_target(pal_100, load_pop_pal(pal)),
  
  ##### Sample assessment
  tar_target(cov_data, "Data/Coverage/covstats.txt"), 
  tar_target(cov_df, read_tsv(cov_data)), 
  tar_target(cov_plot, plot_coverage(cov_df, sample_df, pal_100)), 
  
  ##### Population structure
  #### PCANGSD ####
  tar_target(pca_data, "Data/PCA/PCAngsd_out.cov", format = "file"),
  tar_target(pca_pct, get_pca_pct(pca_data, sample_df)),
  tar_target(pca_df, get_pca_df(pca_data, sample_df)),
  tar_target(pca_1_2, plot_pca(pca_df, pca_pct, c(1,2), pal_100)),
  tar_target(pca_2_3, plot_pca(pca_df, pca_pct, c(2,3), pal_100)),
  tar_target(pca_1_2_save, ggsave(plot = pca_1_2 + theme(legend.position = "none"), filename =  "Plots/Tcon_pca_1_2.pdf", height = 3, width = 3.25)),
  tar_target(pca_2_3_save, ggsave(plot = pca_2_3 + theme(legend.position = "none"), filename = "Plots/Tcon_pca_2_3.pdf", height = 3, width = 3.25)),

  #### admixture ####
  tar_target(admix_pal, "Data/metadata/tcon_admix_pal_100.txt", format = "file"),
  tar_target(admix_pal_100, load_admix_pal(admix_pal)),
  tar_target(admix_files, list.files("Data/NGSadmix", pattern = "qopt", full.names = T)),
  tar_target(admix_df, get_admix_df(admix_files, sample_df)),
  tar_target(k_labels, make_k_labels(c(2,10))),
  tar_target(admix_plot_2_10, plot_admix(admix_df, c(2,10), k_labels, admix_pal_100)),
  tar_target(admix_plot_2_3, plot_admix(admix_df, c(2,3), k_labels, admix_pal_100)),
  tar_target(admix_plot_2_10_save, ggsave(plot = admix_plot_2_10, filename = "Plots/Tcon_admix_plot_2_10.pdf", height = 9.25, width = 8)),
  tar_target(admix_plot_2_3_save, ggsave(plot = admix_plot_2_3, filename = "Plots/Tcon_admix_plot_2_3.pdf", height = 2.25, width = 8)),
  
  #### WG Fst ####
  tar_target(fst_allComps, "Data/metadata/Tcon.allPop.pairwise.txt", format = "file"),
  tar_target(fst_WG_data, "Data/Fst/allPop.norep-sites.saf.wholeGenome.mod.txt", format = "file"),
  tar_target(fst_WG_df, get_wg_fst_df(fst_WG_data, fst_allComps)),
  tar_target(fst_WG_plot_unwt, plot_wg_fst(fst_WG_df, "fst_unwt", c( "goldenrod1", "firebrick1"))),
  tar_target(fst_WG_plot_wt, plot_wg_fst(fst_WG_df, "fst_wt", c( "goldenrod1", "firebrick1"))),

  #### WG Fst ####t
  tar_target(dxy_WG_data, "Data/Dxy/allPops_norep-sites.txt", format = "file"),
  tar_target(dxy_WG_df, get_wg_dxy_df(dxy_WG_data, fst_allComps)),
  tar_target(dxy_WG_plot, plot_wg_fst(dxy_WG_df, "dxy",c("plum", "dodgerblue4"))),

  ##### Demographic structure & Introgression
  #### MSMC ####
  tar_target(msmc_files,  list.files(path= "Data/MSMC", pattern = "8highest_noRepeat_nopbs.final.txt")), 
  tar_target(msmc_df,  get_msmc_df(msmc_files)), 
  tar_target(msmc_plot, plot_msmc(msmc_df, pal_100)),
  
  #### Dstats ####
  tar_target(dstats_comps_file, "Data/Dstats/Tcon_dstats_comps.txt", format = "file"), 
  tar_target(dstats_comps_df, read_tsv(dstats_comps_file)), 
  tar_target(dstats_files, "Data/Dstats/DTparallel_SETS_dsuite_combined_Dmin.txt", format = "file"), 
  tar_target(dstats_df, read_tsv(dstats_files)), 
  tar_target(dstats_plot, plot_dstats(dstats_df, dstats_comps_df, ecotype_pal_100)),
  
  #### Trees #### 
  tar_target(astral_tree_file, "Data/Trees/complete_single_diptera_buscos.astral.nex"), 
  tar_target(astral_tree, read.nexus(astral_tree_file)),
  tar_target(astral_tree_plot, plot_tree(astral_tree, sample_df, pal_100)),

  
  ##### Window Based analyses
  #### General ####
  tar_target(contig_master, read_tsv("Data/refs/pt_042_hifiasm20201214.primary.yahsHiCnoBreak.contigMaster.txt")),
  tar_target(scaffold_master, get_scaffold_master(contig_master)),
  tar_target(sd_cutoff, 3),
  tar_target(GO_path, "Data/refs/braker.adj.UTR.mod.longIso_2line.default.emapper.annotations.GO.mod.tsv"),
  tar_target(gene_names_path, "Data/refs/braker.adj.UTR.mod.longIso_2line.consensus.emapper.annotations.gene_function.gene_name.mod.tsv"),
  
  #### BED #### 
  tar_target(bed_path, "Data/refs/braker.adj.UTR.mod.longIso.mod.slop2kb.50kb.50kb.mod.bed"),
  tar_target(bed, get_bed(bed_path)),
  tar_target(bed_master, bed %>% 
               dplyr::select(chr, slopStart, slopEnd, midPos, geneID) %>% 
               arrange(chr, midPos) %>% 
               left_join(contig_master) %>% unique() %>% 
               arrange(contig_order, midPos) %>% 
               mutate(midPos_scaffold = if_else(contig_dir == "+", midPos, contig_length - midPos),
                      midPos_contigContinuous = midPos_scaffold + contig_start_continuous,
                      uniq = paste(contig_order, midPos, sep = "_"))), 
  
  #### Local PCA ####
  tar_target(local_pca_files, list.files(path = "Data/PCA", pattern = "_1.cov", full.names = T )), 
  tar_target(local_pca_df, get_local_pca(local_pca_files, sample_df)), 
  tar_target(local_pca_plot, plot_local_pca(local_pca_df, contig_master, scaffold_master, pal_100)),
  tar_target(inversion_genes, bed_master %>% 
                                 filter(contig_order >= 611 & contig_order <= 699) %>% 
                                 dplyr::select(geneID) %>% unique()),

  #### Fst ####
  tar_target(fst_files, list.files("Data/Fst", pattern = "windowed", full.names = T)),
  tar_target(fst_list, get_fst_list(fst_files, contig_master)),
  tar_target(fst_lowCov, get_fst_eliminatedSites(fst_files, contig_master)),
  tar_target(fst_plot_list, plot_fst(fst_list, scaffold_master, sd_cutoff)),
  tar_target(fst_comps, read_tsv("Data/metadata/Tcon_FstComps.tsv")),
  tar_target(fst_combined, comb_fst(fst_list, fst_comps)),
  tar_target(fst_dist_plot, plot_fst_dist(fst_combined, ecotype_pal_100)),
  
  #### PBS ####
  tar_target(pbs_files, list.files("Data/PBS", pattern = "pbs.tsv", full.names = T)),
  tar_target(pbs_list, get_pbs_list(pbs_files, contig_master)),
  tar_target(pbs_lowCov, get_pbs_eliminatedSites(pbs_files, contig_master)),
  tar_target(pbs_plot_list, plot_pbs(pbs_list, scaffold_master, sd_cutoff, pal_100)),
  tar_target(pbs_outlier_windows, extract_PBS_outliers(pbs_list, sd_cutoff)),

  #### dxy ####
  tar_target(dxy_files, list.files("Data/Dxy", pattern = "windowed.dxy.txt", full.names = T)),
  tar_target(dxy_list, get_dxy_list(dxy_files, fst_list, contig_master)),
  tar_target(dxy_plot_list, plot_dxy(dxy_list, scaffold_master, sd_cutoff)),
  tar_target(dxy_combined, comb_dxy(dxy_list, fst_comps)),
  tar_target(dxy_dist_plot, plot_dxy_dist(dxy_combined, ecotype_pal_100)),
  
  #### BayPass #### 
  tar_target(baypass_file, "Data/BayPass/all-incl-outgroup_norep-filts.filtered.noOutgroups_ALL_snps.50_50.tsv"),
  tar_target(baypass_df, get_baypass(baypass_file, contig_master)),
  tar_target(baypass_plot_BF_50, plot_baypass(baypass_df, scaffold_master, "BF.db_50", sd_cutoff)),
  tar_target(baypass_plot_M_beta_50, plot_baypass(baypass_df, scaffold_master, "M_beta_abs_50", sd_cutoff)),
  tar_target(baypass_plot_M_XtX_50, plot_baypass(baypass_df, scaffold_master, "M_XtX_50", sd_cutoff)),
  tar_target(baypass_BF_50_outliers, extract_baypass_outliers(baypass_df,"BF.db_50", sd_cutoff)), 
  
  #### RND ####
  tar_target(rnd_files, list.files("Data/RND", pattern = "RND_windowed.tsv", full.names = T)),
  tar_target(rnd_list, get_rnd_list(rnd_files, pbs_list, contig_master)),
  tar_target(rnd_plot_list, plot_rnd(rnd_list, scaffold_master, sd_cutoff, pal_100)),
  # tar_target(pbs_prime_plot_list, plot_PBS_prime(rnd_list, scaffold_master)),
  tar_target(rnd_outlier_windows, extract_rnd_outliers(rnd_list, sd_cutoff)),
  
  
  
  #### Outliers ####
  tar_target(outlier_comb_mat_CO_east, ComplexHeatmap::make_comb_mat(list(`PBS COES`=pbs_outlier_windows$COES,
                                                                          `RND COES`=rnd_outlier_windows$COES.CHES,
                                                                          BayPass = baypass_BF_50_outliers))), 
  tar_target(outlier_comb_mat_CH_east, ComplexHeatmap::make_comb_mat(list(`PBS CHES`=pbs_outlier_windows$CHES,
                                                                          `RND CHES`=rnd_outlier_windows$CHES.COES,
                                                                          BayPass = baypass_BF_50_outliers))), 
  # tar_target(target_combs,  c("11111", "11110", "10111", "11101", "10101", "01011")),
  tar_target(outlier_upset_plot_CH_east, upset_outliers(outlier_comb_mat_CH_east, 18,  pal_100$col[4])),
  tar_target(outlier_upset_plot_CO_east, upset_outliers(outlier_comb_mat_CO_east, 18,  pal_100$col[5])),
  tar_target(outlier_uniq_CH_east, get_outlier_uniq_comb(outlier_comb_mat_CH_east)), 
  tar_target(outlier_uniq_CO_east, get_outlier_uniq_comb(outlier_comb_mat_CO_east)), 
  tar_target(outlier_euler_fit_east, 
             eulerr::euler(list(`CH East` = outlier_uniq_CH_east, `CO East` = outlier_uniq_CO_east))),
  tar_target(outlier_euler_plots_east, 
             plot(outlier_euler_fit_east, quantities = T, fill = pal_100$col[c(4,5)])),
  
  tar_target(outlier_comb_mat_CH_west, ComplexHeatmap::make_comb_mat(list(`PBS CHSK`=pbs_outlier_windows$CHSK, 
                                                                       `RND CHSK`=rnd_outlier_windows$CHSK.COSK, 
                                                                       BayPass = baypass_BF_50_outliers))),
  tar_target(outlier_comb_mat_CO_west, ComplexHeatmap::make_comb_mat(list(`PBS COSK`=pbs_outlier_windows$COSK, 
                                                                       `RND COSK`=rnd_outlier_windows$COSK.CHSK, 
                                                                      BayPass = baypass_BF_50_outliers))),
  tar_target(outlier_upset_plot_CH_west, upset_outliers(outlier_comb_mat_CH_west, 16, pal_100$col[1])),
  tar_target(outlier_upset_plot_CO_west, upset_outliers(outlier_comb_mat_CO_west, 16, pal_100$col[8])),
  tar_target(outlier_uniq_CH_west, get_outlier_uniq_comb(outlier_comb_mat_CH_west)), 
  tar_target(outlier_uniq_CO_west, get_outlier_uniq_comb(outlier_comb_mat_CO_west)), 
  tar_target(outlier_euler_fit_west, eulerr::euler(list(`CH West` = outlier_uniq_CH_west, `CO West` = outlier_uniq_CO_west))),
  tar_target(outlier_euler_plots_west, 
             plot(outlier_euler_fit_west, quantities = T, fill = pal_100$col[c(1,8)])),
 
  
   tar_target(outlier_venn_fit_all, 
             eulerr::venn(list(`CO West` = outlier_uniq_CO_west, `CO East` = outlier_uniq_CO_east, 
                                `CH West` = outlier_uniq_CH_west, `CH East` = outlier_uniq_CH_east ))),
  tar_target(outlier_venn_plots_all, 
             plot(outlier_venn_fit_all, quantities = T, fill = pal_100$col[c(8,5,1,4)])),
  tar_target(outlier_uniq, list(
    intersect(outlier_uniq_CO_west, outlier_uniq_CH_west), 
              intersect(outlier_uniq_CO_east, outlier_uniq_CH_east))),
  tar_target(outlier_genes_intersect, bed_master %>% 
               filter(uniq %in% intersect(outlier_uniq[[1]], outlier_uniq[[2]])) %>% 
               dplyr::select(geneID) %>% unique()),
  tar_target(outlier_genes_save, write(x = outlier_genes_intersect$geneID, file = "Data/Outliers/outliers_intersect_STRING.txt")),
    tar_target(tcon_topgo_outliers, topgo_outliers(outlier_genes_intersect, GO_path)),
  tar_target(tcon_topgo_outlier_plots, plot_topgo_outliers(tcon_topgo_outliers)),
  
  #### Pi & D ####
  tar_target(theta_files, list.files("Data/Theta", pattern = "pestPG", full.names = T)),
  tar_target(theta_list, get_theta_list(theta_files, contig_master)),
  tar_target(pi_plot_list, plot_theta(theta_list, outlier_uniq, scaffold_master, sd_cutoff, pal_100, "nucl_div")),
  tar_target(tajD_plot_list, plot_theta(theta_list, outlier_uniq, scaffold_master, sd_cutoff, pal_100, "Tajima")),
  tar_target(theta_combined, comb_theta(theta_list, outlier_uniq)),
  tar_target(pi_dist_plot, plot_theta_dist(theta_combined, pal_100, "nucl_div" )),
  tar_target(tajD_dist_plot, plot_theta_dist(theta_combined, pal_100, "Tajima" )),
  tar_target(theta_diff_list, get_diff(theta_combined, fst_comps)),
  tar_target(pi_plot_diff_list, plot_diff(theta_diff_list, scaffold_master, ecotype_pal_100, "nucl_div_diff", sd_cutoff)),
  tar_target(tajD_plot_diff_list, plot_diff(theta_diff_list, scaffold_master, ecotype_pal_100, "Tajima_diff", sd_cutoff)),
  tar_target(theta_diff_combined, comb_theta_diff(theta_diff_list, fst_comps, outlier_uniq)),
  tar_target(pi_diff_dist_plot, plot_theta_diff_dist(theta_diff_combined, ecotype_pal_100, "nucl_div_diff")),
  tar_target(tajD_diff_dist_plot, plot_theta_diff_dist(theta_diff_combined, ecotype_pal_100, "Tajima_diff")),
  
  #### XP-EHH ####
  tar_target(xpehh_files, list.files("Data/xpEHH", pattern = ".50kb.windows", full.names = T)),
  tar_target(xpehh_list, get_xpehh_list(xpehh_files, contig_master)),
  tar_target(xpehh_plot_list, plot_xpehh(xpehh_list, outlier_uniq, scaffold_master, sd_cutoff, pal_100)),
  tar_target(xpehh_combined, comb_xpehh(xpehh_list, outlier_uniq)),
  tar_target(xpehh_All_Average_dist_plot, plot_xpehh_dist(xpehh_combined, pal_100, "All_Average" )),
  tar_target(xpehh_average_dist_plot, plot_xpehh_dist(xpehh_combined, pal_100, "average" )),
  tar_target(selected_genes_files, list.files("Data/xpEHH", pattern = "geneID.txt", full.names = T)),
  tar_target(selected_genes_list, get_selected_genes_list(selected_genes_files)),
  tar_target(selected_genes_euler, plot_selected_genes_euler(selected_genes_list, outlier_genes_intersect, inversion_genes, pal_100)),
  
  #### Dstat windows ####
  tar_target(dstat_window_files, list.files(path = "Data/Dstats", pattern = "_unpolarized_fd.csv", full.names = T )),
  tar_target(dstat_window_list, get_dstats(dstat_window_files, dstats_comps_df, contig_master)),
  tar_target(fdm_plot_list, plot_dstats_window(dstat_window_list,outlier_uniq,  scaffold_master, "fdm.", sd_cutoff, pal_100)),
  tar_target(fd_plot_list, plot_dstats_window(dstat_window_list, outlier_uniq, scaffold_master, "fd", sd_cutoff, pal_100)),
  tar_target(dstat_window_combined, comb_dstat_window(dstat_window_list, outlier_uniq)),
  tar_target(dstat_window_mod, lme4::lmer(fdm~P3/is_outlier + (1|uniq), data = dstat_window_combined)),
  tar_target(fdm_window_dist_plot, plot_dstat_diff_dist(dstat_window_combined, pal_100, "fdm")),
  tar_target(fd_window_dist_plot, plot_dstat_diff_dist(dstat_window_combined, pal_100, "fd")), 
  
  #### Outside the inversion #### 
  tar_target(baypass_BF_50_noInversion_outliers, extract_baypass_outliers_noInversion(baypass_df,"BF.db_50", sd_cutoff)),
  tar_target(rnd_outlier_noInversion_windows, extract_rnd_outliers_noInversion(rnd_list, sd_cutoff)),
  tar_target(pbs_outlier_noInversion_windows, extract_PBS_outliers_noInversion(pbs_list, sd_cutoff)), 
  
  
  tar_target(outlier_comb_mat_CO_east_noInversion, ComplexHeatmap::make_comb_mat(list(`PBS COES`=pbs_outlier_noInversion_windows$COES,
                                                                          `RND COES`=rnd_outlier_noInversion_windows$COES.CHES,
                                                                          BayPass = baypass_BF_50_noInversion_outliers))), 
  tar_target(outlier_comb_mat_CH_east_noInversion, ComplexHeatmap::make_comb_mat(list(`PBS CHES`=pbs_outlier_noInversion_windows$CHES,
                                                                          `RND CHES`=rnd_outlier_noInversion_windows$CHES.COES,
                                                                          BayPass = baypass_BF_50_noInversion_outliers))), 

  tar_target(outlier_upset_plot_CH_east_noInversion, upset_outliers(outlier_comb_mat_CH_east_noInversion, 18,  pal_100$col[4])),
  tar_target(outlier_upset_plot_CO_east_noInversion, upset_outliers(outlier_comb_mat_CO_east_noInversion, 18,  pal_100$col[5])),
  tar_target(outlier_uniq_CH_east_noInversion, get_outlier_uniq_comb(outlier_comb_mat_CH_east_noInversion)), 
  tar_target(outlier_uniq_CO_east_noInversion, get_outlier_uniq_comb(outlier_comb_mat_CO_east_noInversion)), 
  
  tar_target(outlier_comb_mat_CH_west_noInversion, ComplexHeatmap::make_comb_mat(list(`PBS CHSK`=pbs_outlier_noInversion_windows$CHSK, 
                                                                          `RND CHSK`=rnd_outlier_noInversion_windows$CHSK.COSK, 
                                                                          BayPass = baypass_BF_50_noInversion_outliers))),
  tar_target(outlier_comb_mat_CO_west_noInversion, ComplexHeatmap::make_comb_mat(list(`PBS COSK`=pbs_outlier_noInversion_windows$COSK, 
                                                                          `RND COSK`=rnd_outlier_noInversion_windows$COSK.CHSK, 
                                                                          BayPass = baypass_BF_50_noInversion_outliers))),
  tar_target(outlier_upset_plot_CH_west_noInversion, upset_outliers(outlier_comb_mat_CH_west_noInversion, 16, pal_100$col[1])),
  tar_target(outlier_upset_plot_CO_west_noInversion, upset_outliers(outlier_comb_mat_CO_west_noInversion, 16, pal_100$col[8])),
  tar_target(outlier_uniq_CH_west_noInversion, get_outlier_uniq_comb(outlier_comb_mat_CH_west_noInversion)), 
  tar_target(outlier_uniq_CO_west_noInversion, get_outlier_uniq_comb(outlier_comb_mat_CO_west_noInversion)), 
  tar_target(outlier_list_noInversion, list( `CH East` = outlier_uniq_CH_east_noInversion, 
                                             `CH West` = outlier_uniq_CH_west_noInversion,
                                             `CO East` = outlier_uniq_CO_east_noInversion,
                                             `CO West` = outlier_uniq_CO_west_noInversion)), 
  
  tar_target(outlier_euler_fit_west_noInversion, eulerr::euler(outlier_list_noInversion[c(2,4)])),
  tar_target(outlier_euler_plots_west_noInversion, 
             plot(outlier_euler_fit_west_noInversion, quantities = T, fill = pal_100$col[c(1,8)])),
  
  tar_target(outlier_euler_fit_east_noInversion, eulerr::euler(outlier_list_noInversion[c(1,3)])),
  tar_target(outlier_euler_plots_east_noInversion, 
             plot(outlier_euler_fit_east_noInversion, quantities = T, fill = pal_100$col[c(4,5)])),
  
  tar_target(outlier_venn_fit_all_noInversion, eulerr::venn(outlier_list_noInversion[c(4, 3, 2, 1)])), 
  tar_target(outlier_venn_plots_all_noInversion, 
             plot(outlier_venn_fit_all_noInversion, quantities = T, fill = pal_100$col[c(8,5,1,4)])),
  
  tar_target(outlier_genes_all_noInversion, bed_master %>% 
               filter(uniq %in% Reduce(intersect,outlier_list_noInversion)) %>% 
               dplyr::select(geneID) %>% unique()),
  
  tar_target(outlier_genes_save_noInversion, write(x = outlier_genes_all_noInversion$geneID, 
                                                   file = "Data/Outliers/outliers_noInversion_intersect_STRING.txt")),
  tar_target(tcon_topgo_outliers_noInversion, topgo_outliers(outlier_genes_all_noInversion %>% 
                                                               filter(str_sub(geneID, 1, 9) != "Tcon_g825" &
                                                                        str_sub(geneID, 1, 9) != "Tcon_g826" & 
                                                                        str_sub(geneID, 1, 10) != "Tcon_g2362" ), GO_path)),
  tar_target(tcon_topgo_outlier_plots_noInversion, plot_topgo_outliers(tcon_topgo_outliers_noInversion)),
  tar_target(gene_loc_noInversion, read_tsv(gene_names_path)%>%
               right_join(outlier_genes_all_noInversion) %>% left_join(bed_master) %>% 
               mutate(slopMidPos = contig_start_continuous + slopStart + (slopEnd - slopStart)/2) %>% 
               dplyr::select(chr, scaffold, geneID, cons_desc, cons_name, slopMidPos) %>% unique()),
  tar_target(baypass_plot_BF_50_noInversion, plot_baypass_noInversion(baypass_df, gene_loc_noInversion, scaffold_master, "BF.db_50", sd_cutoff)),
  tar_target(pbs_plot_list_noInversion, plot_pbs_noInversion(pbs_list, gene_loc_noInversion, scaffold_master, sd_cutoff, pal_100)),
  tar_target(rnd_plot_list_noInversion, plot_rnd_noInversion(rnd_list, gene_loc_noInversion, scaffold_master, sd_cutoff, pal_100))
  
)

