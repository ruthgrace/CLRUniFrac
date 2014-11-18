####Files in [test_script](clr_vs_proportion_unifrac/brazil_study_test) folder

#####[pcoa_test_script.r](clr_vs_proportion_unifrac/brazil_study_test/pcoa_test_script.r)
Test script that runs GUniFrac and CLRUniFrac on sample data.

#####[metrics.r](clr_vs_proportion_unifrac/brazil_study_test/metrics.r)
Script that calculates overlap and average read count between all pairs of samples.

#####[test_plots_with_brazil_study_data.pdf](clr_vs_proportion_unifrac/brazil_study_test/test_plots_with_brazil_study_data.pdf)
Plots that compare proportional abundance weighted GUniFrac, CLRUnifrac, and QIIME output. Plots are:
* overlap between samples vs. CLRUniFrac distances between samples
* overlap between samples vs. GUniFrac distances between samples
* average total read count between samples vs. CLRUniFrac distances between samples
* average total read count between samples vs. GUniFrac distances between samples
* Principle Coordinates of Analysis plot of CLRUniFrac distances
* Colored by Bacterial Vaginosis (as diagnosed by Nugent and Amsel), Intermediate, and Normal. Samples which are composed of at least 50% of the same taxa are colored by that taxa.
* Principle Coordinates of Analysis plot of GUniFrac distances
* Colored by Bacterial Vaginosis (as diagnosed by Nugent and Amsel), Intermediate, and Normal. Samples which are composed of at least 50% of the same taxa are colored by that taxa.
* Principle Coordinates of Analysis plot of weighted UniFrac distances calculated by QIIME
* Principle Coordinates of Analysis plot from QIIME output

####Files in [brazil_study_data](clr_vs_proportion_unifrac/brazil_study_test/brazil_study_data) folder
This study collected 16S rRNA sequencing data from vaginal swabs of women in Brazil.

#####[fasttree_all_seed_OTUs.tre](clr_vs_proportion_unifrac/brazil_study_test/brazil_study_data/fasttree_all_seed_OTUs.tre)
Phylogenetic tree of taxa, based on a sequence alignment.

#####[metadata_BVsamplesonly.txt](clr_vs_proportion_unifrac/brazil_study_test/brazil_study_data/metadata_BVsamplesonly.txt)
Metadata file.

#####[td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt](clr_vs_proportion_unifrac/brazil_study_test/brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt)
Counts per sample per operational taxonomic unit.

#####[weighted_unifrac_dm_from_qiime.txt](clr_vs_proportion_unifrac/brazil_study_test/brazil_study_data/weighted_unifrac_dm_from_qiime.txt)
UniFrac distance matrix output by QIIME.

#####[weighted_unifrac_pc_from_qiime.txt](clr_vs_proportion_unifrac/brazil_study_test/brazil_study_data/weighted_unifrac_pc_from_qiime.txt)
Principle Coordinates of Analysis vectors output by QIIME.
