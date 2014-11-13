CLRUniFrac
==========

Calculate UniFrac difference by weighting by centered log ratio transformed abundances, rather than proportional abundance. Preserves the distance property of UniFrac, while allowing for weighting.


####Files in main folder

#####CLRUniFrac.r
Generalized UniFrac script with weighting replaced by a centered log ratio transformed values based weighting (details below). Variance adjusted UniFrac feature is removed.

#####GUniFrac.r
Generalized UniFrac script ripped straight from the [GUniFrac R package][1]


####Files in test_script folder

#####pcoa_test_script.r
Test script that runs GUniFrac and CLRUniFrac on sample data.

#####metrics.r
Script that calculates overlap and average read count between all pairs of samples.

#####test_plots_with_brazil_study_data.pdf
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

####Files in brazil_study_data folder
This study collected 16S rRNA sequencing data from vaginal swabs of women in Brazil.

#####fasttree_all_seed_OTUs.tre
Phylogenetic tree of taxa, based on a sequence alignment.

#####metadata_BVsamplesonly.txt
Metadata file.

#####td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt
Counts per sample per operational taxonomic unit.

#####weighted_unifrac_dm_from_qiime.txt
UniFrac distance matrix output by QIIME.

#####weighted_unifrac_pc_from_qiime.txt
Principle Coordinates of Analysis vectors output by QIIME.

[1]: http://cran.r-project.org/web/packages/GUniFrac/index.html