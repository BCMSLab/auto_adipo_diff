#!/bin/bash

# Define directory structure; for scripts
FIG_SRC=scripts/figures
TAB_SRC=scripts/tables

# Define directory structure; for output
MANUSCRIPT=manuscript
FIG_DIR=manuscript/figures
TAB_DIR=manuscript/tables

# Define directory structure; for intermediates
DATA=autoreg/data
LOG=log
LOG_FIG=log/figures
LOG_TAB=log/tables

# Define commands
RFIG=R CMD BATCH --vanilla $< $(LOG_FIG)/$(<F).Rout
RTAB=R CMD BATCH --vanilla $< $(LOG_TAB)/$(<F).Rout

# All
all: ## Run all parts of the makefile
all: figures tables clean

# Directories
dir_manuscript: ## Make manuscript directory tree
dir_manuscript:
	test ! -d $(MANUSCRIPT) && mkdir $(MANUSCRIPT) || exit 0
	test ! -d $(TAB_DIR) && mkdir $(TAB_DIR) || exit 0
	test ! -d $(FIG_DIR) && mkdir $(FIG_DIR) || exit 0
dir_logs: ## Make logs directory tree
dir_logs:
	test ! -d $(LOG) && mkdir $(LOG) || exit 0
	test ! -d $(LOG_FIG) && mkdir $(LOG_FIG) || exit 0
	test ! -d $(LOG_TAB) && mkdir $(LOG_TAB) || exit 0
	
figures: ## Generate the figures
figures: dir_manuscript \
	dir_logs \
	$(FIG_DIR)/Autophagy_markers.png \
	$(FIG_DIR)/Lipogenesis_markers.png \
	$(FIG_DIR)/Adipogenesis_markers.png \
	$(FIG_DIR)/mds_gene_group.png \
	$(FIG_DIR)/volcanos_all.png \
	$(FIG_DIR)/volcanos_autophagy.png \
	$(FIG_DIR)/volcanos_binding_all.png \
	$(FIG_DIR)/volcanos_binding_autophagy.png \
	$(FIG_DIR)/correlations_factor_stage.png \
	$(FIG_DIR)/factor_correlations_heatmap.png \
	$(FIG_DIR)/annotation_correlations_heatmap.png \
	$(FIG_DIR)/peak_overlap_heatmap.png \
	$(FIG_DIR)/factor_cofac_stage.png \
	$(FIG_DIR)/factor_cofac_annotation.png \
	$(FIG_DIR)/factor_hm_stage.png \
	$(FIG_DIR)/factor_hm_annotation.png \
	$(FIG_DIR)/signal_tissue_CEBPB_autophagy_genes.png \
	$(FIG_DIR)/signal_tissue_CEBPB_autophagy_tf.png \
	$(FIG_DIR)/signal_tissue_CEBPB_adipogenic_tf.png \
	$(FIG_DIR)/signal_tissue_PPARG_autophagy_genes.png \
	$(FIG_DIR)/signal_tissue_PPARG_autophagy_tf.png \
	$(FIG_DIR)/signal_tissue_PPARG_adipogenic_tf.png \
	$(FIG_DIR)/signal_CEBPB_autophagy_genes.png \
	$(FIG_DIR)/signal_CEBPB_autophagy_tf.png \
	$(FIG_DIR)/signal_CEBPB_adipogenic_tf.png \
	$(FIG_DIR)/signal_PPARG_autophagy_genes.png \
	$(FIG_DIR)/signal_PPARG_autophagy_tf.png \
	$(FIG_DIR)/signal_PPARG_adipogenic_tf.png \
	$(FIG_DIR)/profile_expression_autophagy_genes.png \
	$(FIG_DIR)/profile_expression_autophagy_tf.png \
	$(FIG_DIR)/profile_expression_adipogenic_tf.png \
	$(FIG_DIR)/profile_peaks_autophagy_genes.png \
	$(FIG_DIR)/profile_peaks_autophagy_tf.png \
	$(FIG_DIR)/profile_peaks_adipogenic_tf.png \
	$(FIG_DIR)/coexpres_adipogenic_autophagy_tf.png \
	$(FIG_DIR)/coexpres_adipogenic_autophagy_genes.png \
	$(FIG_DIR)/coexpres_adipogenic_adipogenic.png \
	$(FIG_DIR)/PPARG_occupancy_change.png \
	$(FIG_DIR)/CEBPB_occupancy_change.png \
	$(FIG_DIR)/volcanos_kd_cebpb.png \
	$(FIG_DIR)/volcanos_kd_pparg.png \
	$(FIG_DIR)/kd_heatmap_cebpb.png \
	$(FIG_DIR)/kd_heatmap_pparg.png

tables: ## Generate the tables
tables: dir_manuscript \
	dir_logs \
	$(TAB_DIR)/datasets.tex \
	$(TAB_DIR)/deg_markers.tex \
	$(TAB_DIR)/deg_fc.tex \
	$(TAB_DIR)/dep_fc_genes.tex \
	$(TAB_DIR)/dep_fc_tf_mod.tex \
	$(TAB_DIR)/variance_explained.tex \
	$(TAB_DIR)/kd_fc.tex \
	$(TAB_DIR)/deg_go.tex \
	$(TAB_DIR)/kd_go.tex
	
# Figures
$(FIG_DIR)/%_markers.png: $(FIG_SRC)/markers.R $(DATA)/gene_counts.rds
	$(RFIG)
$(FIG_DIR)/mds_gene_group.png: $(FIG_SRC)/mds_gene_group.R \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds \
	$(DATA)/transformed_counts.rds
	$(RFIG)
$(FIG_DIR)/mds_binding_factors.png: $(FIG_SRC)/mds_binding_factors.R \
	$(DATA)/binding_data.rds
	$(RFIG)
$(FIG_DIR)/volcanos_%.png: $(FIG_SRC)/volcanos.R \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds \
	$(DATA)/deg_res.rds
	$(RFIG)
$(FIG_DIR)/volcanos_binding_%.png: $(FIG_SRC)/volcanos_binding.R \
	$(DATA)/go_annotation.rds \
	$(DATA)/binding_data.rds \
	$(DATA)/dep_res.rds
	$(RFIG)
$(FIG_DIR)/correlations_factor_stage.png: $(FIG_SRC)/correlations_factor_stage.R \
	$(DATA)/deg_res.rds \
	$(DATA)/ddcor.rds
	$(RFIG)
$(FIG_DIR)/occupancy_factor_direction.png: $(FIG_SRC)/occupancy_factor_direction.R \
	$(DATA)/deg_res.rds \
	$(DATA)/factor_occupancy.rds
	$(RFIG)	
$(FIG_DIR)/modification_factor_direction.png: $(FIG_SRC)/modification_factor_direction.R \
	$(DATA)/deg_res.rds \
	$(DATA)/hm_occupancy.rds
	$(RFIG)
$(FIG_DIR)/transcription_factor_direction.png: $(FIG_SRC)/transcription_factor_direction.R \
	$(DATA)/deg_res.rds \
	$(DATA)/factor_occupancy.rds
	$(RFIG)	
#$(FIG_DIR)/string_chip_networks.png: $(FIG_SRC)/string_chip_networks.R \
	$(DATA)/interactions.rds \
	$(DATA)/factor_targets.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
#	$(RFIG)	
#$(FIG_DIR)/profiles_autophagy_genes.png: $(FIG_SRC)/profiles_autophagy_genes.R \
	$(DATA)/gene_counts.rds \
	$(DATA)/factor_occupancy.rds
#	$(RFIG)
#$(FIG_DIR)/profiles_autophagy_tfs.png: $(FIG_SRC)/profiles_autophagy_tfs.R \
	$(DATA)/gene_counts.rds \
	$(DATA)/factor_occupancy.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
#	$(RFIG)
#$(FIG_DIR)/profiles_adipogenic_tfs.png: $(FIG_SRC)/profiles_adipogenic_tfs.R \
	$(DATA)/gene_counts.rds \
	$(DATA)/factor_occupancy.rds
#	$(RFIG)
#$(FIG_DIR)/peaks_autophagy_genes.png: $(FIG_SRC)/peaks_autophagy_genes.R \
	$(DATA)/binding_data.rds
#	$(RFIG)
#$(FIG_DIR)/peaks_autophagy_tfs.png: $(FIG_SRC)/peaks_autophagy_tfs.R \
	$(DATA)/binding_data.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
#	$(RFIG)
#$(FIG_DIR)/peaks_adipogenic_tfs.png: $(FIG_SRC)/peaks_adipogenic_tfs.R \
	$(DATA)/binding_data.rds
#	$(RFIG)
$(FIG_DIR)/coexpres_adipogenic_autophagy_tf.png: $(FIG_SRC)/coexpres_adipogenic_autophagy_tf.R \
	$(DATA)/dgca.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RFIG)
$(FIG_DIR)/coexpres_adipogenic_autophagy_genes.png: $(FIG_SRC)/coexpres_adipogenic_autophagy_genes.R \
	$(DATA)/dgca.rds
	$(RFIG)
$(FIG_DIR)/coexpres_adipogenic_adipogenic.png: $(FIG_SRC)/coexpres_adipogenic_adipogenic.R \
	$(DATA)/dgca.rds
	$(RFIG)
$(FIG_DIR)/factor_correlations_heatmap.png: $(FIG_SRC)/factor_correlations_heatmap.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/go_annotation.rds
	$(RFIG)
$(FIG_DIR)/annotation_correlations_heatmap.png: $(FIG_SRC)/annotation_correlations_heatmap.R \
	$(DATA)/peak_counts.rds $(DATA)/go_annotation.rds
	$(RFIG)
$(FIG_DIR)/factor_%.png: $(FIG_SRC)/factor_cofac_hm.R \
	$(DATA)/peak_counts.rds $(DATA)/go_annotation.rds
	$(RFIG)
$(FIG_DIR)/tf_tracks.png: $(FIG_SRC)/tf_tracks.R \
	$(DATA)/coverage_tracks.rds
	$(RFIG)
$(FIG_DIR)/hm_tracks.png: $(FIG_SRC)/hm_tracks.R \
	$(DATA)/coverage_tracks.rds
	$(RFIG)
$(FIG_DIR)/pparg_signal.png: $(FIG_SRC)/pparg_signal.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/cebpb_signal.png: $(FIG_SRC)/cebpb_signal.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/med1_signal.png: $(FIG_SRC)/med1_signal.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/rxrg_signal.png: $(FIG_SRC)/rxrg_signal.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/pparg_signal_direct.png: $(FIG_SRC)/pparg_signal_direct.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/cebpb_signal_direct.png: $(FIG_SRC)/cebpb_signal_direct.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/cebpb_signal_tissue.png: $(FIG_SRC)/cebpb_signal_tissue.R \
	$(DATA)/data_tracks_tissue.rds
	$(RFIG)
$(FIG_DIR)/pparg_signal_tissue.png: $(FIG_SRC)/pparg_signal_tissue.R \
	$(DATA)/data_tracks_tissue.rds
	$(RFIG)
$(FIG_DIR)/transcription_expression.png: $(FIG_SRC)/transcription_expression.R \
	$(DATA)/transformed_counts.rds $(DATA)/go_annotation.rds \
	$(DATA)/factor_occupancy.rds
	$(RFIG)
$(FIG_DIR)/peak_overlap_heatmap.png: $(FIG_SRC)/peak_overlap_heatmap.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/peak_overlaps_all.rds
	$(RFIG)
$(FIG_DIR)/profile_%.png: $(FIG_SRC)/profiles.R \
	$(DATA)/gene_counts.rds \
	$(DATA)/factor_occupancy.rds \
	$(DATA)/binding_data.rds
	$(RFIG)
$(FIG_DIR)/signal_%.png: $(FIG_SRC)/signal_tracks.R \
	$(DATA)/data_tracks_tidy.rds \
	$(DATA)/data_tracks_tissue_tidy.rds
	$(RFIG)
$(FIG_DIR)/%_occupancy_change.png: $(FIG_SRC)/occupancy_change.R \
	$(DATA)/go_annotation.rds \
	$(DATA)/occupancy_res.rds \
	$(DATA)/factor_targets.rds
	$(RFIG)
$(FIG_DIR)/volcanos_kd_%.png: $(FIG_SRC)/kd_volcanos.R \
	$(DATA)/cebpb_kd_res.rds \
	$(DATA)/pparg_kd_res.rds
	$(RFIG)
$(FIG_DIR)/kd_heatmap_%.png: $(FIG_SRC)/kd_heatmap.R \
	$(DATA)/counts_cebpb_kd.rds \
	$(DATA)/arrays_pparg_kd.rds
	$(RFIG)	
	
# Tables
$(TAB_DIR)/datasets.tex: $(TAB_SRC)/datasets.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/gene_counts.rds
	$(RTAB)
$(TAB_DIR)/deg_markers.tex: $(TAB_SRC)/deg_markers.R \
	$(DATA)/deg_res.rds
	$(RTAB)
$(TAB_DIR)/deg_fc.tex: $(TAB_SRC)/deg_fc.R \
	$(DATA)/deg_res.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RTAB)
$(TAB_DIR)/dep_fc.tex: $(TAB_SRC)/dep_fc.R \
	$(DATA)/binding_data.rds \
	$(DATA)/dep_res.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RTAB)
$(TAB_DIR)/dep_fc_genes.tex: $(TAB_SRC)/dep_fc_genes.R \
	$(DATA)/binding_data.rds \
	$(DATA)/dep_res.rds
	$(RTAB)
$(TAB_DIR)/dep_fc_tf_mod.tex: $(TAB_SRC)/dep_fc_tf_mod.R \
	$(DATA)/binding_data.rds \
	$(DATA)/dep_res.rds
	$(RTAB)
$(TAB_DIR)/variance_explained.tex: $(TAB_SRC)/variance_explained.R \
	$(DATA)/transformed_counts.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RTAB)
$(TAB_DIR)/variance_explained_chip.tex: $(TAB_SRC)/variance_explained_chip.R $(DATA)/binding_data.rds
	$(RTAB)
$(TAB_DIR)/kd_fc.tex: $(TAB_SRC)/kd_fc.R \
	$(DATA)/cebpb_kd_res.rds \
	$(DATA)/pparg_kd_res.rds
	$(RTAB)
$(TAB_DIR)/deg_go.tex: $(TAB_SRC)/deg_go.R \
	$(DATA)/deg_go_res.rds
	$(RTAB)
$(TAB_DIR)/kd_go.tex: $(TAB_SRC)/kd_go.R \
	$(DATA)/kd_go_res.rds
	$(RTAB)

# Clean Up
.PHONY: clean
clean: ## Clean up
clean:
	rm -f *.pdf
	rm -f *.RData

# Source: https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.PHONY: help
help: ## Print the current page
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help
