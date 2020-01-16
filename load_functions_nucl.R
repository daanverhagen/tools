load_test_meth_data <- function(meth_data_file, sample, log_lik_ratio_threshold = 2.5) {
  meth_data <- read.table(meth_data_file, header = TRUE)
  names(meth_data)[c(1, 9, 10)] <- c("ref_seq", "group_size", "group_sequence")
  meth_data <- meth_data[, -8]
  meth_data$sample <- sample
  meth_data$hit_state <- 1*(meth_data$log_lik_ratio >= log_lik_ratio_threshold) - 1*(meth_data$log_lik_ratio <= -log_lik_ratio_threshold)  # 1: methylated, -1: not methylated, 0: ambiguous
  meth_data
}


load_chem_cleav_data <- function(rep_id1, rep_id2, chr, TSS, plot_range_min, plot_range_max, minus_strand) {
  ## load chemical cleavage data(dyad density/occupancy) from paper: Chereji et al. Genome Biology, 2018 19:19
  ## the data is in BigWig format and contains genomic locations and corresponding dyad densities/occupancies stores as scores
  ## 3 replicates
  chem_cleav_list <- list()
  dyad_df_list <- list()
  occu_df_list <- list()
  for(i in 1:length(rep_id1)) {
    dyad_bw <- file.path(paste0("../../../chem_cleavage_data/GSE97290_RAW/GSM256105", rep_id1[i], "_Dyads_H3_CC_rep_", rep_id2[i], ".bw"))
    occu_bw <- file.path(paste0("../../../chem_cleavage_data/GSE97290_RAW/GSM256105", rep_id1[i], "_Occupancy_H3_CC_rep_", rep_id2[i], ".bw"))
    if(minus_strand) {
      which <- GRanges(chr, IRanges(TSS-plot_range_max, TSS+plot_range_min))
      
      dyad = import(dyad_bw, which = which)
      dyad_score_vec <- dyad$score
      dyad_df_list[[i]] = data.frame(position = c(-plot_range_min:plot_range_max), dyad_density = dyad_score_vec[length(dyad_score_vec):1])
      
      occupancy = import(occu_bw, which = which)
      occu_score_vec <- occupancy$score
      occu_df_list[[i]] = data.frame(position = c(-plot_range_min:plot_range_max), occupancy = occu_score_vec[length(occu_score_vec):1])
    } else {
      which <- GRanges(chr, IRanges(TSS-plot_range_min, TSS+plot_range_max))
      
      dyad = import(dyad_bw, which = which)
      dyad_df_list[[i]] = data.frame(position = c(-plot_range_min:plot_range_max), dyad_density = dyad$score)
      
      occupancy = import(occu_bw, which = which)
      occu_df_list[[i]] = data.frame(position = c(-plot_range_min:plot_range_max), occupancy = occupancy$score)
    }
    # dyad_df_list[[i]]$rep <- as.character(paste0('chemi-cleav rep', rep_id2[i]))
    # occu_df_list[[i]]$rep <- as.character(paste0('chemi-cleav rep', rep_id2[i]))
    dyad_df_list[[i]]$rep <- as.character(paste0('rep', rep_id2[i]))
    occu_df_list[[i]]$rep <- as.character(paste0('rep', rep_id2[i]))
  }
  dyad_df <- bind_rows(dyad_df_list)
  occu_df <- bind_rows(occu_df_list)
  
  ## calculate and add the average of three replicates
  dyad_avg <- aggregate(dyad_df[, c('dyad_density')], by = list(dyad_df$position), FUN = mean, na.rm = TRUE)
  colnames(dyad_avg) <- c("position", "dyad_density")
  # dyad_avg <- data.frame(position = dyad_avg$position, dyad_density = dyad_avg$dyad_density, rep = 'chemi-cleav avg')
  dyad_avg <- data.frame(position = dyad_avg$position, dyad_density = dyad_avg$dyad_density, rep = 'avg')
  chem_cleav_list[[1]] <- rbind(dyad_df, dyad_avg)
  
  occu_avg <- aggregate(occu_df[, c('occupancy')], by = list(occu_df$position), FUN = mean, na.rm = TRUE)
  colnames(occu_avg) <- c("position", "occupancy")
  # occu_avg <- data.frame(position = occu_avg$position, occupancy = occu_avg$occupancy, rep = 'chemi-cleav avg')
  occu_avg <- data.frame(position = occu_avg$position, occupancy = occu_avg$occupancy, rep = 'avg')
  chem_cleav_list[[2]] <- rbind(occu_df, occu_avg)
  
  chem_cleav_list
}