# Daan Verhagen
# 20191023
# Functions for aligning different genomic locations


#### align_genes_Np1 function ####



align_genes_Np1 <- function(site_stats, max_promoter_length, max_ORF_length) {
  
  site_stats <- site_stats[!is.na(site_stats$hit_call_fraction_methylated), ]  
  # if hit_call_fraction_methylated is NaN, hits_called is zero, so we don't need to track hits_called into the sum and can omit the whole site
  
  ref_seq_labels <- sort(unique(as.character(site_stats$ref_seq)))
  
  # expected_labels <- c(paste0("chr0", 1:9), paste0("chr", 10:16))
  expected_labels <- sprintf("chr%02d" , c(1:16)) # a shorter way of doing the same thing as the line above
  if(length(ref_seq_labels)==0 || !all(diag(sapply(expected_labels, grepl, ref_seq_labels)))) { # checks whether all chromosomes are present
    warning("Check ref seq labels in input site_stats! Not all 16 chromosomes were found!\n")
    return(NULL)
  }
  
  
  # list of site_stats with ref_seq as index. ie makes a list entry for each chromosome
  site_stats_ref_seq_list <- list()  
  for (ref_seq_label in ref_seq_labels) {
    site_stats_ref_seq_list = c(site_stats_ref_seq_list, list(site_stats[site_stats$ref_seq == ref_seq_label, c('position', 'hit_call_fraction_methylated', 'hits_called')]))
  }
  names(site_stats_ref_seq_list) <- ref_seq_labels
  
  ### old dataset
  # gene_data = read.table(file = '../../external_data/ann_plus1_minus1_filtered.txt', sep = '\t', stringsAsFactors = FALSE, header = TRUE) # old Steinmetz dataset (Xu et al, Nature 2009)
  # colnames(gene_data)[c(2)] <- c("chromosome") # change some names to make compatible with old gene_data version
  ### new(er) dataset
  gene_data = read.delim("../../external_data/Yeast_annotation_Genome_Chereji_2018_Daan2.1.tsv") # Chereji/Henikoff dataset (Chereji et al Genome Biol. 2018)
  colnames(gene_data)[c(2,8,9)] <- c("chromosome", "plus1", "neg1") # change some names to make compatible with old gene_data version
  gene_data$chromosome <- str_remove_all(gene_data$chromosome, "chr") # remove "chr" to make downstream compatible with old gene_data version
  
  chr_names = data.frame(ref_seq_labels, as.character(as.roman(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  
  aligned_genes_list <- list()
  for (i in 1:nrow(gene_data)) {
    chr_temp <- chr_names[chr_names$gene_data == gene_data[i, "chromosome"], "alignment"]
    if (gene_data[i, "strand"] == '1') { # changed + to 1 to reflect new Chereji dataset
      genome_pos_min <- gene_data[i, "plus1"] - max_promoter_length
      genome_pos_max <- min(max_ORF_length + gene_data[i, "plus1"], gene_data[i, "end"])
      new_gene <- subset(site_stats_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_gene$position = new_gene$position - gene_data[i, "plus1"]
    } else {
      genome_pos_min <- max(gene_data[i, "start"], gene_data[i, "plus1"] - max_ORF_length)
      genome_pos_max <- gene_data[i, "plus1"] + max_promoter_length
      new_gene <- subset(site_stats_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_gene$position = gene_data[i, "plus1"] - new_gene$position
      new_gene = new_gene[nrow(new_gene):1,]  # revert position order (not needed)
    }
    colnames(new_gene) = c('position', 'methylation', 'hits_called')
    
    if(nrow(new_gene)>0) {
      new_gene$number = i
      aligned_genes_list[[i]] <- new_gene[, c('number', 'position', 'methylation', 'hits_called')]
    }
  }
  aligned_genes <- bind_rows(aligned_genes_list)  # makes one data.frame from the list
  
  aligned_genes
}


#### NP1 related calc_average_gene_ind_samples function ####

### no major changes were made here compared to original function in calc_functions.R 

calc_average_gene_ind_samples <- function(site_stats_chr, pos_binning, alignment_point = "Np1")  {
  samples <- levels(site_stats_chr$sample)
  aligned_genes_list <- vector(mode = "list", length = length(samples))
  average_gene_list <- vector(mode = "list", length = length(samples))
  
  for(i in 1:length(samples)) {
    site_stats_to_align <- subset(site_stats_chr, sample == samples[i])
    
    if(nrow(site_stats_to_align) > 0) {
      if(alignment_point == "Np1") {
        aligned_genes_list[[i]] <- align_genes_Np1(site_stats_to_align, 2000, 2000)
      } 
      
      aligned_genes_list[[i]]$position_binned = round(aligned_genes_list[[i]]$position/pos_binning,0)*pos_binning
      average_gene_list[[i]] = aggregate(aligned_genes_list[[i]][, c('methylation', 'hits_called')], list(aligned_genes_list[[i]]$position_binned), mean, na.rm = TRUE)
      average_gene_list[[i]]$hits_called_sum = aggregate(aligned_genes_list[[i]][, c('hits_called')], list(aligned_genes_list[[i]]$position_binned), sum, na.rm = TRUE)[, 2]
      
      colnames(average_gene_list[[i]]) = c('position', 'methylation', 'hits_called_mean', 'hits_called_sum')
      average_gene_list[[i]]$sample <- unique(site_stats_to_align[, "sample"])
    }
  }
  
  average_gene_ind_samples <- bind_rows(average_gene_list)  # makes one data.frame from the list
  
  average_gene_ind_samples
}

#### align_genes_GRF function ####


align_genes_GRF <- function(site_stats, max_promoter_length, max_ORF_length, GRF_Subset) {
  
  site_stats <- site_stats[!is.na(site_stats$hit_call_fraction_methylated), ]  
  # if hit_call_fraction_methylated is NaN, hits_called is zero, so we don't need to track hits_called into the sum and can omit the whole site
  
  ref_seq_labels <- sort(unique(as.character(site_stats$ref_seq)))
  
  # expected_labels <- c(paste0("chr0", 1:9), paste0("chr", 10:16))
  expected_labels <- sprintf("chr%02d" , c(1:16)) # a shorter way of doing the same thing as the line above
  if(length(ref_seq_labels)==0 || !all(diag(sapply(expected_labels, grepl, ref_seq_labels)))) { # checks whether all chromosomes are present
    warning("Check ref seq labels in input site_stats! Not all 16 chromosomes were found!\n")
    return(NULL)
  }
  
  # list of site_stats with ref_seq as index. ie makes a list entry for each chromosome
  site_stats_ref_seq_list <- list()  
  for (ref_seq_label in ref_seq_labels) {
    site_stats_ref_seq_list = c(site_stats_ref_seq_list, list(site_stats[site_stats$ref_seq == ref_seq_label, c('position', 'hit_call_fraction_methylated', 'hits_called')]))
  }
  names(site_stats_ref_seq_list) <- ref_seq_labels
  
  gene_data = read.csv("../../external_data/GRF_Coords.csv")
  gene_data = gene_data[!gene_data$seqnames == "chrMT",] # remove MT reads
  
  #optional subset for timing
  gene_data <- gene_data[gene_data$GRF %in% GRF_Subset,]
  colnames(gene_data)[c(1)] <- c("chromosome")
  gene_data$chromosome <- stringr::str_remove_all(gene_data$chromosome, "chr") # remove "chr" to make downstream compatible with old gene_data version
  
  chr_names = data.frame(ref_seq_labels, as.character(as.roman(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  # i = 5
  aligned_genes_list <- list()
  for (i in 1:nrow(gene_data)) {
    chr_temp <- chr_names[chr_names$gene_data == gene_data[i, "chromosome"], "alignment"]
    if (gene_data[i, "strand"] == '*') { # no strand information for GRFs 
      genome_pos_min <- gene_data[i, "start"] - max_promoter_length
      genome_pos_max <- max(max_ORF_length + gene_data[i, "start"], gene_data[i, "end"]) # changed min to max because using start as center genomic coordinate
      new_gene <- subset(site_stats_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max) # subsets coordinates that fall within min to max around genomic coordinate
      new_gene$position = new_gene$position - gene_data[i, "start"] # sets origin to 0 and upstream to negative, downstream to positive
    }
    
    colnames(new_gene) = c('position', 'methylation', 'hits_called')
    
    if(nrow(new_gene)>0) {
      new_gene$number = i # adds number of origin/gene to the created dataframe
      aligned_genes_list[[i]] <- new_gene[, c('number', 'position', 'methylation', 'hits_called')] # adds the gene to the large list
    }
  }
  aligned_genes = bind_rows(aligned_genes_list)  # makes one data.frame from the list
  
  aligned_genes
}


#### GRF related calc_average_gene_ind_samples_origins function ####

calc_average_gene_ind_samples_GRFs <- function(site_stats_chr, pos_binning, alignment_point = "GRF", GRF_subset)  {
  samples <- levels(site_stats_chr$sample)
  aligned_genes_list <- vector(mode = "list", length = length(samples)) # genes
  average_gene_list <- vector(mode = "list", length = length(samples)) # gene!
  
  for(i in 1:length(samples)) {
    site_stats_to_align <- subset(site_stats_chr, sample == samples[i])
    
    if(nrow(site_stats_to_align) > 0) {
      if(alignment_point == "GRF") {
        aligned_genes_list[[i]] <- align_genes_GRF(site_stats_to_align, 1000, 1000, GRF_subset)
      } 
      
      # bins the multiple positions, takes the mean for each bin
      aligned_genes_list[[i]]$position_binned = round(aligned_genes_list[[i]]$position/pos_binning, 0)*pos_binning
      average_gene_list[[i]] = aggregate(aligned_genes_list[[i]][, c('methylation', 'hits_called')], list(aligned_genes_list[[i]]$position_binned), mean, na.rm = TRUE)
      average_gene_list[[i]]$hits_called_sum = aggregate(aligned_genes_list[[i]][, c('hits_called')], list(aligned_genes_list[[i]]$position_binned), sum, na.rm = TRUE)[, 2]
      colnames(average_gene_list[[i]]) = c('position', 'methylation', 'hits_called_mean', 'hits_called_sum')
      average_gene_list[[i]]$sample <- unique(site_stats_to_align[, "sample"])
    }
  }
  
  average_gene_ind_samples_GRFs <<- bind_rows(average_gene_list)  # makes one data.frame from the list
  
  average_gene_ind_samples_GRFs
}


#### align_genes_origins function ####

align_genes_origins <- function(site_stats, max_promoter_length, max_ORF_length) {
  
  site_stats <- site_stats[!is.na(site_stats$hit_call_fraction_methylated), ]  
  # if hit_call_fraction_methylated is NaN, hits_called is zero, so we don't need to track hits_called into the sum and can omit the whole site
  
  ref_seq_labels <- sort(unique(as.character(site_stats$ref_seq)))
  
  # expected_labels <- c(paste0("chr0", 1:9), paste0("chr", 10:16))
  expected_labels <- sprintf("chr%02d" , c(1:16)) # a shorter way of doing the same thing as the line above
  if(length(ref_seq_labels)==0 || !all(diag(sapply(expected_labels, grepl, ref_seq_labels)))) { # checks whether all chromosomes are present
    warning("Check ref seq labels in input site_stats! Not all 16 chromosomes were found!\n")
    return(NULL)
  }
  
  
  # list of site_stats with ref_seq as index. ie makes a list entry for each chromosome
  site_stats_ref_seq_list <- list()  
  for (ref_seq_label in ref_seq_labels) {
    site_stats_ref_seq_list = c(site_stats_ref_seq_list, list(site_stats[site_stats$ref_seq == ref_seq_label, c('position', 'hit_call_fraction_methylated', 'hits_called')]))
  }
  names(site_stats_ref_seq_list) <- ref_seq_labels
  
  gene_data = read.csv("../../external_data/origins_soriano.csv")
  
  #optional subset for timing
  gene_data <- subset(gene_data, gene_data$Timing == "Early")
  
  colnames(gene_data)[c(1,3,4,5,6)] <- c("chromosome", "start", "end", "ori_center", "strand")
  gene_data$chromosome <- str_remove_all(gene_data$chromosome, "chr") # remove "chr" to make downstream compatible with old gene_data version
  
  chr_names = data.frame(ref_seq_labels, as.character(as.roman(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  aligned_genes_list <- list()
  for (i in 1:nrow(gene_data)) {
    chr_temp <- chr_names[chr_names$gene_data == gene_data[i, "chromosome"], "alignment"]
    if (gene_data[i, "strand"] == '+') { 
      genome_pos_min <- gene_data[i, "ori_center"] - max_promoter_length
      genome_pos_max <- min(max_ORF_length + gene_data[i, "ori_center"], gene_data[i, "end"])
      new_gene <- subset(site_stats_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max) # subsets coordinates that fall within min to max around genomic coordinate
      new_gene$position = new_gene$position - gene_data[i, "ori_center"] # sets origin to 0 and upstream to negative, downstream to positive
    } else {
      genome_pos_min <- max(gene_data[i, "start"], gene_data[i, "ori_center"] - max_ORF_length)
      genome_pos_max <- gene_data[i, "ori_center"] + max_promoter_length
      new_gene <- subset(site_stats_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_gene$position = gene_data[i, "ori_center"] - new_gene$position
      new_gene = new_gene[nrow(new_gene):1,]  # revert position order (not needed)
    }
    colnames(new_gene) = c('position', 'methylation', 'hits_called')
    
    if(nrow(new_gene)>0) {
      new_gene$number = i # adds number of origin/gene to the created dataframe
      aligned_genes_list[[i]] <- new_gene[, c('number', 'position', 'methylation', 'hits_called')] # adds the gene to the large list
    }
  }
  aligned_genes <- bind_rows(aligned_genes_list)  # makes one data.frame from the list
  
  aligned_genes
}

#### Origins related calc_average_gene_ind_samples_origins function ####

calc_average_gene_ind_samples_origins <- function(site_stats_chr, pos_binning, alignment_point = "Np1")  {
  samples <- levels(site_stats_chr$sample)
  aligned_genes_list <- vector(mode = "list", length = length(samples)) # genes
  average_gene_list <- vector(mode = "list", length = length(samples)) # gene!
  
  for(i in 1:length(samples)) {
    site_stats_to_align <- subset(site_stats_chr, sample == samples[i])
    
    if(nrow(site_stats_to_align) > 0) {
      if(alignment_point == "Np1") {
        aligned_genes_list[[i]] <- align_genes_origins(site_stats_to_align, 1000, 1000)
      } 
      
      # bins the multiple positions, takes the mean for each bin
      aligned_genes_list[[i]]$position_binned = round(aligned_genes_list[[i]]$position/pos_binning, 0)*pos_binning
      average_gene_list[[i]] = aggregate(aligned_genes_list[[i]][, c('methylation', 'hits_called')], list(aligned_genes_list[[i]]$position_binned), mean, na.rm = TRUE)
      average_gene_list[[i]]$hits_called_sum = aggregate(aligned_genes_list[[i]][, c('hits_called')], list(aligned_genes_list[[i]]$position_binned), sum, na.rm = TRUE)[, 2]
      colnames(average_gene_list[[i]]) = c('position', 'methylation', 'hits_called_mean', 'hits_called_sum')
      average_gene_list[[i]]$sample <- unique(site_stats_to_align[, "sample"])
    }
  }
  
  average_gene_ind_samples_origins <<- bind_rows(average_gene_list)  # makes one data.frame from the list
  
  average_gene_ind_samples_origins
}

