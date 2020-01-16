calc_channel_stats <- function(read_stats) {
  channel_stats <- aggregate(cbind(read_stats[order(read_stats$channel), c("length_incl_hardclipping", "aligned_reference_length", "sites_total", "sites_methylated", "sites_ambiguous")], rep.int(1, nrow(read_stats))), list(read_stats$channel[order(read_stats$channel)]), sum)
  colnames(channel_stats) <- c("channel", "bases_incl_hardclipping", "aligned_bases", "sites_total", "sites_methylated", "sites_ambiguous", "reads")
  channel_stats$reads <- as.integer(channel_stats$reads)
  channel_stats$site_fraction_methylated <- channel_stats$sites_methylated / channel_stats$sites_total
  channel_stats$site_call_fraction_methylated <- channel_stats$sites_methylated / (channel_stats$sites_total - channel_stats$sites_ambiguous)
  channel_stats$site_fraction_ambiguous <- channel_stats$sites_ambiguous / channel_stats$sites_total
  channel_stats$alignment_clipping_loss <- channel_stats$bases_incl_hardclipping - channel_stats$aligned_bases
  channel_stats$log10_aligned_bases <- log10(channel_stats$aligned_bases)
  
  channel_stats
}


# calculates site stats from meth_data with restrictions on basecalling quality, read_length, ...
# alternative to load_site_stats, which uses the output of the calc_site_stats.py script (all reads are considered here)
calc_site_stats <- function(read_names, meth_data) {
  meth_data <- subset(meth_data, read_name %in% read_names)
  meth_data <- meth_data[, c("ref_seq", "start", "end", "group_size", "sample", "hit_state")]
  meth_data$site_ID <- paste(meth_data$sample, meth_data$ref_seq, meth_data$start, sep = "_")  # only hits with the same site_ID will be aggregated
  meth_data$hits_total <- 1
  meth_data$hits_methylated <- as.integer(meth_data$hit_state == 1)
  meth_data$hits_ambiguous <- as.integer(meth_data$hit_state == 0)
  
  site_stats <- meth_data[, c("site_ID", "ref_seq", "start", "end", "group_size", "sample")]
  site_stats <- site_stats[!duplicated(site_stats), ]
  
  site_stats_summed_values <- aggregate(meth_data[, c("hits_total", "hits_methylated", "hits_ambiguous")], list(meth_data$site_ID), sum)
  names(site_stats_summed_values)[1] <- "site_ID"
  
  site_stats <- merge(site_stats, site_stats_summed_values)
  
  site_stats$hits_called <- site_stats$hits_total - site_stats$hits_ambiguous
  site_stats$hit_fraction_methylated <- site_stats$hits_methylated / site_stats$hits_total
  site_stats$hit_call_fraction_methylated <- site_stats$hits_methylated / site_stats$hits_called
  site_stats$hit_fraction_ambiguous <- site_stats$hits_ambiguous / site_stats$hits_total
  
  site_stats$position <- (site_stats$start + site_stats$end) / 2
  
  site_stats
}

calc_sequence_stats <- function(site_stats) {
  sequence_stats <- aggregate(cbind(site_stats[c("hits_total", "hits_methylated", "hits_ambiguous", "hits_called")], 1), list(site_stats$group_sequence), sum)
  
  #ToDo
  #sequence_stats <- merge(sequence_stats, site_stats[, c("group_sequence", "group_size")])
  
  colnames(sequence_stats)[c(1, 6)] <- c("group_sequence", "sites")
  
  sequence_stats$hit_fraction_methylated <- sequence_stats$hits_methylated / sequence_stats$hits_total
  sequence_stats$hit_call_fraction_methylated <- sequence_stats$hits_methylated / sequence_stats$hits_called
  sequence_stats$hit_fraction_ambiguous <- sequence_stats$hits_ambiguous / sequence_stats$hits_total
  sequence_stats <- sequence_stats[order(sequence_stats$hit_fraction_methylated), ]
  
  sequence_stats
}


calc_coverage_for_sample <- function(read_infos, ref_seq_info_file, check_dist) {
  
  ref_seq_infos <- read.table(ref_seq_info_file, sep = "\t", header = TRUE, row.names = 1)
  ref_seqs_present <- unique(as.character(read_infos$ref_seq_detailed))
  sample <- unique(as.character(read_infos$sample))
  if( length(sample) > 1) {
    stop("Input has more than one sample!\n")
  }
  coverage_list <- list()

  # prepare coverage matrix
  for(ref_seq_detailed in ref_seqs_present) {
    ref_seq_length <- ref_seq_infos[ref_seq_detailed, "length"]
    if(check_dist/2 <= ref_seq_length) {
      positions <- seq(round(check_dist/2), ref_seq_length, check_dist)
      coverage_list[[ref_seq_detailed]] <- data.frame(positions, read_count = 0)  # data.frames are slow! list would be better!
      colnames(coverage_list[[ref_seq_detailed]]) <- c("position", "read_count")
    } else {
      warning("check_dist/2 greater than reference sequence length!")
    }
  }
  
  # go through all reads
  for(i in 1:nrow(read_infos)) {
    
    ref_seq_detailed <- as.character(read_infos[i, "ref_seq_detailed"])
    ref_seq_length <- ref_seq_infos[ref_seq_detailed, "length"]
    read_length <- read_infos[i, "aligned_reference_length"]
    read_start <- read_infos[i, "reference_start"]
    positions <- coverage_list[[ref_seq_detailed]][, "position"]
    if(read_start <= max(positions)) {
      first_checkpoint <- min(positions[positions >= read_start])
      read_end <- read_start + read_length
      if(read_end >= first_checkpoint) {
        positions_hit <- seq(first_checkpoint, min(read_end, ref_seq_length), check_dist)
        indeces_hit <- 1+((min(positions_hit)-positions[1]):(max(positions_hit)-positions[1]))/check_dist
        coverage_list[[ref_seq_detailed]][indeces_hit, "read_count"] <- 1 + coverage_list[[ref_seq_detailed]][indeces_hit, "read_count"]
        coverage_list[[ref_seq_detailed]]$ref_seq_detailed <- ref_seq_detailed
      }
    }
  }
  
  # bind matrix entries to data frame
  coverage_df <- bind_rows(coverage_list)
  coverage_df$ref_seq <- sub(".*\\|", "", coverage_df$ref_seq_detailed)
  coverage_df$sample <- sample
  coverage_df
}


calc_coverage <- function(read_infos, ref_seq_info_file, check_dist = 100) {
  coverage_list <- list()
  for(sample in unique(as.character(read_infos$sample))) {
    coverage_list[[sample]] <- calc_coverage_for_sample(read_infos[read_infos$sample == sample, ], ref_seq_info_file, check_dist)
  }
  coverage_df <- bind_rows(coverage_list)
  coverage_df$sample <- factor(coverage_df$sample, levels(read_infos$sample))
  coverage_df$ref_seq_detailed <- factor(coverage_df$ref_seq_detailed, levels(read_infos$ref_seq_detailed))
  coverage_df$ref_seq <- factor(coverage_df$ref_seq, levels(read_infos$ref_seq))
  coverage_df
}


align_genes_TSS <- function(site_stats, max_promoter_length, max_ORF_length) {
  
  site_stats <- site_stats[!is.na(site_stats$hit_call_fraction_methylated), ]  
  # if hit_call_fraction_methylated is NaN, hits_called is zero, so we don't need to track hits_called into the sum and can omit the whole site
  
  ref_seq_labels <- sort(unique(as.character(site_stats$ref_seq)))
  
  expected_labels <- c(paste0("chr0", 1:9), paste0("chr", 10:16))
  if(!all(diag(sapply(expected_labels, grepl, ref_seq_labels)))) stop("Check ref seq labels in input site_stats! Not all 16 chromosomes were found!")
  
  site_stats_ref_seq_list <- list()  # list of site_stats with ref_seq as index
  for (ref_seq_label in ref_seq_labels) {
    site_stats_ref_seq_list = c(site_stats_ref_seq_list, list(site_stats[site_stats$ref_seq == ref_seq_label, c('position', 'hit_call_fraction_methylated', 'hits_called')]))
  }
  names(site_stats_ref_seq_list) <- ref_seq_labels
  
  
  gene_data = read.table(file = '../../external_data/Xu_2009_ORF-Ts_V64.bed', sep = '\t', stringsAsFactors = FALSE)
  colnames(gene_data) <- c("chromosome", "start", "end", "name", "score", "strand")
  chr_names = data.frame(ref_seq_labels, unique(gene_data$chromosome), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  aligned_genes_list <- list()
  for (i in 1:nrow(gene_data)) {
    chr_temp <- chr_names[chr_names$gene_data == gene_data[i, "chromosome"], "alignment"]
    if (gene_data[i, "strand"] == '+') {
      genome_pos_min <- gene_data[i, "start"] - max_promoter_length
      genome_pos_max <- min(max_ORF_length + gene_data[i, "start"], gene_data[i, "end"])
      new_gene <- subset(site_stats_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_gene$position = new_gene$position - gene_data[i, "start"]
    } else {
      genome_pos_min <- max(gene_data[i, "start"], gene_data[i, "end"] - max_ORF_length)
      genome_pos_max <- gene_data[i, "end"] + max_promoter_length
      new_gene <- subset(site_stats_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_gene$position = gene_data[i, "end"] - new_gene$position
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

# align_genes_Np1(site_stats = site_stats_chr, max_promoter_length = 2000, max_ORF_length = 2000)
align_genes_Np1 <- function(site_stats, max_promoter_length, max_ORF_length) {
  
  site_stats <- site_stats[!is.na(site_stats$hit_call_fraction_methylated), ]  
  # if hit_call_fraction_methylated is NaN, hits_called is zero, so we don't need to track hits_called into the sum and can omit the whole site
  
  ref_seq_labels <- sort(unique(as.character(site_stats$ref_seq)))
  
  expected_labels <- c(paste0("chr0", 1:9), paste0("chr", 10:16))
  if(length(ref_seq_labels)==0 || !all(diag(sapply(expected_labels, grepl, ref_seq_labels)))) {
    warning("Check ref seq labels in input site_stats! Not all 16 chromosomes were found!\n")
    return(NULL)
  }
  
  site_stats_ref_seq_list <- list()  # list of site_stats with ref_seq as index
  for (ref_seq_label in ref_seq_labels) {
    site_stats_ref_seq_list = c(site_stats_ref_seq_list, list(site_stats[site_stats$ref_seq == ref_seq_label, c('position', 'hit_call_fraction_methylated', 'hits_called')]))
  }
  names(site_stats_ref_seq_list) <- ref_seq_labels
  
  gene_data = read.table(file = '../../external_data/ann_plus1_minus1_filtered.txt', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  colnames(gene_data)[2] <- c("chromosome")
  chr_names = data.frame(ref_seq_labels, as.character(as.roman(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  aligned_genes_list <- list()
  for (i in 1:nrow(gene_data)) {
    chr_temp <- chr_names[chr_names$gene_data == gene_data[i, "chromosome"], "alignment"]
    if (gene_data[i, "strand"] == '+') {
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


align_meth_data_Np1 <- function(meth_data, chr_ref_seqs, max_promoter_length, max_ORF_length) {
  
  meth_data_ref_seq_list <- list()  # list of meth_data with ref_seq as index
  for (ref_seq_label in chr_ref_seqs) {
    meth_data_ref_seq_list = c(meth_data_ref_seq_list, list(meth_data[meth_data$ref_seq == ref_seq_label, ]))
  }
  names(meth_data_ref_seq_list) <- chr_ref_seqs
  
  # gene_data = load_scer_gene_data() #old line of code
  gene_data = read.delim("../../external_data/Yeast_annotation_Genome_Chereji_2018_Daan2.tsv")
  chr_names = data.frame(chr_ref_seqs, as.character(as.roman(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  aligned_meth_data_list <- list()
  for (i in 1:nrow(gene_data)) {
    chr_temp <- chr_names[paste0('chr', chr_names$gene_data) == gene_data[i, "chr"], "alignment"] #changed "chromosome" to "chr" perhaps it was called this in an older gene_data list
    
    if (gene_data[i, "strand"] == '1') {
      genome_pos_min <- gene_data[i, "plus_1_nucl"] - max_promoter_length
      genome_pos_max <- min(gene_data[i, "plus_1_nucl"] + max_ORF_length, gene_data[i, "TTS"])
      
      new_meth_data <- subset(meth_data_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_meth_data$position = new_meth_data$position - gene_data[i, "plus_1_nucl"]
      new_meth_data$start = new_meth_data$start - gene_data[i, "plus_1_nucl"]
      new_meth_data$end = new_meth_data$end - gene_data[i, "plus_1_nucl"]
      
    } else {
      genome_pos_min <- max(gene_data[i, "TTS"], gene_data[i, "plus_1_nucl"] - max_ORF_length)
      genome_pos_max <- gene_data[i, "plus_1_nucl"] + max_promoter_length
      
      new_meth_data <- subset(meth_data_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_meth_data$position = gene_data[i, "plus_1_nucl"] - new_meth_data$position
      new_meth_data_start <- new_meth_data$start  # swapping start and end on '-' strand, so that start <= end again
      new_meth_data$start = gene_data[i, "plus_1_nucl"] - new_meth_data$end
      new_meth_data$end = gene_data[i, "plus_1_nucl"] - new_meth_data_start

      if(nrow(new_meth_data)>0) {
        new_meth_data = new_meth_data[nrow(new_meth_data):1, ]  # revert position order (not needed)
      }
    }
    
    if(nrow(new_meth_data)>0) {
      new_meth_data$read_name <- paste0(new_meth_data$read_name, "_gene_", i)
      new_meth_data$aligned_gene_number = i
      new_meth_data$aligned_gene_name = gene_data[i, "locus_name"]
      new_meth_data$aligned_gene_strand = gene_data[i, "strand"]
      
      aligned_meth_data_list[[i]] <- new_meth_data
    }
  }
  
  aligned_meth_data <- bind_rows(aligned_meth_data_list)  # makes one data.frame from the list (should be faster than do.call(rbind, ...), but is still the slowest step)
  aligned_meth_data$read_name <- as.factor(aligned_meth_data$read_name)
  
  aligned_meth_data
}

align_meth_data_abf1 <- function(meth_data, chr_ref_seqs, max_promoter_length, max_ORF_length) {
  
  meth_data_ref_seq_list <- list()  # list of meth_data with ref_seq as index
  for (ref_seq_label in chr_ref_seqs) {
    meth_data_ref_seq_list = c(meth_data_ref_seq_list, list(meth_data[meth_data$ref_seq == ref_seq_label, ]))
  }
  names(meth_data_ref_seq_list) <- chr_ref_seqs
  
  # gene_data = load_scer_gene_data() #old line of code
  gene_data = read.csv("../../external_data/abf1_coords.csv")
  chr_names = data.frame(chr_ref_seqs, as.character(as.roman(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  
  aligned_meth_data_list <- list()
  for (i in 1:nrow(gene_data)) {
    
    chr_temp <- chr_names[paste0('chr', chr_names$gene_data) == gene_data[i, "seqnames"], "alignment"] #changed "chromosome" to "chr" perhaps it was called this in an older gene_data list
    
    if (gene_data[i, "strand"] == '*') {
      genome_pos_min <- gene_data[i, "start"] - max_promoter_length
      genome_pos_max <- max(gene_data[i, "end"] + max_ORF_length, gene_data[i, "start"])
      
      new_meth_data <- subset(meth_data_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_meth_data$position = new_meth_data$position - gene_data[i, "start"]
      new_meth_data$start = new_meth_data$start - gene_data[i, "start"]
      new_meth_data$end = new_meth_data$end - gene_data[i, "end"]
      
    } 
    
    if(nrow(new_meth_data)>0) {
      new_meth_data = new_meth_data[nrow(new_meth_data):1, ]  # revert position order (not needed)
    }
    
    
    if(nrow(new_meth_data)>0) {
      new_meth_data$read_name <- paste0(new_meth_data$read_name, "_gene_", i)
      new_meth_data$aligned_gene_number = i
      new_meth_data$aligned_gene_name = gene_data[i, "locus_name"]
      new_meth_data$aligned_gene_strand = gene_data[i, "strand"]
      
      aligned_meth_data_list[[i]] <- new_meth_data
    }
  }
  
  aligned_meth_data <- bind_rows(aligned_meth_data_list)  # makes one data.frame from the list (should be faster than do.call(rbind, ...), but is still the slowest step)
  aligned_meth_data$read_name <- as.factor(aligned_meth_data$read_name)
  
  aligned_meth_data
}

align_meth_data_rap1 <- function(meth_data, chr_ref_seqs, max_promoter_length, max_ORF_length) {
  
  meth_data_ref_seq_list <- list()  # list of meth_data with ref_seq as index
  for (ref_seq_label in chr_ref_seqs) {
    meth_data_ref_seq_list = c(meth_data_ref_seq_list, list(meth_data[meth_data$ref_seq == ref_seq_label, ]))
  }
  names(meth_data_ref_seq_list) <- chr_ref_seqs
  
  # gene_data = load_scer_gene_data() #old line of code
  gene_data = read.csv("../../external_data/rap1_coords.csv")
  chr_names = data.frame(chr_ref_seqs, as.character(as.roman(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  
  aligned_meth_data_list <- list()
  for (i in 1:nrow(gene_data)) {
    
    chr_temp <- chr_names[paste0('chr', chr_names$gene_data) == gene_data[i, "seqnames"], "alignment"] #changed "chromosome" to "chr" perhaps it was called this in an older gene_data list
    
    if (gene_data[i, "strand"] == '*') {
      genome_pos_min <- gene_data[i, "start"] - max_promoter_length
      genome_pos_max <- max(gene_data[i, "end"] + max_ORF_length, gene_data[i, "start"])
      
      new_meth_data <- subset(meth_data_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_meth_data$position = new_meth_data$position - gene_data[i, "start"]
      new_meth_data$start = new_meth_data$start - gene_data[i, "start"]
      new_meth_data$end = new_meth_data$end - gene_data[i, "end"]
      
    } 
    
    if(nrow(new_meth_data)>0) {
      new_meth_data = new_meth_data[nrow(new_meth_data):1, ]  # revert position order (not needed)
    }
    
    
    if(nrow(new_meth_data)>0) {
      new_meth_data$read_name <- paste0(new_meth_data$read_name, "_gene_", i)
      new_meth_data$aligned_gene_number = i
      new_meth_data$aligned_gene_name = gene_data[i, "locus_name"]
      new_meth_data$aligned_gene_strand = gene_data[i, "strand"]
      
      aligned_meth_data_list[[i]] <- new_meth_data
    }
  }
  
  aligned_meth_data <- bind_rows(aligned_meth_data_list)  # makes one data.frame from the list (should be faster than do.call(rbind, ...), but is still the slowest step)
  aligned_meth_data$read_name <- as.factor(aligned_meth_data$read_name)
  
  aligned_meth_data
}


align_meth_data_reb1 <- function(meth_data, chr_ref_seqs, max_promoter_length, max_ORF_length) {
  
  meth_data_ref_seq_list <- list()  # list of meth_data with ref_seq as index
  for (ref_seq_label in chr_ref_seqs) {
    meth_data_ref_seq_list = c(meth_data_ref_seq_list, list(meth_data[meth_data$ref_seq == ref_seq_label, ]))
  }
  names(meth_data_ref_seq_list) <- chr_ref_seqs
  
  # gene_data = load_scer_gene_data() #old line of code
  gene_data = read.csv("../../external_data/reb1_coords.csv")
  chr_names = data.frame(chr_ref_seqs, as.character(as.roman(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  
  aligned_meth_data_list <- list()
  for (i in 1:nrow(gene_data)) {
    
    chr_temp <- chr_names[paste0('chr', chr_names$gene_data) == gene_data[i, "seqnames"], "alignment"] #changed "chromosome" to "chr" perhaps it was called this in an older gene_data list
    
    if (gene_data[i, "strand"] == '*') {
      genome_pos_min <- gene_data[i, "start"] - max_promoter_length
      genome_pos_max <- max(gene_data[i, "end"] + max_ORF_length, gene_data[i, "start"])
      
      new_meth_data <- subset(meth_data_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_meth_data$position = new_meth_data$position - gene_data[i, "start"]
      new_meth_data$start = new_meth_data$start - gene_data[i, "start"]
      new_meth_data$end = new_meth_data$end - gene_data[i, "end"]
      
    } 
    
    if(nrow(new_meth_data)>0) {
      new_meth_data = new_meth_data[nrow(new_meth_data):1, ]  # revert position order (not needed)
    }
    
    
    if(nrow(new_meth_data)>0) {
      new_meth_data$read_name <- paste0(new_meth_data$read_name, "_gene_", i)
      new_meth_data$aligned_gene_number = i
      new_meth_data$aligned_gene_name = gene_data[i, "locus_name"]
      new_meth_data$aligned_gene_strand = gene_data[i, "strand"]
      
      aligned_meth_data_list[[i]] <- new_meth_data
    }
  }
  
  aligned_meth_data <- bind_rows(aligned_meth_data_list)  # makes one data.frame from the list (should be faster than do.call(rbind, ...), but is still the slowest step)
  aligned_meth_data$read_name <- as.factor(aligned_meth_data$read_name)
  
  aligned_meth_data
}

align_meth_data_txn <- function(meth_data, factorname, chr_ref_seqs, max_promoter_length, max_ORF_length) {
  
  meth_data_ref_seq_list <- list()  # list of meth_data with ref_seq as index
  for (ref_seq_label in chr_ref_seqs) {
    meth_data_ref_seq_list = c(meth_data_ref_seq_list, list(meth_data[meth_data$ref_seq == ref_seq_label, ]))
  }
  names(meth_data_ref_seq_list) <- chr_ref_seqs
  
  # gene_data = load_scer_gene_data() #old line of code
  gene_data = read.csv(paste0("../../external_data/", factorname, "_coords.csv"))
  chr_names = data.frame(chr_ref_seqs, as.character(as.roman(1:16)), stringsAsFactors = FALSE)
  colnames(chr_names) <- c('alignment', 'gene_data')
  
  
  aligned_meth_data_list <- list()
  for (i in 1:nrow(gene_data)) {
    
    chr_temp <- chr_names[paste0('chr', chr_names$gene_data) == gene_data[i, "seqnames"], "alignment"] #changed "chromosome" to "chr" perhaps it was called this in an older gene_data list
    
    if (gene_data[i, "strand"] == '*') {
      genome_pos_min <- gene_data[i, "start"] - max_promoter_length
      genome_pos_max <- max(gene_data[i, "end"] + max_ORF_length, gene_data[i, "start"])
      
      new_meth_data <- subset(meth_data_ref_seq_list[[chr_temp]], position >= genome_pos_min & position <= genome_pos_max)
      new_meth_data$position = new_meth_data$position - gene_data[i, "start"]
      new_meth_data$start = new_meth_data$start - gene_data[i, "start"]
      new_meth_data$end = new_meth_data$end - gene_data[i, "end"]
      
    } 
    
    if(nrow(new_meth_data)>0) {
      new_meth_data = new_meth_data[nrow(new_meth_data):1, ]  # revert position order (not needed)
    }
    
    
    if(nrow(new_meth_data)>0) {
      new_meth_data$read_name <- paste0(new_meth_data$read_name, "_gene_", i)
      new_meth_data$aligned_gene_number = i
      new_meth_data$aligned_gene_name = gene_data[i, "locus_name"]
      new_meth_data$aligned_gene_strand = gene_data[i, "strand"]
      
      aligned_meth_data_list[[i]] <- new_meth_data
    }
  }
  
  aligned_meth_data <- bind_rows(aligned_meth_data_list)  # makes one data.frame from the list (should be faster than do.call(rbind, ...), but is still the slowest step)
  aligned_meth_data$read_name <- as.factor(aligned_meth_data$read_name)
  
  aligned_meth_data
}  


filter_aligned_meth_data <- function(aligned_meth_data, pos_min = 0, pos_max = 0) {
  # removes reads that do not cover positions pos_min to pos_max
  meth_data_starts <- aggregate(aligned_meth_data$start, list(as.character(aligned_meth_data$read_name)), min)
  colnames(meth_data_starts) <- c("read_name", "start")
  meth_data_ends <- aggregate(aligned_meth_data$end, list(as.character(aligned_meth_data$read_name)), max)
  colnames(meth_data_ends) <- c("read_name", "end")
  meth_data_starts_ends <- merge(meth_data_starts, meth_data_ends)
  read_names_filtered <- subset(meth_data_starts_ends, start < pos_min & end > pos_max)$read_name
  aligned_meth_data_filtered <- aligned_meth_data[aligned_meth_data$read_name %in% read_names_filtered, ]
  
  aligned_meth_data_filtered
}


pool_site_stats <- function(site_stats, sample_pools, pool_labels) {
  # sample_pools: list of vectors of samples, list entries are pooled together
  # pool_labels: vector of pool names
  
  if(missing(sample_pools)) {
    sample_pools <- list(unique(as.character(site_stats$sample)))
    warning("pool_site_stats: pooling all barcodes.")
  }
  
  if(missing(pool_labels)) {
    pool_labels <- character(0)
    for(i in 1:length(sample_pools)) {
      pool_labels[i] <- paste0("samples ", paste(sample_pools[[i]], collapse = ", "))
    }
  }
  
  site_stats_pooled_list <- vector("list", length = length(sample_pools))
  
  for(i in 1:length(sample_pools)) {
    site_stats_to_pool <- subset(site_stats, sample %in% sample_pools[[i]])
    
    if(nrow(site_stats_to_pool) > 0) {
      
      site_stats_to_pool$ID <- paste(site_stats_to_pool$ref_seq, site_stats_to_pool$start, sep = "_")
      
      site_stats_summed_values <- aggregate(site_stats_to_pool[, c("hits_total", "hits_methylated", "hits_ambiguous", "hits_called")], list(site_stats_to_pool$ID), sum)
      names(site_stats_summed_values)[1] <- "ID"
      site_stats_pooled <- site_stats_to_pool[, c("ID", "ref_seq", "start", "end", "group_size", "position")]
      site_stats_pooled <- site_stats_pooled[!duplicated(site_stats_pooled), ]
      site_stats_pooled <- merge(site_stats_pooled, site_stats_summed_values)
      site_stats_pooled$pool_label <- pool_labels[i]
      
      site_stats_pooled$hit_fraction_methylated <- site_stats_pooled$hits_methylated / site_stats_pooled$hits_total
      site_stats_pooled$hit_call_fraction_methylated <- site_stats_pooled$hits_methylated / site_stats_pooled$hits_called
      site_stats_pooled$hit_fraction_ambiguous <- site_stats_pooled$hits_ambiguous / site_stats_pooled$hits_total
      
      site_stats_pooled_list[[i]] <- site_stats_pooled
    } else {
      warning(paste0("pool_site_stats: no data to pool for samples ", sample_pools[[i]], "!\n"))
    }
  }
  
  site_stats_pooled <- bind_rows(site_stats_pooled_list)
  
  if(!is.null(site_stats_pooled)) {
    site_stats_pooled$pool_label <- factor(site_stats_pooled$pool_label, pool_labels)
    return(site_stats_pooled)
  } else {
    return(data.frame("pool_label" = NA, "position" = NA, "hit_call_fraction_methylated" = NA))
  }
}

add_pool_label_to_read_stats <- function(read_stats, sample_pools, pool_labels) {
  # sample_pools: list of vectors of barcodes, list entries are pooled together
  # pool_labels: vector of pool names
  
  if(missing(sample_pools)) {
    sample_pools <- list(unique(as.character(site_stats$sample)))
    warning("pool_site_stats: pooling all barcodes.")
  }
  
  if(missing(pool_labels)) {
    pool_labels <- character(0)
    for(i in 1:length(sample_pools)) {
      pool_labels[i] <- paste0("samples ", paste(sample_pools[[i]], collapse = ", "))
    }
  }
  
  read_stats$pool_label <- NA
  
  for(i in 1:length(sample_pools)) {
    read_stats[read_stats$sample %in% sample_pools[[i]], "pool_label"] <- pool_labels[[i]]
  }
  
  read_stats$pool_label <- factor(read_stats$pool_label, levels = pool_labels)
  
  return(read_stats)
}


add_pool_label_to_meth_data <- function(meth_data, sample_pools, pool_labels) {
  # sample_pools: list of vectors of barcodes, list entries are pooled together
  # pool_labels: vector of pool names
  
  if(missing(sample_pools)) {
    sample_pools <- list(unique(as.character(site_stats$sample)))
    warning("pool_site_stats: pooling all barcodes.")
  }
  
  if(missing(pool_labels)) {
    pool_labels <- character(0)
    for(i in 1:length(sample_pools)) {
      pool_labels[i] <- paste0("samples ", paste(sample_pools[[i]], collapse = ", "))
    }
  }
  
  meth_data$pool_label <- NA
  
  for(i in 1:length(sample_pools)) {
    meth_data[meth_data$sample %in% sample_pools[[i]], "pool_label"] <- pool_labels[[i]]
  }
  
  meth_data$pool_label <- factor(meth_data$pool_label, levels = pool_labels)
  
  return(meth_data)
}


calc_average_gene_ind_samples <- function(site_stats_chr, pos_binning, alignment_point = "TSS")  {
  samples <- levels(site_stats_chr$sample)
  aligned_genes_list <- vector(mode = "list", length = length(samples))
  average_gene_list <- vector(mode = "list", length = length(samples))
  
  for(i in 1:length(samples)) {
    site_stats_to_align <- subset(site_stats_chr, sample == samples[i])
    
    if(nrow(site_stats_to_align) > 0) {
      if(alignment_point == "Np1") {
        aligned_genes_list[[i]] <- align_genes_Np1(site_stats_to_align, 1000, 2000)
      } else if(alignment_point == "TSS") {
        aligned_genes_list[[i]] <- align_genes_TSS(site_stats_to_align, 1000, 2000)
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


calc_average_gene_pooled_samples <- function(site_stats_chr_pooled, sample_pools, pool_labels, pos_binning, alignment_point = "TSS") {
  aligned_genes_list <- vector(mode = "list", length = length(sample_pools))
  average_gene_list <- vector(mode = "list", length = length(sample_pools))
  
  for(i in 1:length(sample_pools)) {
    if(alignment_point == "Np1") {
      aligned_genes_list[[i]] <- align_genes_Np1(subset(site_stats_chr_pooled, pool_label == pool_labels[i]), 1000, 2000)
    } else if(alignment_point == "TSS") {
      aligned_genes_list[[i]] <- align_genes_TSS(subset(site_stats_chr_pooled, pool_label == pool_labels[i]), 1000, 2000)
    }
    
    aligned_genes_list[[i]]$position_binned = round(aligned_genes_list[[i]]$position/pos_binning,0)*pos_binning
    average_gene_list[[i]] = aggregate(aligned_genes_list[[i]][, c('methylation')], list(aligned_genes_list[[i]]$position_binned), mean, na.rm = TRUE)
    average_gene_list[[i]]$hits_called_sum = aggregate(aligned_genes_list[[i]][, c('hits_called')], list(aligned_genes_list[[i]]$position_binned), sum, na.rm = TRUE)[, 2]
    
    colnames(average_gene_list[[i]]) = c('position', 'methylation', 'hits_called_sum')
    average_gene_list[[i]]$pool_label <- pool_labels[i]
  }
  
  average_gene_pooled_samples <- bind_rows(average_gene_list)  # makes one data.frame from the list
  average_gene_pooled_samples$pool_label <- factor(average_gene_pooled_samples$pool_label, pool_labels)
  
  average_gene_pooled_samples
}


calc_read_perfection <- function(read_name, meth_data_to_plot) {
  agreement_vector <- meth_data_to_plot[meth_data_to_plot$read_name == read_name, ]$hit_state == meth_data_to_plot[meth_data_to_plot$read_name == read_name, ]$perfect_hit_state
  length(agreement_vector[agreement_vector]) / length(meth_data_to_plot[meth_data_to_plot$read_name == read_name, ]$hit_state)
}


calc_meth_data_matrix <- function(meth_data, reads_to_plot, pos_binning, perfect_read) {
  reads_to_plot$read_number <- 1:nrow(reads_to_plot)
  
  meth_data_to_plot <- subset(meth_data, read_name %in% reads_to_plot$read_name)
  meth_data_to_plot <- merge(meth_data_to_plot, reads_to_plot[, c('read_name','read_number')])
  meth_data_to_plot <- meth_data_to_plot[order(meth_data_to_plot$read_number, meth_data_to_plot$start), ]
  if(!missing(perfect_read)) {
    meth_data_to_plot$perfect_hit_state <- sapply(meth_data_to_plot$position_binned, function(pos) perfect_read$hit_state[perfect_read$binned_position == pos])
  }
  position_offset <- min(meth_data_to_plot$position_binned) - pos_binning
  meth_data_matrix <- matrix(NaN, nrow = nrow(reads_to_plot), ncol = (max(meth_data_to_plot$position_binned) - position_offset) / pos_binning)
  for(i in 1:nrow(meth_data_matrix)) {
    values <- rep(NaN, length=ncol(meth_data_matrix))
    sites <- meth_data_to_plot$read_name == as.character(reads_to_plot$read_name[i])
    positions <- (meth_data_to_plot$position_binned[sites] - position_offset) / pos_binning
    values[positions] <- meth_data_to_plot$hit_state[sites]
    meth_data_matrix[i, ] <- values
  }
  
  reads_to_plot$agreement_with_perfect_read <- sapply(reads_to_plot$read_name, calc_read_perfection, meth_data_to_plot)
  
  pos_min = min(meth_data_to_plot$position_binned)/pos_binning;
  pos_max = max(meth_data_to_plot$position_binned)/pos_binning;
  
  calc_meth_data_matrix_output <- list(meth_data_matrix = meth_data_matrix, reads_to_plot = reads_to_plot, pos_min = pos_min, pos_max = pos_max, position_offset = position_offset)
}

calc_agreement_with_perfect_read <- function(meth_data, meth_data_perfect) {
  meth_data <- meth_data[, c("read_name", "position", "hit_state")]
  meth_data$read_name <- as.character(meth_data$read_name)  # split crashes if this is left a factor
  meth_data_read_list <- split(meth_data, meth_data$read_name)
  read_agreement_list <- list()
  read_name_list <- list()
  all_sites <- 1:nrow(meth_data_perfect)
  for(i in 1:length(meth_data_read_list)) {
    common_positions <- intersect(meth_data_read_list[[i]]$position, meth_data_perfect$position)
    if(length(common_positions) == 0) {
      stop("No common site positions found. Check meth_data_perfect. Did you change the 601 25-mer reference genome?\n")
    }
    num_matches <- meth_data_read_list[[i]][meth_data_read_list[[i]]$position %in% common_positions, "hit_state"] == meth_data_perfect[meth_data_perfect$position %in% common_positions, "hit_state"]
    num_matches <- length(num_matches[num_matches])
    common_site_numbers <- all_sites[meth_data_perfect$position %in% common_positions]
    read_agreement_list[[i]] <- num_matches / (max(common_site_numbers) - min(common_site_numbers) + 1)  # there can be missing sites in the meth_data of a read, this way also count missed sites
    read_name_list[[i]] <- meth_data_read_list[[i]][1, "read_name"]
  }
  
  read_agreement <- data.frame(read_name = unlist(read_name_list), agreement_with_perfect_read = unlist(read_agreement_list), stringsAsFactors = FALSE)
}


calc_read_infos_by_sample <- function(read_infos) {
  if(nrow(read_infos)==0) return()
  read_infos_by_sample <- aggregate(read_infos[, c("aligned_reference_length", "quality_average")], list(read_infos$sample), mean, na.rm = TRUE)
  read_infos_by_sample <- merge(read_infos_by_sample, aggregate(cbind(rep(1, length(read_infos$sample)), read_infos$aligned_reference_length/1000000), list(read_infos$sample), sum), sort = FALSE)
  names(read_infos_by_sample) <- c("sample", "av. length", "av. quality", "reads", "Mbases")
  read_infos_by_sample$reference_genomes <- sapply(read_infos_by_sample$sample, function(sample) paste(unique(sub("\\|.*$", "", as.character(read_infos[read_infos$sample == sample, "ref_seq_detailed"]))), collapse=", "))
  
  return(read_infos_by_sample[c("sample", "reference_genomes", "Mbases", "reads", "av. length", "av. quality")])
}


calc_read_stats_by_sample <- function(read_stats) {
  if(nrow(read_stats)==0) return()
  read_stats$over_meth <- with(read_stats, ifelse(aligned_reference_length > 1000 & group_fraction_ambiguous < 0.33, group_call_fraction_methylated > 0.9, NA))
  read_stats$under_meth <- with(read_stats, ifelse(aligned_reference_length > 1000 & group_fraction_ambiguous < 0.33, group_call_fraction_methylated < 0.1, NA))
  read_stats_by_sample <- aggregate(read_stats[, c("group_fraction_ambiguous", "group_call_fraction_methylated", "aligned_reference_length", "quality_average", "over_meth", "under_meth")], list(read_stats$sample), mean, na.rm = TRUE)
  read_stats_by_sample <- merge(read_stats_by_sample, aggregate(cbind(rep(1, length(read_stats$sample)), read_stats$aligned_reference_length/1000000), list(read_stats$sample), sum), sort = FALSE)
  names(read_stats_by_sample) <- c("sample", "av. ambiguity", "av. methylation", "av. length", "av. quality", "over-meth. frac.", "under-meth. frac.", "reads", "Mbases")
  read_stats_by_sample[c("sample", "av. methylation", "over-meth. frac.", "under-meth. frac.", "av. ambiguity")]
}


calc_site_stats_by_sample <- function(site_stats) {
  if(nrow(site_stats) > 0) {
    site_stats_by_sample <- aggregate(site_stats[, c("hits_called", "hit_fraction_ambiguous", "hit_call_fraction_methylated", "group_size")], list(site_stats$sample), mean, na.rm = TRUE)
    names(site_stats_by_sample) <- c("sample", "av. hits_called", "av. ambiguity", "av. methylation", "av. group_size")
    return(site_stats_by_sample)
  } else {
    return()
  }
}


calc_site_stats_by_pool_label <- function(site_stats) {
  if(!is.null(nrow(site_stats)) && nrow(site_stats) > 0 && !is.na(site_stats$pool_label[1])) {
    site_stats_by_sample <- aggregate(site_stats[, c("hits_called", "hit_fraction_ambiguous", "hit_call_fraction_methylated", "group_size")], list(site_stats$pool_label), mean, na.rm = TRUE)
    names(site_stats_by_sample) <- c("pool_label", "av. hits_called", "av. ambiguity", "av. methylation", "av. group_size")
    return(site_stats_by_sample)
  } else {
    warning("Input has zero rows or is NA!\n")
    return(data.frame("pool_label" = NA, "av. hits_called" = NA, "av. ambiguity" = NA, "av. methylation" = NA, "av. group_size" = NA))
  }
}


calc_perfect_601_25mer_meth_data <- function(motif = "CG") {
  
  sequence_file <- "../../../reference_genomes/Daan_pFMP233/601_25mer_601_linkers.fa"
  sequence <- gsub("[\r\n]", "", readChar(sequence_file, file.info(sequence_file)$size))
  site_positions <- as.integer(gregexpr(motif, sequence)[[1]])
  
  site_starts <- site_positions[1]
  site_ends <- integer(0)
  group_sizes <- integer(0)
  group_size <- 1
  
  # find CG groups with distance <= 10 bases
  for(i in 1:(length(site_positions)-1)) {
    if(group_size == 1 && site_positions[i+1] - site_positions[i] <= 10) {
      group_size <- group_size + 1
    } else if(group_size > 1 && site_positions[i+1] - site_positions[i] <= 10) {
      group_size <- group_size + 1
    } else if(group_size > 1 && site_positions[i+1] - site_positions[i] > 10) {
      site_ends <- c(site_ends, site_positions[i])
      group_sizes <- c(group_sizes, group_size)
      site_starts <- c(site_starts,  site_positions[i+1])
      group_size <- 1
    } else if(group_size == 1 && site_positions[i+1] - site_positions[i] > 10) {
      site_ends <- c(site_ends, site_positions[i])
      group_sizes <- c(group_sizes, group_size)
      site_starts <- c(site_starts,  site_positions[i+1])
    }
  }
  site_ends <- c(site_ends, site_positions[length(site_positions)])
  group_sizes <- c(group_sizes, group_size)
  
  first_601_start <- 31
  length_601 <- 147
  length_linker <- 50
  linker_positions <- integer(0)
  nucl_positions <- integer(0)
  for(i in 0:24) {
    linker_positions <- c(linker_positions, i*(length_601+length_linker)+length_601+(1:length_linker))
    nucl_positions <- c(nucl_positions, i*(length_601+length_linker)+(1:length_601))
  }
  linker_positions <- c(1:(first_601_start-1), linker_positions + first_601_start-1)
  nucl_positions <- nucl_positions + first_601_start-1
  
  perfect_meth_data <- data.frame(start = site_starts, end = site_ends, position = (site_ends + site_starts)/2, group_size = group_sizes)
  
  perfect_meth_data$hit_state <- 1*(perfect_meth_data$start %in% linker_positions & perfect_meth_data$end %in% linker_positions) -
                            1*(perfect_meth_data$start %in% nucl_positions & perfect_meth_data$end %in% nucl_positions)

  if(any(perfect_meth_data$hit_state == 0)) {
    warning("Perfect read has groups combining linker and nucleosome regions")
  }
  perfect_meth_data$read_name = "perfect_601_25mer_read"
  perfect_meth_data$ref_seq = "pFMP233|601_25mer"
  
  perfect_meth_data
}


calc_autocorrelation_from_meth_data_for_sample <- function(meth_data, max_lag = 500) {

  meth_data <- meth_data[, c("read_name", "position", "hit_state")]
  meth_data$read_name <- as.character(meth_data$read_name)  # split crashes if this is left a factor
  meth_data_read_list <- split(meth_data, meth_data$read_name)
  
  meth_vec_list <- list()
  for(read_name in names(meth_data_read_list)) {
    meth_data_read <- meth_data_read_list[[read_name]]
    meth_vec_list[[read_name]] <- rep(NA, max(meth_data_read$position) - min(meth_data_read$position))
    
    # do not use the first site, but everything directly after it, so that last postition of current read and first position of next read won't come directly next to each other after unlisting
    meth_vec_list[[read_name]][round(meth_data_read$position[-1]) - min(meth_data_read$position)] <- meth_data_read$hit_state[-1]
  }
  meth_vec <- unlist(meth_vec_list, use.names = FALSE)
  
  return(acf(meth_vec, na.action=na.pass, lag.max=max_lag, plot=FALSE)[["acf"]])
}

calc_autocorrelation_from_meth_data <- function(meth_data, max_lag = 500) {
  
  autocorrelation_list <- list()
  for(sample in unique(meth_data$sample)) {
    autocorrelation_temp <- calc_autocorrelation_from_meth_data_for_sample(meth_data[meth_data$sample == sample, ], max_lag = max_lag)
    autocorrelation_list[[sample]] <- data.frame(lag = 0:max_lag, autocorrelation = autocorrelation_temp)
    autocorrelation_list[[sample]]$sample <- sample
  }
  
  meth_data_autocorrelation <- bind_rows(autocorrelation_list)
  rownames(meth_data_autocorrelation) <- NULL
  if(nrow(meth_data_autocorrelation) > 0) {
    meth_data_autocorrelation$sample <- factor(meth_data_autocorrelation$sample, levels=levels(meth_data$sample))
  }
  
  return(meth_data_autocorrelation)
}


calc_conditional_methylation_probabilities <- function(meth_data, max_lag=500, dist_binning=20) {
  count_array_list <- count_state_distances(meth_data, max_lag, dist_binning)
  prob_array_list <- list()
  
  for(s in names(count_array_list)) {
    prob_array_list[[s]] <- count_array_list[[s]] / aperm(apply(count_array_list[[s]], c(1,3), function(x) rep(sum(x), 3)), c(2,1,3))
  }
  
  return(prob_array_list)
}


count_state_distances <- function(meth_data, max_lag=500, dist_binning=20) {
  # for each sample: counts binned distances of all state combinations for each read in meth_data and sums them
  
  count_array_list <- list()
  seq_temp <- seq(dist_binning, max_lag, dist_binning)

  for(s in unique(meth_data$sample)) {

    meth_data_s <- subset(meth_data, sample == s)
    
    meth_data_s$read_name <- as.character(meth_data_s$read_name)  # split crashes if this is left a factor
    meth_data_read_list <- split(meth_data_s, meth_data_s$read_name)
    
    count_array <- array(0, dim=c(2*length(seq_temp),3,3), dimnames=list("distance" = c(-rev(seq_temp), seq_temp), "state" = c(-1,0,1), "given_state" = c(-1,0,1)))    
    for(read_name in names(meth_data_read_list)) {
      count_array <- count_array + count_state_distances_for_read(meth_data_read_list[[read_name]], max_lag, dist_binning) 
    }
    
    count_array_list[[s]] <- count_array
  }
  
  return(count_array_list)
}


count_state_distances_for_read <- function(read_meth_data, max_lag, dist_binning) {
  # counts binned distances of all state combinations in given read_meth_data
  
  seq_temp <- seq(dist_binning, max_lag, dist_binning)
  count_array <- array(0, dim=c(2*length(seq_temp),3,3), dimnames=list("distance" = c(-rev(seq_temp), seq_temp), "state" = c(-1,0,1), "given_state" = c(-1,0,1)))    
  
  for(state in dimnames(count_array)[[2]]) {
    for(given_state in dimnames(count_array)[[3]]) {
      x <- read_meth_data[read_meth_data$hit_state == given_state, 'position']
      y <- read_meth_data[read_meth_data$hit_state == state, 'position']
      if(length(x)>0 && length(y)>0) {
        count_array[, state, given_state] <- count_distances(x, y, max_lag, dist_binning)
      }
    }
  }
  
  return(count_array)
}

count_distances <- function(x, y, max_lag, dist_binning) {  
  # counts binned distances from points in x to points in y between -max_lag to +max_lag excluding 0
  # example: dist_binning=10, max_lag=20: bins are [-20, -10), [-10, 0), (0, 10], (10, 20] with bin names -20, -10, 10, 20
  
  A <- matrix(rep(x, length(y)), nrow = length(x))
  B <- t(matrix(rep(y, length(x)), nrow = length(y)))
  D <- B-A
  D <- ceiling(sign(D) * D / dist_binning) * dist_binning * sign(D)
  d <- D[D != 0 & abs(D) <= max_lag]
  
  seq_temp <- seq(dist_binning, max_lag, dist_binning)
  counts <- rep(0, 2*length(seq_temp))
  names(counts) <- c(-rev(seq_temp), seq_temp)
  
  t <- table(d)
  counts[names(t)] <- t
  
  return(counts)
}


calc_object_sizes <- function() {  # gives object sizes in MB of all object in workspace
  object_sizes <- vector()
  for (var in ls(envir=globalenv())) {
    object_sizes[var] <- object.size(get(var))
  }
  object_sizes <- sort(round(object_sizes / 1000000), decreasing = TRUE)
  object_sizes[object_sizes>1]
}


calc_local_density_maxima <- function(vec, xlim) {
  if(!missing(xlim)) {
    d <- xlim[2] - xlim[1]
    vec <- vec[vec >= xlim[1] - 0.25 * d & vec <= xlim[2] + 0.25 * d]
  }
  dens <- density(vec)  # include points outside of xlim into density calculation (it shouldn't immediately go to zero outside)
  maxima <- dens$x[which(diff(sign(diff(dens$y))) == -2) + 1]
  return(maxima[maxima >= xlim[1] & maxima <= xlim[2]])
}


calc_local_density_minima <- function(vec, xlim) {
  if(!missing(xlim)) {
    d <- xlim[2] - xlim[1]
    vec <- vec[vec >= xlim[1] - 0.25 * d & vec <= xlim[2] + 0.25 * d]
  }
  dens <- density(vec)  # include points outside of xlim into density calculation (it shouldn't immediately go to zero outside)
  minima <- dens$x[which(diff(sign(diff(dens$y))) == 2) + 1]
  return(minima[minima >= xlim[1] & minima <= xlim[2]])
}
